"""
Bootstrap uncertainty for sediment-tracer transport metrics.

This routine estimates uncertainty in tracer centre-of-mass positions and
derived transport metrics using non-parametric bootstrap resampling of sampling
stations. The public version uses only Cartesian x/y coordinates. Any
along-body or cross-body coordinate calculations have been removed.

Expected workbook structure
---------------------------
The input workbook must contain:

- IP: initial injection/reference position
- one or more campaign sheets, for example C1, C2, C3, C4, C5

Required columns in IP:

- x
- y
- DATE_TIME or End_date

Required columns in campaign sheets:

- x
- y
- Ac or Ac with trailing whitespace
- m
- DATE_TIME or End_date

Optional columns in campaign sheets:

- sample or SAMPLE_ID
- delta_mix

Outputs
-------
- Bootstrap summary workbook
- Ellipse-vertex workbook for GIS import
- Multi-panel diagnostic figures
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Config:
    input_file: Path
    output_workbook: Path
    vertices_workbook: Path
    plots_dir: Path
    campaigns: tuple[str, ...]
    n_bootstrap: int
    seed: int
    sigma_x_m: float
    sigma_y_m: float
    n_vertices: int


CAMPAIGN_COLORS = {
    "C1": (255 / 255, 168 / 255, 0 / 255),
    "C2": (102 / 255, 51 / 255, 0 / 255),
    "C3": (22 / 255, 182 / 255, 78 / 255),
    "C4": (0 / 255, 150 / 255, 180 / 255),
    "C5": (220 / 255, 30 / 255, 30 / 255),
}


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(c).strip() for c in out.columns]
    return out


def get_ac_column(df: pd.DataFrame) -> str:
    for col in df.columns:
        if str(col).strip() == "Ac":
            return col
    raise KeyError("Missing 'Ac' column.")


def get_date_column(df: pd.DataFrame, sheet_name: str) -> str:
    for col in ["End_date", "DATE_TIME", "Date", "Datetime", "date/time"]:
        if col in df.columns:
            return col
    raise KeyError(f"Missing DATE_TIME/End_date column in sheet '{sheet_name}'.")


def standardize_campaign_df(df: pd.DataFrame, sheet_name: str, require_weights: bool) -> pd.DataFrame:
    out = normalize_columns(df)
    out = out.rename(columns={"SAMPLE_ID": "sample", "DATE_TIME": "End_date"})

    required = ["x", "y"]
    if require_weights:
        required.extend([get_ac_column(out), "m"])

    missing = [col for col in required if col not in out.columns]
    if missing:
        raise KeyError(f"Missing required columns {missing} in sheet '{sheet_name}'.")

    date_col = get_date_column(out, sheet_name)
    if date_col != "End_date":
        out = out.rename(columns={date_col: "End_date"})

    out["x"] = pd.to_numeric(out["x"], errors="coerce")
    out["y"] = pd.to_numeric(out["y"], errors="coerce")
    out["End_date"] = pd.to_datetime(out["End_date"], errors="coerce")

    if require_weights:
        ac_col = get_ac_column(out)
        out[ac_col] = pd.to_numeric(out[ac_col], errors="coerce")
        out["m"] = pd.to_numeric(out["m"], errors="coerce")
        out = out.dropna(subset=["x", "y", "End_date", ac_col, "m"]).reset_index(drop=True)
    else:
        out = out.dropna(subset=["x", "y", "End_date"]).reset_index(drop=True)

    if out.empty:
        raise ValueError(f"Sheet '{sheet_name}' has no valid rows after quality control.")
    return out


def get_first_valid_datetime(df: pd.DataFrame, sheet_name: str) -> pd.Timestamp:
    dt = df["End_date"].dropna()
    if dt.empty:
        raise ValueError(f"Could not parse any valid date in sheet '{sheet_name}'.")
    return pd.Timestamp(dt.iloc[0])


def get_delta_mix_value(df: pd.DataFrame) -> float:
    for col in ["delta_mix", "delta_mix ", "mixing_layer_thickness"]:
        if col in df.columns:
            values = pd.to_numeric(df[col], errors="coerce").dropna()
            if not values.empty:
                return float(values.iloc[0])
    return float("nan")


def tracer_weights(df: pd.DataFrame) -> np.ndarray:
    ac_col = get_ac_column(df)
    ac = pd.to_numeric(df[ac_col], errors="coerce").to_numpy(dtype=float)
    mass = pd.to_numeric(df["m"], errors="coerce").to_numpy(dtype=float)
    return np.where(mass > 0, ac / mass, 0.0)


def kish_effective_sample_size(weights: np.ndarray) -> float:
    weights = np.asarray(weights, dtype=float)
    total = float(np.nansum(weights))
    if total <= 0 or np.isclose(total, 0.0):
        return float("nan")
    p = weights / total
    denom = float(np.nansum(p**2))
    return 1.0 / denom if denom > 0 else float("nan")


def weighted_centre_of_mass(df: pd.DataFrame, weights: np.ndarray) -> dict[str, float]:
    total = float(np.nansum(weights))
    if total <= 0 or np.isclose(total, 0.0):
        return {"CMx": float("nan"), "CMy": float("nan"), "TM_sum": total}

    x = df["x"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=float)
    return {
        "CMx": float(np.nansum(weights * x) / total),
        "CMy": float(np.nansum(weights * y) / total),
        "TM_sum": total,
    }


def bootstrap_metrics_cloud(
    df: pd.DataFrame,
    n_bootstrap: int,
    seed: int,
    sigma_x_m: float,
    sigma_y_m: float,
    previous_cmx: float,
    previous_cmy: float,
    delta_t_seconds: float,
    delta_mix: float,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    base_weights = tracer_weights(df)
    x0 = df["x"].to_numpy(dtype=float)
    y0 = df["y"].to_numpy(dtype=float)
    n = len(df)
    rows: list[dict[str, float]] = []

    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        weights = base_weights[idx]
        total = float(np.nansum(weights))

        if total <= 0 or np.isclose(total, 0.0):
            rows.append(
                {
                    "CMx": np.nan,
                    "CMy": np.nan,
                    "TCOMd_xy": np.nan,
                    "D_x": np.nan,
                    "D_y": np.nan,
                    "D_total_xy": np.nan,
                    "Vm_x": np.nan,
                    "Vm_y": np.nan,
                    "Q_x": np.nan,
                    "Q_y": np.nan,
                    "Q_xy_resultant": np.nan,
                    "TM_sum": total,
                }
            )
            continue

        x = x0[idx] + rng.normal(0, sigma_x_m, size=n)
        y = y0[idx] + rng.normal(0, sigma_y_m, size=n)
        cmx = float(np.nansum(weights * x) / total)
        cmy = float(np.nansum(weights * y) / total)

        dist_x = x - previous_cmx
        dist_y = y - previous_cmy
        weighted_square_sum = float(np.nansum(weights**2))

        if delta_t_seconds > 0 and weighted_square_sum > 0:
            d_x = float(np.nansum((weights * dist_x) ** 2) / (weighted_square_sum * delta_t_seconds))
            d_y = float(np.nansum((weights * dist_y) ** 2) / (weighted_square_sum * delta_t_seconds))
        else:
            d_x = d_y = np.nan

        dcm_x = cmx - previous_cmx
        dcm_y = cmy - previous_cmy
        if delta_t_seconds > 0:
            vm_x = dcm_x / delta_t_seconds
            vm_y = dcm_y / delta_t_seconds
        else:
            vm_x = vm_y = np.nan

        if np.isfinite(delta_mix):
            q_x = vm_x * delta_mix
            q_y = vm_y * delta_mix
        else:
            q_x = q_y = np.nan

        rows.append(
            {
                "CMx": cmx,
                "CMy": cmy,
                "TCOMd_xy": math.sqrt(dcm_x**2 + dcm_y**2),
                "D_x": d_x,
                "D_y": d_y,
                "D_total_xy": (d_x + d_y) if np.isfinite(d_x) and np.isfinite(d_y) else np.nan,
                "Vm_x": vm_x,
                "Vm_y": vm_y,
                "Q_x": q_x,
                "Q_y": q_y,
                "Q_xy_resultant": math.sqrt(q_x**2 + q_y**2) if np.isfinite(q_x) and np.isfinite(q_y) else np.nan,
                "TM_sum": total,
            }
        )

    return pd.DataFrame(rows)


def summary_stats(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    return {
        "Median": float(np.nanmedian(values)),
        "Mean": float(np.nanmean(values)),
        "Std": float(np.nanstd(values, ddof=1)),
        "CI2p5": float(np.nanpercentile(values, 2.5)),
        "CI97p5": float(np.nanpercentile(values, 97.5)),
    }


def ellipse_params(x: np.ndarray, y: np.ndarray, level: float = 0.95) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if len(x) < 3:
        return {"a": np.nan, "b": np.nan, "angle_deg": np.nan}

    cov = np.cov(x, y)
    values, vectors = np.linalg.eigh(cov)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]

    chi2_value = 5.991 if abs(level - 0.95) < 1e-12 else 2.279
    return {
        "a": float(np.sqrt(values[0] * chi2_value)),
        "b": float(np.sqrt(values[1] * chi2_value)),
        "angle_deg": float(np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0]))),
    }


def ellipse_vertices(cx: float, cy: float, a: float, b: float, angle_deg: float, n_vertices: int) -> np.ndarray:
    t = np.linspace(0, 2 * np.pi, n_vertices, endpoint=False)
    x = a * np.cos(t)
    y = b * np.sin(t)
    angle = np.deg2rad(angle_deg)
    xr = x * np.cos(angle) - y * np.sin(angle)
    yr = x * np.sin(angle) + y * np.cos(angle)
    return np.column_stack([cx + xr, cy + yr])


def ellipse_points(cx: float, cy: float, a: float, b: float, angle_deg: float, n: int = 240):
    if not np.isfinite(a) or not np.isfinite(b) or a <= 0 or b <= 0:
        return None
    vertices = ellipse_vertices(cx, cy, a, b, angle_deg, n)
    closed = np.vstack([vertices, vertices[0]])
    return closed[:, 0], closed[:, 1]


def compute_bootstrap_workbook(cfg: Config) -> tuple[pd.ExcelFile, pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    if not cfg.input_file.exists():
        raise FileNotFoundError(f"Input workbook not found: {cfg.input_file}")

    cfg.output_workbook.parent.mkdir(parents=True, exist_ok=True)
    xls = pd.ExcelFile(cfg.input_file)
    if "IP" not in xls.sheet_names:
        raise ValueError("Sheet 'IP' not found. It is required as the reference position.")

    ip_df = standardize_campaign_df(pd.read_excel(xls, sheet_name="IP"), "IP", require_weights=False)
    reference_time = get_first_valid_datetime(ip_df, "IP")
    previous_cmx = float(ip_df.iloc[0]["x"])
    previous_cmy = float(ip_df.iloc[0]["y"])

    original_rows: list[dict[str, float | str]] = []
    summary_rows: list[dict[str, float | str]] = []
    ellipse_rows: list[dict[str, float | str]] = []
    invalid_rows: list[dict[str, float | str]] = []
    boot_clouds: dict[str, pd.DataFrame] = {}

    rng = np.random.default_rng(cfg.seed)
    metrics = ["CMx", "CMy", "TCOMd_xy", "D_x", "D_y", "D_total_xy", "Vm_x", "Vm_y", "Q_x", "Q_y", "Q_xy_resultant"]

    for campaign in cfg.campaigns:
        if campaign not in xls.sheet_names:
            logging.warning("Campaign sheet '%s' not found. Skipping.", campaign)
            continue

        df = standardize_campaign_df(pd.read_excel(xls, sheet_name=campaign), campaign, require_weights=True)
        weights = tracer_weights(df)
        n_eff = kish_effective_sample_size(weights)
        cm = weighted_centre_of_mass(df, weights)
        campaign_time = get_first_valid_datetime(df, campaign)
        delta_t_seconds = float((campaign_time - reference_time).total_seconds())
        delta_mix = get_delta_mix_value(df)

        dcm_x = cm["CMx"] - previous_cmx if np.isfinite(cm["CMx"]) else np.nan
        dcm_y = cm["CMy"] - previous_cmy if np.isfinite(cm["CMy"]) else np.nan
        original_rows.append(
            {
                "Campaign": campaign,
                "n_eff": n_eff,
                "delta_t_s": delta_t_seconds,
                "delta_mix": delta_mix,
                **cm,
                "TCOMd_xy": math.sqrt(dcm_x**2 + dcm_y**2) if np.isfinite(dcm_x) and np.isfinite(dcm_y) else np.nan,
            }
        )

        seed = int(rng.integers(0, 2**31 - 1))
        boot_raw = bootstrap_metrics_cloud(
            df=df,
            n_bootstrap=cfg.n_bootstrap,
            seed=seed,
            sigma_x_m=cfg.sigma_x_m,
            sigma_y_m=cfg.sigma_y_m,
            previous_cmx=previous_cmx,
            previous_cmy=previous_cmy,
            delta_t_seconds=delta_t_seconds,
            delta_mix=delta_mix,
        )
        n_total = len(boot_raw)
        n_invalid = int(boot_raw["CMx"].isna().sum())
        boot = boot_raw.dropna(subset=["CMx", "CMy"]).reset_index(drop=True)
        boot_clouds[campaign] = boot

        invalid_rows.append(
            {
                "Campaign": campaign,
                "B_requested": cfg.n_bootstrap,
                "B_total_rows": n_total,
                "B_invalid": n_invalid,
                "B_valid": len(boot),
                "p_invalid": (n_invalid / n_total) if n_total > 0 else np.nan,
            }
        )

        for metric in metrics:
            summary_rows.append(
                {
                    "Campaign": campaign,
                    "Metric": metric,
                    "n_eff": n_eff,
                    "delta_t_s": delta_t_seconds,
                    "delta_mix": delta_mix,
                    "B_valid": len(boot),
                    **summary_stats(boot[metric].to_numpy(dtype=float)),
                }
            )

        ellipse_rows.append(
            {
                "Campaign": campaign,
                "level": 0.95,
                "B_valid": len(boot),
                **ellipse_params(boot["CMx"], boot["CMy"], level=0.95),
            }
        )

        if np.isfinite(cm["CMx"]) and np.isfinite(cm["CMy"]):
            reference_time = campaign_time
            previous_cmx = cm["CMx"]
            previous_cmy = cm["CMy"]

        logging.info("%s: %s valid bootstrap replicates from %s requested.", campaign, len(boot), cfg.n_bootstrap)

    original_df = pd.DataFrame(original_rows)
    summary_df = pd.DataFrame(summary_rows)
    ellipse_df = pd.DataFrame(ellipse_rows)
    invalid_df = pd.DataFrame(invalid_rows)

    params_df = pd.DataFrame(
        [
            {
                "input_file": str(cfg.input_file),
                "B_requested": cfg.n_bootstrap,
                "seed": cfg.seed,
                "sigma_x_m": cfg.sigma_x_m,
                "sigma_y_m": cfg.sigma_y_m,
                "coordinate_system": "Cartesian x/y only",
                "note": "Weights are Ac/m. Bootstrap replicates with zero total weight are removed before summaries and ellipse fitting.",
            }
        ]
    )

    wide_summaries = {
        "CM_TCOM_SUMMARY": ["CMx", "CMy", "TCOMd_xy"],
        "DIFFUSION_SUMMARY": ["D_x", "D_y", "D_total_xy"],
        "VELOCITY_SUMMARY": ["Vm_x", "Vm_y"],
        "Q_TRANSPORT_SUMMARY": ["Q_x", "Q_y", "Q_xy_resultant"],
    }

    with pd.ExcelWriter(cfg.output_workbook, engine="openpyxl") as writer:
        params_df.to_excel(writer, sheet_name="PARAMS", index=False)
        invalid_df.to_excel(writer, sheet_name="BOOT_INVALID_COUNTS", index=False)
        original_df.to_excel(writer, sheet_name="ORIGINAL_METRICS", index=False)
        summary_df.to_excel(writer, sheet_name="BOOT_SUMMARY_LONG", index=False)
        ellipse_df.to_excel(writer, sheet_name="ELLIPSE95_XY", index=False)

        for sheet_name, sheet_metrics in wide_summaries.items():
            subset = summary_df[summary_df["Metric"].isin(sheet_metrics)]
            if subset.empty:
                continue
            wide = subset.pivot(index="Campaign", columns="Metric", values=["Median", "Mean", "Std", "CI2p5", "CI97p5"])
            wide.columns = [f"{stat}_{metric}" for stat, metric in wide.columns]
            wide.reset_index().to_excel(writer, sheet_name=sheet_name, index=False)

        for campaign, boot in boot_clouds.items():
            boot.to_excel(writer, sheet_name=f"BOOT_{campaign}"[:31], index=False)

    logging.info("Saved bootstrap workbook: %s", cfg.output_workbook)
    return xls, original_df, ellipse_df, boot_clouds


def export_vertices(cfg: Config, original_df: pd.DataFrame, ellipse_df: pd.DataFrame) -> None:
    cfg.vertices_workbook.parent.mkdir(parents=True, exist_ok=True)
    merged = original_df.merge(ellipse_df, on="Campaign", how="inner").set_index("Campaign")
    with pd.ExcelWriter(cfg.vertices_workbook, engine="openpyxl") as writer:
        for campaign in cfg.campaigns:
            if campaign not in merged.index:
                continue
            row = merged.loc[campaign]
            if not all(np.isfinite([row["CMx"], row["CMy"], row["a"], row["b"], row["angle_deg"]])):
                continue
            vertices = ellipse_vertices(row["CMx"], row["CMy"], row["a"], row["b"], row["angle_deg"], cfg.n_vertices)
            pd.DataFrame(
                {
                    "Order": np.arange(len(vertices), dtype=int),
                    "X": vertices[:, 0].astype(float),
                    "Y": vertices[:, 1].astype(float),
                }
            ).to_excel(writer, sheet_name=f"{campaign}_XY"[:31], index=False)

        pd.DataFrame(
            [
                {
                    "source_bootstrap_excel": str(cfg.output_workbook),
                    "n_vertices_per_ellipse": cfg.n_vertices,
                    "coordinate_system": "Cartesian x/y",
                    "notes": "Ellipse vertices computed from campaign centre-of-mass coordinates and 95% covariance-ellipse parameters.",
                }
            ]
        ).to_excel(writer, sheet_name="META", index=False)
    logging.info("Saved ellipse vertices workbook: %s", cfg.vertices_workbook)


def compute_plot_limits(dfs: list[pd.DataFrame]) -> tuple[float, float, float, float, float, float]:
    all_x = np.concatenate([df["x"].to_numpy(dtype=float) for df in dfs])
    all_y = np.concatenate([df["y"].to_numpy(dtype=float) for df in dfs])
    xmin, xmax = float(np.nanmin(all_x) - 5), float(np.nanmax(all_x) + 5)
    ymin, ymax = float(np.nanmin(all_y) - 5), float(np.nanmax(all_y) + 5)
    side = max(xmax - xmin, ymax - ymin)
    xmid, ymid = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
    xmin, xmax = xmid - side / 2, xmid + side / 2
    ymin, ymax = ymid - side / 2, ymid + side / 2
    origin_x = math.floor(xmin)
    origin_y = math.floor(ymin)
    return origin_x, origin_y, xmin, xmax, ymin, ymax


def make_plot_page(
    cfg: Config,
    xls: pd.ExcelFile,
    original_df: pd.DataFrame,
    ellipse_df: pd.DataFrame,
    campaigns: tuple[str, ...],
    page_name: str,
) -> None:
    cfg.plots_dir.mkdir(parents=True, exist_ok=True)
    dfs = [standardize_campaign_df(pd.read_excel(xls, sheet_name=c), c, require_weights=True) for c in campaigns]
    origin_x, origin_y, xmin, xmax, ymin, ymax = compute_plot_limits(dfs)
    fig, axes = plt.subplots(3, len(campaigns), figsize=(4.6 * len(campaigns), 10), gridspec_kw={"height_ratios": [1, 1, 0.6]})

    if len(campaigns) == 1:
        axes = np.asarray(axes).reshape(3, 1)

    for index, (campaign, df) in enumerate(zip(campaigns, dfs)):
        color = CAMPAIGN_COLORS.get(campaign, (0.2, 0.2, 0.2))
        ac_col = get_ac_column(df)
        ac = pd.to_numeric(df[ac_col], errors="coerce").to_numpy(dtype=float)
        has_tracer = ac > 0
        orig = original_df[original_df["Campaign"] == campaign].iloc[0]
        ellipse = ellipse_df[ellipse_df["Campaign"] == campaign].iloc[0]

        cmx = float(orig["CMx"])
        cmy = float(orig["CMy"])
        n_eff = float(orig["n_eff"])

        ax = axes[0, index]
        ax.scatter(df.loc[~has_tracer, "x"] - origin_x, df.loc[~has_tracer, "y"] - origin_y, facecolors="none", edgecolors="0.6", s=26, linewidths=0.8)
        ax.scatter(df.loc[has_tracer, "x"] - origin_x, df.loc[has_tracer, "y"] - origin_y, c="black", s=30)
        ax.scatter(cmx - origin_x, cmy - origin_y, c=[color], s=75, edgecolors="black", linewidths=0.7)
        ax.set_title(campaign)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(xmin - origin_x, xmax - origin_x)
        ax.set_ylim(ymin - origin_y, ymax - origin_y)
        ax.grid(alpha=0.25)

        ax = axes[1, index]
        ax.scatter(df["x"] - origin_x, df["y"] - origin_y, facecolors="white", edgecolors="black", s=18, linewidths=0.55)
        points = ellipse_points(cmx - origin_x, cmy - origin_y, float(ellipse["a"]), float(ellipse["b"]), float(ellipse["angle_deg"]))
        if points is not None:
            ax.plot(points[0], points[1], color="black", lw=1.6)
        ax.scatter(cmx - origin_x, cmy - origin_y, c=[color], s=65, edgecolors="black", linewidths=0.7)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(xmin - origin_x, xmax - origin_x)
        ax.set_ylim(ymin - origin_y, ymax - origin_y)
        ax.grid(alpha=0.25)

        ax = axes[2, index]
        ax.bar([0], [n_eff], color=color, width=0.6)
        ax.set_ylim(0, max(len(df), 1))
        ax.set_xticks([])
        if index == 0:
            ax.set_ylabel("n_eff")

    axes[0, 0].set_ylabel("Sampling geometry")
    axes[1, 0].set_ylabel("95% CM uncertainty")
    fig.supxlabel(f"Delta X (m) from {int(origin_x)}")
    fig.supylabel(f"Delta Y (m) from {int(origin_y)}")
    fig.suptitle("Tracer sampling geometry, centre-of-mass uncertainty, and effective sample size", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = cfg.plots_dir / f"tracer_metrics_uncertainty_{page_name}.png"
    fig.savefig(out, dpi=300, facecolor="white")
    plt.close(fig)
    logging.info("Saved plot: %s", out)


def make_all_plots(cfg: Config, xls: pd.ExcelFile, original_df: pd.DataFrame, ellipse_df: pd.DataFrame) -> None:
    available = tuple(c for c in cfg.campaigns if c in xls.sheet_names)
    for start in range(0, len(available), 3):
        page = available[start : start + 3]
        if not page:
            continue
        make_plot_page(cfg, xls, original_df, ellipse_df, page, f"{page[0]}_{page[-1]}")


def parse_args(argv: list[str] | None = None) -> Config:
    parser = argparse.ArgumentParser(description="Bootstrap uncertainty analysis for sediment-tracer transport metrics.")
    parser.add_argument("-i", "--input", default="examples/example_tracer_uncertainty_input.xlsx", help="Input Excel workbook.")
    parser.add_argument("-o", "--output", default="results/tracer_metrics_uncertainty_results.xlsx", help="Output bootstrap workbook.")
    parser.add_argument("--vertices-output", default="results/tracer_metrics_uncertainty_ellipse_vertices.xlsx", help="Output ellipse vertices workbook.")
    parser.add_argument("--plots-dir", default="results/plots", help="Directory for diagnostic plots.")
    parser.add_argument("--campaigns", default="C1,C2,C3,C4,C5", help="Comma-separated campaign sheet names.")
    parser.add_argument("--n-bootstrap", type=int, default=20000, help="Number of bootstrap realizations.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--sigma-x", type=float, default=2.5, help="Station-position uncertainty in x, in metres.")
    parser.add_argument("--sigma-y", type=float, default=2.5, help="Station-position uncertainty in y, in metres.")
    parser.add_argument("--n-vertices", type=int, default=180, help="Number of vertices exported per uncertainty ellipse.")
    args = parser.parse_args(argv)

    campaigns = tuple(c.strip() for c in args.campaigns.split(",") if c.strip())
    return Config(
        input_file=Path(args.input),
        output_workbook=Path(args.output),
        vertices_workbook=Path(args.vertices_output),
        plots_dir=Path(args.plots_dir),
        campaigns=campaigns,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
        sigma_x_m=args.sigma_x,
        sigma_y_m=args.sigma_y,
        n_vertices=args.n_vertices,
    )


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    cfg = parse_args(argv)
    xls, original_df, ellipse_df, _boot_clouds = compute_bootstrap_workbook(cfg)
    export_vertices(cfg, original_df, ellipse_df)
    make_all_plots(cfg, xls, original_df, ellipse_df)
    logging.info("Completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

