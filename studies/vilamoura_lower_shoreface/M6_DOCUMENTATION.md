# M6 — Tracer Centre-of-Mass Uncertainty Analysis

**Vilamoura lower-shoreface tracer experiment — reproducible, study-specific bootstrap workflow**

This case-study module adapts the author's original `M_7 - tcom_uncertainty_statistics.py` (January 2026) for publication as **M6**. The bootstrap estimators, treatment of invalid replicates, confidence-ellipse calculations, and figure construction remain unchanged. The revised entry point accepts command-line arguments for input, output folder, bootstrap count, and random seed. This is a *renumbering and packaging*, not a scientific-method change.

## Scientific approach

For each campaign **C2–C7**, tracer-weighted centre-of-mass coordinates are calculated from the sample stations using `TMcal = Ac / m`. The script processes four **supplied** coordinate columns: `x`, `y`, `as`, and `cs`. The rotated coordinates must already be in the input; M6 does not derive them from `x` and `y`.

- Sample stations are resampled **with replacement**, using the original number of stations per bootstrap replicate.
- Each resampled station receives Gaussian positional perturbations. Defaults are **σ = 2.5 m** for each of `x`, `y`, `as`, `cs`. These represent *assumed* coordinate uncertainty, not measured GNSS errors.
- The original workflow generates **20,000 replicates per campaign**, with base seed **42**.
- Replicates with zero total tracer weight have an undefined centre of mass and are dropped from percentile summaries, saved bootstrap clouds and covariance calculations. Requested, valid and invalid replicate counts are written to `BOOT_INVALID_COUNTS`.
- Percentile intervals (2.5% and 97.5%) are computed separately for each coordinate. Covariance-based **95% ellipses** are derived from the valid bootstrap cloud using the χ²₂,0.95 factor of 5.991; these are **descriptive ellipses of the bootstrap distribution**, rather than direct hypothesis tests for significant displacement between campaigns.
- The Kish effective sample size (`n_eff`) is reported as a measure of concentration of normalized tracer weights. It is not equivalent to the number of sampling stations.

**Interpretation:** The script quantifies positional uncertainty of individual campaign TCOM estimates. It does not automatically calculate a joint bootstrap confidence interval for *differences* between two campaigns. Overlapping ellipses are suggestive but not a formal displacement-significance test. Gaussian coordinate perturbations in `x/y` and `as/cs` are drawn independently in the original code; the two simulated coordinate systems are therefore not strictly coupled by a rotation at each bootstrap replicate.

## Public files

- `M6_tcom_uncertainty_bootstrap.py`: scientific analysis and plotting routine.
- `M6_DOCUMENTATION.md`: this methodological README.
- `examples/INPUT_TCOM_SYNTHETIC.xlsx`: invented sample geometry and tracer weights for all six campaigns.

The synthetic data use arbitrary coordinates and are **not observational Vilamoura records or manuscript results**. The full private experimental workbook, real bootstrap output, original figures and ellipse vertices are **not** distributed in this repository.

## Input format

An `.xlsx` file with sheets **C2, C3, C4, C5, C6, C7**, each containing:

| Column | Meaning |
| --- | --- |
| `x`, `y` | Projected east / north positions (m) |
| `as`, `cs` | Along-shelf / cross-shelf sample coordinates (m), already computed externally |
| `Ac` (or `Ac `) | Image-derived fluorescent tracer area fraction |
| `m` | Calibration slope; must be positive for nonzero weights |

Each campaign needs sufficient positive tracer weights for valid bootstrap centres and covariance estimation. No campaign source coordinates or published results have been embedded in this case-study release.

## Installation and usage

Dependencies: `numpy`, `pandas`, `matplotlib`, `openpyxl`.

Run from the case-study directory:

```bash
python -m pip install numpy pandas matplotlib openpyxl
python M6_tcom_uncertainty_bootstrap.py
```

By default, this reads `examples/INPUT_TCOM_SYNTHETIC.xlsx`, uses 20,000 bootstrap replicates and seed 42, and writes into `outputs/`.

For a **rapid synthetic test**, use fewer replicates (this is not the manuscript configuration):

```bash
python M6_tcom_uncertainty_bootstrap.py --bootstrap 120 --output-dir outputs/test
```

For **private scientific analysis**, supply the actual workbook and use the original bootstrap configuration:

```bash
python M6_tcom_uncertainty_bootstrap.py --input "C:/private/input_TRACER_DATA_delta_mix_p95.xlsx" --bootstrap 20000 --seed 42 --output-dir "C:/private/M6_RESULTS"
```

In Spyder, the script can be run without arguments using its synthetic default. For the observational workbook, configure command-line arguments in the Run settings or call the file from the console with arguments. Do **not** place a copy of the experimental workbook in the public `examples` folder.

## Output products

1. `TRACER_CM_uncertainty_bootstrap.xlsx`: `PARAMS`, `BOOT_INVALID_COUNTS`, `ORIGINAL_CM`, `BOOT_MEDIAN_CM`, `COM_SUMMARY`, `ELLIPSE95_XY`, `ELLIPSE95_ASCS`, and `BOOT_C2`–`BOOT_C7`.
2. `Combined_geometry_uncertainty_neff_C2-C4_relativeAxes.png` and `Combined_geometry_uncertainty_neff_C5-C7_relativeAxes.png`: sampling geometry, bootstrap ellipse and effective-sample-size graphics.
3. `CM_uncertainty_ellipses_95_vertices_by_campaign.xlsx`: 180 ellipse vertices per campaign, intended for GIS import, plus a metadata worksheet.

For the map-ready vertex export the source routine labels the reference system `ETRS89 / UTM zone 29N (EPSG:25829)`. Verify the coordinate reference system of any alternative workbook before using these outputs in GIS.

## Validation and availability

The case-study code was tested using its synthetic workbook. In an equivalence test with **120 replicates**, both the original M7 and renamed M6 produced **identical** numerical/worksheet outputs (13 bootstrap workbook sheets and 7 GIS-vertices workbook sheets), excluding expected differences in absolute file-path metadata. The manuscript setting of **20,000** replicates remains unchanged as the default and should be used for formal analysis.

The complete observational input and resulting scientific output are intended for release separately after journal publication, subject to the journal's requirements.

**Author:** João Cascalho, FCUL / Instituto Dom Luiz, University of Lisbon, Portugal.
