# Tracer Metrics Uncertainty

Python routine for estimating bootstrap uncertainty in sediment-tracer centre-of-mass positions and derived transport metrics.

The workflow uses non-parametric bootstrap resampling of sampling stations, optional station-position uncertainty in Cartesian `x` and `y` coordinates, and tracer weights calculated as `Ac / m`. It exports confidence intervals, 95% centre-of-mass uncertainty ellipses, GIS-ready ellipse vertices, and diagnostic figures.

This public version uses only Cartesian `x/y` coordinates. Along-body and cross-body coordinate calculations are not included.

## Features

- Reads an Excel workbook with an `IP` reference sheet and campaign sheets such as `C1` to `C5`.
- Calculates tracer centre-of-mass uncertainty from bootstrap resampling.
- Estimates confidence intervals for centre-of-mass displacement, diffusion, mean velocity, and transport proxy metrics.
- Exports 95% covariance ellipse parameters and ellipse vertices for GIS import.
- Produces diagnostic figures for sampling geometry, centre-of-mass uncertainty, and effective sample size.
- Includes a small artificial example workbook for testing.

## Installation

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

If you are using Anaconda or Spyder, install the packages in your active environment.

## Input Format

The input workbook must include:

- `IP`: initial injection/reference position
- campaign sheets such as `C1`, `C2`, `C3`, `C4`, and `C5`

Required columns in `IP`:

- `x`
- `y`
- `DATE_TIME` or `End_date`

Required columns in campaign sheets:

- `x`
- `y`
- `Ac`
- `m`
- `DATE_TIME` or `End_date`

Optional columns:

- `sample` or `SAMPLE_ID`
- `delta_mix`

The example file in `examples/example_tracer_uncertainty_input.xlsx` is artificial and is provided only to demonstrate the expected workbook structure.

## Usage

Run the included artificial example:

```bash
python tracer_metrics_uncertainty.py
```

Use fewer bootstrap iterations for a quick test:

```bash
python tracer_metrics_uncertainty.py --n-bootstrap 500
```

Run with your own workbook:

```bash
python tracer_metrics_uncertainty.py --input path/to/your_input.xlsx --output results/tracer_metrics_uncertainty_results.xlsx
```

Select campaign sheets:

```bash
python tracer_metrics_uncertainty.py --campaigns C1,C2,C3
```

## Outputs

Default outputs:

- `results/tracer_metrics_uncertainty_results.xlsx`
- `results/tracer_metrics_uncertainty_ellipse_vertices.xlsx`
- `results/plots/tracer_metrics_uncertainty_C1_C3.png`
- `results/plots/tracer_metrics_uncertainty_C4_C5.png`

The main workbook includes:

- `PARAMS`
- `BOOT_INVALID_COUNTS`
- `ORIGINAL_METRICS`
- `BOOT_SUMMARY_LONG`
- `ELLIPSE95_XY`
- summary sheets for centre-of-mass, diffusion, velocity, and transport proxy metrics
- one bootstrap cloud sheet per campaign

## Notes

Bootstrap realizations with zero total tracer weight are removed before summary statistics and ellipse fitting.

Do not upload unpublished or confidential field datasets to this repository. Use the `examples` folder only for small artificial demonstration files or data explicitly approved for public release.

## Citation

If you use this routine, please cite the archived Zenodo release and the GitHub repository. Citation metadata are provided in `CITATION.cff`.

