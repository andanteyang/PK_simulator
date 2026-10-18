# PK Simulator

A browser-based pharmacokinetic workbench built with native HTML, CSS, and JavaScript ES modules. No build, database, backend, account, CDN, or runtime dependency is required. All calculations run locally. This is a simulation and educational tool, **not a validated clinical dosing application**.

## Run locally

From this repository, run `python3 -m http.server 8000` (or `npm start`). Open the address printed by the server in your browser. Serve the files over HTTP; opening `index.html` via `file://` does not reliably load ES modules. Python is only a static development server, not an application backend.

Select a compartment model, route, and dosing mode. Edit numeric parameters or move the sliders; results update automatically. Numeric fields permit values outside slider ranges. Change the simulation end to include the desired post-dose decline. Invalid settings display an explanation and disable exports. Switch between linear and logarithmic concentration axes, hover for exact model concentrations, and export CSV data or a PNG plot. Mobile devices support the same configuration and exports.

## Models and units

All eight combinations of one/two compartments, IV bolus/oral (PO), and single/multiple dosing are supported.

- Dose: mg/kg; CL and Q: L/h/kg; V, Vc, Vp: L/kg; time: h.
- Concentrations are calculated in mg/L, then multiplied by 1,000 to give ng/mL. AUC is ng·h/mL.
- Linear time-invariant PK; central first-order elimination; PO first-order absorption with constant F and no lag time. Parameters stay constant across doses.

### Equations

One-compartment IV: `C(t) = Dose/V × exp(-kel × t)`, with `kel = CL/V`.

One-compartment PO: `C(t) = F × Dose/V × ka × (exp(-kel*t) - exp(-ka*t))/(ka-kel)`.

Two compartments: `k10 = CL/Vc`, `k12 = Q/Vc`, `k21 = Q/Vp`. Alpha and beta solve `λ² - (k10+k12+k21)λ + k10*k21 = 0`. Beta is calculated from the root product to avoid subtractive cancellation. The eigenvalue discriminant is evaluated as `sqrt((k10+k12-k21)² + 4*k12*k21)`.

IV: `C(t) = Dose/Vc × [a exp(-alpha*t) + b exp(-beta*t)]`, where `a=(alpha-k21)/(alpha-beta)` and `b=(k21-beta)/(alpha-beta)`. Thus A = Dose/Vc × a and B = Dose/Vc × b.

PO convolves each disposition exponential with `ka*exp(-ka*t)` and multiplies by F. The convolution is evaluated using `expm1` and the slower exponential to avoid cancellation and overflow. Equal rates use the exact limit `ka*t*exp(-ka*t)`.

Multiple doses superpose the single-dose functions at `0, τ, …, (N−1)τ`. IV event rows contain distinct pre/post values at the exact same time; the graph draws vertical jumps. PO dosing produces no instantaneous concentration jump. The time axis is continuous and numerical. True zero concentrations are omitted on log plots, never replaced with artificial positive values.

### Results

- Single-dose Cmax and Tmax are numerical estimates of the full single-dose profile, even if the plotting window ends before the peak. Peak search includes 40 slowest time constants.
- Multiple-dose global Cmax and occurrence time apply to the simulated window only.
- First-interval metrics cover `[0, τ]`; final-dose interval metrics cover `[(N−1)τ, Nτ]`. The right endpoint excludes the next dose. Tmax is relative to the relevant dose. If the full requested interval is outside the simulation, its result is N/A; a partial interval is never presented as complete. For N = 1, first and final intervals coincide.
- Peak searches use linear and logarithmic age grids, followed by bracket refinement, rather than simply reporting a plotting-grid maximum.
- Simulated AUC₀–end uses adaptive Simpson integration over a segmented logarithmic/linear grid for each dose. Theoretical cumulative AUC₀–∞ is `1,000 × Dose × N / CL` for IV and additionally multiplied by F for PO; N = 1 for single dosing.
- Disposition terminal half-life is `ln(2)/kel` or `ln(2)/beta`. Apparent PO terminal half-life is `ln(2)/min(ka, kel or beta)`, exposing flip-flop kinetics. These are asymptotic rate descriptors. At zero Dose or F there is no measurable terminal concentration slope.
- CSV preserves numeric precision and records time, concentration, event state, model, relevant parameters, dose, regimen, and end time. PNG includes axes, units, model, route, dosing, dose, and CL.

## Validation

Requires Node.js 18 or newer. No `npm install` is necessary:

```sh
npm test
```

Tests independently import `pk-model.js`; they cover reference IV results, PO Tmax, exact/near-equal absorption and disposition rates, integrated AUC, two-compartment eigenvalues, an independent RK4 three-amount ODE reference, absorption-limited terminal behavior, all eight combinations, superposition, IV pre/post jumps, PO continuity, interval availability, CSV, and invalid inputs. Tests use absolute/relative floating-point tolerances. No UI or chart library is needed for mathematical tests.

An optional browser smoke test is provided in `tests/browser-check.py`. It requires Python Playwright and a Chromium installation, uses a running local static server, and checks all eight UI combinations, scale switches, tooltip, exports, input errors, reset, and mobile overflow. Install Playwright separately if desired; it is not an application dependency. Set `PK_CHROMIUM` to override the Chromium executable path and `PK_TEST_URL` to override the static server URL:

```sh
python3 tests/browser-check.py
```

## GitHub Pages

1. Push these files to the repository's `main` branch.
2. In GitHub repository **Settings → Pages**, choose **Deploy from a branch**.
3. Select **main** and **/ (root)**, then save.

All asset references are relative, including the module import, so deployment under a repository subpath works without changes. No build workflow, environment secrets, or external service is required. Publication must be performed in GitHub; this project does not automatically publish itself.

## Files

- `index.html`: accessible English controls and dashboard structure.
- `style.css`: responsive layout.
- `main.js`: input handling, Canvas chart, results, tooltips, and downloads.
- `pk-model.js`: independent equations, validation, peak search, integration, and CSV serialization.
- `tests/`: mathematical validation and optional browser smoke test.

## V1 limitations

No nonlinear or time-varying PK, lag time, infusion, loading/maintenance regimens, steady-state accumulation ratio, fitting, variability, or clinical recommendations. Repeated regimens are limited to 100 doses to keep browser calculations responsive. Very extreme parameters beyond floating-point range are rejected; numerical results are estimates, not arbitrary-precision computations. The chart uses an adaptive finite grid; exported samples are not a continuous curve. Tooltip values are evaluated from the analytical model at the pointer time (snapping to dosing events). There is no keyboard-driven plot cursor; CSV provides exact numeric access. PNG export captures the plot and its caption, rather than the entire dashboard. No clinical or regulatory validation has been performed.
