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

## Optional NCA Sampling Mode (single dose)

Enable **NCA Sampling Mode** below the existing parameter controls. This extension leaves `pk-model.js` and its existing tests unchanged. With the mode off, the original curve, results and simulation CSV are unchanged; multiple-dose simulations remain available with NCA disabled. The mode adds `nca.js` for NCA and `nca-ui.js` for its optional controls and chart layer. There is still no build step, external chart dependency, or backend.

Enter up to 1,000 nonnegative, unique sampling times separated by commas or newlines. Numeric equivalents such as `1` and `1.0` are duplicates. Times are sorted numerically; points beyond the simulation end produce an error instead of silently being omitted. LLOQ must be finite and positive. A concentration **equal to LLOQ is quantifiable**; lower values are BLQ.

Samples reuse the existing single-dose analytical concentration function. The adapter retains true concentrations for the export, but `reportedSamples()` removes them before `analyzeNCA()` runs. BLQ is passed as `c: null`, never a hidden numeric value or a substituted LLOQ. Observed Cmax/Tmax use only quantified post-dose samples (earliest time wins a tie). If all samples are BLQ, parameters are N/A with reasons.

### Zero time and integration policy

- **IV 0 h defaults to predose**, excluded from NCA exposure and terminal regression. Selecting an explicit **postdose (0+)** sample evaluates the simulated sampled concentration at 0+. It is an idealized sampling choice, not an inferred measurement. Without a quantifiable sampled IV 0+, AUC0-last and AUC0-inf are **N/A**; the available partial AUCfirst-last is labeled separately. Theoretical C0 is never used to fill the gap or connect predose to the first post-dose sample.
- **PO** uses a structural zero origin. Leading BLQ samples are explicitly assigned zero **only for AUC integration**; they remain BLQ/null in reported data and cannot enter Cmax or regression. If there is no origin sample, a zero origin is inserted for integration. This exploratory convention may bias exposure and is displayed in the analysis warnings.
- **Intermediate BLQ** (between quantifiable samples) blocks AUC by default. The optional **explicit bridging** policy skips BLQ rows and integrates directly between the quantifiable endpoints. Warnings and the CSV list every affected bracket and its BLQ times. No hidden BLQ concentration is used. This is a selected approximation, not a reconstruction of the missing values.
- **Trailing BLQ** samples are excluded, stopping integration at the last quantifiable concentration (`Clast`). There is no LLOQ/2 substitution.
- AUC uses linear-up/log-down trapezoids: rising, equal, or zero-endpoint intervals use `Δt*(C1+C2)/2`; falling positive intervals use `Δt*(C1-C2)/ln(C1/C2)`, with stable logarithms near equality. Units are ng·h/mL.

### Terminal selection, extrapolation, and diagnostics

Automatic selection uses an **unweighted OLS regression of ln(reported concentration) against actual sampling time**. It considers consecutive suffixes of the final uninterrupted quantifiable run, always ending at the last quantifiable sample. PO points must be strictly after observed Tmax; IV may include its observed peak, including an explicit 0+ sample. No regression crosses a BLQ gap. At least three eligible points and a finite negative slope are required. Among all valid suffixes, the highest adjusted R² is preferred; candidates within 0.0001 of the best prefer more points. This is a transparent Best Fit-inspired rule, **not Phoenix WinNonlin validation or an assertion of numerical equivalence to Phoenix**.

Manual mode exposes checkboxes in the sampling table and uses exactly those eligible points. Invalid retained selections are explained and can be unchecked; they are not silently replaced. Manual selection may omit the last quantifiable sample. Regardless of selection, extrapolation uses **observed Clast at the last quantifiable sample**, not the regression-predicted value.

`λz = -slope`; terminal half-life is `ln(2)/λz`; `R² = 1 - SSE/SST` and `adjusted R² = 1 - (1-R²)*(n-1)/(n-2)` (intercept plus one slope); `AUC0-inf = AUC0-last + Clast/λz`; extrapolated percentage is `100*(Clast/λz)/AUC0-inf`. Terminal span is `(last regression time - first regression time)/estimated half-life`. No theoretical terminal rate fills missing NCA results. An invalid regression or unavailable AUC origin/gap yields N/A and a reason for dependent quantities.

The panel shows selected times, count, R², adjusted R², terminal span, and extrapolated percentage. Default configurable diagnostic thresholds are adjusted R² < 0.9, span < 2 half-lives, and extrapolated AUC > 20%. These are **exploratory flags, not universal regulatory acceptance criteria**. In addition, a half-life difference exceeding 20% from the known simulated observable asymptote triggers a model-comparison caution. High R² alone cannot show that a true terminal phase was sampled: early two-compartment distribution may be almost log-linear. Sparse PO samples can miss the true peak or remain in an absorption/distribution transition.

The comparison table presents sampled Cmax, Tmax, AUC0-inf and terminal half-life alongside the full single-dose model values. Relative bias is `100*(NCA-model)/model`; it is N/A if a value is unavailable or the model benchmark is zero (notably IV Tmax). For PO the terminal benchmark is the **apparent** asymptotic half-life `ln(2)/min(ka, terminal disposition rate)`, including flip-flop. The disposition half-life is displayed separately. Finite sampling windows, especially near-equal rates, need not have reached this asymptote. Rate descriptors are not measurable slopes for zero dose or zero bioavailability.

### Chart and exports

The existing continuous curve remains visible. Quantifiable samples are blue hollow circles; λz points are purple circles and the fitted segment is purple dashed. LLOQ is a labeled horizontal dashed line. **BLQ is an orange cross on a separate top rail with no numeric Y value**, never a point at LLOQ or at its hidden true concentration. The legend and exported PNG caption explain this convention. Sample hover shows exact time, quantitation state, and reported concentration only if quantifiable. Linear and log modes use real numeric times; true zero model concentrations remain omitted on log axes.

**NCA CSV** is a separate export; the original simulation CSV retains its original format. Each sample row includes sampling time, underlying theoretical concentration, reported concentration or BLQ, quantitation/event state, actual terminal selection, LLOQ, model/regimen parameters, NCA results, regression diagnostics, thresholds, requested manual times, integration policy/brackets, N/A reasons, warnings, and model comparisons. Underlying BLQ values in this explicitly labeled comparison export must not be interpreted as assayed concentrations. PNG export includes the optional sampling overlay when active.

### NCA tests

`npm test` runs both existing PK and new UI-independent NCA tests. NCA coverage includes ideal one-compartment IV regression and exact exponential AUC, missed PO Tmax, early two-compartment sampling with misleadingly high R², LLOQ equality and insufficient points, hidden-value isolation, synthetic intermediate BLQ with both integration policies, manual selection, predose/postdose distinction, leading/trailing BLQ, flip-flop comparisons, units, extrapolation, CSV, input validation, and non-mutation of all eight original PK combinations. Intermediate BLQ tests use synthetic reported observations because the noise-free V1 single-dose profiles generally do not generate such gaps after a descending phase.

Optional Chromium tests (same Python Playwright prerequisites as the original browser smoke test):

```sh
python3 tests/browser-check.py
python3 tests/nca-browser-check.py
```

The additional browser test checks that NCA activation leaves original PK results and CSV unchanged and that deactivation restores the same chart pixels, all four single-dose combinations, NCA unavailability for all multiple-dose combinations, manual regression, zero-time state, CSV/PNG downloads, malformed times/LLOQ, sample and BLQ tooltips, both scales, dosing transitions, and mobile layout. It uses `PK_TEST_URL` and `PK_CHROMIUM` overrides.

NCA V1 does not add assay noise, replicate samples at the same time, automatic IV C0 back-extrapolation, multiple-dose NCA, Phoenix equivalence, or clinical/regulatory validation. It intentionally makes missing information explicit instead of replacing sampled estimates with model truth.

### IV bolus C0 and AUC back-extrapolation

For single-dose IV bolus, a 0 h predose BLQ is **not** a postdose C0 and is excluded from IV AUC integration. If no observed 0+ concentration exists, the NCA module estimates C0 by log-linear back-extrapolation from the first two quantifiable postdose samples when their log slope is negative; otherwise it falls back to the first quantifiable postdose concentration. With fewer than two quantifiable postdose points and no observed 0+, C0 and AUClast are unavailable. The inserted C0 is used only for the AUC interval from 0 to the first postdose sample, not for Cmax/Tmax or terminal lambda-z regression. The NCA panel displays C0 method, the back-extrapolated area, and its percentage of AUCinf (when estimable). IV AUClast integrates from dosing time, not from observed Tmax. The first-two-point rule is adapted from Certara Phoenix NCA computation rules; this tool does not replicate Phoenix's outlier flags or regression weighting. Source: https://onlinehelp.certara.com/phoenix/8.2/topics/ncacomprules.htm
