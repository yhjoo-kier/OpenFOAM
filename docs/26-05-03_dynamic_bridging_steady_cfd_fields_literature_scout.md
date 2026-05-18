# Literature scout: dynamic bridging of steady CFD fields for transient flow surrogate modeling

## Question

Can transient CFD step-response information be used to bridge or dynamically interpolate a library of steady CFD fields at multiple flow-rate conditions?

This is relevant to the paint-booth/HVAC surrogate direction: use steady OpenFOAM fields at discrete supply flow rates, then add transient response/settling-time information so the surrogate can represent blower setpoint changes without running full transient CFD for every control trajectory.

## Short answer

Yes. The exact phrase "bridging steady CFD fields" is not the most common wording, but the idea overlaps strongly with several established research areas:

1. Quasi-steady / frozen-equilibrium modeling with dynamic lag correction.
2. Reduced-order modeling (ROM) using POD/Galerkin or non-intrusive regression.
3. Parametric ROM and interpolation on reduced spaces/manifolds.
4. Step-response / system-identification models for control-oriented CFD surrogates.
5. Dynamic Mode Decomposition (DMD), Koopman/operator-based models.
6. Neural-operator and sequence surrogate models for parametric/transient PDEs.
7. HVAC/indoor-airflow control-oriented CFD surrogate modeling.

For our project, the most defensible near-term formulation is:

```text
steady CFD library: F_ss(Q)
transient correction: dF/dt = [F_ss(Q_command(t)) - F] / tau(region, Q_i -> Q_j)
```

or, in discrete time:

```text
F_{k+1} = F_k + alpha [Interp_Q(F_ss, Q_k) - F_k]
alpha = 1 - exp(-dt / tau)
```

This is a first-order dynamic quasi-steady surrogate. It can be extended from scalar metrics to region-wise or field-wise lag models, and later to neural operators or operator-inference ROMs.

## Literature categories and representative references

### 1. Dynamic Mode Decomposition and data-driven flow ROMs

DMD is a standard method for extracting linear dynamical models from snapshot sequences of fluid flows. It is directly relevant once we have transient snapshots at `0, 2, ..., 60 s`.

Representative references found:

- Schmid, P. J. (2010). **Dynamic mode decomposition of numerical and experimental data**. *Journal of Fluid Mechanics*. DOI: `10.1017/S0022112010001217`. Crossref citation count observed: ~4974.
- Schmid, P. J. (2022). **Dynamic Mode Decomposition and Its Variants**. *Annual Review of Fluid Mechanics*. DOI: `10.1146/annurev-fluid-030121-015835`. Crossref citation count observed: ~565.
- Parametric DMD example:
  - **Parametric dynamic mode decomposition for reduced order modeling**. *Journal of Computational Physics* (2023). DOI: `10.1016/j.jcp.2022.111852`. Crossref citation count observed: ~42.

Relevance to us:

- DMD needs time-resolved snapshots, so it is not merely steady-field interpolation.
- Our 60 s step-response case can become the first DMD/Koopman-style training sequence.
- Parametric DMD becomes relevant when we run multiple steps, e.g. `2.18 -> 4.36`, `4.36 -> 5.45`, `6.54 -> 4.36`.

### 2. Parametric reduced-order models and interpolation of ROMs

This area is close to the idea of having a family of CFD solutions over parameters such as flow rate and interpolating reduced representations.

Representative references found:

- Amsallem, D., Farhat, C. (2008). **Interpolation Method for Adapting Reduced-Order Models and Application to Aeroelasticity**. *AIAA Journal*. DOI: `10.2514/1.35374`. Crossref citation count observed: ~558.
- Bui-Thanh, T., Willcox, K., Ghattas, O. (2008). **Parametric Reduced-Order Models for Probabilistic Analysis of Unsteady Aerodynamic Applications**. *AIAA Journal*. DOI: `10.2514/1.35850`. Crossref citation count observed: ~151.
- Xiao, D. et al. (2015). **Non-intrusive reduced order modelling of the Navier-Stokes equations based on RBF interpolation**. *International Journal for Numerical Methods in Fluids*. DOI: `10.1002/fld.4066`. Crossref citation count observed: ~165.
- Xiao, D. et al. (2015). **Non-intrusive reduced order modelling of the Navier-Stokes equations**. *Computer Methods in Applied Mechanics and Engineering*. DOI: `10.1016/j.cma.2015.05.015`. Crossref citation count observed: ~147.

Relevance to us:

- Our steady U-sweep library is a parametric CFD dataset with parameter `Q` or `U_supply`.
- RBF/interpolation-based non-intrusive ROMs are close to interpolating fields or reduced coefficients between steady cases.
- Interpolation of reduced bases/manifolds is more robust than naive cell-wise interpolation of full CFD fields, especially if topology or recirculation changes.

### 3. Operator inference and dynamic data-driven ROMs

Operator inference learns reduced dynamical systems from high-dimensional data. It is a principled alternative to an ad hoc first-order lag model once more transient data exist.

Representative references found:

- Peherstorfer, B., Willcox, K. (2016). **Data-driven operator inference for nonintrusive projection-based model reduction**. *Computer Methods in Applied Mechanics and Engineering*. DOI: `10.1016/j.cma.2016.03.025`. Crossref citation count observed: ~376.
- Peherstorfer, B., Willcox, K. (2015). **Dynamic data-driven reduced-order models**. *Computer Methods in Applied Mechanics and Engineering*. DOI: `10.1016/j.cma.2015.03.018`. Crossref citation count observed: ~164.
- Kramer, B., Peherstorfer, B., Willcox, K. (2024). **Learning Nonlinear Reduced Models from Data with Operator Inference**. *Annual Review of Fluid Mechanics*. DOI: `10.1146/annurev-fluid-121021-025220`. Crossref citation count observed: ~103.

Relevance to us:

- This provides a research-backed path from scalar response times to learned reduced dynamics.
- A possible pipeline is:
  1. collect steady fields and transient step snapshots;
  2. compute POD basis;
  3. project fields to low-dimensional coefficients;
  4. learn coefficient dynamics with flow-rate input;
  5. reconstruct fields/metrics for control.

### 4. POD-based ROMs for flow control

POD/Galerkin and POD-regression methods are classic for flow-field compression and control-oriented flow surrogates.

Representative references found:

- **Optimal flow control using a POD-based reduced-order model** (2016). *Numerical Heat Transfer, Part B*. DOI: `10.1080/10407790.2016.1173472`.
- **A POD reduced-order model for wake steering control** (2018). *Journal of Physics: Conference Series*. DOI: `10.1088/1742-6596/1037/3/032014`.
- Recent example found via Semantic Scholar: **Investigation on the Reduced-Order Model for the Hydrofoil of the Blended-Wing-Body Underwater Glider Flow Control with Steady-Stream Suction and Jets Based on the POD Method** (2024), *Actuators*. DOI: `10.3390/act13060194`.

Relevance to us:

- POD modes can be learned from our steady + transient snapshot fields.
- Control variables can be supply velocity or blower frequency.
- The first-order lag model can be a baseline against which POD/OpInf/DMD models are compared.

### 5. HVAC / indoor airflow control-oriented CFD surrogates

There is precedent for using CFD-derived airflow information for HVAC control, although many papers focus on room airflow rather than paint booths or neural operators.

Representative reference found:

- Pichurov, G., Stankov, P., Markov, D. (2006). **HVAC Control Based on CFD Analysis of Room Airflow**. *IFAC Proceedings Volumes*. DOI: `10.3182/20061002-4-BG-4905.00036`.

Other related titles found:

- **Indoor HVAC Prediction in Multi-Use Building Using Reduced Order Model** (2024). DOI: `10.3795/ksme-c.2024.12.2.135`.
- **Data-Driven Reduced-Order Model for Urban Airflow Prediction** (2023), conference chapter.

Relevance to us:

- The control-oriented HVAC literature supports replacing repeated CFD calls with reduced/surrogate models.
- Our contribution would be more specific: paint-booth/plenum/filter CFD, OpenFOAM-generated dataset, transient response identification, and neural-operator/control integration.

### 6. Neural operators and field-level surrogates

Neural operators learn mappings between functions/fields and are commonly used for parametric PDE surrogates. They become relevant for field-level prediction over geometry/parameter/time.

Representative known areas/papers to review next:

- Fourier Neural Operator (FNO) for parametric PDEs.
- Graph Neural Operator / GINO-like models for irregular meshes and CFD point clouds.
- Local neural operators for fluid dynamics.

A Crossref query found:

- **On the locality of local neural operator in learning fluid dynamics** (2024). *Computer Methods in Applied Mechanics and Engineering*. DOI: `10.1016/j.cma.2024.117035`.
- **Surrogate convolutional neural network models for steady computational fluid dynamics simulations** (2022). DOI: `10.1553/etna_vol56s235`.

Relevance to us:

- Neural operators can learn `F_ss(Q)` from steady samples, but dynamic control needs temporal input/state.
- A hybrid model is plausible:
  - neural operator learns the steady manifold `F_ss(Q, geometry)`;
  - transient CFD identifies dynamic lag/response kernels;
  - sequence neural operator or recurrent correction learns residual dynamics.

## How this maps to our paint-booth project

### Immediate baseline model

Use steady sweep fields as target manifold:

```text
F_target(t) = Interp_Q(F_ss(Q_command(t)))
```

Use transient step-response to estimate dynamic lag:

```text
F_dyn(t+dt) = F_dyn(t) + alpha_region [F_target(t) - F_dyn(t)]
alpha_region = 1 - exp(-dt / tau_region)
```

Start with scalar metrics:

- `work_zone_Uz_mean(t)`
- `near_car_Uz_mean(t)`
- `near_car_reverse_fraction(t)`
- `filter_layer_Uz_mean(t)`

Then extend to regions or field coefficients.

### Required CFD data

Minimum useful transient matrix:

- upward steps:
  - `2.18 -> 4.36`
  - `4.36 -> 5.45`
  - `4.36 -> 6.54`
- downward steps:
  - `6.54 -> 4.36`
  - `5.45 -> 4.36`
- optional local low-flow kink investigation:
  - `2.72 -> 3.27`
  - `3.27 -> 2.72`

For each case, compute time-series metrics and response times: 63%, 90%, 95%, overshoot, and settling time.

### Suggested paper framing

Possible framing for a future paper:

```text
A dynamic quasi-steady surrogate for paint-booth airflow control:
combining steady CFD libraries with transient step-response identification
```

Potential novelty relative to generic ROM literature:

- application-specific paint-booth/plenum/filter domain;
- CFD-to-control workflow using OpenFOAM transient step responses;
- field-level steady library plus region-wise dynamic lag;
- eventual GINO/neural-operator coupling for real-time HVAC/booth-flow control;
- reproducible HPC/closed-network automation workflow.

## Risks / cautions

1. Naive full-field linear interpolation can fail when recirculation topology changes.
2. A single time constant may be insufficient; region-wise or modal time constants are more defensible.
3. Upward and downward steps may be asymmetric.
4. Scalar settling does not guarantee field-level settling near car wakes.
5. Real blower dynamics should be separated from CFD boundary-flow dynamics:

```text
frequency command -> actual supply flow -> booth flow-field response
```

Our OpenFOAM step case identifies the second part.

## Recommended next actions

1. Implement time-series metric extraction for the completed 60 s case.
2. Estimate response times for key scalar metrics.
3. Build the first dynamic quasi-steady scalar surrogate and compare against the transient CFD time series.
4. Add at least one downward step case to test asymmetry/hysteresis.
5. After 3--5 step cases, test POD/DMD/OpInf or neural-operator sequence models.
