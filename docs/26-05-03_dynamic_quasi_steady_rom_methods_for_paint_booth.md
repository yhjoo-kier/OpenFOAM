# Dynamic quasi-steady / ROM methods for paint-booth airflow surrogate

## Purpose

This note converts the literature scout into directly usable methodology for the paint-booth/HVAC control surrogate project.

The immediate research objective is pragmatic:

> Use existing steady CFD libraries well first. Add transient response information only to the extent needed for control/surrogate use. If existing methods expose limitations for paint-booth/plenum/filter flows, those gaps become later research/paper topics.

## Problem setting

Let the CFD state be

```text
F(x, t) = { U(x,t), p(x,t), k(x,t), omega(x,t), ... }
```

and let the control/input parameter be supply flow or supply velocity:

```text
q(t) = U_supply(t) or Q_supply(t)
```

We already have a steady CFD library:

```text
D_ss = { q_i, F_ss(x; q_i), metrics_ss(q_i) }_{i=1..N}
```

For the current OpenFOAM paint-booth case, `q_i` corresponds to supply velocities such as:

```text
2.18, 2.72, 3.27, 3.81, 4.36, 4.91, 5.45, 6.00, 6.54 m/s
```

We also now have one transient step-response sequence:

```text
q: 4.36 -> 5.45 m/s at t = 0
F(x,t), t = 0, 2, 4, ..., 60 s
```

## Method family A: dynamic quasi-steady scalar surrogate

### Concept

For control, we often need a few scalar outputs rather than the full field:

```text
y(t) = [
  work_zone_Uz_mean(t),
  near_car_Uz_mean(t),
  near_car_reverse_fraction(t),
  filter_layer_Uz_mean(t),
  mass_imbalance(t),
  ...
]
```

The simplest useful model is:

```text
y_target(t) = Interp_q[y_ss(q(t))]
```

and then the actual output follows the target with a finite response time:

```text
dy/dt = (y_target(t) - y(t)) / tau_y
```

Discrete form:

```text
y_{k+1} = y_k + alpha_y [ y_target,k - y_k ]
alpha_y = 1 - exp(-Delta_t / tau_y)
```

### Step-response identification

For a step from `q0` to `q1`, if the response is monotonic and approximately first order:

```text
y(t) = y_final + [y_initial - y_final] exp(-t / tau_y)
```

Normalized response fraction:

```text
r(t) = (y(t) - y_initial) / (y_final - y_initial)
```

For an increasing response, the ideal first-order milestones are:

```text
r(tau)        = 1 - exp(-1)      = 0.632
r(2.3026 tau) = 0.90
r(2.9957 tau) = 0.95
r(4.6052 tau) = 0.99
```

Therefore:

```text
tau_y ≈ t_63
settling_90 ≈ 2.3 tau_y
settling_95 ≈ 3.0 tau_y
```

For decreasing responses, use the same normalized formula; the denominator handles sign:

```text
r(t) = (y(t) - y_initial) / (y_final - y_initial)
```

### Overshoot and settling

If the transient is not monotonic, use more robust metrics:

```text
overshoot = max_t | y(t) - y_final | / | y_final - y_initial |
```

A practical settling time with tolerance epsilon:

```text
t_settle,epsilon = min T such that
  |y(t) - y_final| <= epsilon |y_final - y_initial|
  for all t >= T
```

Common tolerances:

```text
epsilon = 0.10, 0.05, 0.02
```

### How to use in our current project

1. Build `y_ss(q)` from the 9-point steady sweep.
2. Extract `y(t)` from transient snapshots.
3. Estimate `tau_y` and/or settling times for each metric.
4. Use the dynamic model for fast control simulation.

The first target outputs should be:

```text
work_zone_Uz_mean
near_car_Uz_mean
near_car_reverse_fraction_Uz_positive
filter_layer_panel_Uz_mean
floor_outflow_m3s
relative_mass_imbalance
```

## Method family B: dynamic quasi-steady full-field surrogate

### Concept

The field-level extension replaces scalar `y` with CFD state `F`:

```text
F_target(x,t) = Interp_q[F_ss(x; q(t))]
```

Then:

```text
partial F / partial t = [F_target(x,t) - F(x,t)] / tau_F(x)
```

or with one global time constant:

```text
F_{k+1}(x) = F_k(x) + alpha [F_target,k(x) - F_k(x)]
```

### Region-wise version

A more defensible middle ground uses different time constants per physical region:

```text
F_{k+1}(x in R_j)
= F_k(x) + alpha_j [F_target,k(x) - F_k(x)]

alpha_j = 1 - exp(-Delta_t / tau_j)
```

Candidate regions:

```text
R_1: plenum
R_2: filter layer
R_3: work zone
R_4: near-car region
R_5: wake/reverse-flow-sensitive region
```

### Risk

Cell-wise field interpolation may be poor if:

- recirculation zones move;
- wake topology changes;
- low-flow kink/non-monotonic behavior appears;
- different variables need different lags.

Therefore, use full-field dynamic quasi-steady interpolation as a baseline, not as the final claim.

## Method family C: POD / reduced-basis model

### Snapshot matrix

Flatten fields or selected variables into vectors:

```text
f_j = vec(F(x, t_j; q_case)) in R^n
```

Build snapshot matrix:

```text
X = [ f_1, f_2, ..., f_m ] in R^{n x m}
```

Subtract mean if needed:

```text
X_c = X - mean(X)
```

Compute SVD:

```text
X_c = Phi Sigma V^T
```

Retain `r` modes:

```text
Phi_r = [phi_1, ..., phi_r]
```

Approximate field:

```text
f(t) ≈ f_mean + Phi_r a(t)
```

where `a(t)` are reduced coefficients:

```text
a(t) = Phi_r^T [f(t) - f_mean]
```

### Parametric interpolation of coefficients

For steady fields:

```text
a_ss(q_i) = Phi_r^T [f_ss(q_i) - f_mean]
```

Interpolate in parameter space:

```text
a_target(q) = Interp_q[a_ss(q_i)]
```

Then apply dynamic lag in coefficient space:

```text
a_{k+1} = a_k + Alpha [a_target(q_k) - a_k]
```

where `Alpha` can be scalar, diagonal, or learned.

### Why this is useful

This is safer and more compact than full-field interpolation:

```text
field interpolation -> interpolate millions of values
POD coefficient interpolation -> interpolate r coefficients
```

It also creates a bridge toward DMD/operator inference.

## Method family D: DMD / Koopman-style linear dynamics

### Basic DMD formulation

Given time-resolved snapshots:

```text
X  = [ f_1, f_2, ..., f_{m-1} ]
X' = [ f_2, f_3, ..., f_m     ]
```

Assume linear map:

```text
X' ≈ A X
```

DMD seeks a low-rank approximation of `A`.

SVD:

```text
X = U Sigma V^T
```

Reduced operator:

```text
A_tilde = U_r^T X' V_r Sigma_r^{-1}
```

Eigen-decomposition:

```text
A_tilde W = W Lambda
```

DMD modes:

```text
Phi_DMD = X' V_r Sigma_r^{-1} W
```

Time evolution:

```text
f_k ≈ Phi_DMD Lambda^{k-1} b
```

### With input/control

For controlled systems, use DMD with control (DMDc):

```text
x_{k+1} = A x_k + B u_k
```

where:

```text
x_k = reduced flow state
u_k = supply velocity or flow-rate command
```

### Relevance to us

DMD/DMDc is attractive once we have several transient step sequences. For now, one sequence is enough only for a first diagnostic, not a robust parametric model.

## Method family E: operator inference in POD space

### General model

Operator inference learns the reduced dynamics:

```text
da/dt = c + A a + H(a ⊗ a) + B u
```

or discrete:

```text
a_{k+1} = c + A a_k + H(a_k ⊗ a_k) + B u_k
```

where:

```text
a_k: POD coefficients
u_k: supply flow input
A, H, B: learned operators
```

### Relevance

This is more research-grade than the first-order lag model. It can represent nonlinear transient dynamics and input effects, but it needs more transient data.

Suggested progression:

```text
1. scalar first-order dynamic model
2. region-wise first-order dynamic model
3. POD coefficient lag model
4. DMDc / operator inference model
5. neural operator / sequence model
```

## Method family F: neural-operator / sequence surrogate

### Steady operator

Learn a steady mapping:

```text
N_ss: (geometry, mesh/cell centers, q) -> F_ss(x; q)
```

### Dynamic operator

Learn transient update:

```text
N_dyn: (F_k, q_k, q_{k-1}, ..., x, geometry) -> F_{k+1}
```

or residual form:

```text
F_{k+1} = F_k + N_res(F_k, F_target(q_k), q_history, x)
```

### Hybrid practical form

Use steady CFD/neural operator as the target manifold and learn only the residual dynamics:

```text
F_target = N_ss(q_k)
F_{k+1} = F_k + alpha(x,q) [F_target - F_k] + R_theta(...)
```

This is a possible future contribution if simple lag/POD/DMD methods fail in paint-booth wake regions.

## Recommended implementation sequence for our project

### Stage 1: scalar time-series extraction

From the completed 60 s step case, compute at each written time:

```text
t
supply_inflow_m3s
floor_outflow_m3s
relative_imbalance
work_zone_Uz_mean
near_car_Uz_mean
near_car_reverse_fraction_Uz_positive
filter_layer_panel_Uz_mean
```

Then compare with initial and target steady values.

### Stage 2: scalar first-order model

Fit per-output response times:

```text
tau_work_zone_Uz
tau_near_car_Uz
tau_reverse_fraction
tau_filter_layer_Uz
```

Then build a baseline model:

```text
y_dyn(k+1) = y_dyn(k) + alpha_y [Interp_q(y_ss, q_k) - y_dyn(k)]
```

### Stage 3: add more step cases

Minimum next cases:

```text
2.18 -> 4.36
4.36 -> 6.54
5.45 -> 4.36
6.54 -> 4.36
```

Reason:

- upward/downward asymmetry;
- low/high flow dependence;
- possible hysteresis or recirculation topology changes.

### Stage 4: reduced-basis field model

Export consistent `.npz` snapshots and build:

```text
POD basis + coefficient interpolation + dynamic lag
```

This becomes the field-level baseline surrogate before GINO/neural-operator training.

## Candidate validation metrics

For scalar model:

```text
MAE/RMSE over y(t)
error in t63/t90/t95
overshoot prediction error
settling-time prediction error
```

For field model:

```text
relative L2 error of U
relative L2 error of Uz
region-wise mean error
reverse-flow mask IoU or F1
near-car/wake-local error
```

For control use:

```text
prediction speed
stability under long command sequence
error accumulation over repeated setpoint changes
```

## Literature anchors

The following references support the methodology categories:

- Schmid, P. J. (2010). Dynamic mode decomposition of numerical and experimental data. *Journal of Fluid Mechanics*. DOI: `10.1017/S0022112010001217`.
- Schmid, P. J. (2022). Dynamic Mode Decomposition and Its Variants. *Annual Review of Fluid Mechanics*. DOI: `10.1146/annurev-fluid-030121-015835`.
- Amsallem, D., Farhat, C. (2008). Interpolation Method for Adapting Reduced-Order Models and Application to Aeroelasticity. *AIAA Journal*. DOI: `10.2514/1.35374`.
- Bui-Thanh, T., Willcox, K., Ghattas, O. (2008). Parametric Reduced-Order Models for Probabilistic Analysis of Unsteady Aerodynamic Applications. *AIAA Journal*. DOI: `10.2514/1.35850`.
- Xiao, D. et al. (2015). Non-intrusive reduced order modelling of the Navier-Stokes equations based on RBF interpolation. *International Journal for Numerical Methods in Fluids*. DOI: `10.1002/fld.4066`.
- Xiao, D. et al. (2015). Non-intrusive reduced order modelling of the Navier-Stokes equations. *Computer Methods in Applied Mechanics and Engineering*. DOI: `10.1016/j.cma.2015.05.015`.
- Peherstorfer, B., Willcox, K. (2016). Data-driven operator inference for nonintrusive projection-based model reduction. *Computer Methods in Applied Mechanics and Engineering*. DOI: `10.1016/j.cma.2016.03.025`.
- Peherstorfer, B., Willcox, K. (2015). Dynamic data-driven reduced-order models. *Computer Methods in Applied Mechanics and Engineering*. DOI: `10.1016/j.cma.2015.03.018`.
- Kramer, B., Peherstorfer, B., Willcox, K. (2024). Learning Nonlinear Reduced Models from Data with Operator Inference. *Annual Review of Fluid Mechanics*. DOI: `10.1146/annurev-fluid-121021-025220`.
- Pichurov, G., Stankov, P., Markov, D. (2006). HVAC Control Based on CFD Analysis of Room Airflow. *IFAC Proceedings Volumes*. DOI: `10.3182/20061002-4-BG-4905.00036`.

## Practical conclusion

For our immediate goal, the most useful existing methodology is not a heavy neural operator first. It is:

```text
steady CFD library
+ scalar/region-wise transient step-response identification
+ dynamic quasi-steady surrogate baseline
```

This gives a fast, interpretable, control-ready baseline. If this baseline fails for near-car wake, reverse-flow topology, or repeated setpoint changes, then those failures define the need for POD/DMD/OpInf/neural-operator methods and may become the research novelty.
