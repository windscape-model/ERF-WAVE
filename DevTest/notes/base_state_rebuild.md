# Moving-mesh base state: incremental evolution vs rebuild (2026-09-01)

## Finding

On `MovingFittedMesh` the base state is evolved with the mesh every stage
(`ERF_TI_slow_rhs_pre.H`): (rho theta)_0 is advanced in flux form with the grid velocity
and p_0 is re-created from it through the equation of state. The pow() round trip
p_0 -> (rho theta)_0 -> p_0 accumulates one-sided round-off of about 5e-11 Pa per model
step (3 ulp of 1e5 Pa), uniform in the column, present even with the mesh at rest
(a = 0 wave: pert_pres drifts to +3.7e-8 Pa after 1000 steps, +5.5e-8 Pa after 2000).
The drift enters the perturbation pressure p - p_0 written to plotfiles one to one and is
independent of the wave amplitude. Two runs with different step counts at the same acoustic
step therefore differ by a constant that floors fixed-dtau pressure comparisons once the
physical differences fall below ~1e-7 Pa (reflection-free wave, quiescent start: the finest
dtau = 2.5e-4 s pair gave 1.36 on pert_pres and 1.98 on the state pressure). The dynamics
are essentially unaffected: the vertical pressure gradient is not in perturbational form.
Round-off floor of the pressure diagnostic itself (twin runs): 1.1e-9 relative.

## Option

`erf.moving_base_state_rebuild = 1` (default 0) rebuilds the base state on the stage-target
geometry from the initial hydrostatic profile at every stage, exactly as at initialization
(`prob->erf_init_dens_hse` at the new cell heights, then `erf_enforce_hse`), instead of
evolving it. It keeps dp_0/dz = -rho_0 g. Restricted to the custom isentropic base state
(`erf.init_type` none or Isentropic); a sounding-based base state is not re-derivable from
`erf_init_dens_hse` and the code aborts with a message.

## Validation (local, serial, reflection-free wave 40x4x200, quiescent start)

- flag off: gates g1-g6 PASS (bit-identical to the stored references, trap-off runs on the
  Apple-Silicon machine); the finest quiescent run reproduces the jcp-writing-pass binary
  bitwise.
- flag on, wave-free column, 2000 steps: mean pert_pres -1.9e-11 Pa (evolved: +5.5e-8 Pa).
- flag on, dtau = 2.5e-4 s finest pair (dt 1e-3 / 5e-4): pert_pres column rate 2.00
  (evolved: 1.36); state-pressure rate 2.00 (evolved: 1.98).
- mass and rho*theta drift unchanged (+2.024e-12 over 1 s for both treatments).
- the two treatments give solutions that differ by 1.3e-4 of the perturbation-pressure norm
  in the pressure field (structured, hydrostatic-like), 7e-7 in u and 2e-9 in w: a
  discretization-level difference of the base-state/perturbation split, not a defect of
  either.

Raw runs: `ERF-Research/rerun_quiescent_2026-09-01/runs/{rebuild_amp0,rebuild_fixdtau,
dev_default}`; runner `run_rebuild_check.sh` there.

## Open observation (not addressed)

With beta_s = beta_d = 0 the fast integration's error is second order in dtau for the
pressure perturbation but first order in dtau for u, w and theta (M-arm at dt = 0.004
against M = 32: the differences follow (1/M - 1/32) with a constant coefficient to within
4 percent; pure-first-order rates 1.10/1.22/1.58 against observed 1.13/1.24/1.59). It
cancels exactly at fixed dtau (coupling rates 2.00 in every norm). It is not a half-sub-step
time offset of the field (no correlation with du/dt) and a static-mesh test from the
potential-flow field was transient-dominated and inconclusive. Mechanism undiagnosed.

Resolved 2026-10-29: the fork's review of 2026-09-19 (`review_substep_mt_2026-09-19/REPORT.md`)
traced the w part to the missing volume-change tendency of the implicit predictor (D1; with it
the sub-step order of w rises from 1.0 to 1.6-1.9) and the u and theta part to the dtau/2 lag of
the forward-backward sub-step scheme. D1 and the other corrections are in this branch since the
merge of the tag `jcp-manuscript-2026-10-02` (commit `ecaf57f4`).
