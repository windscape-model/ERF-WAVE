# Capability roadmap: DNS / wall-resolved LES / wall-modeled LES x compressible / anelastic

Written 2026-09-01 on branch `dev-wave-capabilities`. The JCP manuscript is frozen;
nothing here enters it. Target: after the manuscript is submitted, release the full
code with both equation sets and all three fidelity levels over the moving wave.

## Capability matrix (state at branch point, with evidence)

| | DNS | Wall-resolved LES | Wall-modeled LES |
|---|---|---|---|
| Compressible, moving wave | WORKS laminar (JCP Stokes case, TV wall, 0.95 percent amplitude error) | WIRING WORKS (smoke test below); physics validation open | WORKS (JCP turbulent cases, MOST wave drag) |
| Anelastic, flat | WORKS upstream (Channel_DNS deck, Re_tau 395) | same knobs as compressible; untested | expected (upstream MOST); untested |
| Anelastic, static wave | runs clean (stage A test below) | untested | untested |
| Anelastic, moving wave | WORKS: Stokes 0.99 percent amplitude / 0.9 degree phase vs analytic (stage B + C) | WIRING WORKS (smoke, SGS + molecular + orbital wall) | WIRING WORKS (mini, MOST wave drag) |

Evidence for the shared viscosity knobs:

- Molecular and eddy viscosity are ADDITIVE in momentum stress:
  `Source/Diffusion/ERF_ComputeStress_T.cpp` (`mu_tot = rhoAlpha(i,j,k) + 2*mu_turb(i,j,k,EddyDiff::Mom_v)`)
  and in scalar diffusion: `Source/Diffusion/ERF_GetRhoAlpha.H`.
  So DNS = `molec_diff_type=Constant, les_type=None`; WRLES = both on; WMLES = LES + MOST.
  The JCP turbulent deck already ran Smagorinsky + Constant(1e-2) together.
- The orbital no-slip wall for resolved runs is the TERRAIN_VELOCITY build
  (`ERF_BoundaryConditionsXvel.cpp`), verified by the JCP Stokes case.

## CI state (checked 2026-09-01)

- GitHub Actions are DISABLED at the repository level on the fork
  (`gh api repos/tli-fluid/ERF-mrigark/actions/permissions` returns enabled=false);
  the 16 upstream workflow files have never run there. Only Dependabot ran.
- The in-tree CTest harness (`Tests/CTestList.cmake` + `Tests/ERFGoldFiles`, including
  MovingTerrain_nosub/sub) exists but needs a CMake build and its gold files predate
  the fork's moving-mesh default (midpoint metric sampling), so its currency on this
  branch is unverified.
- Working protection is the local gate suite `DevTest/gates/run_gates.sh` (five gates,
  bit-level references, all PASS at the branch point). Enabling cloud CI on the private
  fork is the user's call (Actions minutes billing); the gates do not depend on it.

## Direction 1: compressible DNS / WRLES (mostly enablement + validation)

Done here:

- `wrles_wave/inputs_smoke`: Smagorinsky + Constant(0.02) + TV orbital wall + moving
  wave, 200 steps, 4 ranks. Runs clean; eddy viscosity active (max |Kmv| 2.5e-3 m2/s
  vs molecular nu 1.7e-2); fields finite; spanwise symmetry preserved.

Open items, in order:

1. Physics validation of WRLES: a resolved channel-over-wave at moderate Re_tau against
   the literature (or against a DNS of the same case). Needs an HPC-scale run; design
   the deck from `Exec/CanonicalTests/Channel_DNS` plus the wave bottom.
2. RESOLVED (2026-10-29). The May crash and the July 14-percent error of `erf.vert_implicit=1`
   on the moving wall were the wall term and the doubled coefficient of the implicit solve,
   fixed in the fork on 2026-05-27 and 2026-07-25. On 2026-10-29 the terrain kernel took the
   face metric of the strain kernel and the row closure at grid ends (erf-model/ERF, October
   2026); the flat and stretched kernels took the one-sided wall gradient. Measured on the
   manuscript Stokes layer at full duration: implicit and explicit agree to 1.07e-5 m/s at
   dt 0.004 s and 1.7e-7 m/s at dt 5e-4 s (second order in dt, the splitting error), with the
   same error against the analytic solution at every period. Implicit vertical diffusion is
   available for the wall-resolved cases.
3. RESOLVED (2026-09-01, full read): the tau13/tau23 metric cross-term in
   `ERF_ComputeStress_T.cpp` averages `rhoAlpha_bar` (molecular) alongside `mu_bar`
   (eddy) and sums both into `mu_tot` -- the molecular part IS included. The earlier
   flag came from a truncated grep window; no change needed.
4. Near-wall SGS behavior (no wall damping in plain Smagorinsky): acceptable for
   capability release; note in the docs, consider Deardorff for production WRLES.

## Direction 2: anelastic moving wave (the real development)

### What exists (all read and confirmed at file level)

- Anelastic integrator: RK2, no acoustic substepping (`ERF_MRI.H`), stage update via
  `ERF_TI_no_substep_fun.H`, projection at the end of each stage
  (`project_momenta`, `Source/LinearSolvers/ERF_PoissonSolve.cpp`).
- `no_substep_fun` ALREADY has the MovingFittedMesh branch: J-weighted cons update,
  h_zeta-weighted x/y momenta, dJ-face-weighted z momentum, bottom zmom set from
  WFromOmega(rho * z_t). This is what upstream MovingTerrain_nosub exercises.
- Terrain Poisson solve: GMRES with FFT preconditioner
  (`ERF_SolveWithGMRES.cpp`, `ERF_TerrainPoisson.cpp`); REQUIRES a USE_FFT build
  (hard abort otherwise). Local FFT build works (see Build variants below).
- Divergence operator: `ERF_ComputeDivergence.cpp` (ax, ay, detJ metrics; z-flux is
  Omega). Momentum correction adds solver fluxes directly to (rho0 u, rho0 v, Omega).
- Stage A tests (2026-09-01, FFT build, 20 steps each, `anelastic_wave/`):
  - static wave: runs clean, projected divergence 3e-13.
  - moving wave AS-IS: runs, no assert stops it, solver reports divergence 5e-14 —
    but of the WRONG (static-mesh) constraint. Silent inconsistency, see below.

### Stage B status (2026-09-01): IMPLEMENTED and validated at the 1 percent level

All three defects below are fixed behind `erf.anelastic_moving_mesh = 1` (default:
the combination aborts). A FOURTH defect was found during validation and fixed:
the anelastic second-RK-stage Heun refold in `ERF_SlowRhsPre.cpp` reconstructed the
first-stage tendency as (S1-S0)/dt, a static-mesh identity; on the moving mesh that
difference contains the volume-dilution rate -S d(lnJ)/dt, which corrupted theta by
~0.9 K (dt-independent, phase-locked, surface-peaked) and the momenta likewise.
The refold is now J-weighted on the moving mesh (`ERF_TI_slow_rhs_pre.H`):
dt Jsrc F0 = J1 S1 - J0 S0 for cons, h_zeta-weighted for x/y momenta.

Validation (reflection-free wave geometry, a=0.16, lambda=100, rigid lid, ramped
start, fit of the oscillating mode over t = 16..20 s at dt = 0.0125):

- w amplitude vs the analytic lid-bounded potential flow: 0.07 percent at the
  surface, max 1.2 percent, mean 0.7 percent over the column; u max 0.7 percent;
  phase 0.3-0.45 degrees; quasi-steady residual mode 4e-4 m/s (0.3 percent).
- density pin exact (max|rho - rho0| = 0); projected divergence at solver
  tolerance with the d(J rho0)/dt source included (the only remainder is the
  spatially constant solvability mean, the base-state inventory term).
- The compressible twin (same grid and duration, no sponge) agrees near the
  surface but is contaminated aloft by lid-trapped acoustic modes (two-mode fit
  residual 100x larger) -- the anelastic solver needs no acoustic damping
  apparatus at all, which is the point of the capability.

PROJECTION-LAG REMEDY: IMPLEMENTED. The step-end projection satisfies the interval
just completed, so the step-start momenta lagged the mesh by one interval and theta
carried a phase-locked first-order error ~ 0.5 * 300 * z_t * dt / dz (0.47 K at
dt=0.05). The fix re-projects the stage-input momenta against the UPCOMING
interval at the top of stage 1 (`ERF_TI_slow_rhs_pre.H`, right after
update_terrain_stage; one extra Poisson solve per step; only the stage-evaluation
copy S_data is modified, S_old stays the conservative anchor; base_state_new still
holds the step-start base state there, which is the correct source to within the
base-state transport term ~1e-5 relative). Measured after the fix: theta error
0.0070 K at dt=0.05 and 0.0058 K at dt=0.005 -- both at the mK physical floor,
no longer dt-limited. Wave validation improved to w max 0.95 percent / mean 0.56
percent, u max 0.52 percent, phase < 0.44 degrees.
Still inherited upstream: slow scalars (RhoKE etc.) get no Heun average on the
moving mesh (`ERF_SlowRhsPost.cpp` moving branch). `vert_implicit` +
`anelastic_moving_mesh` now ABORTS with a clear message (untested combination;
use explicit vertical diffusion).

TEMPORAL ORDER (measured 2026-09-01): the anelastic integration is FIRST order in
time, and that is an UPSTREAM property, not a moving-mesh one. Evidence: dt-halving
Richardson on the moving wave gives order 1.0 in u, w, theta, pert_pres; the same
study on the STATIC-terrain anelastic density current (pure upstream path) gives
0.84-0.93. The stage-projection Heun with incremental pressure carries the classic
projection-splitting order limit. The JCP compressible scheme is the second-order
one; a second-order anelastic would need a time-centered pressure treatment in the
upstream integrator -- out of scope here. At practical dt the first-order error is
small (the 0.95 percent wave validation at dt = 0.0125 already contains it).

RESTART (measured 2026-09-01): checkpoint at step 100, restart to step 200,
compare with the uninterrupted run -- bit-for-bit identical in every field.

### The three defects fixed (stage B history)

In `project_momenta` and `compute_divergence`, on MovingFittedMesh:

1. GEOMETRY CURRENCY: the projection uses `z_phys_nd[lev]`, `detJ_cc[lev]`,
   `ax/ay/az[lev]` (old step geometry). It must use the stage-target geometry
   (the `*_new` fields that `no_substep_fun` already uses), both in the divergence
   and in the GMRES/TerrainPoisson coefficients, rebuilt every stage.
2. GRID VELOCITY IN OMEGA: `OmegaFromW` is the static formula. The projected flux
   must be the flux relative to the moving grid: Omega_rel = OmegaFromW(rho0 w) -
   rho0 z_t. With the mesh following the material wall, Omega_rel(k=0) = 0 is then
   exactly the current bottom condition (the existing k=0 zero is correct once the
   subtraction is in place; convert back after the solve).
3. MESH-MOTION SOURCE: the anelastic constraint on a moving mesh is
   d(J rho0)/dt + div(rho0 U_contravariant,rel) = 0,
   so the Poisson RHS needs the discrete GCL-consistent source
   (dJ_new*r0_new - dJ_old*r0_old) / (dt_stage * dJ_new)
   added to the computed divergence. r0_new is the base state on the new mesh
   (`base_state_new`, maintained by the existing moving-terrain machinery).
   With this discrete form the projected field satisfies the same discrete GCL the
   compressible scheme satisfies (the JCP construction carries over).

Plus one guard: until stage B is complete, make anelastic + MovingFittedMesh abort
with a clear message instead of running silently wrong; the development flag
(suggestion: `erf.anelastic_moving_mesh`, default false = abort) lifts the guard.
Base-state currency in `ConvertForProjection` (r_hse alias) must move to the
new-time base state as part of item 1.

### Stage B validation ladder (each step gated)

1. Freestream/GCL: uniform horizontal flow, no gravity, moving wave, anelastic.
   Exact answer: the flow stays uniform. The (J rho0)_t source must cancel the
   mesh-motion divergence to solver tolerance. Sharpest single test of items 1-3.
2. Inviscid moving wave vs the potential-flow solution (the JCP reflection-free
   configuration, same geometry): compare wave-induced w and pert_pres against the
   compressible result and the analytic solution. The anelastic solution has no
   acoustic mode, so agreement should be at least as good.
3. Viscous Stokes layer under anelastic (TV wall + Constant viscosity): amplitude
   and phase against the analytic Stokes solution, mirroring the JCP 0.95 percent /
   0.910 degree numbers.
4. Miniature WMLES (the g4 configuration with anelastic=1): stability + sanity of
   surface stress; then an A/B against the compressible g4 statistics window.

### Stage C: the fidelity ladder under anelastic

DNS / WRLES / WMLES knobs are equation-set independent (same diffusion and MOST
code paths), so after stage B this is configuration + testing, not development.
Status:

- WMLES mini (`inputs_wmles_anelastic`, the g4 configuration with anelastic=1,
  no substepping, no sponge): runs 100 steps clean on the FFT build; rho pin
  exact, theta perturbation 0.03 K (bounded), log-law u range preserved, fields
  finite. Ready for a longer statistics A/B against compressible g4 when needed.
- Viscous Stokes DNS (`inputs_stokes_anelastic`, the JCP in.uniform grid on the
  TV+FFT build, no sponge, no substepping, t = 12 s = 2T): amplitude L2 error
  0.99 percent, mean phase lag 0.9 degrees against the analytic Stokes solution
  (`Exec/visc_test/stokes_compare.py`) -- matching the JCP compressible
  verification (0.95 percent, +0.910 degrees). Theta stays within 7 mK over the
  two periods; rho pin exact.
- WRLES smoke under anelastic (wrles_wave deck + anelastic overrides on the
  TV+FFT build, `erf.rayleigh_damp_W=0` to switch the deck's sponge off): runs
  clean; eddy viscosity active (Kmv 2e-3 m2/s), theta within 2 mK, rho pin
  exact. NOTE: the "off" switch for Rayleigh damping is the per-variable damp
  flags; `rayleigh_damping_type` has no None token.

## Build variants (GNU make, from Exec/)

    make -j8 USE_NETCDF=FALSE                          -> ERF3d.gnu.TEST.MPI.ex      (gates g1,g2,g4,g5)
    make -j8 USE_NETCDF=FALSE USE_TERRAIN_VELOCITY=TRUE    -> ERF3d.gnu.TEST.MPI.TV.ex   (gate g3, WRLES smoke)
    make -j8 USE_NETCDF=FALSE USE_FFT=TRUE FFTW_HOME=$FFTW \
                                                       -> ERF3d.gnu.TEST.MPI.FFT.ex  (anelastic wave work)
    make -j8 USE_NETCDF=FALSE USE_FFT=TRUE FFTW_HOME=$FFTW USE_TERRAIN_VELOCITY=TRUE \
                                                       -> ERF3d.gnu.TEST.MPI.TV.FFT.ex (anelastic Stokes)

FFTW is the FFTW prefix: `$HOME/local` on the laptop of September 2026, `/opt/homebrew` (Homebrew
fftw) since October 2026; `run_gates.sh` finds it. USE_NETCDF=FALSE: no wave or gate case reads
NetCDF, and neither laptop has a parallel NetCDF. Under Apple clang add
`XTRA_CXXFLAGS=-ffp-exception-behavior=maytrap` to every build, after removing `tmp_build_dir`
(make does not track flag changes): the decks set `amrex.fpe_trap_invalid = 1`, and the default
exception behavior lets a speculated floating-point operation of the moving-terrain initializer
trip the trap (SIGILL on Apple Silicon); the results are unchanged bit for bit. Each variant has
its own object directory, so the executables coexist.

## Standing protections

- Every source change runs `DevTest/gates/run_gates.sh` before commit; the five gates
  must stay bit-identical with new flags off.
- New features default off; behavior changes only behind named input flags.
- Test runs write only under `DevTest/`; `Exec/` run directories, `verification_runs/`,
  and `mrigark_paper/` stay untouched.
- `jcp-writing-pass` receives nothing from this branch before submission.

## Status 2026-10-29: the manuscript code merged in

- Merge commit `ecaf57f4` brings the tag `jcp-manuscript-2026-10-02` (= `eff4101c` on
  `jcp-writing-pass`) into this branch: the acoustic sub-step corrections D1 (volume-change
  tendency of the implicit predictor), D9 (perturbation flux with each state on its own mesh),
  D4 (conservative w reconstruction), D2 (grid-velocity flux in the predictor) and D5 (slow
  tendency with the stage volume); the stress-divergence correction D10 (bottom Omega = z_t on
  the moving surface); the input `prob.wave_c`; the diagnostic options `erf.slow_geometry_time`
  and `erf.substep_volume_sampling`. Fork commits `5b8a38ad`, `827054ac`, `eff4101c`; records
  `review_substep_mt_2026-09-19/` and `review_devbranch_fixes_2026-10-29.md` in the fork's top
  folder.
- Gates: before the merge all six PASS on the laptop of October 2026 (Apple clang 21, Homebrew
  open-mpi and FFTW, the trap flag above). After the merge g5 and g6 are unchanged to the last bit
  (no acoustic sub-step, no stress on a moving surface), and g1 to g4 change. Their references were
  regenerated from the merged binaries (FORCE_REFS=1), and the binaries built from the tag itself
  reproduce the new references bit for bit. Size of the change in the final plotfile, relative:
  inviscid wave (g1, g2) u and w 3e-3, p' 2e-2; Stokes layer (g3) u and w 9e-3, p' 0.46 (D10 adds
  the isotropic first-cell stress to the small p' of this case); wall-modeled LES (g4) u 1e-6,
  w 3e-5, p' 8e-6 after 30 steps.
- Capability matrix: the compressible DNS, wall-resolved LES and wall-modeled LES over the moving
  wave now run the corrected sub-step and stress. The anelastic moving-mesh path has no acoustic
  sub-step, and its slow tendency already enters with the stage volume (the J_src weight of the
  refold), so only D10 reaches its viscous cases (Stokes DNS, WRLES smoke, WMLES mini): the
  first-cell pressure changes, the velocities at second order only. Their validation numbers above
  were measured before D10 and are re-measured when those cases run again.
- The open observation of `base_state_rebuild.md` (first-order fast error of u, w and theta at
  beta_s = 0) is explained: the w part was D1, the u and theta part is the dtau/2 lag of the
  forward-backward sub-step scheme (fork record `review_substep_mt_2026-09-19/REPORT.md`).
