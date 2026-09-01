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
| Anelastic, moving wave | BLOCKED: projection not moving-mesh consistent (stage B is this development) | after stage B | after stage B |

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
2. `erf.vert_implicit=1` with the moving wall: the implicit vertical-diffusion path is
   inconsistent with the moving-terrain metric update (documented in
   `Exec/visc_test/SETUP_NOTES.md`: fast crash; earlier measurement: ~20 percent Stokes
   amplitude error). Wall-resolved grids want implicit vertical diffusion, so this is
   worth fixing behind a flag. Entry points: `Source/Diffusion/ERF_ImplicitDiff_T.cpp`,
   `Source/TimeIntegration/ERF_ImplicitPre.H` / `ERF_ImplicitPost.H`.
3. Check one stress detail found while reading: in `ERF_ComputeStress_T.cpp` the
   metric cross-term correction to tau13 uses an average built from `mu_turb` only
   (search `mu_bar`); confirm whether the molecular part belongs there too and add it
   if so (matters only for resolved runs over steep waves).
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
- Stage A tests (this session, FFT build, 20 steps each, `anelastic_wave/`):
  - static wave: runs clean, projected divergence 3e-13.
  - moving wave AS-IS: runs, no assert stops it, solver reports divergence 5e-14 —
    but of the WRONG (static-mesh) constraint. Silent inconsistency, see below.

### The three defects to fix (stage B)

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
code paths), so after stage B this is configuration + testing, not development:
anelastic Stokes DNS (stage B item 3 doubles as this), anelastic WRLES smoke
(wrles_wave deck + anelastic=1), anelastic WMLES mini (g4 + anelastic=1).

## Build variants (GNU make, from Exec/)

    make -j8 USE_NETCDF=FALSE                          -> ERF3d.gnu.TEST.MPI.ex      (gates g1,g2,g4,g5)
    make -j8 USE_NETCDF=FALSE TERRAIN_VELOCITY=TRUE    -> ERF3d.gnu.TEST.MPI.TV.ex   (gate g3, WRLES smoke)
    CPATH=$HOME/local/include LIBRARY_PATH=$HOME/local/lib \
    make -j8 USE_NETCDF=FALSE USE_FFT=TRUE             -> ERF3d.gnu.TEST.MPI.FFT.ex  (anelastic wave work)

USE_NETCDF=FALSE because the local NetCDF (4.10.1 in ~/local) lost parallel support
(`netcdf_par.h` gone); no wave or gate case reads NetCDF. The July NetCDF-enabled
executables are preserved as `Exec/ERF3d.jul25-nc-keep.ex` / `.TV.ex`.
The FFT variant links FFTW from ~/local via CPATH/LIBRARY_PATH; each variant has its
own object directory, so the three executables coexist.

## Standing protections

- Every source change runs `DevTest/gates/run_gates.sh` before commit; the five gates
  must stay bit-identical with new flags off.
- New features default off; behavior changes only behind named input flags.
- Test runs write only under `DevTest/`; `Exec/` run directories, `verification_runs/`,
  and `mrigark_paper/` stay untouched.
- `jcp-writing-pass` receives nothing from this branch before submission.
