# Wind over waves

These decks simulate air over a prescribed monochromatic water wave on a moving terrain-following
mesh, with the split-explicit (acoustic sub-step) time integration of ERF. They are the cases of the
manuscript "A split-explicit moving-mesh method for wind over ocean waves in the Energy Research and
Forecasting (ERF) model", submitted to the Journal of Computational Physics.

## Build

From `Exec/`:

    make -j8                              # all cases except the Stokes layer
    make -j8 USE_TERRAIN_VELOCITY=TRUE    # Stokes layer: executable ERF3d.*.TV.ex

`USE_TERRAIN_VELOCITY = TRUE` (CMake: `-DERF_ENABLE_TERRAIN_VELOCITY=ON`) sets the velocity of a
no-slip lower wall to the orbital velocity of the wave. Without it, a no-slip lower wall is at rest.

## Cases

| Deck | Case |
|---|---|
| `inputs_inviscid_wave` | traveling wave under a slip lid, no viscosity |
| `inputs_reflection_free_wave` | standing wave started from rest in a deep domain, no off-centering and no divergence damping |
| `inputs_stratified_wave` | the reflection-free wave over a stably stratified base state (`input_sounding_stratified`) |
| `inputs_stokes_layer` | viscous Stokes layer on the no-slip wave surface (TV build) |
| `inputs_wave_les` | turbulent boundary layer over the wave, driven by a pressure gradient, with the wave drag |
| `inputs_geostrophic_stable` | geostrophic boundary layer with a cooled surface (`input_sounding_geo`) |
| `inputs_geostrophic_neutral` | geostrophic boundary layer with zero surface heat flux (`input_sounding_geo_1400`) |
| `inputs_geostrophic_convective` | geostrophic boundary layer with a heated surface (`input_sounding_geo_1400`) |

Run a deck from a folder that contains its sounding file, if it uses one.

Studies of the manuscript:

- Time step and acoustic sub-step: change `erf.fixed_dt` and `erf.fixed_mri_dt_ratio` (M); the
  sub-step is `erf.fixed_dt / M`.
- Wave age of the turbulent boundary layer: for the phase speed c = sqrt(g lambda / 2 pi) = 8.84 m/s
  of the 50-m wave and the friction velocity u* = c / (c/u*), set `prob.init_ustar = u*` and
  `erf.abl_pressure_grad = -1.16 u*^2 / 200  0.  0.` (density 1.16 kg/m3, domain height 200 m).

## Inputs of the moving wave

| Input | Default | Meaning |
|---|---|---|
| `prob.custom_terrain_type` | | `MovingSineWave`: h = a sin(kx - omega t); `MovingSineWaveRest`: h = a sin(kx) (1 - cos(omega t)); `MovingSineWaveRamp`: h = a R(t) sin(kx - omega t), with R(t) = (1 - cos(omega t))/2 for omega t < pi and 1 after |
| `prob.Ampl` | 0 | amplitude a [m] |
| `prob.wavelength` | 100 | wavelength [m], k = 2 pi / wavelength |
| `prob.wave_c` | 0 | phase speed c [m/s]; omega = c k if c > 0, else sqrt(g k) |
| `prob.init_wave_field` | 1 | 1: initial velocity and pressure of the potential flow over the wave; 0: air at rest |
| `prob.init_ustar`, `prob.init_z0` | 0, 2e-4 | initial wind u = (u*/kappa) ln((z + z0)/z0) above the surface |
| `prob.init_pert_amp` | 0 | amplitude [m/s] of a smooth 3-D initial perturbation |
| `erf.most.wave_drag` | 0 | surface stress from the wind relative to the orbital velocity (Zhu et al. 2023, Mon. Wea. Rev. 151, Eqs. 26-29); needs `MovingSineWaveRamp` and `zlo.type = "surface_layer"`, and uses `erf.most.z0` and `erf.most.zref` |
| `erf.most.wave_drag_stability` | 0 | stability functions and surface heat flux at each point along the wave (Sullivan et al. 2014, J. Atmos. Sci. 71) |
| `erf.most.wave_drag_theta_mode` | `flux` | `flux`: `erf.most.wave_drag_surf_flux` [K m/s]; `surface_temp`: `erf.most.wave_drag_surf_temp` [K] and `erf.most.wave_drag_surf_temp_rate` [K/s] |
| `erf.most.wave_drag_z0h` | -1 | thermal roughness [m]; a negative value selects `erf.most.z0` |
| `erf.face_sync_int` | 1 | steps between the synchronization of the two copies of a periodic face velocity; 0: never |

The time integration on the moving mesh uses these evaluation times in the acoustic sub-steps. The
defaults are the scheme of the manuscript; the other values only show its conditions.

| Input | Default | Other value |
|---|---|---|
| `erf.substep_metric_sampling` | `midpoint` (moving mesh): metric terms of the pressure gradient and operator Jacobians at the sub-step midpoint | `endpoint`: at the sub-step start |
| `erf.grid_velocity_sampling` | `DiscreteDisplacement`: displacement over the sub-step | `IntervalStart`: rate over the preceding interval |
| `erf.slow_geometry_time` | `StageTime`: slow tendency on the mesh at the stage time | `StepStart`: on the mesh at the step start |
| `erf.substep_volume_sampling` | `Endpoint`: cell volume at the sub-step start on the sub-step-start state | `Midpoint`: mean of the start and end volumes |
| `erf.beta_d` | -1: 0.1, or 0 with midpoint metric terms and `erf.beta_s = 0` | divergence damping coefficient |
