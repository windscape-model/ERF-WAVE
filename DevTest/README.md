# DevTest

New-capability development area for the wind-over-wave code.
Everything new (tests, references, development notes) lives under this directory.
Nothing here touches `Exec/` run directories, `verification_runs/`, or `mrigark_paper/`.

Branch: `dev-wave-capabilities` (branched from `jcp-writing-pass`).
The JCP manuscript does not present any capability developed here.

## Protection rules

1. Every new feature sits behind an input flag whose default reproduces current behavior.
2. The gates in `gates/` must stay bit-identical to their stored references whenever the
   new flags are off. Run them before and after every source change.
3. Test runs write only under `DevTest/` (the `_work/` scratch areas).
4. The `jcp-writing-pass` branch stays untouched; merges go one way, from it into this branch.
5. Reference plotfiles are never regenerated silently. `run_gates.sh --make-refs` refuses to
   overwrite; a deliberate regeneration (`FORCE_REFS=1`) must be recorded in the git log.

## Gates (gates/run_gates.sh)

| Gate | Configuration | Protects |
|------|---------------|----------|
| g1_inviscid_wave_sub | inviscid moving wave, implicit substepping, midpoint metrics | moving-mesh metrics, GCL sampling, WS-RK3 (the JCP default) |
| g2_inviscid_wave_endpoint | same with endpoint metric sampling | the legacy sampling path and its input branch |
| g3_stokes_viscous_tv | viscous Stokes layer, orbital no-slip wall (TV build) | TERRAIN_VELOCITY wall BC, constant viscosity, FastImplicit sponge |
| g4_wmles_most_wavedrag | miniature wall-modeled LES over the wave | MOST wave drag, Smagorinsky + molecular viscosity, ABL forcing |
| g5_anelastic_densitycurrent | anelastic density current, flat static grid | the anelastic projection machinery the wave work extends |

Usage, from `DevTest/gates/`:

    ./run_gates.sh                 # build both executables, run all gates, compare
    ./run_gates.sh g1 g3           # selected gates only
    ./run_gates.sh --make-refs     # generate missing references (first time only)

PASS means the final plotfile matches the stored reference within 1e-12 on every field.
Builds here use `USE_NETCDF=FALSE` because the local NetCDF (4.10.1, `~/local`) has no
parallel support since the August environment change; no gate or wave case reads NetCDF.

## Layout

    gates/          regression gates + references (references are committed)
    notes/          development and design notes (ROADMAP.md is the plan of record)
    wrles_wave/     wall-resolved LES over the moving wave (compressible, TV build)
    anelastic_wave/ anelastic wave cases: static baseline + moving as-is record (FFT build)
