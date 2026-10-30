#!/bin/bash
# Regression gates for the wind-over-wave code.
#
# The gates lock the behavior shown in the JCP manuscript (and the upstream anelastic
# baseline) to stored reference plotfiles. Run them before and after every source change
# on the development branch. New capabilities must keep every gate bit-identical when
# their input flags are off.
#
# Usage:
#   ./run_gates.sh              build both binaries, run all gates, compare against refs
#   ./run_gates.sh g1 g3        run only the named gates
#   ./run_gates.sh --make-refs  generate missing reference plotfiles (refuses to overwrite;
#                               set FORCE_REFS=1 to regenerate on purpose)
#   SKIP_BUILD=1 ./run_gates.sh skip the make step (binaries already current)
#
# PASS = final plotfile matches the stored reference within 1e-12 (abs and rel).
# With an unchanged code path this is bitwise identity; the tolerance only absorbs
# last-digit noise if a compiler reassociates.

set -u
cd "$(dirname "$0")"
GATES_DIR="$(pwd)"
EXEC_DIR="$GATES_DIR/../../Exec"
FCOMPARE="$GATES_DIR/../../Submodules/AMReX/Tools/Plotfile/fcompare.gnu.ex"
PLAIN_EXE="$EXEC_DIR/ERF3d.gnu.TEST.MPI.ex"
TV_EXE="$EXEC_DIR/ERF3d.gnu.TEST.MPI.TV.ex"
FFT_EXE="$EXEC_DIR/ERF3d.gnu.TEST.MPI.FFT.ex"
NRANKS=4

MAKE_REFS=0
GATE_LIST=()
for arg in "$@"; do
    if [ "$arg" = "--make-refs" ]; then MAKE_REFS=1; else GATE_LIST+=("$arg"); fi
done
if [ ${#GATE_LIST[@]} -eq 0 ]; then
    GATE_LIST=(g1_inviscid_wave_sub g2_inviscid_wave_endpoint g3_stokes_viscous_tv \
               g4_wmles_most_wavedrag g5_anelastic_densitycurrent g6_anelastic_moving_wave)
fi

# Expand short names (g1 -> g1_inviscid_wave_sub)
expand_gate () {
    case "$1" in
        g1) echo g1_inviscid_wave_sub ;;
        g2) echo g2_inviscid_wave_endpoint ;;
        g3) echo g3_stokes_viscous_tv ;;
        g4) echo g4_wmles_most_wavedrag ;;
        g5) echo g5_anelastic_densitycurrent ;;
        g6) echo g6_anelastic_moving_wave ;;
        *)  echo "$1" ;;
    esac
}

gate_binary () {
    case "$1" in
        g3_stokes_viscous_tv)      echo "$TV_EXE" ;;
        g6_anelastic_moving_wave)  echo "$FFT_EXE" ;;
        *)                         echo "$PLAIN_EXE" ;;
    esac
}

if [ ! -x "$FCOMPARE" ]; then
    echo "ERROR: fcompare not found at $FCOMPARE"
    echo "Build it with: make -C $(dirname "$FCOMPARE") programs=fcompare"
    exit 2
fi

if [ "${SKIP_BUILD:-0}" != "1" ]; then
    # USE_NETCDF=FALSE: the local NetCDF (4.10.1) has no parallel support; no gate needs it.
    echo "== Building the plain executable =="
    ( cd "$EXEC_DIR" && make -j8 USE_NETCDF=FALSE ) > "$GATES_DIR/_build_plain.log" 2>&1 \
        || { echo "BUILD FAILED (plain); see _build_plain.log"; exit 2; }
    echo "== Building the TERRAIN_VELOCITY executable =="
    ( cd "$EXEC_DIR" && make -j8 USE_NETCDF=FALSE TERRAIN_VELOCITY=TRUE ) > "$GATES_DIR/_build_tv.log" 2>&1 \
        || { echo "BUILD FAILED (TV); see _build_tv.log"; exit 2; }
    # FFTW for the FFT build (anelastic gate g6): FFTW_HOME from the environment, else ~/local
    # (the laptop of September 2026), else Homebrew.
    if [ -z "${FFTW_HOME:-}" ]; then
        if   [ -f "$HOME/local/include/fftw3.h" ]; then FFTW_HOME="$HOME/local";
        elif [ -f /opt/homebrew/include/fftw3.h ]; then FFTW_HOME=/opt/homebrew;
        else echo "BUILD FAILED (FFT): fftw3.h not found; set FFTW_HOME"; exit 2; fi
    fi
    echo "== Building the FFT executable (FFTW_HOME=$FFTW_HOME) =="
    ( cd "$EXEC_DIR" && make -j8 USE_NETCDF=FALSE USE_FFT=TRUE FFTW_HOME="$FFTW_HOME" ) > "$GATES_DIR/_build_fft.log" 2>&1 \
        || { echo "BUILD FAILED (FFT); see _build_fft.log"; exit 2; }
fi

FAIL=0
for g in "${GATE_LIST[@]}"; do
    gate="$(expand_gate "$g")"
    if [ ! -f "$GATES_DIR/$gate/inputs" ]; then
        echo "?? $gate: no such gate"; FAIL=1; continue
    fi
    exe="$(gate_binary "$gate")"
    if [ ! -x "$exe" ]; then
        echo "?? $gate: executable missing: $exe"; FAIL=1; continue
    fi

    work="$GATES_DIR/_work/$gate"
    rm -rf "$work"; mkdir -p "$work"
    ( cd "$work" && mpirun -np $NRANKS "$exe" "$GATES_DIR/$gate/inputs" ) \
        > "$work/run.log" 2>&1
    rc=$?
    last_plt="$(ls -d "$work"/plt[0-9]* 2>/dev/null | sort | tail -1)"
    if [ $rc -ne 0 ] || [ -z "$last_plt" ] || [ "$(basename "$last_plt")" = "plt00000" ]; then
        echo "FAIL $gate: run failed (rc=$rc); see $work/run.log"; FAIL=1; continue
    fi

    ref="$GATES_DIR/refs/$gate/$(basename "$last_plt")"
    if [ "$MAKE_REFS" = "1" ]; then
        if [ -d "$ref" ] && [ "${FORCE_REFS:-0}" != "1" ]; then
            echo "SKIP $gate: reference exists (set FORCE_REFS=1 to regenerate)"
        else
            mkdir -p "$GATES_DIR/refs/$gate"
            rm -rf "$ref"
            cp -r "$last_plt" "$ref"
            echo "REF  $gate: stored $(basename "$last_plt")"
        fi
        continue
    fi

    if [ ! -d "$ref" ]; then
        echo "?? $gate: no reference; run with --make-refs first"; FAIL=1; continue
    fi
    if "$FCOMPARE" --abs_tol 1e-12 --rel_tol 1e-12 "$ref" "$last_plt" \
        > "$work/fcompare.log" 2>&1; then
        echo "PASS $gate"
    else
        echo "FAIL $gate: differs from reference; see $work/fcompare.log"; FAIL=1
    fi
done

exit $FAIL
