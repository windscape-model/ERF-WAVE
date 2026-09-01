#!/usr/bin/env python3
"""
Validate the anelastic moving-wave run against the linear potential-flow solution
under a rigid lid.

Bottom h = a sin(kx - w t) (after the amplitude ramp), lid at z = H. The linear
incompressible solution is
    w(x,z,t) = Re{ W(z) e^{i(kx - w t)} },   W(z) = -a w sinh(k(H-z))/sinh(kH)
    u'(x,z,t) = Re{ U(z) e^{i(kx - w t)} },  U(z) = +i a w cosh(k(H-z))/sinh(kH)
The script extracts the streamwise fundamental of the spanwise-averaged w and u at
every level, pairs it with the level's mean physical height, and reports amplitude
and phase against the analytic profile. It also reports max|rho - rho0| (the
anelastic density pin) and the fundamental of the lowest-level w against the wall
kinematics.

Plotfile reader copied from Exec/visc_test/stokes_compare.py (native, numpy only).
Usage: check_wave_solution.py <plotfile>
"""
import cmath, math, os, re, sys
import numpy as np

A_WAVE  = 0.16
LAMBDA  = 100.0
H_LID   = 100.0
GRAV    = 9.81
KP      = 2.0 * math.pi / LAMBDA
OMEGA   = math.sqrt(GRAV * KP)

def _boxes(cell_h):
    with open(cell_h) as f:
        txt = f.read()
    ncomp = int(txt.splitlines()[2])
    boxes = re.findall(r'\(\((\d+),(\d+),(\d+)\) \((\d+),(\d+),(\d+)\) '
                       r'\((\d+),(\d+),(\d+)\)\)', txt)
    fabs = re.findall(r'FabOnDisk: (\S+) (\d+)', txt)
    if len(boxes) != len(fabs):
        boxes = boxes[-len(fabs):]
    return ncomp, [(tuple(int(x) for x in b[0:3]),
                    tuple(int(x) for x in b[3:6]), fn, int(off))
                   for b, (fn, off) in zip(boxes, fabs)]

def read_plotfile(plotfile, varnames):
    with open(os.path.join(plotfile, 'Header')) as f:
        lines = f.read().splitlines()
    nvar = int(lines[1])
    names = lines[2:2 + nvar]
    time = float(lines[3 + nvar])
    idx = 5 + nvar
    mm = re.match(r'\(\((\d+),(\d+),(\d+)\) \((\d+),(\d+),(\d+)\)',
                  lines[idx + 3].strip())
    dims = tuple(int(mm.group(i + 4)) - int(mm.group(i + 1)) + 1
                 for i in range(3))
    dx = [float(x) for x in lines[idx + 5].split()]
    ncomp, fabs = _boxes(os.path.join(plotfile, 'Level_0', 'Cell_H'))
    assert ncomp == nvar
    data = {v: np.empty(dims) for v in varnames}
    lvl = os.path.join(plotfile, 'Level_0')
    for lo, hi, fn, off in fabs:
        n = tuple(hi[d] - lo[d] + 1 for d in range(3))
        with open(os.path.join(lvl, fn), 'rb') as f:
            f.seek(off)
            f.readline()
            raw = np.fromfile(f, dtype='<f8', count=n[0]*n[1]*n[2]*ncomp)
        raw = raw.reshape(ncomp, n[2], n[1], n[0])
        for v in varnames:
            c = names.index(v)
            data[v][lo[0]:hi[0]+1, lo[1]:hi[1]+1, lo[2]:hi[2]+1] = \
                raw[c].transpose(2, 1, 0)
    return time, dx, data

def fundamentals(plt):
    """Return (time, zlev, w_hat, u_hat, rho_pin): per-level streamwise fundamentals
    of the spanwise-averaged w and u (no time rotation applied here)."""
    time, dx, d = read_plotfile(plt, ['x_velocity', 'z_velocity', 'density',
                                      'dens_hse', 'z_phys'])
    nx, ny, nz = d['z_velocity'].shape
    x = (np.arange(nx) + 0.5) * dx[0]
    phase_x = np.exp(-1j * KP * x)
    wbar = d['z_velocity'].mean(axis=1)
    ubar = d['x_velocity'].mean(axis=1)
    zlev = d['z_phys'].mean(axis=1).mean(axis=0)
    w_hat = 2.0 / nx * (wbar * phase_x[:, None]).sum(axis=0)
    u_hat = 2.0 / nx * (ubar * phase_x[:, None]).sum(axis=0)
    rho_pin = np.max(np.abs(d['density'] - d['dens_hse']))
    return time, zlev, w_hat, u_hat, rho_pin

def fit_main(plts):
    """Multi-snapshot mode: fit m_j(z) = W_pot(z) + W_frozen(z) e^{i omega t_j}
    per level, where m_j = w_hat(t_j) e^{i omega t_j}. The oscillating part W_pot
    is compared against the analytic profile; W_frozen collects any quasi-steady
    residual k-mode (startup vorticity), which does not oscillate at omega."""
    times, zlev, mw, mu = [], None, [], []
    rho_pin = 0.0
    for p in plts:
        t, z, w_hat, u_hat, rp = fundamentals(p)
        rot = cmath.exp(1j * OMEGA * t)
        times.append(t); zlev = z
        mw.append(w_hat * rot); mu.append(u_hat * rot)
        rho_pin = max(rho_pin, rp)
    times = np.array(times)
    mw = np.array(mw); mu = np.array(mu)           # (nt, nz)
    E = np.exp(1j * OMEGA * times)                  # (nt,)
    A = np.stack([np.ones_like(E), E], axis=1)      # (nt, 2)
    coefw, _, _, _ = np.linalg.lstsq(A, mw, rcond=None)
    coefu, _, _, _ = np.linalg.lstsq(A, mu, rcond=None)
    resw = np.abs(mw - A @ coefw).max(axis=0)
    denom = math.sinh(KP * H_LID)
    print(f"# fit over {len(plts)} snapshots, t = {times.min():.2f} .. {times.max():.2f} s")
    print(f"# density pin: max|rho - rho0| = {rho_pin:.3e}")
    print("#  z[m]    |w_pot|   |w_ana|   rel.err  ph.err[deg]  |w_frozen|  fit.res   |u_pot|   |u_ana|   rel.err")
    werrs, uerrs, perrs = [], [], []
    for k in range(len(zlev)):
        z = zlev[k]
        w_ana = -A_WAVE * OMEGA * math.sinh(KP * (H_LID - z)) / denom
        u_ana = 1j * A_WAVE * OMEGA * math.cosh(KP * (H_LID - z)) / denom
        if abs(w_ana) < 0.02 * A_WAVE * OMEGA:
            continue
        w_pot, w_fro = coefw[0, k], coefw[1, k]
        u_pot = coefu[0, k]
        werr = abs(abs(w_pot) - abs(w_ana)) / abs(w_ana)
        uerr = abs(abs(u_pot) - abs(u_ana)) / abs(u_ana)
        perr = math.degrees(cmath.phase(w_pot / w_ana))
        werrs.append(werr); uerrs.append(uerr); perrs.append(abs(perr))
        print(f"{z:8.3f}  {abs(w_pot):.5f}  {abs(w_ana):.5f}  {werr:8.4f}  {perr:9.3f}   "
              f"{abs(w_fro):.6f}  {resw[k]:.2e}  {abs(u_pot):.5f}  {abs(u_ana):.5f}  {uerr:8.4f}")
    print(f"# oscillating w amplitude rel.err: max {max(werrs):.4f}  mean {np.mean(werrs):.4f}")
    print(f"# oscillating u amplitude rel.err: max {max(uerrs):.4f}  mean {np.mean(uerrs):.4f}")
    print(f"# oscillating w phase err [deg]  : max {max(perrs):.3f}  mean {np.mean(perrs):.3f}")

def main():
    if len(sys.argv) > 2:
        fit_main(sys.argv[1:])
        return
    plt = sys.argv[1]
    time, dx, d = read_plotfile(plt, ['x_velocity', 'z_velocity', 'density',
                                      'dens_hse', 'z_phys'])
    nx, ny, nz = d['z_velocity'].shape
    x = (np.arange(nx) + 0.5) * dx[0]
    phase_x = np.exp(-1j * KP * x)
    rot_t   = cmath.exp(1j * OMEGA * time)

    # anelastic density pin
    rho_pin = np.max(np.abs(d['density'] - d['dens_hse']))
    rho_ref = np.max(np.abs(d['dens_hse']))

    # spanwise average, then streamwise fundamental per level
    wbar = d['z_velocity'].mean(axis=1)   # (nx, nz)
    ubar = d['x_velocity'].mean(axis=1)
    zbar = d['z_phys'].mean(axis=(0, 1) if d['z_phys'].ndim == 3 else 0)
    zlev = d['z_phys'].mean(axis=1).mean(axis=0)  # (nz,) mean height per level

    denom = math.sinh(KP * H_LID)
    print(f"# {plt}: t = {time:.6f} s  (omega t / 2pi = {OMEGA*time/2/math.pi:.4f} periods)")
    print(f"# density pin: max|rho - rho0| = {rho_pin:.3e}  (rho0 max {rho_ref:.4f})")
    print("#  z[m]    |w_num|    |w_ana|   rel.err    phase.err[deg]  |u_num|   |u_ana|   rel.err")
    werrs, uerrs, perrs = [], [], []
    for k in range(nz):
        z = zlev[k]
        w_hat = 2.0 / nx * np.sum(wbar[:, k] * phase_x) * rot_t
        u_hat = 2.0 / nx * np.sum(ubar[:, k] * phase_x) * rot_t
        w_ana = -A_WAVE * OMEGA * math.sinh(KP * (H_LID - z)) / denom
        u_ana = 1j * A_WAVE * OMEGA * math.cosh(KP * (H_LID - z)) / denom
        if abs(w_ana) < 0.02 * A_WAVE * OMEGA:
            continue
        werr = abs(abs(w_hat) - abs(w_ana)) / abs(w_ana)
        uerr = abs(abs(u_hat) - abs(u_ana)) / abs(u_ana)
        perr = math.degrees(cmath.phase(w_hat / w_ana))
        werrs.append(werr); uerrs.append(uerr); perrs.append(abs(perr))
        print(f"{z:8.3f}  {abs(w_hat):.5f}  {abs(w_ana):.5f}  {werr:8.4f}   {perr:10.3f}   "
              f"{abs(u_hat):.5f}  {abs(u_ana):.5f}  {uerr:8.4f}")
    print(f"# w amplitude rel.err: max {max(werrs):.4f}  mean {np.mean(werrs):.4f}")
    print(f"# u amplitude rel.err: max {max(uerrs):.4f}  mean {np.mean(uerrs):.4f}")
    print(f"# w phase err [deg]  : max {max(perrs):.3f}  mean {np.mean(perrs):.3f}")

if __name__ == '__main__':
    main()
