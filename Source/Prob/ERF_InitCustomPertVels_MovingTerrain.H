    ParmParse pp_prob("prob");
    Real Ampl = zero;        pp_prob.query("Ampl", Ampl);
    Real wavelength = one; pp_prob.get("wavelength", wavelength);

    // --- Equilibrium turbulent start (opt-in; default 0 => unchanged rest-start, so the inviscid
    //     and viscous moving-wave cases that share this prob are bit-identical). ---
    //   init_ustar    : friction velocity of the initial log-law mean wind u(z)=(u*/kappa) ln((z+z0)/z0),
    //                   z = height above the wavy surface. Set to the target u*=c/(c/u*) for the wave age
    //                   so the ABL starts near equilibrium instead of spinning up from rest (~F*t, slow).
    //   init_z0       : roughness for the initial profile (match erf.most.z0).
    //   init_pert_amp : amplitude [m/s] of a small, smooth, 3-D velocity perturbation that breaks the
    //                   spanwise symmetry and seeds transition to 3-D turbulence (a 2-D-symmetric start
    //                   over a 2-D wave would stay 2-D). Default 0 => no perturbation.
    Real init_ustar    = zero;          pp_prob.query("init_ustar",    init_ustar);
    Real init_z0       = Real(2.0e-4);  pp_prob.query("init_z0",       init_z0);
    Real init_pert_amp = zero;          pp_prob.query("init_pert_amp", init_pert_amp);
    // init_wave_field : 1 (default) imposes the potential-flow field of the traveling wave
    //                   (u_wave, w_wave below and the matching pressure perturbation in the
    //                   state initializer); 0 starts the atmosphere at rest, the consistent
    //                   choice when the surface itself starts flat and at rest.
    int  init_wave_field = 1;           pp_prob.query("init_wave_field", init_wave_field);
    const Real Ampl_field = (init_wave_field != 0) ? Ampl : zero;

    Real H           = geomdata.ProbHi()[2];
    Real Lx          = geomdata.ProbHi()[0] - geomdata.ProbLo()[0];
    Real Ly          = geomdata.ProbHi()[1] - geomdata.ProbLo()[1];
    Real kp          = two * PI / wavelength;
    Real g           = CONST_GRAV;
    Real omega       = std::sqrt(g * kp);
    Real wave_c      = zero;            pp_prob.query("wave_c", wave_c);   // prescribed phase speed (optional)
    if (wave_c > zero) { omega = wave_c * kp; }
    Real l_kappa     = KAPPA;

    // Set the x-velocity:  analytic wave orbital + log-law mean wind + spanwise-breaking seed.
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();

        const Real x = prob_lo[0] + i * dx[0];
        const Real y = prob_lo[1] + (j + myhalf) * dx[1];
        Real z = fourth * (z_nd(i,j,k) + z_nd(i,j+1,k) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1));

        Real z_base = Ampl * std::sin(kp * x);
        z -= z_base;                                        // height above the wavy surface

        Real fac     = std::cosh( kp * (z - H) ) / std::sinh(kp * H);
        Real u_wave  = -Ampl_field * omega * fac * std::sin(kp * x);

        Real u_mean  = (init_ustar > zero && z > zero) ?
                       (init_ustar / l_kappa) * std::log((z + init_z0) / init_z0) : zero;

        // smooth multi-mode 3-D seed, zero at surface and top (sin(pi z/H))
        Real fz   = (z > zero && z < H) ? std::sin(PI * z / H) : zero;
        Real u_pert = init_pert_amp * fz *
                      ( std::sin(two*PI*Real(3.0)*x/Lx) * std::cos(two*PI*Real(2.0)*y/Ly)
                      + std::sin(two*PI*Real(5.0)*x/Lx) * std::cos(two*PI*Real(4.0)*y/Ly) );

        x_vel_pert(i, j, k) = u_wave + u_mean + u_pert;
    });

    // Set the y-velocity: spanwise-breaking seed only (zero spanwise mean).
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];
        const Real y = prob_lo[1] + j * dx[1];
        Real z = fourth * (z_nd(i,j,k) + z_nd(i+1,j,k) + z_nd(i,j,k+1) + z_nd(i+1,j,k+1));
        Real z_base = Ampl * std::sin(kp * x);
        z -= z_base;

        Real fz   = (z > zero && z < H) ? std::sin(PI * z / H) : zero;
        y_vel_pert(i, j, k) = init_pert_amp * fz *
                      ( std::cos(two*PI*Real(3.0)*x/Lx) * std::sin(two*PI*Real(2.0)*y/Ly)
                      + std::cos(two*PI*Real(5.0)*x/Lx) * std::sin(two*PI*Real(4.0)*y/Ly) );
    });

    // Set the z-velocity from the wave (impenetrable) condition + a small seed away from boundaries.
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx         = geomdata.CellSize();

        Real x   = prob_lo[0] + (i + myhalf) * dx[0];
        Real y   = prob_lo[1] + (j + myhalf) * dx[1];
        Real z   = fourth * ( z_nd(i,j,k) + z_nd(i+1,j,k) + z_nd(i,j+1,k) + z_nd(i+1,j+1,k) );
        Real z_base = Ampl * std::sin(kp * x);

        z -= z_base;
        Real fac    = std::sinh( kp * (z - H) ) / std::sinh(kp * H);
        Real w_wave = Ampl_field * omega * fac * std::cos(kp * x);

        Real fz     = (z > zero && z < H) ? std::sin(PI * z / H) : zero;
        Real w_pert = init_pert_amp * Real(0.5) * fz *
                      std::sin(two*PI*Real(4.0)*x/Lx) * std::sin(two*PI*Real(3.0)*y/Ly);

        z_vel_pert(i, j, k) = w_wave + w_pert;
    });
