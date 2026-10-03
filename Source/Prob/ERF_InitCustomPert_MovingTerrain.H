
    // Parse params
    ParmParse pp("prob");
    Real Ampl = zero       ; pp.query("Ampl", Ampl);
    Real wavelength = one; pp.get("wavelength", wavelength);
    // Same switch as the velocity initializer: 0 = no potential-flow pressure perturbation
    int  init_wave_field = 1; pp.query("init_wave_field", init_wave_field);
    const Real Ampl_field = (init_wave_field != 0) ? Ampl : zero;

    Real H           = geomdata.ProbHi()[2];
    Real kp          = two * PI / wavelength;
    Real g           = CONST_GRAV;
    Real omega       = std::sqrt(g * kp);
    Real wave_c      = zero;            pp.query("wave_c", wave_c);        // prescribed phase speed (optional)
    if (wave_c > zero) { omega = wave_c * kp; }

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + myhalf) * dx[0];

        Real z = z_cc(i,j,k);
        Real z_base = Ampl * std::sin(kp * x);
        z -= z_base;

        Real fac     = std::cosh( kp * (z - H) ) / std::sinh(kp * H);
        Real p_prime = -(Ampl_field * omega * omega / kp) * fac * std::sin(kp * x) * r_hse(i,j,k);
        Real p_total = p_hse(i,j,k) + p_prime;

        // Define (rho theta) given pprime
        state_pert(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p_total) - getRhoThetagivenP(p_hse(i,j,k));

        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = state_pert(i,j,k,Rho_comp);
    });
