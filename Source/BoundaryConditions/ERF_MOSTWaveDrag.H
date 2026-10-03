#ifndef ERF_MOST_WAVE_DRAG_H_
#define ERF_MOST_WAVE_DRAG_H_

#include <cmath>
#include <AMReX_Array4.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include "ERF_Constants.H"
#include "ERF_IndexDefines.H"

/**
 * Relative-velocity log-law wall drag for a turbulent ABL over a MOVING monochromatic wave.
 *
 * Implements Zhu et al. (2023, Mon. Wea. Rev., 151, 2877-2896), Eqs. 26-29 -- the WRF-LES
 * quadratic bottom drag modified so the stress is built from the wind RELATIVE to the water
 * (orbital) velocity:
 *      u_r = u - u_w                                                         (Eq 26)
 *      u_par = [ u_r + w_r (dh/dx) ] / sqrt( (dh/dx)^2 + 1 )                 (Eq 27a)
 *      v_par = [ v_r + w_r (dh/dy) ] / sqrt( (dh/dy)^2 + 1 )                 (Eq 27b)
 *      |u_r| = sqrt( u_r^2 + v_r^2 + w_r^2 )                                  (full 3-component rel. speed)
 *      tau13 = Cd |u_r| u_par ,   tau23 = Cd |u_r| v_par                      (Eq 28)
 *      Cd    = [ kappa / ln( zref / z0 ) ]^2 ,  zref = first cell-center ht   (Eq 29)
 *
 * This is OPT-IN (erf.most.wave_drag=1) and used ONLY with MovingFittedMesh + surface_layer and
 * with the MOST rotate path OFF (Eq 27 already projects onto the surface tangent; rotating again
 * would double-count the slope -- asserted in the SurfaceLayer ctor). When the gate is off, the
 * standard MOST stress (compute_u_flux/compute_v_flux) runs unchanged and all other cases are
 * byte-identical.
 *
 * The wave (a, k, omega, ramp) MUST match the moving terrain and the orbital-u_w BC, i.e.
 *   ERF_InitCustomTerrain.cpp  "MovingSineWaveRamp" :  h = a r(t) sin(kx - omega t)
 *   ERF_BoundaryConditionsXvel.cpp                   :  orbital u_w on a no-slip wall
 * with k = 2 pi / wavelength,  omega = sqrt(g k),  r(t) = 1/2 (1 - cos omega t) for omega t < pi.
 * (NB: WRF used a tanh(2 t / tau) ramp; here we deliberately use ERF's cosine ramp so the orbital
 *  drag and the moving-mesh motion z_t ramp up IN PHASE from the shared rest state.)
 *
 * Staggering (Arakawa C-grid), matching ERF_MOSTStress.H::compute_u_flux/compute_v_flux:
 *   tau13 lives on x-faces (co-located with u);  velx is local, vely/velz 4-pt-averaged to it.
 *   tau23 lives on y-faces (co-located with v);  vely is local, velx/velz 4-pt-averaged to it.
 * Sign: returns rho*<u'w'> = -rho Cd |u_r| u_par (negative -> a drag opposing the relative wind),
 * matching compute_u_flux which returns -rho u*^2 (...).
 */
namespace WaveDrag {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real ramp (const amrex::Real omega, const amrex::Real t)
    {
        // r(t) = 1/2 (1 - cos(omega t)) for omega t < pi, else 1  (== MovingSineWaveRamp)
        return (omega * t < amrex::Real(PI)) ?
               amrex::Real(0.5) * (amrex::Real(1.0) - std::cos(omega * t)) : amrex::Real(1.0);
    }

    // tau13 at the x-face (i,j,k). xu = x-coordinate of this u-node = prob_lo_x + i*dx.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    compute_wave_drag_u (const int& i, const int& j, const int& k,
                         const amrex::Array4<const amrex::Real>& cons_arr,
                         const amrex::Array4<const amrex::Real>& velx_arr,
                         const amrex::Array4<const amrex::Real>& vely_arr,
                         const amrex::Array4<const amrex::Real>& velz_arr,
                         const amrex::Real& xu,   const amrex::Real& time,
                         const amrex::Real& a,    const amrex::Real& kw,
                         const amrex::Real& omega,const amrex::Real& Cd)
    {
        const amrex::Real r    = ramp(omega, time);
        const amrex::Real ph   = kw * xu - omega * time;
        const amrex::Real u_w  =  a * r * omega * std::sin(ph);   // orbital u at surface
        const amrex::Real w_w  = -a * r * omega * std::cos(ph);   // orbital w at surface  (= dh/dt)
        const amrex::Real dhdx =  a * r * kw    * std::cos(ph);    // surface slope at the u-face

        // velx at the face; vely, velz averaged to the same x-face (mirror compute_u_flux)
        const amrex::Real uloc = velx_arr(i,j,k);
        const amrex::Real vloc = amrex::Real(0.25) * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                                     + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        const amrex::Real wloc = amrex::Real(0.25) * ( velz_arr(i  ,j,k) + velz_arr(i  ,j,k+1)
                                                     + velz_arr(i-1,j,k) + velz_arr(i-1,j,k+1) );

        const amrex::Real u_r = uloc - u_w;
        const amrex::Real v_r = vloc;            // v_wave = 0 for a wave travelling in x
        const amrex::Real w_r = wloc - w_w;

        const amrex::Real u_par = (u_r + w_r * dhdx) / std::sqrt(dhdx*dhdx + amrex::Real(1.0)); // Eq 27a
        // Eq 26/28: |u_r| is the FULL 3-component relative speed (u_r, v_r, w_r)
        const amrex::Real V0    = std::sqrt(u_r*u_r + v_r*v_r + w_r*w_r) + amrex::Real(1.0e-15);
        const amrex::Real rho   = amrex::Real(0.5) * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        return -rho * Cd * V0 * u_par;           // = rho<u'w'> (drag opposing the relative wind)
    }

    // tau23 at the y-face (i,j,k). xv = x-coordinate of this v-node = prob_lo_x + (i+1/2)*dx.
    // NOTE: uses +dh/dy in u_par (Eq 27b), correcting the -dh/dy sign in the WRF v_para (a typo that
    // is harmless for a wave travelling in x, dh/dy=0, but wrong for oblique/3-D waves).
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    compute_wave_drag_v (const int& i, const int& j, const int& k,
                         const amrex::Array4<const amrex::Real>& cons_arr,
                         const amrex::Array4<const amrex::Real>& velx_arr,
                         const amrex::Array4<const amrex::Real>& vely_arr,
                         const amrex::Array4<const amrex::Real>& velz_arr,
                         const amrex::Real& xv,   const amrex::Real& time,
                         const amrex::Real& a,    const amrex::Real& kw,
                         const amrex::Real& omega,const amrex::Real& Cd)
    {
        const amrex::Real r    = ramp(omega, time);
        const amrex::Real ph   = kw * xv - omega * time;
        const amrex::Real u_w  =  a * r * omega * std::sin(ph);
        const amrex::Real w_w  = -a * r * omega * std::cos(ph);
        // 2-D wave in x: v_wave = 0, dh/dy = 0  -> v_par = v_r (the projection term vanishes).
        const amrex::Real dhdy = amrex::Real(0.0);

        const amrex::Real vloc = vely_arr(i,j,k);
        const amrex::Real uloc = amrex::Real(0.25) * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                                     + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        const amrex::Real wloc = amrex::Real(0.25) * ( velz_arr(i,j  ,k) + velz_arr(i,j  ,k+1)
                                                     + velz_arr(i,j-1,k) + velz_arr(i,j-1,k+1) );

        const amrex::Real u_r = uloc - u_w;
        const amrex::Real v_r = vloc;            // v_wave = 0
        const amrex::Real w_r = wloc - w_w;

        const amrex::Real v_par = (v_r + w_r * dhdy) / std::sqrt(dhdy*dhdy + amrex::Real(1.0)); // Eq 27b (+dh/dy)
        // Eq 26/28: |u_r| is the FULL 3-component relative speed (u_r, v_r, w_r)
        const amrex::Real V0    = std::sqrt(u_r*u_r + v_r*v_r + w_r*w_r) + amrex::Real(1.0e-15);
        const amrex::Real rho   = amrex::Real(0.5) * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        return -rho * Cd * V0 * v_par;
    }

    // ==================================================================================
    // Gated stability extension (erf.most.wave_drag_stability=1, default OFF).
    //
    // Adds Monin-Obukhov stability corrections and a surface heat flux to the
    // relative-velocity drag, applied POINT BY POINT along the wave following the
    // wave-resolving LES practice of Sullivan, McWilliams & Patton (2014, JAS, 71,
    // 4001-4027, appendix): bulk formulas at the first grid level with the stability
    // functions of Large et al. (1994), no plane averaging.  The neutral functions
    // above are untouched; with the gate off none of the code below executes and all
    // existing cases are bit-identical.
    //
    //   u_*  = kappa |u_r| / [ ln(zref/z0)  - psi_m(zref/L) ]
    //   th_* = kappa (th1 - th_s) / [ ln(zref/z0h) - psi_h(zref/L) ]   (surface_temp mode)
    //   th_* = -Q / u_*                                                 (flux mode)
    //   L    = u_*^2 th1 / (kappa g th_*),  iterated to a fixed point per cell/face.
    //
    // Stability functions (Businger-Dyer family as given by Large et al. 1994):
    //   zeta < 0 : x = (1 - 16 zeta)^(1/4)
    //              psi_m = 2 ln((1+x)/2) + ln((1+x^2)/2) - 2 atan(x) + pi/2
    //              psi_h = 2 ln((1+x^2)/2)
    //   zeta >= 0: psi_m = psi_h = -5 zeta
    //
    // Momentum keeps Zhu's tangential projection and the full 3-component relative
    // speed; only Cd = [kappa/(ln(zref/z0) - psi_m)]^2 replaces the neutral Cd.  The
    // heat flux returned is rho*<th'w'> = -rho u_* th_* (the hfx3 convention of
    // ERF_MOSTStress.H).  In flux mode this reduces to +rho*Q identically.
    // ==================================================================================

    struct StabParams {
        int         mode;     // 0 = prescribed kinematic heat flux Q; 1 = prescribed surface temperature
        amrex::Real zref;     // first cell-center height above the surface
        amrex::Real z0;       // momentum roughness
        amrex::Real z0h;      // thermal roughness
        amrex::Real Q;        // kinematic surface heat flux <w'th'> [K m/s], mode 0
        amrex::Real theta_s;  // surface potential temperature [K], mode 1
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real psi_m_stab (const amrex::Real zeta)
    {
        if (zeta >= amrex::Real(0.0)) { return -amrex::Real(5.0) * zeta; }
        const amrex::Real x = std::pow(amrex::Real(1.0) - amrex::Real(16.0)*zeta, amrex::Real(0.25));
        return amrex::Real(2.0)*std::log(amrex::Real(0.5)*(amrex::Real(1.0)+x))
             + std::log(amrex::Real(0.5)*(amrex::Real(1.0)+x*x))
             - amrex::Real(2.0)*std::atan(x) + amrex::Real(0.5)*amrex::Real(PI);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real psi_h_stab (const amrex::Real zeta)
    {
        if (zeta >= amrex::Real(0.0)) { return -amrex::Real(5.0) * zeta; }
        const amrex::Real x = std::pow(amrex::Real(1.0) - amrex::Real(16.0)*zeta, amrex::Real(0.25));
        return amrex::Real(2.0)*std::log(amrex::Real(0.5)*(amrex::Real(1.0)+x*x));
    }

    // Fixed-point iteration for (u_*, th_*) from the local relative speed V0 and the
    // local first-cell theta.  zeta is clamped to [-10, 2] for robustness; eight
    // iterations converge to round-off at these parameters.  Returns u_*, th_*, and
    // psi_m at the converged zeta.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void iterate_ustar_tstar (const amrex::Real& V0, const amrex::Real& theta1,
                              const StabParams& sp,
                              amrex::Real& ustar, amrex::Real& tstar, amrex::Real& psim)
    {
        const amrex::Real ln_z0  = std::log(sp.zref / sp.z0);
        const amrex::Real ln_z0h = std::log(sp.zref / sp.z0h);
        amrex::Real zeta = amrex::Real(0.0);
        psim  = amrex::Real(0.0);
        ustar = KAPPA * V0 / ln_z0;
        tstar = amrex::Real(0.0);
        for (int it = 0; it < 8; ++it) {
            psim = psi_m_stab(zeta);
            const amrex::Real psih  = psi_h_stab(zeta);
            const amrex::Real den_m = amrex::max(ln_z0  - psim, amrex::Real(1.0));
            const amrex::Real den_h = amrex::max(ln_z0h - psih, amrex::Real(1.0));
            ustar = KAPPA * V0 / den_m;
            tstar = (sp.mode == 0) ? -sp.Q / amrex::max(ustar, amrex::Real(1.0e-6))
                                   : KAPPA * (theta1 - sp.theta_s) / den_h;
            if (std::abs(tstar) < amrex::Real(1.0e-10)) { zeta = amrex::Real(0.0); continue; }
            const amrex::Real L = ustar*ustar*theta1 / (KAPPA * CONST_GRAV * tstar);
            zeta = sp.zref / L;
            zeta = amrex::min(amrex::max(zeta, amrex::Real(-10.0)), amrex::Real(2.0));
        }
    }

    // Stability-corrected tau13 at the x-face (i,j,k).  Mirrors compute_wave_drag_u
    // with Cd built from the local Obukhov iteration.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    compute_wave_drag_u_stab (const int& i, const int& j, const int& k,
                              const amrex::Array4<const amrex::Real>& cons_arr,
                              const amrex::Array4<const amrex::Real>& velx_arr,
                              const amrex::Array4<const amrex::Real>& vely_arr,
                              const amrex::Array4<const amrex::Real>& velz_arr,
                              const amrex::Real& xu,   const amrex::Real& time,
                              const amrex::Real& a,    const amrex::Real& kw,
                              const amrex::Real& omega,const StabParams& sp)
    {
        const amrex::Real r    = ramp(omega, time);
        const amrex::Real ph   = kw * xu - omega * time;
        const amrex::Real u_w  =  a * r * omega * std::sin(ph);
        const amrex::Real w_w  = -a * r * omega * std::cos(ph);
        const amrex::Real dhdx =  a * r * kw    * std::cos(ph);

        const amrex::Real uloc = velx_arr(i,j,k);
        const amrex::Real vloc = amrex::Real(0.25) * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                                     + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        const amrex::Real wloc = amrex::Real(0.25) * ( velz_arr(i  ,j,k) + velz_arr(i  ,j,k+1)
                                                     + velz_arr(i-1,j,k) + velz_arr(i-1,j,k+1) );

        const amrex::Real u_r = uloc - u_w;
        const amrex::Real v_r = vloc;
        const amrex::Real w_r = wloc - w_w;

        const amrex::Real u_par = (u_r + w_r * dhdx) / std::sqrt(dhdx*dhdx + amrex::Real(1.0));
        const amrex::Real V0    = std::sqrt(u_r*u_r + v_r*v_r + w_r*w_r) + amrex::Real(1.0e-15);
        const amrex::Real rho    = amrex::Real(0.5) * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );
        const amrex::Real theta1 = amrex::Real(0.5) * ( cons_arr(i-1,j,k,RhoTheta_comp)/cons_arr(i-1,j,k,Rho_comp)
                                                      + cons_arr(i  ,j,k,RhoTheta_comp)/cons_arr(i  ,j,k,Rho_comp) );

        amrex::Real ustar, tstar, psim;
        iterate_ustar_tstar(V0, theta1, sp, ustar, tstar, psim);
        const amrex::Real dvis = amrex::max(std::log(sp.zref/sp.z0) - psim, amrex::Real(1.0));
        const amrex::Real Cd   = (KAPPA/dvis) * (KAPPA/dvis);

        return -rho * Cd * V0 * u_par;
    }

    // Stability-corrected tau23 at the y-face (i,j,k).  Mirrors compute_wave_drag_v.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    compute_wave_drag_v_stab (const int& i, const int& j, const int& k,
                              const amrex::Array4<const amrex::Real>& cons_arr,
                              const amrex::Array4<const amrex::Real>& velx_arr,
                              const amrex::Array4<const amrex::Real>& vely_arr,
                              const amrex::Array4<const amrex::Real>& velz_arr,
                              const amrex::Real& xv,   const amrex::Real& time,
                              const amrex::Real& a,    const amrex::Real& kw,
                              const amrex::Real& omega,const StabParams& sp)
    {
        const amrex::Real r    = ramp(omega, time);
        const amrex::Real ph   = kw * xv - omega * time;
        const amrex::Real u_w  =  a * r * omega * std::sin(ph);
        const amrex::Real w_w  = -a * r * omega * std::cos(ph);
        const amrex::Real dhdy = amrex::Real(0.0);

        const amrex::Real vloc = vely_arr(i,j,k);
        const amrex::Real uloc = amrex::Real(0.25) * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                                     + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        const amrex::Real wloc = amrex::Real(0.25) * ( velz_arr(i,j  ,k) + velz_arr(i,j  ,k+1)
                                                     + velz_arr(i,j-1,k) + velz_arr(i,j-1,k+1) );

        const amrex::Real u_r = uloc - u_w;
        const amrex::Real v_r = vloc;
        const amrex::Real w_r = wloc - w_w;

        const amrex::Real v_par = (v_r + w_r * dhdy) / std::sqrt(dhdy*dhdy + amrex::Real(1.0));
        const amrex::Real V0    = std::sqrt(u_r*u_r + v_r*v_r + w_r*w_r) + amrex::Real(1.0e-15);
        const amrex::Real rho    = amrex::Real(0.5) * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );
        const amrex::Real theta1 = amrex::Real(0.5) * ( cons_arr(i,j-1,k,RhoTheta_comp)/cons_arr(i,j-1,k,Rho_comp)
                                                      + cons_arr(i,j  ,k,RhoTheta_comp)/cons_arr(i,j  ,k,Rho_comp) );

        amrex::Real ustar, tstar, psim;
        iterate_ustar_tstar(V0, theta1, sp, ustar, tstar, psim);
        const amrex::Real dvis = amrex::max(std::log(sp.zref/sp.z0) - psim, amrex::Real(1.0));
        const amrex::Real Cd   = (KAPPA/dvis) * (KAPPA/dvis);

        return -rho * Cd * V0 * v_par;
    }

    // Surface heat flux rho*<th'w'> = -rho u_* th_* at the cell center (i,j,k).
    // xc = x-coordinate of the cell center.  In flux mode this returns +rho*Q exactly.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    compute_wave_t_flux (const int& i, const int& j, const int& k,
                         const amrex::Array4<const amrex::Real>& cons_arr,
                         const amrex::Array4<const amrex::Real>& velx_arr,
                         const amrex::Array4<const amrex::Real>& vely_arr,
                         const amrex::Array4<const amrex::Real>& velz_arr,
                         const amrex::Real& xc,   const amrex::Real& time,
                         const amrex::Real& a,    const amrex::Real& kw,
                         const amrex::Real& omega,const StabParams& sp)
    {
        const amrex::Real r    = ramp(omega, time);
        const amrex::Real ph   = kw * xc - omega * time;
        const amrex::Real u_w  =  a * r * omega * std::sin(ph);
        const amrex::Real w_w  = -a * r * omega * std::cos(ph);

        const amrex::Real uloc = amrex::Real(0.5) * ( velx_arr(i,j,k) + velx_arr(i+1,j,k) );
        const amrex::Real vloc = amrex::Real(0.5) * ( vely_arr(i,j,k) + vely_arr(i,j+1,k) );
        const amrex::Real wloc = amrex::Real(0.5) * ( velz_arr(i,j,k) + velz_arr(i,j,k+1) );

        const amrex::Real u_r = uloc - u_w;
        const amrex::Real v_r = vloc;
        const amrex::Real w_r = wloc - w_w;
        const amrex::Real V0  = std::sqrt(u_r*u_r + v_r*v_r + w_r*w_r) + amrex::Real(1.0e-15);

        const amrex::Real rho    = cons_arr(i,j,k,Rho_comp);
        const amrex::Real theta1 = cons_arr(i,j,k,RhoTheta_comp) / rho;

        amrex::Real ustar, tstar, psim;
        iterate_ustar_tstar(V0, theta1, sp, ustar, tstar, psim);

        return -rho * ustar * tstar;
    }

} // namespace WaveDrag
#endif
