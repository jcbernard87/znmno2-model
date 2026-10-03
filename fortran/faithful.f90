!> Faithful ports of the original research programs (docs/model.md section 9, docs/deviations.md):
!>   line = 'charge'   ZnMn02_v3 (N = 5, R2 and R3, ZHS-solubility pH)
!>   line = 'phcell'   ZnMn02_v2.1_GITT_no_probe (N = 6, H+ transported, R5)
!> Every operation is done in the original's order and precision, defects included (single-precision
!> constants and spline, the original's limits and explicit solid updates), so the output file is the
!> original's byte for byte when the program is built without fused multiply-adds (-ffp-contract=off).
!> The block solve is bandsolver's legacy pivot.
!>
!> SPDX-License-Identifier: BSD-3-Clause
module zn_faithful
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use bandsolver_kernel, only: band_solve, BAND_OK, BAND_NON_FINITE, PIVOT_LEGACY, SINGULAR_EXACT
    implicit none
    private
    public :: run_faithful

    ! ------------------------------------------------------------------ the line and its inputs
    logical :: ph                            ! .true. for the pH-cell line
    integer :: N, NJ, S
    ! template values of the original run generator (single-precision literals)
    real(sp) :: rxnk_2_, rxnk_3_, frac_zmcx_, frac_zmcmax_, xmax_, applied_current_, porosity_, volfrac_mno2_, &
                stated_mass_loading_, zhs_ksp_, rxnk_5_, fraction_kmno2_

    ! ------------------------------------------------------------------ constants (as the original stored them)
    real(dp) :: Rigc, Temp, Fconst, diff_Zn, diff_Mn, diff_SO4, diff_H, pH_init
    real(dp) :: molar_mass_KMn8O16, density_KMn8O16, molar_mass_ZHS, density_ZHS, density_ZMC
    real(dp) :: molar_mass_MnO2, density_MnO2, xmax_c, V_at_Zmin, V_at_Zmax, BL_thickness, sigma_sep
    real(dp) :: Rxn2_K, Rxn3_K, Rxn5_K, K_sp, Zmin, Zmax, ratio_initial, cbulk_Zn, cbulk_Mn, cbulk_H
    real(dp) :: sigma, xmax, Area_CS, len_sep, eps_sep, eps, volfrac_inert, AM_Grams, AM_Grams_sim
    real(dp) :: KMn8O16_init, ZHS_init, ZMCx_init, ZMC_max_init, MnO2_init, Phi_1_init, molar_mass_ZMC_max
    real(dp) :: resevoir_scaling, current_target
    real(dp) :: z_ion(6), c_initial(6), diff_ion(6)
    integer, parameter :: RAMP_ITERS = 1000
    real(dp), parameter :: RAMP_DELT = 1.0e-9_dp, DELT_NOMINAL = 1.0_dp, WRITE_DENSITY = 10.0_dp

    ! ------------------------------------------------------------------ spline tables (float32)
    real(sp) :: zn_pts(153), mn_pts(153), coef(149,149), ocp_t(51), ocp_c(50,4)

    ! ------------------------------------------------------------------ state
    real(dp), allocatable :: c(:,:), delC(:,:), xx(:), delx(:), diff_term(:,:), mig_term(:,:)
    real(dp), allocatable :: KMn(:), ZHS(:), ZMCx(:), ZMCm(:), MnO2(:), por(:), tort(:)
    real(dp), allocatable :: a_K(:), a_ZHS(:), a_ZMCx(:), a_ZMCm(:), a_MnO2(:), znm(:), mnm(:), ratio(:), MW(:)
    real(dp), allocatable :: pH_phreeqc(:), pH_ZHS(:), pH_standard(:)
    real(dp), allocatable :: Ab(:,:,:), Bb(:,:,:), Db(:,:,:), Gb(:,:)
    real(dp) :: time, delT, current, c_density, c_specific, mAhg, anode_pot, ramp_initial
    integer  :: ramp_count, last_write_time
    logical  :: ramp_on
    character(len=1) :: state
    logical  :: last_nan                     ! a NaN in the last assembled coefficients (check_good_vals)
    integer  :: ou

contains

    ! ================================================================== entry point
    subroutine run_faithful(line, nml_unit, data_dir, out_file, reason)
        character(len=*), intent(in) :: line, data_dir, out_file
        integer, intent(in) :: nml_unit
        character(len=*), intent(out) :: reason
        real(sp) :: rxnk_2, rxnk_3, frac_zmcx, frac_zmcmax, xmax_t, applied_current, porosity, volfrac_mno2, &
                    stated_mass_loading, zhs_ksp, rxnk_5, fraction_kmno2
        integer :: ios
        namelist /faithful/ rxnk_2, rxnk_3, frac_zmcx, frac_zmcmax, xmax_t, applied_current, porosity, volfrac_mno2, &
                            stated_mass_loading, zhs_ksp, rxnk_5, fraction_kmno2
        ! illustrative defaults (not fitted to a cell)
        rxnk_2 = -8.5; rxnk_3 = -7.5; frac_zmcx = 0.4; frac_zmcmax = 0.03; xmax_t = 0.022
        applied_current = 0.000121; porosity = 0.8; volfrac_mno2 = 0.07; stated_mass_loading = 0.00121; zhs_ksp = 30.0
        rxnk_5 = -9.0; fraction_kmno2 = 0.7
        rewind(nml_unit)
        read(nml_unit, nml=faithful, iostat=ios)
        if (ios > 0) stop 'error reading &faithful'
        rxnk_2_ = rxnk_2; rxnk_3_ = rxnk_3; frac_zmcx_ = frac_zmcx; frac_zmcmax_ = frac_zmcmax; xmax_ = xmax_t
        applied_current_ = applied_current; porosity_ = porosity; volfrac_mno2_ = volfrac_mno2
        stated_mass_loading_ = stated_mass_loading; zhs_ksp_ = zhs_ksp; rxnk_5_ = rxnk_5; fraction_kmno2_ = fraction_kmno2
        ph = (line == 'phcell')
        if (.not. ph .and. line /= 'charge') stop 'faithful line must be charge or phcell'
        call read_tables(data_dir)
        call constants()
        call initial_condition()
        open(newunit=ou, file=out_file, status='replace', action='write')
        call main_loop(reason)
        close(ou)
    end subroutine run_faithful

    real(sp) function pow10_sp(x)
        real(sp), intent(in) :: x
        ! 10.0**(x) with a single-precision literal x: gfortran folds it to the correctly rounded float32
        if (x == aint(x)) then
            pow10_sp = real(10.0_dp**nint(x), sp)
        else
            pow10_sp = real(10.0_dp**real(x, dp), sp)
        end if
    end function pow10_sp

    subroutine read_tables(dir)
        character(len=*), intent(in) :: dir
        integer :: u, i, n1, n2, n3, n4
        open(newunit=u, file=trim(dir)//'/ph_spline_phreeqc.txt', status='old', action='read')
        call skip_comments(u)
        read(u, *) n1, n2, n3, n4
        read(u, *) zn_pts, mn_pts
        read(u, *) (coef(i,:), i = 1, 149)
        close(u)
        open(newunit=u, file=trim(dir)//'/r3_ocp_spline.txt', status='old', action='read')
        call skip_comments(u)
        read(u, *) n1
        read(u, *) ocp_t
        read(u, *) (ocp_c(i,:), i = 1, 50)
        close(u)
    end subroutine read_tables

    subroutine skip_comments(u)
        integer, intent(in) :: u
        character(len=512) :: line
        do
            read(u, '(A)') line
            if (adjustl(line(1:1)) /= '#') exit
        end do
        backspace(u)
    end subroutine skip_comments

    ! ================================================================== constants (sampling_variables)
    subroutine constants()
        real(dp) :: tmp
        Rigc = 8.314; Temp = 298; Fconst = 96485
        diff_Zn = 7.15d-6; diff_Mn = 6.88d-6; diff_SO4 = 1.07d-5; diff_H = 9.0d-5
        pH_init = 5.5
        molar_mass_KMn8O16 = 734.59; density_KMn8O16 = 5.03
        molar_mass_ZHS = 549.819; density_ZHS = 2.67; density_ZMC = 5.0
        molar_mass_MnO2 = 86.9368; density_MnO2 = 5.03
        xmax_c = 200.0d-5
        BL_thickness = 0.5*1.0d-4
        sigma_sep = 1.0d-20
        Area_CS = 0.178134094
        len_sep = 600.0*1.0d-4
        eps_sep = 0.9
        z_ion = 0; c_initial = 0; diff_ion = 0
        if (.not. ph) then
            N = 5; NJ = 122; S = 51
            V_at_Zmin = 1.75; V_at_Zmax = 1.45
            Rxn2_K = pow10_sp(rxnk_2_); Rxn3_K = pow10_sp(rxnk_3_)
            if (zhs_ksp_ == aint(zhs_ksp_)) then
                K_sp = real(10.0_dp**nint(zhs_ksp_), sp)
            else
                K_sp = pow10_sp(zhs_ksp_)
            end if
            Zmin = 0.2; Zmax = 0.5
            ratio_initial = Zmin*1.001
            cbulk_Zn = 0.002; cbulk_Mn = 0.00005
            sigma = 0.1
            xmax = xmax_
            current_target = applied_current_
            eps = porosity_
            ZMCx_init = real(volfrac_mno2_*frac_zmcx_, dp)
            ZMC_max_init = real(volfrac_mno2_*frac_zmcmax_, dp)
            tmp = 0.0001
            volfrac_inert = 1.0 - eps - tmp - real(0.0001, dp) - ZMCx_init - ZMC_max_init
            AM_Grams = stated_mass_loading_
            AM_Grams_sim = ZMCx_init*Area_CS*xmax*density_ZMC
            KMn8O16_init = tmp*density_KMn8O16
            ZHS_init = real(0.0001, dp)*density_ZHS
            ZMCx_init = ZMCx_init*density_ZMC
            ZMC_max_init = ZMC_max_init*density_ZMC
            MnO2_init = 0
            Phi_1_init = 1.7
            molar_mass_ZMC_max = 65.38*Zmax + 54.93 + (15.999*2)
            z_ion(3) = 2.0; z_ion(4) = 0.0; z_ion(5) = -2.0          ! z_ion(4) = z_Mn0, zero (M-1)
            c_initial(3) = cbulk_Zn; c_initial(4) = cbulk_Mn
            c_initial(5) = -1*((z_ion(3)*c_initial(3))+(z_ion(4)*c_initial(4)))/z_ion(5)
            diff_ion(3) = diff_Zn; diff_ion(4) = diff_Mn; diff_ion(5) = diff_SO4
            resevoir_scaling = (real(0.0001, dp)+(len_sep*Area_CS*eps_sep)+(xmax*Area_CS*eps))/(len_sep*Area_CS*eps_sep)
            len_sep = len_sep*resevoir_scaling
        else
            N = 6; NJ = 134; S = 63
            Rxn5_K = pow10_sp(rxnk_5_)
            K_sp = 7.0d-26
            Zmin = 0.35; Zmax = 0.65
            ratio_initial = Zmax*0.9999
            cbulk_Zn = 0.002; cbulk_Mn = 0.00005; cbulk_H = 0.000543
            if ((-1.0*log10(1000*cbulk_H)) <= pH_init) pH_init = -1.0*log10(1000*cbulk_H)
            sigma = 0.01
            xmax = 0.0218
            current_target = 0.107/1000.0
            eps = 0.815
            MnO2_init = real(0.0616*fraction_kmno2_, dp)
            tmp = 0.00001                                              ! volfrac_ZMC_max
            volfrac_inert = 1.0 - eps - real(0.00001, dp) - real(0.000001, dp) - tmp*0.01 - tmp - MnO2_init
            AM_Grams = 0.00108
            AM_Grams_sim = tmp*0.01*Area_CS*xmax*density_ZMC
            KMn8O16_init = real(0.00001, dp)*density_KMn8O16
            MnO2_init = MnO2_init*density_MnO2
            ZHS_init = real(0.000001, dp)*density_ZHS
            ZMCx_init = tmp*0.01*density_ZMC
            ZMC_max_init = tmp*density_ZMC
            Phi_1_init = 1.1
            molar_mass_ZMC_max = 65.38*Zmax + 54.93 + (15.999*2)
            z_ion(3) = 2.0; z_ion(4) = 2.0; z_ion(5) = -2.0; z_ion(6) = 1.0
            c_initial(3) = cbulk_Zn; c_initial(4) = cbulk_Mn; c_initial(6) = cbulk_H
            c_initial(5) = -1*((z_ion(3)*c_initial(3))+(z_ion(4)*c_initial(4))+(z_ion(6)*c_initial(6)))/z_ion(5)
            diff_ion(3) = diff_Zn; diff_ion(4) = diff_Mn; diff_ion(5) = diff_SO4; diff_ion(6) = diff_H
            resevoir_scaling = 1.0                                      ! negative scaling, reset to 1 by the original
        end if
    end subroutine constants

    ! ================================================================== initial condition
    subroutine initial_condition()
        real(dp) :: h_sep, h_cath
        integer :: j, ic
        allocate(c(N,NJ), delC(N,NJ), xx(NJ), delx(NJ), diff_term(N,NJ), mig_term(N,NJ))
        allocate(KMn(NJ), ZHS(NJ), ZMCx(NJ), ZMCm(NJ), MnO2(NJ), por(NJ), tort(NJ), a_K(NJ), a_ZHS(NJ), a_ZMCx(NJ), &
                 a_ZMCm(NJ), a_MnO2(NJ), znm(NJ), mnm(NJ), ratio(NJ), MW(NJ), pH_phreeqc(NJ), pH_ZHS(NJ), pH_standard(NJ))
        allocate(Ab(N,N,NJ), Bb(N,N,NJ), Db(N,N,NJ), Gb(N,NJ))
        h_sep = len_sep/float(S-2)
        h_cath = xmax/float(NJ-S-1)
        xx = 0; delx = 0
        do j = 1, NJ
            if (j == 1) then
                xx(j) = 0.0
            else if (j < S) then
                xx(j) = h_sep*float(j-1) - h_sep/2.0
            else if (j == S) then
                xx(j) = len_sep
            else if (j == NJ) then
                xx(j) = xmax + len_sep
            else
                xx(j) = len_sep + h_cath*float(j-S) - h_cath/2.0
            end if
        end do
        do j = 2, NJ-1
            if (j < S) then
                delx(j) = h_sep
            else if (j > S) then
                delx(j) = h_cath
            end if
        end do
        delx(1) = 0; delx(S) = 0; delx(NJ) = 0
        KMn = 0; ZHS = 0; ZMCx = 0; ZMCm = 0; MnO2 = 0; a_K = 0; a_ZHS = 0; a_ZMCx = 0; a_ZMCm = 0; a_MnO2 = 0
        znm = 0; mnm = 0; ratio = 0; MW = 0
        pH_phreeqc = pH_init; pH_ZHS = pH_init; pH_standard = pH_init
        do j = 1, NJ
            c(1,j) = Phi_1_init
            c(2,j) = 0
            do ic = 3, N
                c(ic,j) = c_initial(ic)
            end do
            if (j < S) then
                por(j) = eps_sep
                tort(j) = 2*por(j)**(-0.5)
            else
                KMn(j) = KMn8O16_init
                MnO2(j) = MnO2_init
                ZHS(j) = ZHS_init
                ZMCx(j) = ZMCx_init
                ZMCm(j) = ZMC_max_init
                por(j) = eps
                tort(j) = 2*por(j)**(-0.5)
                a_K(j) = 3.0*KMn(j)/(density_KMn8O16*(xmax_c))
                if (ph) a_MnO2(j) = 3.0*MnO2(j)/(density_MnO2*(xmax_c))
                a_ZHS(j) = 3.0*ZHS(j)/(density_ZHS*(xmax_c))
                a_ZMCx(j) = 3.0*ZMCx(j)/(density_ZMC*(xmax_c))
                a_ZMCm(j) = 3.0*ZMCm(j)/(density_ZMC*(xmax_c))
                ratio(j) = ratio_initial
                MW(j) = 65.38*ratio(j) + 54.938 + (15.999*2)
                znm(j) = ZMCx(j)*ratio(j)/MW(j)
                mnm(j) = ZMCx(j)/MW(j)
            end if
            call transport_terms(j)
        end do
        delC = 0
        time = 0; delT = DELT_NOMINAL; current = 0; c_density = 0; c_specific = 0; mAhg = 0
        anode_pot = 0; ramp_on = .true.; ramp_count = 0; ramp_initial = 0; last_write_time = 0
        if (current_target >= 0.0) then
            state = 'D'
        else
            state = 'C'
        end if
    end subroutine initial_condition

    subroutine transport_terms(j)
        integer, intent(in) :: j
        integer :: ic
        do ic = 3, N
            if (j < S .and. .not. ph) then
                diff_term(ic,j) = (resevoir_scaling**2)*por(j)*diff_ion(ic)/tort(j)
                mig_term(ic,j) = (resevoir_scaling**2)*por(j)*z_ion(ic)*diff_ion(ic)*Fconst/(Rigc*Temp*tort(j))
            else
                diff_term(ic,j) = por(j)*diff_ion(ic)/tort(j)
                mig_term(ic,j) = por(j)*z_ion(ic)*diff_ion(ic)*Fconst/(Rigc*Temp*tort(j))
            end if
        end do
    end subroutine transport_terms

    ! ================================================================== pH
    recursive function basis(i, k, x, t) result(bf)
        integer, intent(in) :: i, k
        real(sp), intent(in) :: x
        real(sp), intent(in) :: t(:)
        real(sp) :: bf, a, b
        if (i + k > size(t)) then                ! the original reads past the knot array here
            if (ieee_is_nan(x)) then
                bf = ieee_value(bf, ieee_quiet_nan)
                return
            end if
            stop 'pH spline knot read out of bounds: not reproducible'
        end if
        if (k == 1) then
            if (t(i) <= x .and. x < t(i+1)) then
                bf = 1.0
            else
                bf = 0.0
            end if
        else
            if (t(i+k-1) /= t(i)) then
                a = (x - t(i))/(t(i+k-1) - t(i))
            else
                a = 0.0
            end if
            if (t(i+k) /= t(i+1)) then
                b = (t(i+k) - x)/(t(i+k) - t(i+1))
            else
                b = 0.0
            end if
            bf = a*basis(i, k-1, x, t) + b*basis(i+1, k-1, x, t)
        end if
    end function basis

    integer function find_interval(x, knots)
        real(sp), intent(in) :: x, knots(:)
        integer :: low, high, mid
        low = 1; high = size(knots)
        do while (low < high)
            mid = (low + high)/2
            if (x < knots(mid)) then
                high = mid
            else
                low = mid + 1
            end if
        end do
        find_interval = max(min(low - 1, size(knots) - 2), 0)
    end function find_interval

    real(sp) function manual_spline(conc1, conc2)
        real(sp), intent(in) :: conc1, conc2
        real(sp) :: zn, mn, lg_zn, lg_mn, bx(4), by(4)
        integer :: ix, iy, i, j, ci, cj
        integer, parameter :: kx = 3, ky = 3
        zn = conc2; mn = conc1                   ! swapped in the original (undoes a transposed read)
        lg_zn = log10(zn); lg_mn = log10(mn)
        ix = find_interval(lg_zn, zn_pts); iy = find_interval(lg_mn, mn_pts)
        do i = 0, kx
            if ((ix - kx + i) >= 1 .and. (ix - kx + i) <= 153 - kx) then
                bx(i+1) = basis(ix - kx + i, kx + 1, lg_zn, zn_pts)
            else
                bx(i+1) = 0.0
            end if
        end do
        do i = 0, ky
            if ((iy - ky + i) >= 1 .and. (iy - ky + i) <= 153 - ky) then
                by(i+1) = basis(iy - ky + i, ky + 1, lg_mn, mn_pts)
            else
                by(i+1) = 0.0
            end if
        end do
        manual_spline = 0.0
        do i = 1, 4
            do j = 1, 4
                ci = ix - kx + i - 1; cj = iy - ky + j - 1
                if (ci >= 1 .and. ci <= 149 .and. cj >= 1 .and. cj <= 149) then
                    manual_spline = manual_spline + coef(ci, cj)*bx(i)*by(j)
                end if
            end do
        end do
    end function manual_spline

    real(sp) function eval_pH(zn, mn, j)
        real(sp), intent(in) :: zn, mn
        integer, intent(in) :: j
        real(sp) :: zn_m, mn_m, c_so4, ph_ksp
        zn_m = zn*1000.0; mn_m = mn*1000.0
        c_so4 = zn_m + mn_m
        ph_ksp = -log10((zn_m**4*c_so4/K_sp)**(1.0/6.0))
        eval_pH = manual_spline(zn_m, mn_m)
        if (j > S) then
            if (ph_ksp < eval_pH) eval_pH = ph_ksp
        end if
    end function eval_pH

    real(dp) function zhs_pH(c_zn, c_so4)
        real(dp), intent(in) :: c_zn, c_so4
        real(dp) :: zn, so4, c_h
        zn = c_zn*1000.0d0; so4 = c_so4*1000.0d0
        if (ph) then
            c_h = ((zn**4)*so4/K_sp)**(-1.0d0/6.0d0)
        else
            c_h = ((zn**4)*so4/K_sp)**(1.0d0/6.0d0)
        end if
        zhs_pH = -log10(c_h)
    end function zhs_pH

    real(dp) function zn_anode_pot(c_zn)
        real(dp), intent(in) :: c_zn
        real(dp) :: zn_ref, n_e
        zn_ref = 0.001; n_e = 2.0
        zn_anode_pot = -0.762 + (Rigc*Temp/(n_e*Fconst))*log(c_zn/zn_ref)
    end function zn_anode_pot

    real(dp) function r3_ocp_spline(theta)
        real(dp), intent(in) :: theta
        integer :: ind
        real(dp) :: d, p1, p2, p3, p4
        ind = 1
        do while (ind <= 51)
            if (theta <= ocp_t(ind)) exit
            ind = ind + 1
        end do
        ind = ind - 1
        if (ieee_is_nan(theta)) then
            r3_ocp_spline = theta
            return
        end if
        if (ind < 1 .or. ind > 50) stop 'R3 OCP spline read out of bounds (M-17): not reproducible'
        p1 = ocp_c(ind,1); p2 = ocp_c(ind,2); p3 = ocp_c(ind,3); p4 = ocp_c(ind,4)
        d = theta - ocp_t(ind)
        r3_ocp_spline = p1*d**3 + p2*d**2 + p3*d + p4
        r3_ocp_spline = ((V_at_Zmin - V_at_Zmax)*r3_ocp_spline) + V_at_Zmax
    end function r3_ocp_spline

    ! ================================================================== reactions (charge line)
    subroutine reaction_2(p1, p2, c3, c4, c5, j, out)
        real(dp), intent(in) :: p1, p2, c3, c4, c5
        integer, intent(in) :: j
        real(dp), intent(out) :: out(13)
        real(dp) :: xreact, pH_func, alpha_a, alpha_c, rk, n_e, OCP, exi, eta, area, irxn
        real(dp) :: c03_ref, c04_ref, c05_ref, Solid_consumption_limit, M
        real(dp) :: dZn, dMn, dH, dZMC, dSO4, dZHS
        c03_ref = 0.002; c04_ref = 0.0001; c05_ref = 0.0021; Solid_consumption_limit = 0.1
        out = 0
        xreact = Zmax
        pH_func = zhs_pH(c3, c5)
        alpha_a = 0.5
        alpha_c = 1.0 - alpha_a
        rk = Rxn2_K
        n_e = 2.0*(1.0 - xreact)
        OCP = (1.78+0.76) &
            + ((Rigc*Temp/(2*Fconst))*(1.0d0*(-1)*log(c3/c03_ref))) &
            + ((Rigc*Temp/(2*Fconst))*(1.0d0*(-2.0)*log(c4/c04_ref))) &
            - 0.0592*(8.0/2.0)*pH_func
        exi = Fconst*rk &
            *((c3/c03_ref)**(1.0d0*((8.0/3)-xreact)*alpha_a/n_e)) &
            *((c4/c04_ref)**(1.0d0*(-1.0)*alpha_c/n_e)) &
            *((c5/c05_ref)**(1.0d0*(2.0/3.0)*alpha_a/n_e))
        eta = p1 - p2 - OCP
        if (state == 'C') then
            area = a_ZHS(j)
        else
            area = a_ZMCm(j)
        end if
        irxn = area*exi*(exp(alpha_a*Fconst*eta/(Rigc*Temp)) - exp(-alpha_c*Fconst*eta/(Rigc*Temp)))
        M = molar_mass_ZMC_max
        call rates(.true.)
        if ((-1.0*dMn/area) > (c4*diff_Mn/BL_thickness)) then
            irxn = ((c4*diff_Mn/BL_thickness)/((-1.0*dMn)/area))*irxn
            call rates(.true.)
        end if
        if ((-1.0*dZn/area) > (c3*diff_Zn/BL_thickness)) then
            irxn = ((c3*diff_Zn/BL_thickness)/((-1.0*dZn)/area))*irxn
            call rates(.true.)
        end if
        if ((-1.0*dSO4/area) > (c5*diff_SO4/BL_thickness)) then
            irxn = ((c5*diff_SO4/BL_thickness)/((-1.0*dSO4)/area))*irxn
            call rates(.true.)
        end if
        if ((-1.0*(dZMC*M*delT)/ZMCm(j)) >= Solid_consumption_limit) then
            irxn = irxn*abs(ZMCm(j)/(dZMC*M*delT))*Solid_consumption_limit
            call rates(.false.)
        end if
        if ((-1.0*(dZHS*M*delT)/ZHS(j)) >= Solid_consumption_limit) then
            irxn = irxn*abs(ZHS(j)/(dZHS*M*delT))*Solid_consumption_limit
            call rates(.false.)
        end if
        if ((ZMCm(j) + (dZMC*M*delT)) <= 0.0) then
            irxn = irxn*abs(ZMCm(j)/(dZMC*M*delT))*0.1
            call rates(.true.)
        end if
        if ((ZHS(j) + (dZHS*M*delT)) <= 0.0) then
            irxn = irxn*abs(ZHS(j)/(dZHS*M*delT))*0.1
            call rates(.true.)
        end if
        out(1) = OCP; out(2) = eta; out(3) = exi; out(4) = irxn
        out(5) = dZn; out(6) = dMn; out(7) = dSO4; out(8) = dH; out(10) = dZHS; out(11) = dZMC
    contains
        subroutine rates(with_h)
            logical, intent(in) :: with_h
            dZn = ((8.0/3.0) - xreact)*irxn/(n_e*Fconst)
            dMn = -irxn/(n_e*Fconst)
            if (with_h) dH = 4*irxn/(n_e*Fconst)
            dZMC = irxn/(n_e*Fconst)
            dSO4 = (2.0/3.0)*irxn/(n_e*Fconst)
            dZHS = (-2.0/3.0)*irxn/(n_e*Fconst)
        end subroutine rates
    end subroutine reaction_2

    subroutine reaction_3(p1, p2, c3, c4, c5, j, out)
        real(dp), intent(in) :: p1, p2, c3, c4, c5
        integer, intent(in) :: j
        real(dp), intent(out) :: out(13)
        real(dp) :: pH_func, alpha_a, alpha_c, rk, n_e, theta, nernst, OCP, exi, eta, area, irxn, dZn, c03_ref
        real(sp) :: Vint                         ! implicitly REAL(4) in the original
        c03_ref = 0.002
        out = 0
        pH_func = zhs_pH(c3, c5)                 ! computed and unused in the original
        alpha_a = 0.5
        alpha_c = 1.0 - alpha_a
        rk = Rxn3_K
        n_e = 2.0
        theta = (ratio(j) - Zmin)/(Zmax - Zmin)
        if (theta <= (-0.1)) then
            Vint = 1.8
        else
            Vint = r3_ocp_spline(real(theta, 8)) - 0.762
        end if
        nernst = ((Rigc*Temp/(n_e*Fconst))*(1.0*log(c3/c03_ref)))
        OCP = nernst + Vint - anode_pot
        if (ratio(j) < Zmax) then
            exi = Fconst*rk*((c3/c03_ref)**(0.5*alpha_c))*(ratio(j)**(alpha_a))*((Zmax - ratio(j))**(alpha_c))
        else
            exi = 1.0d-15
        end if
        eta = p1 - p2 - OCP
        area = a_ZMCx(j)
        irxn = area*exi*(exp(alpha_a*Fconst*eta/(Rigc*Temp)) - exp(-alpha_c*Fconst*eta/(Rigc*Temp)))
        dZn = irxn/(n_e*Fconst)
        if ((-1.0*dZn/area) > (c3*diff_Zn/BL_thickness)) then
            irxn = ((c5*diff_Zn/BL_thickness)/abs(dZn))*irxn       ! c5: the original's typo (M-6)
            dZn = irxn/(n_e*Fconst)
        end if
        if ((znm(j) + (-1.0*dZn*delT)) <= 0.0) then
            irxn = irxn*abs(znm(j)/(-1.0*dZn*delT))*0.1
            dZn = irxn/(n_e*Fconst)
        end if
        out(1) = OCP; out(2) = eta; out(3) = exi; out(4) = irxn; out(5) = dZn
        out(12) = -1.0*dZn
    end subroutine reaction_3

    ! ================================================================== reaction 5 (pH-cell line)
    subroutine reaction_5(p1, p2, c3, c4, c5, c6, j, out)
        real(dp), intent(in) :: p1, p2, c3, c4, c5, c6
        integer, intent(in) :: j
        real(dp), intent(out) :: out(13)
        real(dp), parameter :: c04_ref = 0.0001d0, Solid_consumption_limit = 0.1d0
        real(dp) :: pH_std, pH_zhs_r, alpha_a, alpha_c, Uref, Uanode, rk, n_e, area
        real(dp) :: pH_func, OCP, exi, eta, irxn, dMn, dH, dZn, dSO4, dZHS, dMnO2
        real(dp) :: H_post, ZHS_H, H_needed, H_rate, r_ZHS
        out = 0
        pH_std = -log10(1000.0d0*c6)
        pH_zhs_r = zhs_pH(c3, c5)
        alpha_a = 0.5d0
        alpha_c = 1.0d0 - alpha_a
        Uref = 1.2225d0
        Uanode = anode_pot
        rk = Rxn5_K
        n_e = 2.0d0
        area = a_MnO2(j)
        if (pH_std < pH_zhs_r) then
            pH_func = pH_std
            OCP = (Uref - Uanode) + (Rigc*Temp/(n_e*Fconst))*(-log(c4/c04_ref)) - 0.0592d0*(4.0d0/n_e)*pH_func
            exi = Fconst*rk*((c4/c04_ref)**(-alpha_c*1.0d0/n_e))
            eta = p1 - p2 - OCP
            irxn = area*exi*(exp(alpha_a*Fconst*eta/(Rigc*Temp)) - exp(-alpha_c*Fconst*eta/(Rigc*Temp)))
            dMn = -irxn/(n_e*Fconst)
            dH = 4.0d0*irxn/(n_e*Fconst)
            dMnO2 = irxn/(n_e*Fconst)
            dZHS = 0.0d0; dZn = 0.0d0; dSO4 = 0.0d0
            if ((-dMn/area) > (c4*diff_Mn/BL_thickness)) then
                irxn = ((c4*diff_Mn/BL_thickness)/((-dMn)/area))*irxn
                dMn = -irxn/(n_e*Fconst); dH = 4.0d0*irxn/(n_e*Fconst); dMnO2 = irxn/(n_e*Fconst)
            end if
            if ((-dMnO2*molar_mass_MnO2*delT)/MnO2(j) >= Solid_consumption_limit) then
                irxn = irxn*abs(MnO2(j)/(dMnO2*molar_mass_MnO2*delT))*Solid_consumption_limit
                dMn = -irxn/(n_e*Fconst); dH = 4.0d0*irxn/(n_e*Fconst); dMnO2 = irxn/(n_e*Fconst)
            end if
            H_post = (c6*por(j)*delx(j)*Area_CS) + (dH*delT)
            ZHS_H = 10.0d0**(-pH_zhs_r)/1000.0d0
            if ((H_post < ZHS_H) .and. (1.1*pH_std >= pH_zhs_r)) then
                H_needed = ZHS_H - H_post
                H_rate = H_needed/delT
                r_ZHS = H_rate/6.0d0
                dH = dH + H_rate
                dZn = dZn - 4.0d0*r_ZHS
                dSO4 = dSO4 - 1.0d0*r_ZHS
                dZHS = dZHS + r_ZHS
            end if
        else
            pH_func = pH_zhs_r
            OCP = (Uref - Uanode) + (Rigc*Temp/(n_e*Fconst))*(-log(c4/c04_ref)) - 0.0592d0*(4.0d0/n_e)*pH_func
            exi = Fconst*rk*((c4/c04_ref)**(-alpha_c*1.0d0/n_e))
            eta = p1 - p2 - OCP
            irxn = area*exi*(exp(alpha_a*Fconst*eta/(Rigc*Temp)) - exp(-alpha_c*Fconst*eta/(Rigc*Temp)))
            call rates5()
            if ((-dZn/area) > (c3*diff_Zn/BL_thickness)) then
                irxn = ((c3*diff_Zn/BL_thickness)/((-dZn)/area))*irxn
                call rates5()
            end if
            if ((-dSO4/area) > (c5*diff_SO4/BL_thickness)) then
                irxn = ((c5*diff_SO4/BL_thickness)/((-dSO4)/area))*irxn
                call rates5()
            end if
            if ((-dMnO2*molar_mass_MnO2*delT)/MnO2(j) >= Solid_consumption_limit) then
                irxn = irxn*abs(MnO2(j)/(dMnO2*molar_mass_MnO2*delT))*Solid_consumption_limit
                call rates5()
            end if
            if ((-dZHS*molar_mass_MnO2*delT)/ZHS(j) >= Solid_consumption_limit) then
                irxn = irxn*abs(ZHS(j)/(dZHS*molar_mass_MnO2*delT))*Solid_consumption_limit
                call rates5()
            end if
        end if
        out(1) = OCP; out(2) = eta; out(3) = exi; out(4) = irxn
        out(5) = dZn; out(6) = dMn; out(7) = dSO4; out(8) = dH; out(10) = dZHS; out(13) = dMnO2
    contains
        subroutine rates5()
            dMn = -irxn/(n_e*Fconst)
            dMnO2 = irxn/(n_e*Fconst)
            dH = 0.0d0
            dZHS = -(4.0d0*irxn/(n_e*Fconst))/6.0d0
            dZn = (8.0d0/3.0d0)*irxn/(n_e*Fconst)
            dSO4 = (2.0d0/3.0d0)*irxn/(n_e*Fconst)
        end subroutine rates5
    end subroutine reaction_5

    subroutine react_tot(v, j, tot)
        real(dp), intent(in) :: v(:)
        integer, intent(in) :: j
        real(dp), intent(out) :: tot(13)
        real(dp) :: r2(13), r3(13), r5(13), acc
        integer :: n, k
        tot = 0
        if (.not. ph) then
            call reaction_2(v(1), v(2), v(3), v(4), v(5), j, r2)
            call reaction_3(v(1), v(2), v(3), v(4), v(5), j, r3)
            do n = 4, 12
                acc = 0.0
                acc = 0.0 + acc                  ! R1 (off)
                acc = r2(n) + acc
                acc = r3(n) + acc
                tot(n) = acc
            end do
        else
            call reaction_5(v(1), v(2), v(3), v(4), v(5), v(6), j, r5)
            do n = 4, 13
                acc = 0.0
                do k = 1, 4
                    acc = 0.0 + acc
                end do
                acc = r5(n) + acc
                do k = 1, 4
                    acc = 0.0 + acc
                end do
                tot(n) = acc
            end do
        end if
    end subroutine react_tot

    subroutine drx_dc(v, j, d)
        real(dp), intent(in) :: v(:)
        integer, intent(in) :: j
        real(dp), intent(out) :: d(:,:)                ! (N variables, N-1 outputs: irxn and the species)
        real(dp) :: st(6), a(6), b(6), t1(13), t2(13)
        integer :: m, nout
        nout = N - 1
        st(1:N) = v(1:N)*0.001
        if (.not. ph) then
            st(2) = v(1)*0.001                           ! the phi2 step uses phi1 (M-7)
        else
            where (st(1:N) == 0.0d0) st(1:N) = 1.0d-6
        end if
        do m = 1, N
            a(1:N) = v(1:N); a(m) = v(m) + st(m)
            if (m >= 3 .and. v(m) <= st(m)) then
                call react_tot(a(1:N), j, t1)
                call react_tot(v(1:N), j, t2)
                d(m,1:nout) = (t1(4:3+nout) - t2(4:3+nout))/(st(m))
            else
                b(1:N) = v(1:N); b(m) = v(m) - st(m)
                call react_tot(a(1:N), j, t1)
                call react_tot(b(1:N), j, t2)
                d(m,1:nout) = (t1(4:3+nout) - t2(4:3+nout))/(2.0*st(m))
            end if
        end do
    end subroutine drx_dc

    ! ================================================================== fillmat + ABDGXY
    subroutine assemble()
        real(dp) :: dE(N,N), dW(N,N), fE(N,N), fW(N,N), rj(N,N), smG(N), rxn(13), drx(N,N-1)
        real(dp) :: cE(N), cW(N), dcdxE(N), dcdxW(N), alphaE, alphaW, betaE, betaW, pW, pE, p
        real(dp) :: mW(N), mE(N), dfW(N), dfE(N), v(N)
        integer :: j, ic, iq, iv
        alphaE = 0; alphaW = 0; betaE = 0; betaW = 0
        last_nan = .false.
        Ab = 0; Bb = 0; Db = 0; Gb = 0
        do j = 1, NJ
            dE = 0; dW = 0; fE = 0; fW = 0; rj = 0; smG = 0; cE = 0; cW = 0; dcdxE = 0; dcdxW = 0
            v = c(:,j)
            if (j == 1) then
                alphaE = delx(j)/(delx(j+1) + delx(j))
                betaE = 2.0/(delx(j) + delx(j+1))
                do ic = 1, N
                    cE(ic) = alphaE*c(ic,j+1) + (1.d0 - alphaE)*c(ic,j)
                    dcdxE(ic) = betaE*(c(ic,j+1) - c(ic,j))
                end do
                dE(1,1) = -1.0
                smG(1) = -(dE(1,1)*dcdxE(1))
                smG(2) = en_row(v)
                do ic = 3, N
                    rj(2,ic) = z_ion(ic)
                end do
                dE(3,3) = -1.0*diff_term(3,j)
                dE(3,2) = -1.0*mig_term(3,j)*cE(3)
                fE(3,3) = -1.0*mig_term(3,j)*dcdxE(2)
                smG(3) = -1.0*c_density/(z_ion(3)*Fconst) + (dE(3,3)*dcdxE(3) + fE(3,3)*cE(3))
                do ic = 4, N-1
                    dE(ic,ic) = -1.0*diff_term(ic,j)
                    dE(ic,2) = -1.0*mig_term(ic,j)*cE(ic)
                    fE(ic,ic) = -1.0*mig_term(ic,j)*dcdxE(2)
                    smG(ic) = 0.0 + (dE(ic,ic)*dcdxE(ic) + fE(ic,ic)*cE(ic))
                end do
                smG(N) = 0.0 - c(2,j)
                rj(N,2) = 1.0
                Bb(:,:,j) = rj - (1.d0 - alphaE)*fE + betaE*dE
                Db(:,:,j) = -alphaE*fE - betaE*dE
                Gb(:,j) = smG
                cycle
            end if
            alphaW = delx(j-1)/(delx(j-1) + delx(j))
            betaW = 2.0/(delx(j-1) + delx(j))
            if (j < NJ) then
                alphaE = delx(j)/(delx(j+1) + delx(j))
                betaE = 2.0/(delx(j) + delx(j+1))
            end if
            do ic = 1, N
                cW(ic) = alphaW*c(ic,j) + (1.d0 - alphaW)*c(ic,j-1)
                dcdxW(ic) = betaW*(c(ic,j) - c(ic,j-1))
                if (j < NJ) then
                    cE(ic) = alphaE*c(ic,j+1) + (1.d0 - alphaE)*c(ic,j)
                    dcdxE(ic) = betaE*(c(ic,j+1) - c(ic,j))
                end if
            end do
            if (j == S) then
                dE(1,1) = -1.0
                smG(1) = -(dE(1,1)*dcdxE(1))
                smG(2) = en_row(v)
                do ic = 3, N
                    rj(2,ic) = z_ion(ic)
                end do
                do ic = 3, N
                    dW(ic,ic) = -1*diff_term(ic,j-1)
                    dE(ic,ic) = -1*diff_term(ic,j+1)
                    fW(ic,ic) = -1.0*mig_term(ic,j-1)*dcdxW(2)
                    fE(ic,ic) = -1*mig_term(ic,j+1)*dcdxE(2)
                    dW(ic,2) = -1*mig_term(ic,j-1)*cW(ic)
                    dE(ic,2) = -1*mig_term(ic,j+1)*cE(ic)
                    smG(ic) = 0.0 - (fW(ic,ic)*cW(ic) + dW(ic,ic)*dcdxW(ic)) + (fE(ic,ic)*cE(ic) + dE(ic,ic)*dcdxE(ic))
                end do
            else if (j == NJ) then
                dW(1,1) = -(1.0 - por(j))*sigma
                smG(1) = (c_density) - dW(1,1)*dcdxW(1)
                smG(2) = en_row(v)
                do ic = 3, N
                    rj(2,ic) = z_ion(ic)
                end do
                do ic = 3, N
                    dW(ic,ic) = -1.0*diff_term(ic,j)
                    dW(ic,2) = -1.0*mig_term(ic,j)*cW(ic)
                    fW(ic,ic) = -1.0*mig_term(ic,j)*dcdxW(2)
                    smG(ic) = 0.0 - (dW(ic,ic)*dcdxW(ic) + fW(ic,ic)*cW(ic))
                end do
                Ab(:,:,j) = (1.d0 - alphaW)*fW - betaW*dW
                Bb(:,:,j) = rj + betaW*dW + alphaW*fW
                Gb(:,j) = smG
                call nan_check(dE, dW, fE, fW, rj, smG)
                cycle
            else if (j < S) then
                p = por(j)
                dE(1,1) = -(1.0 - p)*sigma_sep
                dW(1,1) = -(1.0 - p)*sigma_sep
                smG(1) = 0.0 - (fW(1,1)*cW(1) + dW(1,1)*dcdxW(1)) + (fE(1,1)*cE(1) + dE(1,1)*dcdxE(1))
                smG(2) = en_row(v)
                do ic = 3, N
                    rj(2,ic) = z_ion(ic)
                end do
                do ic = 3, N
                    dW(ic,ic) = -1*diff_term(ic,j)
                    dE(ic,ic) = -1*diff_term(ic,j)
                    fW(ic,ic) = -1*mig_term(ic,j)*dcdxW(2)
                    fE(ic,ic) = -1*mig_term(ic,j)*dcdxE(2)
                    dW(ic,2) = -1*mig_term(ic,j)*cW(ic)
                    dE(ic,2) = -1*mig_term(ic,j)*cE(ic)
                    smG(ic) = -(fW(ic,ic)*cW(ic) + dW(ic,ic)*dcdxW(ic)) + (fE(ic,ic)*cE(ic) + dE(ic,ic)*dcdxE(ic))
                end do
                do iq = 3, N
                    do iv = 1, N
                        if (iq == iv) then
                            rj(iq,iv) = -(p/delT)*delx(j)
                        else
                            rj(iq,iv) = 0.0
                        end if
                    end do
                end do
            else                                                ! cathode interior
                call react_tot(v, j, rxn)
                call drx_dc(v, j, drx)
                if (ph) then                                    ! POROSITY_UPDATE_ON: face values interpolated
                    pW = alphaW*por(j) + (1.d0 - alphaW)*por(j-1)
                    pE = alphaE*por(j+1) + (1.d0 - alphaE)*por(j)
                    mW = alphaW*mig_term(:,j) + (1.d0 - alphaW)*mig_term(:,j-1)
                    mE = alphaE*mig_term(:,j+1) + (1.d0 - alphaE)*mig_term(:,j)
                    dfW = alphaW*diff_term(:,j) + (1.d0 - alphaW)*diff_term(:,j-1)
                    dfE = alphaE*diff_term(:,j+1) + (1.d0 - alphaE)*diff_term(:,j)
                else
                    pW = por(j); pE = por(j)
                    mW = mig_term(:,j); mE = mig_term(:,j)
                    dfW = diff_term(:,j); dfE = diff_term(:,j)
                end if
                dE(1,1) = -(1.0 - pE)*sigma
                dW(1,1) = -(1.0 - pW)*sigma
                do ic = 1, N
                    rj(1,ic) = -drx(ic,1)*delx(j)
                end do
                smG(1) = rxn(4)*delx(j) - (fW(1,1)*cW(1) + dW(1,1)*dcdxW(1)) + (fE(1,1)*cE(1) + dE(1,1)*dcdxE(1))
                smG(2) = en_row(v)
                do ic = 3, N
                    rj(2,ic) = z_ion(ic)
                end do
                do ic = 3, N
                    dW(ic,ic) = -1*dfW(ic)
                    dE(ic,ic) = -1*dfE(ic)
                    fW(ic,ic) = -1*mW(ic)*dcdxW(2)
                    fE(ic,ic) = -1*mE(ic)*dcdxE(2)
                    dW(ic,2) = -1*mW(ic)*cW(ic)
                    dE(ic,2) = -1*mE(ic)*cE(ic)
                    smG(ic) = -(rxn(ic+2))*delx(j) - (fW(ic,ic)*cW(ic) + dW(ic,ic)*dcdxW(ic)) &
                              + (fE(ic,ic)*cE(ic) + dE(ic,ic)*dcdxE(ic))
                end do
                p = por(j)
                do iq = 3, N
                    do iv = 1, N
                        if (iq == iv) then
                            rj(iq,iv) = drx(iv,iq-1)*delx(j) - (p/delT)*delx(j)
                        else
                            rj(iq,iv) = drx(iv,iq-1)*delx(j)
                        end if
                    end do
                end do
            end if
            Ab(:,:,j) = (1.d0 - alphaW)*fW - betaW*dW
            Bb(:,:,j) = rj + betaW*dW + alphaW*fW - (1.d0 - alphaE)*fE + betaE*dE
            Db(:,:,j) = -alphaE*fE - betaE*dE
            Gb(:,j) = smG
        end do
    end subroutine assemble

    real(dp) function en_row(v)
        real(dp), intent(in) :: v(:)
        if (ph) then
            en_row = -z_ion(3)*v(3) - z_ion(4)*v(4) - z_ion(5)*v(5) - z_ion(6)*v(6)
        else
            en_row = -z_ion(3)*v(3) - z_ion(4)*v(4) - z_ion(5)*v(5)
        end if
    end function en_row

    subroutine nan_check(dE, dW, fE, fW, rj, smG)
        real(dp), intent(in) :: dE(:,:), dW(:,:), fE(:,:), fW(:,:), rj(:,:), smG(:)
        last_nan = any(ieee_is_nan(dE)) .or. any(ieee_is_nan(dW)) .or. any(ieee_is_nan(fE)) .or. &
                   any(ieee_is_nan(fW)) .or. any(ieee_is_nan(rj)) .or. any(ieee_is_nan(smG))
    end subroutine nan_check

    subroutine solve()
        integer :: status
        call band_solve(N, NJ, Ab, Bb, Db, Gb, delC, status, pivot=PIVOT_LEGACY, singular=SINGULAR_EXACT)
        if (status /= BAND_OK) delC = ieee_value(1.0_dp, ieee_quiet_nan)    ! the original propagates NaN
    end subroutine solve

    ! ================================================================== after the solve
    subroutine update_band_variables()
        integer :: j, k
        do j = 1, NJ
            do k = 1, N
                c(k,j) = c(k,j) + delC(k,j)
            end do
            if (ph) then
                do k = 3, N
                    c(k,j) = max(c(k,j), 1.0d-20)
                end do
            end if
        end do
    end subroutine update_band_variables

    subroutine update_other_variables(stopped)
        logical, intent(out) :: stopped
        real(dp) :: r(13), v(6), cutoff_theta, zn_t, mn_t, mn1, mn2, zn1, por_past
        integer :: j, q
        stopped = .false.
        cutoff_theta = 0.985
        do j = S+1, NJ-1
            do q = 1, N
                v(q) = c(q,j) - (delC(q,j)/2)
            end do
            call react_tot(v(1:N), j, r)
            KMn(j) = KMn(j) + (molar_mass_KMn8O16*r(9)*delT)
            ZHS(j) = ZHS(j) + (molar_mass_ZHS*r(10)*delT)
            ZMCm(j) = ZMCm(j) + (molar_mass_ZMC_max*r(11)*delT)
            if (ph) MnO2(j) = MnO2(j) + (molar_mass_MnO2*r(13)*delT)
            znm(j) = znm(j) + (r(12)*delT)
            ratio(j) = znm(j)/mnm(j)
            if (ratio(j) >= (cutoff_theta*Zmax)) then
                zn_t = znm(j); mn_t = mnm(j)
                mn1 = (zn_t - Zmax*mn_t)/(cutoff_theta*Zmax - Zmax)
                mn2 = mn_t - mn1
                zn1 = cutoff_theta*Zmax*mn1
                znm(j) = zn1; mnm(j) = mn1
                ratio(j) = znm(j)/mnm(j)
                ZMCm(j) = ZMCm(j) + (mn2*molar_mass_ZMC_max)
            end if
            MW(j) = 65.38*ratio(j) + 54.93 + (15.999*2)
            ZMCx(j) = MW(j)*mnm(j)
            if (KMn(j) < 0.0) then
                stopped = .true.                      ! 'KMn8O16 is negative': the original stops
                return
            end if
            if (ZHS(j) < 0.0) ZHS(j) = 1.0d-50
            if (ZMCx(j) < 0.0) ZMCx(j) = 1.0d-50
            if (ZMCm(j) < 0.0) ZMCm(j) = 1.0d-50
            if (ph .and. MnO2(j) < 0.0) MnO2(j) = 1.0d-50
            if (znm(j) < 0.0) znm(j) = 1.7d-50
            if (mnm(j) < 0.0) mnm(j) = 1.0d-50
            a_K(j) = 3.0*KMn(j)/(density_KMn8O16*(xmax_c))
            a_ZMCx(j) = 3.0*ZMCx(j)/(density_ZMC*(xmax_c))
            a_ZMCm(j) = 3.0*ZMCm(j)/(density_ZMC*(xmax_c))
            a_ZHS(j) = 3.0*ZHS(j)/(density_ZHS*(xmax_c))
            if (ph) a_MnO2(j) = 3.0*MnO2(j)/(density_MnO2*(xmax_c))
            por_past = por(j)
            if (ph) then
                por(j) = 1.0 - volfrac_inert - (KMn(j)/density_KMn8O16) - (ZMCx(j)/density_ZMC) - (ZHS(j)/density_ZHS) &
                         - (ZMCm(j)/density_ZMC) - (MnO2(j)/density_MnO2)
                tort(j) = 2.0*por(j)**(-0.5)
                call transport_terms(j)
            end if
            do q = 3, N
                c(q,j) = c(q,j)*por_past/por(j)
            end do
            pH_phreeqc(j) = eval_pH(real(c(3,j), 4), real(c(4,j), 4), j)
            pH_ZHS(j) = zhs_pH(real(c(3,j), 8), real(c(5,j), 8))
            if (ph) then
                pH_standard(j) = -1.0*log10(real(1000*c(6,j), 4))
            else
                pH_standard(j) = -1.0*log10(real(1000*c(1,j+1), 4))        ! cprev(6,j) with N = 5 (M-14)
            end if
        end do
        call copy_boundary(S, S+1)
        pH_phreeqc(NJ) = pH_phreeqc(S+1); pH_ZHS(NJ) = pH_ZHS(S+1); pH_standard(NJ) = pH_standard(S+1)   ! M-16
        call copy_boundary(NJ, NJ-1)
        pH_phreeqc(NJ) = pH_phreeqc(NJ-1); pH_ZHS(NJ) = pH_ZHS(NJ-1); pH_standard(NJ) = pH_standard(NJ-1)
        anode_pot = zn_anode_pot(c(3,1))
    end subroutine update_other_variables

    subroutine copy_boundary(to, from)
        integer, intent(in) :: to, from
        KMn(to) = KMn(from); ZMCx(to) = ZMCx(from); ZMCm(to) = ZMCm(from); ZHS(to) = ZHS(from)
        MnO2(to) = MnO2(from); znm(to) = znm(from); mnm(to) = mnm(from); ratio(to) = ratio(from)
        MW(to) = MW(from); por(to) = por(from); tort(to) = tort(from)
        a_K(to) = a_K(from); a_ZMCx(to) = a_ZMCx(from); a_ZMCm(to) = a_ZMCm(from); a_ZHS(to) = a_ZHS(from)
        a_MnO2(to) = a_MnO2(from)
    end subroutine copy_boundary

    ! ================================================================== output
    subroutine write_header()
        if (ph) then
            write(ou, '(A)') ' time Voltage Current mAh/g State  C_Zn(avg) C_Mn(avg) C_SO4(avg) C_H(avg)  pH_phreeq(avg) '// &
                'pH_ZHS(avg) pH_standard(avg)  KMn8O16(avg) ZMC(avg) ZMC_max(avg) ZHS(avg) MnO2(avg) '// &
                'intercalation_Zn_to_Mn_ratio(avg)  R1_OCP(avg) R2_OCP(avg) R3_OCP(avg) R4_OCP(avg) R5_OCP(avg)  '// &
                'R1_Eta(avg) R2_Eta(avg) R3_Eta(avg) R4_Eta(avg) R5_Eta(avg)  R1_Exi(avg) R2_Exi(avg) R3_Exi(avg) '// &
                'R4_Exi(avg) R5_Exi(avg)  R1_I(avg) R2_I(avg) R3_I(Avg) R4_I(avg) R5_I(avg)  Area_KMn8O16 Area_ZMC '// &
                'Area_ZMC_max Area_ZHS Area_MnO2  mAh/g_sim dom_react'
        else
            write(ou, '(A)') ' time Voltage Current mAh/g State C_Zn(avg) C_Mn(avg) C_SO4(avg) pH_phreeq(avg) pH_ZHS(avg) '// &
                'pH_standard(avg)  KMn8O16(avg) ZMC(avg) ZMC_max(avg) ZHS(avg) intercalation_Zn_to_Mn_ratio(avg)  '// &
                'R1_OCP(avg) R2_OCP(avg) R3_OCP(avg)  R1_Eta(avg) R2_Eta(avg) R3_Eta(avg)  R1_Exi(avg) R2_Exi(avg) '// &
                'R3_Exi(avg)  R1_I(avg) R2_I(avg) R3_I(Avg)  Area_KMn8O16 Area_ZMC Area_ZMC_max Area_ZHS  '// &
                'mAh/g_sim dom_react'
        end if
    end subroutine write_header

    subroutine write_row()
        real(dp) :: rall(4,2,NJ), out(13), zr(4), r2(4), r3(4), r5(4), temp(4)
        integer :: j, nn, q, dom
        nn = NJ - S - 1
        rall = 0
        do j = S+1, NJ-1
            if (.not. ph) then
                call reaction_2(c(1,j), c(2,j), c(3,j), c(4,j), c(5,j), j, out)
                rall(:,1,j) = out(1:4)
                call reaction_3(c(1,j), c(2,j), c(3,j), c(4,j), c(5,j), j, out)
                rall(:,2,j) = out(1:4)
            else
                call reaction_5(c(1,j), c(2,j), c(3,j), c(4,j), c(5,j), c(6,j), j, out)
                rall(:,1,j) = out(1:4)
            end if
        end do
        do q = 1, 4
            zr(q) = sum_div0()
            r2(q) = sum(rall(q,1,S+1:NJ-1))/nn
            r3(q) = sum(rall(q,2,S+1:NJ-1))/nn
        end do
        r5 = r2
        if (.not. ph) then
            temp = [abs(zr(4)), abs(r2(4)), abs(r3(4)), abs(zr(4))]
            dom = maxloc(temp, 1)
            write(ou, *) time, c(1,NJ), current, mAhg, state, &
                sum(c(3,S+1:NJ-1)/nn), sum(c(4,S+1:NJ-1)/nn), sum(c(5,S+1:NJ-1)/nn), &
                sum(pH_phreeqc(S+1:NJ-1)/nn), sum(pH_ZHS(S+1:NJ-1)/nn), sum(pH_standard(S+1:NJ-1)/nn), &
                sum(KMn(S+1:NJ-1)/nn), sum(ZMCx(S+1:NJ-1)/nn), sum(ZMCm(S+1:NJ-1)/nn), sum(ZHS(S+1:NJ-1)/nn), &
                sum(ratio(S+1:NJ-1)/nn), &
                zr(1), r2(1), r3(1), zr(2), r2(2), r3(2), zr(3), r2(3), r3(3), zr(4), r2(4), r3(4), &
                sum(a_K(S+1:NJ-1))/nn, sum(a_ZMCx(S+1:NJ-1))/nn, sum(a_ZMCm(S+1:NJ-1))/nn, sum(a_ZHS(S+1:NJ-1))/nn, &
                mAhg*AM_Grams/AM_Grams_sim, dom
        else
            dom = 1                                         ! maxloc over |R1..R4| currents, all zero
            write(ou, *) time, c(1,NJ), current, mAhg, state, &
                sum(c(3,S+1:NJ-1)/nn), sum(c(4,S+1:NJ-1)/nn), sum(c(5,S+1:NJ-1)/nn), sum(c(6,S+1:NJ-1)/nn), &
                sum(pH_phreeqc(S+1:NJ-1)/nn), sum(pH_ZHS(S+1:NJ-1)/nn), sum(pH_standard(S+1:NJ-1)/nn), &
                sum(KMn(S+1:NJ-1)/nn), sum(ZMCx(S+1:NJ-1)/nn), sum(ZMCm(S+1:NJ-1)/nn), sum(ZHS(S+1:NJ-1)/nn), &
                sum(MnO2(S+1:NJ-1)/nn), sum(ratio(S+1:NJ-1)/nn), &
                zr(1), zr(1), zr(1), zr(1), r5(1), zr(2), zr(2), zr(2), zr(2), r5(2), &
                zr(3), zr(3), zr(3), zr(3), r5(3), zr(4), zr(4), zr(4), zr(4), r5(4), &
                sum(a_K(S+1:NJ-1))/nn, sum(a_ZMCx(S+1:NJ-1))/nn, sum(a_ZMCm(S+1:NJ-1))/nn, sum(a_ZHS(S+1:NJ-1))/nn, &
                sum(a_MnO2(S+1:NJ-1))/nn, mAhg*AM_Grams/AM_Grams_sim, dom
        end if
    contains
        real(dp) function sum_div0()
            real(dp) :: zero(NJ)
            zero = 0
            sum_div0 = sum(zero(S+1:NJ-1))/nn
        end function sum_div0
    end subroutine write_row

    ! ================================================================== main loop
    subroutine current_ramp()
        if (ramp_on) then
            if (ramp_count == 0) ramp_initial = current
            delT = RAMP_DELT
            current = (real(ramp_count, kind=8)*(current_target - ramp_initial)/real(RAMP_ITERS, kind=8)) + ramp_initial
            ramp_count = ramp_count + 1
            if (ramp_count == RAMP_ITERS) then
                current = current_target
                ramp_on = .false.
                ramp_count = 0
                delT = DELT_NOMINAL
            end if
        end if
        c_density = current/Area_CS
        c_specific = current/AM_Grams
    end subroutine current_ramp

    subroutine main_loop(reason)
        character(len=*), intent(out) :: reason
        integer :: it
        logical :: stopped
        call write_header()
        it = 0
        do
            it = it + 1
            if ((c(2,NJ) >= 99.0) .and. state == 'C') then
                reason = 'EXIT BECAUSE END OF CHARGE'; return
            else if (ieee_is_nan(delC(1,1))) then
                reason = 'EXIT BECAUSE delC ISNAN'; return
            else if (time >= 99.0*3600.0) then
                reason = 'EXIT BECAUSE END OF SIMULATION TIME'; return
            end if
            if (.not. ph) then
                if ((c(1,NJ) - c(2,NJ)) <= 0.8) then
                    reason = 'EXIT BECAUSE LOWER VOLTAGE CUTOFF'; return
                else if (c(1,NJ) >= 1.8) then
                    reason = 'EXIT BECAUSE Upper VOLTAGE CUTOFF'; return
                end if
            else
                if (c(1,NJ) <= 1.0) then
                    reason = 'EXIT BECAUSE LOWER VOLTAGE CUTOFF'; return
                else if (c(1,NJ) >= 2.0) then
                    reason = 'EXIT BECAUSE Upper VOLTAGE CUTOFF'; return
                end if
            end if
            call current_ramp()
            if (current > 1.0d-10) then
                state = 'D'
            else if (current <= -1.0d-10) then
                state = 'C'
            else
                state = 'R'
            end if
            if (ph) then
                if (mAhg >= 89.0) delT = 0.1
                if (mAhg >= 100.0) delT = 1.0
            end if
            if (it >= RAMP_ITERS) then
                if ((time - WRITE_DENSITY) >= last_write_time) then
                    call write_row()
                    last_write_time = int(time)
                end if
            else if (it == 1) then
                call write_row()
                last_write_time = int(time)
            end if
            call assemble()
            call solve()
            call update_band_variables()
            call update_other_variables(stopped)
            if (stopped) then
                reason = 'KMn8O16 is negative'; return
            end if
            if (any(ieee_is_nan(c))) then
                reason = 'NaN in cprev'; return
            end if
            if (any(ieee_is_nan(delC))) then
                reason = 'NaN in delC'; return
            end if
            if (last_nan) then
                reason = 'NaN in coefficients'; return
            end if
            if ((time - WRITE_DENSITY) >= last_write_time) then
                call write_row()
                last_write_time = int(time)
            end if
            time = time + delT
            if (state == 'D') then
                mAhg = mAhg + 1000.0*c_specific*delT/3600.0
            else if (state == 'C') then
                mAhg = mAhg - 1000.0*c_specific*delT/3600.0
            end if
        end do
    end subroutine main_loop

end module zn_faithful
