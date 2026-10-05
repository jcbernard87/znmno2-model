!> Corrected model (docs/model.md): finite volumes, implicit chemistry, Newton with bandsolver, and the
!> protocol driver. A port of python/znmno2_model (model.py, simulate.py, driver.py, protocol.py), with
!> the same equations, options, numerics and output table; tests/test_ports.py compares the three.
!>
!> Unknowns per cell: 1 phi1, 2 phi2, 3 Zn_T, 4 Mn_T, 5 S_T, 6 H_T, 7 n_MnO2, 8 n_ZMO, 9 s (insertion
!> log-odds), 10 n_ZHS, 11 n_ZnO, 12 n_Zn(OH)2. Cells run from the Zn anode through the probe region,
!> the separator and the cathode to the current collector.
!>
!> SPDX-License-Identifier: BSD-3-Clause
module zn_corrected
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf, ieee_negative_inf
    use, intrinsic :: iso_c_binding, only: c_double
    use zn_params, only: dp, params_t, read_params, M_ZMO, M_host, active_mass
    use zn_equilibria, only: eq_t, eq_init, eq_solve, eq_species, eq_sensitivity, LN10
    use bandsolver_kernel, only: band_solve, BAND_OK
    implicit none
    private
    public :: run_corrected

    interface
        pure real(c_double) function c_log1p(x) bind(C, name='log1p')
            import :: c_double
            real(c_double), value, intent(in) :: x
        end function c_log1p
        pure real(c_double) function c_expm1(x) bind(C, name='expm1')
            import :: c_double
            real(c_double), value, intent(in) :: x
        end function c_expm1
    end interface

    integer, parameter :: N = 12
    integer, parameter :: P1 = 1, P2 = 2, ZN = 3, MN = 4, SO = 5, H = 6, MO = 7, ZM = 8, TH = 9, ZH = 10, &
                          ZO = 11, ZX = 12
    integer, parameter :: SPECIES(4) = [ZN, MN, SO, H]
    integer, parameter :: SOLIDS(5) = [MO, ZM, ZH, ZO, ZX]
    integer, parameter :: PRECIP(3) = [ZH, ZO, ZX]
    integer, parameter :: TOT_COLS(4) = [H, ZN, MN, SO]      ! column of each speciation total
    integer, parameter :: PROBE = 0, SEP = 1, CATH = 2
    real(dp), parameter :: CHARGE(N) = [0.0_dp, 0.0_dp, 2.0_dp, 2.0_dp, -2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: W_SMOOTH = 0.01_dp, S_MAX = 200.0_dp, S_STEP = 20.0_dp, TH_FLOOR = 1e-8_dp
    real(dp), parameter :: RES_TOL = 1e-6_dp, TRACE = 1.0e-9_dp, FB_REF = 1.0e-6_dp
    real(dp), parameter :: TYP(N) = [1.0_dp, 1.0_dp, 1e-3_dp, 1e-3_dp, 1e-3_dp, 1e-6_dp, 1e-3_dp, 1e-3_dp, &
                                     1.0_dp, 1e-3_dp, 1e-3_dp, 1e-3_dp]
    real(dp), parameter :: STEP_FLOOR(N) = [1.0_dp, 1.0_dp, 1e-15_dp, 1e-15_dp, 1e-15_dp, 1e-8_dp, 1e-15_dp, &
                                            1e-15_dp, 1.0_dp, 1e-15_dp, 1e-15_dp, 1e-15_dp]
    ! driver
    real(dp), parameter :: EVENT_DV = 1.0e-4_dp, MIN_SUBSTEP = 1.0e-10_dp, EVENT_MIN_DT = 1.0e-12_dp, &
                           CV_TOL = 1.0e-9_dp, CV_ACCEPT = 1.0e-6_dp, LIMIT = 1.0e-3_dp, THETA_TOL = 1.0e-6_dp
    integer, parameter :: MAX_FAILURES = 200

    ! ------------------------------------------------------------------ model state
    type(params_t) :: p
    type(eq_t) :: eq
    integer :: nc, c0                          ! cells; first cathode cell
    real(dp), allocatable :: xc(:), dxc(:), area_c(:), vol(:)
    integer, allocatable :: region(:)
    logical, allocatable :: incath(:)
    real(dp) :: fRT, Vm(N), eps_fixed, n_host, a_host, Dk(N), mass
    logical :: transported(N), prec_on(N)
    real(dp) :: prec_logk(N), prec_zn(N), prec_s(N), prec_h(N), prec_k(N)
    real(dp) :: ocp_t(51), ocp_c(50,4)
    real(dp) :: ph_zn(153), ph_mn(153), ph_coef(149,149)
    real(dp), allocatable :: D_sp(:)
    real(dp), allocatable :: H0(:), spec_x(:,:), cache_T(:,:), cache_lx(:,:), x0_state(:,:)
    logical :: have_H0 = .false., spec_valid = .false., cache_valid = .false.
    real(dp) :: res_scale(N)
    integer :: i_probe

    ! ------------------------------------------------------------------ protocol
    type :: step_t
        character(len=4) :: kind = ''
        real(dp) :: I = 0, V = 0, t = -1, Vmin = 0, Vmax = 0, Imin = -1   ! t, Imin < 0: not given
    end type step_t

    ! ------------------------------------------------------------------ output
    integer :: ou
    integer :: n_rows = 0

contains

    ! ================================================================== set-up
    subroutine setup(data_dir)
        character(len=*), intent(in) :: data_dir
        integer :: j, k, nparts, reg(3), nn(3), q, iq
        real(dp) :: lens(3), start, hh, i_ref, v_cath
        nparts = 0
        if (p%L_probe > 0) then
            nparts = nparts + 1; reg(nparts) = PROBE; lens(nparts) = p%L_probe; nn(nparts) = p%n_probe
        end if
        nparts = nparts + 1; reg(nparts) = SEP; lens(nparts) = p%L_sep; nn(nparts) = p%n_sep
        nparts = nparts + 1; reg(nparts) = CATH; lens(nparts) = p%L_cath; nn(nparts) = p%n_cath
        nc = sum(nn(1:nparts))
        allocate(xc(nc), dxc(nc), area_c(nc), vol(nc), region(nc), incath(nc))
        j = 0
        start = 0.0_dp
        do k = 1, nparts
            hh = lens(k)/nn(k)
            do q = 1, nn(k)
                j = j + 1
                xc(j) = start + hh*(real(q - 1, dp) + 0.5_dp)
                dxc(j) = hh
                region(j) = reg(k)
            end do
            start = start + lens(k)
        end do
        area_c = p%A_cell
        vol = area_c*dxc
        incath = region == CATH
        c0 = findloc(incath, .true., 1)
        i_probe = 1
        if (p%L_probe > 0) i_probe = p%n_probe/2 + 1
        call eq_init(eq, trim(p%equilibria_db), trim(p%logk_file))
        fRT = p%F/(p%R*p%T)
        Vm = 0
        Vm(MO) = p%M_MnO2/p%rho_MnO2; Vm(ZM) = M_ZMO(p)/p%rho_ZMO; Vm(ZH) = p%M_ZHS/p%rho_ZHS
        Vm(ZO) = p%M_ZnO/p%rho_ZnO; Vm(ZX) = p%M_ZnOH2/p%rho_ZnOH2
        eps_fixed = 1.0_dp - p%eps_cath - p%vf_MnO2 - p%vf_ZMO - p%vf_ZHS
        n_host = p%vf_host*p%rho_host/M_host(p)
        a_host = 3.0_dp*p%vf_host/p%r_host
        Dk = 0
        Dk(ZN) = p%D_Zn; Dk(MN) = p%D_Mn; Dk(SO) = p%D_SO4; Dk(H) = p%D_H
        prec_on = .false.
        prec_on(ZH) = p%zhs == 'kinetic'; prec_on(ZO) = p%ZnO_on; prec_on(ZX) = p%ZnOH2_on
        prec_logk(ZH) = p%logK_ZHS; prec_zn(ZH) = 4.0_dp; prec_s(ZH) = 1.0_dp; prec_h(ZH) = 6.0_dp; prec_k(ZH) = p%k_ZHS
        prec_logk(ZO) = p%logK_ZnO; prec_zn(ZO) = 1.0_dp; prec_s(ZO) = 0.0_dp; prec_h(ZO) = 2.0_dp; prec_k(ZO) = p%k_ZnO
        prec_logk(ZX) = p%logK_ZnOH2; prec_zn(ZX) = 1.0_dp; prec_s(ZX) = 0.0_dp; prec_h(ZX) = 2.0_dp
        prec_k(ZX) = p%k_ZnOH2
        transported = .false.
        transported(ZN) = .true.; transported(MN) = .true.; transported(SO) = .true.
        transported(H) = p%species /= 'no_H'
        if (p%r3_ocp /= 'nernst') call read_ocp(data_dir)
        if (p%ph_mode == 'spline') call read_ph_spline(data_dir)
        if (p%transport == 'quasi') then
            allocate(D_sp(eq%ns))
            do iq = 1, eq%ns
                select case (trim(eq%names(iq)))
                case ('H+'); D_sp(iq) = p%D_H
                case ('Zn+2'); D_sp(iq) = p%D_Zn
                case ('Mn+2'); D_sp(iq) = p%D_Mn
                case ('SO4-2'); D_sp(iq) = p%D_SO4
                case ('OH-'); D_sp(iq) = p%D_OH
                case ('HSO4-'); D_sp(iq) = p%D_HSO4
                case default; D_sp(iq) = p%D_complex
                end select
            end do
        end if
        allocate(H0(nc), spec_x(4,nc), cache_T(4,nc), cache_lx(4,nc), x0_state(N,nc))
        mass = active_mass(p)
        i_ref = 1.0e-3_dp*mass
        v_cath = 0
        do j = 1, nc
            if (incath(j)) v_cath = v_cath + area_c(j)*dxc(j)
        end do
        res_scale = i_ref/(p%F*v_cath)
        res_scale(P1) = i_ref
        res_scale(P2) = 1.0_dp
        do k = 1, 4
            res_scale(SPECIES(k)) = i_ref/p%F
        end do
        if (p%zhs == 'equilibrium') res_scale(ZH) = 1.0_dp
    end subroutine setup

    subroutine skip_comments(u)
        integer, intent(in) :: u
        character(len=512) :: line
        do
            read(u, '(A)') line
            if (adjustl(line(1:1)) /= '#') exit
        end do
        backspace(u)
    end subroutine skip_comments

    subroutine read_ocp(dir)
        character(len=*), intent(in) :: dir
        real :: t4(51), c4(50,4)
        integer :: u, i, n1, ios
        open(newunit=u, file=trim(dir)//'/r3_ocp_spline.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'r3_ocp_spline.txt not found: set data_dir in &run'
        call skip_comments(u)
        read(u, *) n1
        read(u, *) t4
        read(u, *) (c4(i,:), i = 1, 50)
        close(u)
        ocp_t = real(t4, dp)
        ocp_c = real(c4, dp)
    end subroutine read_ocp

    subroutine read_ph_spline(dir)
        character(len=*), intent(in) :: dir
        real :: z4(153), m4(153)
        real, allocatable :: c4(:,:)
        integer :: u, i, n1, n2, n3, n4, ios
        allocate(c4(149,149))
        open(newunit=u, file=trim(dir)//'/ph_spline_phreeqc.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'ph_spline_phreeqc.txt not found: set data_dir in &run'
        call skip_comments(u)
        read(u, *) n1, n2, n3, n4
        read(u, *) z4, m4
        read(u, *) (c4(i,:), i = 1, 149)
        close(u)
        ph_zn = real(z4, dp); ph_mn = real(m4, dp); ph_coef = real(c4, dp)
    end subroutine read_ph_spline

    ! ================================================================== small functions
    elemental real(dp) function bernoulli(a)
        real(dp), intent(in) :: a
        if (abs(a) < 1e-6_dp) then
            bernoulli = 1.0_dp - a/2.0_dp + a*a/12.0_dp
        else
            bernoulli = a/c_expm1(a)
        end if
    end function bernoulli

    elemental real(dp) function dbernoulli(a)
        real(dp), intent(in) :: a
        real(dp) :: em
        if (abs(a) < 1e-4_dp) then
            dbernoulli = -0.5_dp + a/6.0_dp
        else
            em = c_expm1(a)
            dbernoulli = (em - a*(em + 1.0_dp))/(em*em)
        end if
    end function dbernoulli

    elemental subroutine theta_pair(s, thv, omv)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: thv, omv
        real(dp) :: e, big, small
        e = exp(-abs(s))
        big = 1.0_dp/(1.0_dp + e)
        small = e/(1.0_dp + e)
        if (s >= 0) then
            thv = big; omv = small
        else
            thv = small; omv = big
        end if
    end subroutine theta_pair

    elemental real(dp) function sigmoid(s)
        real(dp), intent(in) :: s
        real(dp) :: omv
        call theta_pair(s, sigmoid, omv)
    end function sigmoid

    !> A C1 max(z, 0), exactly zero for z <= 0: z^2 / (z + width) for z > 0.
    elemental real(dp) function pos(z, width)
        real(dp), intent(in) :: z, width
        real(dp) :: zp
        zp = max(z, 0.0_dp)
        pos = zp*zp/(zp + width)
    end function pos

    !> numpy.logaddexp(x, y)
    elemental real(dp) function logaddexp(x, y)
        real(dp), intent(in) :: x, y
        real(dp) :: tmp
        if (x == y) then
            logaddexp = x + log(2.0_dp)
        else
            tmp = x - y
            if (tmp > 0) then
                logaddexp = x + c_log1p(exp(-tmp))
            else if (tmp <= 0) then
                logaddexp = y + c_log1p(exp(tmp))
            else
                logaddexp = tmp
            end if
        end if
    end function logaddexp

    real(dp) function onoff(flag)
        logical, intent(in) :: flag
        onoff = 0.0_dp
        if (flag) onoff = 1.0_dp
    end function onoff

    ! ================================================================== pH spline (ph_mode = spline)
    subroutine bspline_basis(t, i, x, bas)
        real(dp), intent(in) :: t(:), x
        integer, intent(in) :: i                ! 0-based knot span
        real(dp), intent(out) :: bas(0:3)
        real(dp) :: left(0:3), right(0:3), saved, temp
        integer :: j, r
        bas = 0; left = 0; right = 0
        bas(0) = 1.0_dp
        do j = 1, 3
            left(j) = x - t(i + 1 - j + 1)
            right(j) = t(i + j + 1) - x
            saved = 0.0_dp
            do r = 0, j - 1
                temp = bas(r)/(right(r + 1) + left(j - r))
                bas(r) = saved + right(r + 1)*temp
                saved = left(j - r)*temp
            end do
            bas(j) = saved
        end do
    end subroutine bspline_basis

    integer function span(t, x)
        real(dp), intent(in) :: t(:), x
        integer :: i, cnt
        cnt = 0
        do i = 1, size(t)
            if (t(i) <= x) cnt = cnt + 1
        end do
        span = min(max(cnt - 1, 3), size(t) - 3 - 2)
    end function span

    real(dp) function ph_spline(zn_molar, mn_molar)
        real(dp), intent(in) :: zn_molar, mn_molar
        real(dp) :: a, b, bx(0:3), by(0:3), row(0:3)
        integer :: i, j, r, q
        ! the original evaluates its first knot axis at log10(Mn) and the second at log10(Zn)
        a = min(max(log10(mn_molar), ph_zn(4)), ph_zn(153 - 3) - 1e-12_dp)
        b = min(max(log10(zn_molar), ph_mn(4)), ph_mn(153 - 3) - 1e-12_dp)
        i = span(ph_zn, a); j = span(ph_mn, b)
        call bspline_basis(ph_zn, i, a, bx)
        call bspline_basis(ph_mn, j, b, by)
        do q = 0, 3
            row(q) = 0
            do r = 0, 3
                row(q) = row(q) + bx(r)*ph_coef(i - 3 + r + 1, j - 3 + q + 1)
            end do
        end do
        ph_spline = 0
        do q = 0, 3
            ph_spline = ph_spline + row(q)*by(q)
        end do
    end function ph_spline

    ! ================================================================== speciation and constitutive laws
    !> log10 free concentrations [mol/L] (4, nc) of H+, Zn2+, Mn2+, SO4 2- from the totals.
    subroutine lx_of(x, warm, lx, fail)
        real(dp), intent(in) :: x(N,nc)
        logical, intent(in) :: warm
        real(dp), intent(out) :: lx(4,nc)
        logical, intent(out) :: fail
        real(dp) :: T(4,nc)
        integer :: j, q
        logical :: ok
        fail = .false.
        do j = 1, nc
            do q = 1, 4
                T(q,j) = x(TOT_COLS(q),j)*1e3_dp
            end do
        end do
        if (cache_valid) then
            if (all(T == cache_T)) then
                lx = cache_lx
                return
            end if
        end if
        do j = 1, nc
            if (spec_valid) lx(:,j) = spec_x(:,j)
            call eq_solve(eq, T(:,j), lx(:,j), spec_valid, ok)
            if (.not. ok) call eq_solve(eq, T(:,j), lx(:,j), .false., ok)
            if (.not. ok) then
                fail = .true.
                return
            end if
        end do
        if (warm) then
            spec_x = lx; spec_valid = .true.
            cache_T = T; cache_lx = lx; cache_valid = .true.
        end if
    end subroutine lx_of

    subroutine free_of(x, warm, fr, fail)
        real(dp), intent(in) :: x(N,nc)
        logical, intent(in) :: warm
        real(dp), intent(out) :: fr(4,nc)
        logical, intent(out) :: fail
        real(dp) :: lx(4,nc)
        call lx_of(x, warm, lx, fail)
        fr = 10.0_dp**lx
    end subroutine free_of

    !> Concentrations [mol/L] of H+, Zn, Mn, SO4 in the Nernst and rate terms (options basis, ph_mode).
    subroutine basis_of(x, fr, cb)
        real(dp), intent(in) :: x(N,nc), fr(4,nc)
        real(dp), intent(out) :: cb(4,nc)
        integer :: j
        cb = fr
        do j = 1, nc
            if (p%basis == 'totals') then
                cb(2,j) = x(ZN,j)*1e3_dp; cb(3,j) = x(MN,j)*1e3_dp; cb(4,j) = x(SO,j)*1e3_dp
            end if
            if (p%ph_mode == 'zhs_equilibrium') then
                cb(1,j) = exp((4.0_dp*log(cb(2,j)) + log(cb(4,j)) - p%logK_ZHS*LN10)/6.0_dp)
            else if (p%ph_mode == 'fixed') then
                cb(1,j) = 10.0_dp**(-p%pH_fixed)
            else if (p%ph_mode == 'spline') then
                cb(1,j) = 10.0_dp**(-ph_spline(x(ZN,j)*1e3_dp, x(MN,j)*1e3_dp))
            end if
        end do
    end subroutine basis_of

    subroutine porosity(x, eps)
        real(dp), intent(in) :: x(N,nc)
        real(dp), intent(out) :: eps(nc)
        integer :: j, q
        real(dp) :: sv
        do j = 1, nc
            if (region(j) == PROBE) then
                eps(j) = p%eps_probe
            else
                eps(j) = p%eps_sep
            end if
            if (incath(j)) then
                sv = 0
                do q = 1, 5
                    sv = sv + Vm(SOLIDS(q))*x(SOLIDS(q),j)
                end do
                eps(j) = 1.0_dp - eps_fixed - sv
            end if
        end do
    end subroutine porosity

    subroutine tortuosity(eps, tau, fail)
        real(dp), intent(in) :: eps(nc)
        real(dp), intent(out) :: tau(nc)
        logical, intent(out) :: fail
        integer :: j
        fail = any(eps <= 0.0_dp)
        if (fail) return
        do j = 1, nc
            tau(j) = p%tau_probe
            if (region(j) == SEP) tau(j) = p%tau_factor_sep*eps(j)**p%bruggeman_sep
            if (region(j) == CATH) tau(j) = p%tau_factor_cath*eps(j)**p%bruggeman_cath
        end do
    end subroutine tortuosity

    !> R3 equilibrium potential from the insertion log-odds s and Zn2+ [mol/L].
    real(dp) function u3(s, c_zn)
        real(dp), intent(in) :: s, c_zn
        real(dp) :: thv, omv, endt, epsw, w, lt, lo, d, vs
        integer :: k, i, cnt
        if (p%r3_ocp == 'nernst') then
            u3 = p%U3_nernst + 0.5_dp/fRT*(log(c_zn) - s)
            return
        end if
        call theta_pair(s, thv, omv)
        endt = 0.0_dp
        if (p%r3_ocp == 'spline_nernst') then
            epsw = p%r3_end_width; w = 0.5_dp
            lt = min(s, 0.0_dp) - c_log1p(exp(-abs(s)))
            lo = min(-s, 0.0_dp) - c_log1p(exp(-abs(s)))
            endt = -(0.5_dp/fRT)*(soft(lt - log(epsw)) - soft(lo - log(epsw)))
        end if
        cnt = 0
        do i = 1, 51
            if (ocp_t(i) < thv) cnt = cnt + 1
        end do
        k = min(max(cnt - 1, 0), 49) + 1
        d = thv - ocp_t(k)
        vs = ((ocp_c(k,1)*d + ocp_c(k,2))*d + ocp_c(k,3))*d + ocp_c(k,4)
        u3 = (p%V_at_zmin - p%V_at_zmax)*vs + p%V_at_zmax + 0.5_dp/fRT*log(c_zn/p%c_ref3) + endt
    contains
        real(dp) function soft(a)                ! smooth min(a, 0)
            real(dp), intent(in) :: a
            soft = -w*logaddexp(0.0_dp, -a/w)
        end function soft
    end function u3

    real(dp) function area(x, k, j)
        real(dp), intent(in) :: x(N,nc)
        integer, intent(in) :: k, j
        real(dp) :: r
        select case (k)
        case (MO); r = p%r_MnO2
        case (ZM); r = p%r_ZMO
        case default; r = p%r_ZHS
        end select
        area = 3.0_dp*Vm(k)*max(x(k,j), 0.0_dp)/r
    end function area

    !> Kinetic precipitation rate of solid k [mol/cm3/s] in cathode cell j (positive = forms).
    real(dp) function precipitation(x, fr, k, j)
        real(dp), intent(in) :: x(N,nc), fr(4,nc)
        integer, intent(in) :: k, j
        real(dp) :: lq, w, w_nuc, sterm
        precipitation = 0.0_dp
        if (.not. prec_on(k)) return
        sterm = 0.0_dp
        if (prec_s(k) /= 0) sterm = prec_s(k)*log(fr(4,j))
        lq = prec_zn(k)*log(fr(2,j)) + sterm - prec_logk(k)*LN10
        w = exp(lq/prec_h(k))/fr(1,j)
        w_nuc = p%zhs_nucleation**(prec_zn(k)/prec_h(k))
        precipitation = prec_k(k)*(area(x, k, j)*(w - 1.0_dp) + p%a_seed_ZHS*pos(w - w_nuc, W_SMOOTH))
    end function precipitation

    real(dp) function zhs_w(fr, j)
        real(dp), intent(in) :: fr(4,nc)
        integer, intent(in) :: j
        zhs_w = exp((4.0_dp*log(fr(2,j)) + log(fr(4,j)) - p%logK_ZHS*LN10)/6.0_dp)/fr(1,j)
    end function zhs_w

    !> Reaction currents (A/cm3, positive anodic) and the kinetic ZHS rate (mol/cm3/s) in cathode cell j.
    subroutine reactions(x, fr, cb, j, i1, i2, i3, r)
        real(dp), intent(in) :: x(N,nc), fr(4,nc), cb(4,nc)
        integer, intent(in) :: j
        real(dp), intent(out) :: i1, i2, i3, r
        real(dp) :: cH, cZn, cMn, ph1, ph2, a_MO, a_ZM, a_ZH, u1, e1, bv1, n2, u2, e2, bv2, thv, omv, e3, i0
        cH = cb(1,j); cZn = cb(2,j); cMn = cb(3,j)
        ph1 = x(P1,j); ph2 = x(P2,j)
        a_MO = area(x, MO, j); a_ZM = area(x, ZM, j); a_ZH = area(x, ZH, j)
        ! R1: irreversible dissolution, only below its equilibrium potential
        u1 = p%U1 - 0.5_dp/fRT*(log(cMn) - 4.0_dp*log(cH))
        e1 = ph1 - ph2 - u1
        bv1 = exp(p%alpha1*2.0_dp*fRT*e1) - exp(-(1.0_dp - p%alpha1)*2.0_dp*fRT*e1)
        i1 = (-a_MO)*p%F*p%k1*pos(-bv1, 1.0e-2_dp)*onoff(p%R1_on)
        ! R2: oxidation on ZMO + ZHS + seed, reduction on ZMO (lumped: deposition on ZHS)
        n2 = 2.0_dp - 2.0_dp*p%z_ZMO
        u2 = p%U2 - (1.0_dp/(n2*fRT))*(p%z_ZMO*log(cZn) + log(cMn) - 4.0_dp*log(cH))
        e2 = ph1 - ph2 - u2
        bv2 = exp(p%alpha2*n2*fRT*e2) - exp(-(1.0_dp - p%alpha2)*n2*fRT*e2)
        if (p%zhs == 'lumped') then
            i2 = p%F*p%k2*(a_ZH*pos(bv2, 1.0e-2_dp) - a_ZM*pos(-bv2, 1.0e-2_dp))*onoff(p%R2_on)
        else
            i2 = p%F*p%k2*(a_ZM*bv2 + (a_ZH + p%a_seed_R2)*pos(bv2, 1.0e-2_dp))*onoff(p%R2_on)
        end if
        ! R3: insertion
        call theta_pair(x(TH,j), thv, omv)
        e3 = ph1 - ph2 - u3(x(TH,j), cZn)
        i0 = p%F*p%k3*sqrt(cZn*thv*omv)
        i3 = a_host*i0*(exp(p%alpha3*2.0_dp*fRT*e3) - exp(-(1.0_dp - p%alpha3)*2.0_dp*fRT*e3))*onoff(p%R3_on)
        r = precipitation(x, fr, ZH, j)
    end subroutine reactions

    !> Zn anode current density [A/cm2] (positive = Zn dissolves) on the first cell's face.
    real(dp) function anode_current(x, fr)
        real(dp), intent(in) :: x(N,nc), fr(4,nc)
        real(dp) :: czn, eta
        if (p%basis == 'free') then
            czn = fr(2,1)
        else
            czn = x(ZN,1)*1e3_dp
        end if
        eta = 0.0_dp - x(P2,1) - 0.5_dp/fRT*log(czn)
        anode_current = p%F*p%k_an*sqrt(czn)*(exp(p%alpha_an*2.0_dp*fRT*eta) &
                                              - exp(-(1.0_dp - p%alpha_an)*2.0_dp*fRT*eta))
    end function anode_current

    real(dp) function d_theta(s, s_old)
        real(dp), intent(in) :: s, s_old
        real(dp) :: t, o, t0, o0
        call theta_pair(s, t, o)
        call theta_pair(s_old, t0, o0)
        if (s > 0 .and. s_old > 0) then
            d_theta = o0 - o
        else
            d_theta = t - t0
        end if
    end function d_theta

    ! ================================================================== residual
    !> Linear terms and their Jacobian: storage, electroneutrality, trivial rows, solid storage, the
    !> collector current, and (zhs = equilibrium) the ZHS sinks.
    subroutine linear(x, old, dt, I, R, J)
        real(dp), intent(in) :: x(N,nc), old(N,nc), dt, I
        real(dp), intent(out) :: R(N,nc), J(N,N,nc)
        real(dp) :: eps(nc), eps_old(nc), cap, thv, omv, nu
        integer :: jj, q, k, s_, kk
        R = 0; J = 0
        call porosity(x, eps); call porosity(old, eps_old)
        do jj = 1, nc
            do q = 1, 4
                k = SPECIES(q)
                if (.not. transported(k)) cycle
                R(k,jj) = vol(jj)*(eps(jj)*x(k,jj) - eps_old(jj)*old(k,jj))/dt
                J(k,k,jj) = vol(jj)*eps(jj)/dt
                if (incath(jj)) then
                    do s_ = 1, 5
                        J(k,SOLIDS(s_),jj) = vol(jj)*x(k,jj)*(-Vm(SOLIDS(s_)))/dt
                    end do
                end if
            end do
            if (p%species == 'no_H') then
                R(H,jj) = x(H,jj) - H0(jj)
                J(H,H,jj) = 1.0_dp
            end if
            R(P2,jj) = (2.0_dp*x(ZN,jj) + 2.0_dp*x(MN,jj) + x(H,jj) - 2.0_dp*x(SO,jj))*1e3_dp
            J(P2,ZN,jj) = 2e3_dp; J(P2,MN,jj) = 2e3_dp; J(P2,H,jj) = 1e3_dp; J(P2,SO,jj) = -2e3_dp
            if (.not. incath(jj)) then
                R(P1,jj) = x(P1,jj); J(P1,P1,jj) = 1.0_dp
                do s_ = 1, 5
                    k = SOLIDS(s_)
                    R(k,jj) = x(k,jj); J(k,k,jj) = 1.0_dp
                end do
                R(TH,jj) = x(TH,jj) - old(TH,jj)
                J(TH,TH,jj) = 1.0_dp
            else
                do s_ = 1, 5
                    k = SOLIDS(s_)
                    if (k == ZH .and. p%zhs == 'equilibrium') cycle
                    R(k,jj) = (x(k,jj) - old(k,jj))/dt
                    J(k,k,jj) = 1.0_dp/dt
                end do
                if (p%zhs == 'equilibrium') then         ! ZHS formed this step is taken from the solution
                    do q = 1, 3
                        if (q == 1) then
                            kk = ZN; nu = 4.0_dp
                        else if (q == 2) then
                            kk = SO; nu = 1.0_dp
                        else
                            kk = H; nu = -6.0_dp
                        end if
                        if (.not. transported(kk)) cycle
                        R(kk,jj) = R(kk,jj) + nu*vol(jj)*(x(ZH,jj) - old(ZH,jj))/dt
                        J(kk,ZH,jj) = J(kk,ZH,jj) + nu*vol(jj)/dt
                    end do
                end if
                cap = n_host*(p%zmax - p%zmin)
                call theta_pair(x(TH,jj), thv, omv)
                R(TH,jj) = cap*d_theta(x(TH,jj), old(TH,jj))/dt
                J(TH,TH,jj) = cap*thv*omv/dt
            end if
        end do
        R(P1,nc) = R(P1,nc) + I                      ! the solid current I leaves the last cell
    end subroutine linear

    !> Local nonlinear terms: reactions and precipitation in the cathode, the anode flux into cell 1.
    subroutine sources(x, warm, R, fail)
        real(dp), intent(in) :: x(N,nc)
        logical, intent(in) :: warm
        real(dp), intent(out) :: R(N,nc)
        logical, intent(out) :: fail
        real(dp) :: fr(4,nc), cb(4,nc), i1, i2, i3, r_zhs, n2, xi1, xi2, xi3, rk(N), dzn, dso, dh, vc, a, b
        integer :: j, q, k
        R = 0
        call free_of(x, warm, fr, fail)
        if (fail) return
        call basis_of(x, fr, cb)
        n2 = 2.0_dp - 2.0_dp*p%z_ZMO
        do j = c0, nc
            vc = vol(j)
            call reactions(x, fr, cb, j, i1, i2, i3, r_zhs)
            xi1 = -i1/(2.0_dp*p%F); xi2 = -i2/(n2*p%F); xi3 = -i3/(2.0_dp*p%F)
            if (p%zhs == 'lumped') then                  ! the 4 H+ per Mn come from ZHS forming (2024 paper)
                r_zhs = (2.0_dp/3.0_dp)*(xi1 + xi2)
            else if (p%zhs == 'off' .or. p%zhs == 'equilibrium') then
                r_zhs = 0.0_dp
            end if
            rk = 0
            rk(ZH) = r_zhs
            rk(ZO) = precipitation(x, fr, ZO, j)
            rk(ZX) = precipitation(x, fr, ZX, j)
            dzn = p%z_ZMO*xi2 - xi3
            dso = 0.0_dp*xi2
            dh = -4.0_dp*xi1 - 4.0_dp*xi2
            do q = 1, 3
                k = PRECIP(q)
                dzn = dzn - prec_zn(k)*rk(k)
                dso = dso - prec_s(k)*rk(k)
                dh = dh + prec_h(k)*rk(k)
            end do
            R(ZN,j) = -vc*dzn
            R(MN,j) = -vc*(xi1 + xi2)
            R(SO,j) = -vc*dso
            if (p%species == 'with_H') R(H,j) = -vc*dh
            R(P1,j) = vc*(i1 + i2 + i3)
            R(MO,j) = xi1
            R(ZM,j) = xi2
            R(TH,j) = -xi3
            do q = 1, 3
                k = PRECIP(q)
                R(k,j) = -rk(k)
            end do
            if (p%zhs == 'equilibrium') then             ! n_ZHS >= 0, w <= 1, n_ZHS (1 - w) = 0 (Fischer-Burmeister)
                a = x(ZH,j)/FB_REF
                b = 1.0_dp - zhs_w(fr, j)
                R(ZH,j) = a + b - sqrt(a*a + b*b + 1e-20_dp)
            end if
        end do
        R(ZN,1) = R(ZN,1) - anode_current(x, fr)*area_c(1)/(2.0_dp*p%F)
    end subroutine sources

    !> Add the derivatives of an east-face outflow of row `row` (nc-1 faces) to the blocks.
    subroutine add_face(A, B, D, row, dL, dR)
        real(dp), intent(inout) :: A(N,N,nc), B(N,N,nc), D(N,N,nc)
        integer, intent(in) :: row
        real(dp), intent(in) :: dL(nc-1,N), dR(nc-1,N)
        integer :: j, k
        do j = 1, nc - 1
            do k = 1, N
                B(row,k,j) = B(row,k,j) + dL(j,k)
                D(row,k,j) = D(row,k,j) + dR(j,k)
            end do
        end do
        do j = 2, nc
            do k = 1, N
                A(row,k,j) = A(row,k,j) - dL(j-1,k)
                B(row,k,j) = B(row,k,j) - dR(j-1,k)
            end do
        end do
    end subroutine add_face

    !> Face fluxes between neighbouring cells (outflow east minus inflow west) and, with jac, the blocks
    !> (A, B, D) of their derivatives.
    subroutine transport(x, jac, R, A, B, D, fail)
        real(dp), intent(in) :: x(N,nc)
        logical, intent(in) :: jac
        real(dp), intent(out) :: R(N,nc), A(N,N,nc), B(N,N,nc), D(N,N,nc)
        logical, intent(out) :: fail
        real(dp) :: eps(nc), tau(nc), half(nc), G(nc-1), dhalf(nc), dphi(nc-1), bexp
        real(dp) :: dL(nc-1,N), dR(nc-1,N), flow(nc-1), zf, arg, Bp, Bm, dfa, sig, hs(nc), Gs(nc-1), dhs(nc)
        real(dp), allocatable :: lx(:,:), cs(:,:), fl(:,:), dfas(:,:), dc(:,:,:), Bps(:,:), Bms(:,:)
        real(dp) :: T(4), S(4,4), e_, wL, wR, dsum, fsum
        integer :: j, q, k, s_, ss, qq, kk, m_
        R = 0; A = 0; B = 0; D = 0
        call porosity(x, eps)
        call tortuosity(eps, tau, fail)
        if (fail) return
        do j = 1, nc
            half(j) = 0.5_dp*dxc(j)*tau(j)/(eps(j)*area_c(j))
            bexp = 0.0_dp
            if (incath(j)) bexp = p%bruggeman_cath
            dhalf(j) = 0.0_dp
            if (incath(j)) dhalf(j) = (bexp - 1.0_dp)*half(j)/eps(j)
        end do
        do j = 1, nc - 1
            G(j) = 1.0_dp/(half(j) + half(j+1))
            dphi(j) = x(P2,j+1) - x(P2,j)
        end do
        if (p%transport == 'ions') then
            do q = 1, 4
                k = SPECIES(q)
                if (.not. transported(k)) cycle
                zf = CHARGE(k)*fRT
                dL = 0; dR = 0
                do j = 1, nc - 1
                    arg = zf*dphi(j)
                    Bp = bernoulli(arg); Bm = bernoulli(-arg)
                    flow(j) = Dk(k)*G(j)*(Bp*x(k,j) - Bm*x(k,j+1))
                    if (jac) then
                        dL(j,k) = Dk(k)*G(j)*Bp
                        dR(j,k) = -Dk(k)*G(j)*Bm
                        dfa = Dk(k)*G(j)*(dbernoulli(arg)*x(k,j) + dbernoulli(-arg)*x(k,j+1))*zf
                        dR(j,P2) = dfa; dL(j,P2) = -dfa
                    end if
                end do
                call add_flow(k)
                if (jac) then
                    call solid_terms()
                    call add_face(A, B, D, k, dL, dR)
                end if
            end do
        else                                             ! quasi-particles: each species' flux, summed per total
            allocate(lx(4,nc), cs(eq%ns,nc), fl(eq%ns,nc-1), dfas(eq%ns,nc-1), dc(eq%ns,4,nc), &
                     Bps(eq%ns,nc-1), Bms(eq%ns,nc-1))
            call lx_of(x, .true., lx, fail)
            if (fail) return
            do j = 1, nc
                call eq_species(eq, lx(:,j), cs(:,j))
                cs(:,j) = cs(:,j)*1e-3_dp
            end do
            do j = 1, nc - 1
                do ss = 1, eq%ns
                    arg = (eq%z(ss)*fRT)*dphi(j)
                    Bps(ss,j) = bernoulli(arg); Bms(ss,j) = bernoulli(-arg)
                    fl(ss,j) = D_sp(ss)*G(j)*(Bps(ss,j)*cs(ss,j) - Bms(ss,j)*cs(ss,j+1))
                    if (jac) dfas(ss,j) = D_sp(ss)*G(j)*(dbernoulli(arg)*cs(ss,j) + dbernoulli(-arg)*cs(ss,j+1)) &
                                          *eq%z(ss)*fRT
                end do
            end do
            if (jac) then
                do j = 1, nc
                    do qq = 1, 4
                        T(qq) = x(TOT_COLS(qq),j)*1e3_dp
                    end do
                    call eq_sensitivity(eq, lx(:,j), T, S)
                    do ss = 1, eq%ns
                        do qq = 1, 4
                            e_ = 0
                            do m_ = 1, 4
                                e_ = e_ + eq%nu(ss,m_)*S(m_,qq)
                            end do
                            dc(ss,qq,j) = cs(ss,j)*LN10*e_*1e3_dp
                        end do
                    end do
                end do
            end if
            do q = 1, 4
                k = TOT_COLS(q)
                if (.not. transported(k)) cycle
                dL = 0; dR = 0
                do j = 1, nc - 1
                    fsum = 0
                    do ss = 1, eq%ns
                        fsum = fsum + fl(ss,j)*eq%nu(ss,q)
                    end do
                    flow(j) = fsum
                    if (jac) then
                        do qq = 1, 4
                            kk = TOT_COLS(qq)
                            dsum = 0
                            do ss = 1, eq%ns
                                wL = D_sp(ss)*G(j)*Bps(ss,j)*eq%nu(ss,q)
                                dsum = dsum + wL*dc(ss,qq,j)
                            end do
                            dL(j,kk) = dsum
                            dsum = 0
                            do ss = 1, eq%ns
                                wR = -D_sp(ss)*G(j)*Bms(ss,j)*eq%nu(ss,q)
                                dsum = dsum + wR*dc(ss,qq,j+1)
                            end do
                            dR(j,kk) = dsum
                        end do
                        dsum = 0
                        do ss = 1, eq%ns
                            dsum = dsum + dfas(ss,j)*eq%nu(ss,q)
                        end do
                        dR(j,P2) = dsum; dL(j,P2) = -dsum
                    end if
                end do
                call add_flow(k)
                if (jac) then
                    call solid_terms()
                    call add_face(A, B, D, k, dL, dR)
                end if
            end do
        end if
        ! solid current (cathode faces only)
        do j = 1, nc
            sig = 1.0_dp
            if (incath(j)) sig = p%sigma*(1.0_dp - eps(j))
            hs(j) = 0.5_dp*dxc(j)/(sig*area_c(j))
            dhs(j) = 0.0_dp
            if (incath(j)) dhs(j) = hs(j)/(1.0_dp - eps(j))
        end do
        do j = 1, nc - 1
            Gs(j) = 0.0_dp
            if (incath(j) .and. incath(j+1)) Gs(j) = 1.0_dp/(hs(j) + hs(j+1))
            flow(j) = Gs(j)*(x(P1,j) - x(P1,j+1))
        end do
        call add_flow(P1)
        if (jac) then
            dL = 0; dR = 0
            do j = 1, nc - 1
                dL(j,P1) = Gs(j); dR(j,P1) = -Gs(j)
                do s_ = 1, 5
                    ss = SOLIDS(s_)
                    dL(j,ss) = -flow(j)*Gs(j)*dhs(j)*(-Vm(ss))
                    dR(j,ss) = -flow(j)*Gs(j)*dhs(j+1)*(-Vm(ss))
                end do
            end do
            call add_face(A, B, D, P1, dL, dR)
        end if
    contains
        subroutine add_flow(row)
            integer, intent(in) :: row
            integer :: jf
            do jf = 1, nc - 1
                R(row,jf) = R(row,jf) + flow(jf)
            end do
            do jf = 2, nc
                R(row,jf) = R(row,jf) - flow(jf-1)
            end do
        end subroutine add_flow
        subroutine solid_terms()
            integer :: jf, sf, kf
            do sf = 1, 5
                kf = SOLIDS(sf)
                do jf = 1, nc - 1
                    dL(jf,kf) = dL(jf,kf) + (-flow(jf))*G(jf)*dhalf(jf)*(-Vm(kf))
                    dR(jf,kf) = dR(jf,kf) + (-flow(jf))*G(jf)*dhalf(jf+1)*(-Vm(kf))
                end do
            end do
        end subroutine solid_terms
    end subroutine transport

    ! ================================================================== Jacobian and Newton
    subroutine blocks(x, old, dt, I, R, A, B, D, fail)
        real(dp), intent(in) :: x(N,nc), old(N,nc), dt, I
        real(dp), intent(out) :: R(N,nc), A(N,N,nc), B(N,N,nc), D(N,N,nc)
        logical, intent(out) :: fail
        real(dp) :: RL(N,nc), S0(N,nc), RT(N,nc), BT(N,N,nc), hstep(N,nc), xp(N,nc), Sp(N,nc)
        real(dp) :: spec_base(4,nc), cT(4,nc), clx(4,nc)
        logical :: sv_base, cv_base
        integer :: j, k, ii
        call linear(x, old, dt, I, RL, B)
        call sources(x, .true., S0, fail)
        if (fail) return
        call transport(x, .true., RT, A, BT, D, fail)
        if (fail) return
        B = B + BT
        do j = 1, nc
            do k = 1, N
                hstep(k,j) = 1e-7_dp*max(abs(x(k,j)), STEP_FLOOR(k))
            end do
        end do
        spec_base = spec_x; sv_base = spec_valid
        cT = cache_T; clx = cache_lx; cv_base = cache_valid
        do k = 1, N
            xp = x
            xp(k,:) = xp(k,:) + hstep(k,:)
            spec_x = spec_base; spec_valid = sv_base
            call sources(xp, .false., Sp, fail)
            if (fail) return
            do j = 1, nc
                do ii = 1, N
                    B(ii,k,j) = B(ii,k,j) + (Sp(ii,j) - S0(ii,j))/hstep(k,j)
                end do
            end do
        end do
        spec_x = spec_base; spec_valid = sv_base
        cache_T = cT; cache_lx = clx; cache_valid = cv_base
        R = RL + S0 + RT
    end subroutine blocks

    !> Largest lam <= 1 keeping the Zn, Mn and S totals positive and limiting potential changes to 0.2 V.
    real(dp) function damping(x, dx)
        real(dp), intent(in) :: x(N,nc), dx(N,nc)
        integer :: q, k, j
        real(dp) :: big
        damping = 1.0_dp
        do q = 1, 3
            k = SPECIES(q)
            do j = 1, nc
                if (x(k,j) + dx(k,j) < 0 .and. x(k,j) > 0) damping = min(damping, 0.9_dp*x(k,j)/(-dx(k,j)))
            end do
        end do
        big = max(maxval(abs(dx(P1,:))), maxval(abs(dx(P2,:))))
        if (big > 0.2_dp) damping = min(damping, 0.2_dp/big)
    end function damping

    !> One backward-Euler step of length dt at current I, solved by Newton. fail = .true. if it does not converge.
    subroutine newton_step(old, dt, I, xout, fail)
        real(dp), intent(in) :: old(N,nc), dt, I
        real(dp), intent(out) :: xout(N,nc)
        logical, intent(out) :: fail
        real(dp) :: x(N,nc), R(N,nc), A(N,N,nc), B(N,N,nc), D(N,N,nc), G(N,nc), dx(N,nc), fr(4,nc)
        real(dp) :: mag(N,nc), rtol(N,nc), thv(nc), omv(nc), rs, lam, dth(nc), dstep, thn, omn, s_theta, s_s, upd, sc
        logical :: in_theta(nc), conv, nan_upd
        integer :: it, j, ii, k, status, q
        if (.not. have_H0) then
            H0 = old(H,:); have_H0 = .true.
        end if
        x = old
        cache_valid = .false.; spec_valid = .false.
        call free_of(old, .true., fr, fail)           ! warm start from the step's start
        if (fail) return
        do it = 1, p%newton_max_iter
            call blocks(x, old, dt, I, R, A, B, D, fail)
            if (fail) return
            if (.not. all(ieee_is_finite(R))) then
                fail = .true.; return
            end if
            ! round-off floor of each residual: ~1000 eps times the size of its terms (|J| |x| over the band)
            do j = 1, nc
                do ii = 1, N
                    mag(ii,j) = 0
                    do k = 1, N
                        mag(ii,j) = mag(ii,j) + abs(B(ii,k,j))*abs(x(k,j))
                    end do
                    if (j > 1) then
                        do k = 1, N
                            mag(ii,j) = mag(ii,j) + abs(A(ii,k,j))*abs(x(k,j-1))
                        end do
                    end if
                    if (j < nc) then
                        do k = 1, N
                            mag(ii,j) = mag(ii,j) + abs(D(ii,k,j))*abs(x(k,j+1))
                        end do
                    end if
                    rtol(ii,j) = max(RES_TOL*res_scale(ii), 1e3_dp*epsilon(1.0_dp)*mag(ii,j))
                end do
            end do
            ! the insertion unknown: solved for in theta where theta > 1/2, in s elsewhere
            call theta_pair(x(TH,:), thv, omv)
            in_theta = x(TH,:) > 0.0_dp
            do j = 1, nc
                if (in_theta(j)) B(:,TH,j) = B(:,TH,j)/(thv(j)*omv(j))
            end do
            ! scale the columns by the typical size of each unknown and every row to unit maximum
            do j = 1, nc
                do k = 1, N
                    A(:,k,j) = A(:,k,j)*TYP(k); B(:,k,j) = B(:,k,j)*TYP(k); D(:,k,j) = D(:,k,j)*TYP(k)
                end do
                do ii = 1, N
                    rs = max(max(maxval(abs(A(ii,:,j))), maxval(abs(B(ii,:,j))), maxval(abs(D(ii,:,j)))), 1e-300_dp)
                    A(ii,:,j) = A(ii,:,j)/rs; B(ii,:,j) = B(ii,:,j)/rs; D(ii,:,j) = D(ii,:,j)/rs
                    G(ii,j) = -R(ii,j)/rs
                end do
            end do
            call band_solve(N, nc, A, B, D, G, dx, status)
            if (status /= BAND_OK) then
                fail = .true.; return
            end if
            do j = 1, nc
                dx(:,j) = dx(:,j)*TYP
            end do
            dth = dx(TH,:)
            dx(TH,:) = 0.0_dp
            lam = damping(x, dx)
            xout = x + lam*dx
            ! solids are local: a step that would make one negative takes it to a tenth of its value instead
            do q = 1, 5
                k = SOLIDS(q)
                do j = 1, nc
                    if (xout(k,j) < 0) xout(k,j) = 0.1_dp*x(k,j)
                end do
            end do
            upd = 0; nan_upd = .false.
            do j = 1, nc
                dstep = lam*dth(j)
                thn = thv(j) + dstep; omn = omv(j) - dstep
                if (thn < TH_FLOOR*thv(j)) then
                    thn = TH_FLOOR*thv(j); omn = 1.0_dp - thn
                end if
                if (omn < TH_FLOOR*omv(j)) then
                    omn = TH_FLOOR*omv(j); thn = 1.0_dp - omn
                end if
                s_theta = log(max(thn, 1e-300_dp)) - log(max(omn, 1e-300_dp))
                s_s = x(TH,j) + min(max(dstep, -S_STEP), S_STEP)
                if (in_theta(j)) then
                    xout(TH,j) = min(max(s_theta, -S_MAX), S_MAX)
                else
                    xout(TH,j) = min(max(s_s, -S_MAX), S_MAX)
                end if
                do k = 1, N
                    if (k == TH) then
                        if (in_theta(j)) then
                            sc = abs(dth(j))
                        else
                            sc = abs(dth(j))*thv(j)*omv(j)
                        end if
                    else
                        sc = abs(dx(k,j))/TYP(k)
                    end if
                    if (sc /= sc) nan_upd = .true.
                    upd = max(upd, sc)
                end do
            end do
            x = xout
            if (nan_upd .or. .not. ieee_is_finite(upd)) exit
            conv = lam == 1.0_dp .and. upd < p%newton_tol .and. all(abs(R) <= rtol)
            if (conv) then
                call free_of(x, .true., fr, fail)       ! refresh the warm start at the converged state
                if (fail) return
                xout = x
                return
            end if
        end do
        fail = .true.
    end subroutine newton_step

    ! ================================================================== initial state and outputs
    subroutine initial_state(x, fail)
        real(dp), intent(out) :: x(N,nc)
        logical, intent(out) :: fail
        real(dp) :: c_zn, c_mn, fr(4,nc), cb(4,nc), u_zn, xs(N,nc)
        integer :: j
        x = 0
        c_zn = max(p%c_ZnSO4, TRACE); c_mn = max(p%c_MnSO4, TRACE)
        x(ZN,:) = c_zn*1e-3_dp
        x(MN,:) = c_mn*1e-3_dp
        x(SO,:) = (c_zn + c_mn + p%c_H2SO4)*1e-3_dp
        x(H,:) = 2.0_dp*p%c_H2SO4*1e-3_dp
        H0 = x(H,:); have_H0 = .true.
        do j = c0, nc
            x(MO,j) = p%vf_MnO2/Vm(MO)
            x(ZM,j) = p%vf_ZMO/Vm(ZM)
            x(ZH,j) = p%vf_ZHS/Vm(ZH)
        end do
        x(TH,:) = log(p%theta0/(1.0_dp - p%theta0))
        call free_of(x, .true., fr, fail)
        if (fail) return
        u_zn = 0.5_dp/fRT*log(fr(2,1))
        x(P2,:) = -u_zn                                 ! anode at equilibrium
        call basis_of(x, fr, cb)
        do j = c0, nc
            x(P1,j) = x(P2,j) + u3(x(TH,j), cb(2,j))    ! insertion at equilibrium (a guess)
        end do
        ! settle the potentials at open circuit: a step so short that the compositions barely change
        call newton_step(x, 1.0e-4_dp, 0.0_dp, xs, fail)
        x = xs
    end subroutine initial_state

    real(dp) function voltage(x, I)
        real(dp), intent(in) :: x(N,nc), I
        real(dp) :: eps(nc)
        call porosity(x, eps)
        voltage = x(P1,nc) - I*0.5_dp*dxc(nc)/(p%sigma*(1.0_dp - eps(nc))*area_c(nc))
    end function voltage

    !> The physical limit a state has reached when a step at current I cannot be solved, or ''.
    function limit_reason(x, I) result(why)
        real(dp), intent(in) :: x(N,nc), I
        character(len=32) :: why
        real(dp) :: thv(nc), eps(nc), mx, mx0
        integer :: j
        why = ''
        thv = sigmoid(x(TH,:))
        call porosity(x, eps)
        mx = -huge(1.0_dp); mx0 = -huge(1.0_dp)
        do j = c0, nc
            mx = max(mx, x(MO,j) + x(ZM,j))
            mx0 = max(mx0, x0_state(MO,j) + x0_state(ZM,j))
        end do
        if (minval(x(ZN,:)) < LIMIT*minval(x0_state(ZN,:))) then
            why = 'zinc_depleted'
        else if (p%c_MnSO4 > 0 .and. minval(x(MN,:)) < LIMIT*minval(x0_state(MN,:))) then
            why = 'manganese_depleted'
        else if (p%R3_on .and. I > 0 .and. maxval(thv(c0:nc)) > 1.0_dp - LIMIT) then
            why = 'insertion_full'
        else if (p%R3_on .and. I < 0 .and. minval(thv(c0:nc)) < LIMIT) then
            why = 'insertion_empty'
        else if (mx < LIMIT*max(mx0, 1e-300_dp)) then
            why = 'dissolvable_mno2_exhausted'
        else if (minval(eps(c0:nc)) < LIMIT) then
            why = 'pores_clogged'
        end if
    end function limit_reason

    character(len=15) function fmt_e(v)
        real(dp), intent(in) :: v
        if (v /= v) then
            fmt_e = '            NAN'
        else if (.not. ieee_is_finite(v)) then
            if (v > 0) then
                fmt_e = '            INF'
            else
                fmt_e = '           -INF'
            end if
        else if (v /= 0 .and. (abs(v) >= 9.999999995e99_dp .or. abs(v) < 9.999999995e-100_dp)) then
            write(fmt_e, '(ES15.7E3)') v
        else
            write(fmt_e, '(ES15.7E2)') v
        end if
    end function fmt_e

    subroutine write_header()
        character(len=15), parameter :: cols(19) = [character(len=15) :: 't_h', 'V', 'I_mAg', 'mAhg', 'step', &
            'pH_cath', 'pH_probe', 'pH_anode', 'Zn_cath_M', 'Mn_cath_M', 'S_cath_M', 'vf_MnO2', 'vf_ZMO', &
            'vf_ZHS', 'theta', 'i_R1_mAg', 'i_R2_mAg', 'i_R3_mAg', 'eps_cath']
        character(len=16*19) :: line
        integer :: i
        line = ''
        do i = 1, 19
            line(16*(i-1)+1:16*i-1) = adjustr(cols(i))
        end do
        write(ou, '(A)') line(1:16*19-1)
    end subroutine write_header

    subroutine write_row(t, x, mAhg, I, k)
        real(dp), intent(in) :: t, x(N,nc), mAhg, I
        integer, intent(in) :: k
        real(dp) :: fr(4,nc), cb(4,nc), eps(nc), i1, i2, i3, rz, s1, s2, s3, vsum, row(19)
        character(len=16*19) :: line
        integer :: j, q
        logical :: fail
        call free_of(x, .true., fr, fail)
        call basis_of(x, fr, cb)
        call porosity(x, eps)
        s1 = 0; s2 = 0; s3 = 0
        do j = c0, nc
            call reactions(x, fr, cb, j, i1, i2, i3, rz)
            s1 = s1 + vol(j)*i1; s2 = s2 + vol(j)*i2; s3 = s3 + vol(j)*i3
        end do
        vsum = sum(vol(c0:nc))
        row(1) = t/3600.0_dp
        row(2) = voltage(x, I)
        row(3) = I*1.0e3_dp/mass
        row(4) = mAhg
        row(5) = k
        row(6) = -log10(cmean(fr(1,:)))
        row(7) = -log10(fr(1,i_probe))
        row(8) = -log10(fr(1,1))
        row(9) = cmean(x(ZN,:))*1e3_dp
        row(10) = cmean(x(MN,:))*1e3_dp
        row(11) = cmean(x(SO,:))*1e3_dp
        row(12) = cmean(x(MO,:)*Vm(MO))
        row(13) = cmean(x(ZM,:)*Vm(ZM))
        row(14) = cmean(x(ZH,:)*Vm(ZH))
        row(15) = cmean(sigmoid(x(TH,:)))
        row(16) = s1*(1.0e3_dp/mass)
        row(17) = s2*(1.0e3_dp/mass)
        row(18) = s3*(1.0e3_dp/mass)
        row(19) = cmean(eps)
        line = ''
        do q = 1, 19
            if (q == 5) then
                write(line(16*(q-1)+1:16*q-1), '(I15)') k
            else
                line(16*(q-1)+1:16*q-1) = fmt_e(row(q))
            end if
        end do
        write(ou, '(A)') line(1:16*19-1)
        n_rows = n_rows + 1
    contains
        real(dp) function cmean(v)                   ! volume average over the cathode
            real(dp), intent(in) :: v(nc)
            integer :: jj
            cmean = 0
            do jj = c0, nc
                cmean = cmean + vol(jj)*v(jj)
            end do
            cmean = cmean/vsum
        end function cmean
    end subroutine write_row

    ! ================================================================== protocol
    subroutine parse_protocol(text, steps, nsteps)
        character(len=*), intent(in) :: text
        type(step_t), allocatable, intent(out) :: steps(:)
        integer, intent(out) :: nsteps
        type(step_t) :: tmp(200), st
        character(len=512) :: rest, part, word, key
        real(dp) :: val
        integer :: sc, sp, eqp, n
        logical :: has_i, has_v
        rest = text
        n = 0
        if (len_trim(rest) == 0) call die('empty protocol')
        do
            sc = index(rest, ';')
            if (sc > 0) then
                part = adjustl(rest(1:sc-1)); rest = rest(sc+1:)
            else
                part = adjustl(rest); rest = ''
            end if
            if (len_trim(part) > 0) then
                n = n + 1
                sp = index(part, ' ')
                word = part(1:sp-1)
                call lower(word)
                part = adjustl(part(sp:))
                st = step_t()
                st%kind = trim(word)
                if (word /= 'cc' .and. word /= 'cv' .and. word /= 'rest') &
                    call die('unknown step type '''//trim(word)//''' (expected cc, cv or rest)')
                has_i = .false.; has_v = .false.
                st%Vmin = p%V_min; st%Vmax = p%V_max
                do while (len_trim(part) > 0)
                    sp = index(part, ' ')
                    word = part(1:sp-1)
                    part = adjustl(part(sp:))
                    eqp = index(word, '=')
                    if (eqp == 0) call die('expected key=value, got '''//trim(word)//'''')
                    key = word(1:eqp-1)
                    call lower(key)
                    call read_real(word(eqp+1:), val)
                    select case (trim(st%kind)//':'//trim(key))
                    case ('cc:i'); st%I = val; has_i = .true.
                    case ('cc:t', 'cv:t', 'rest:t'); st%t = val
                    case ('cc:vmin'); st%Vmin = val
                    case ('cc:vmax'); st%Vmax = val
                    case ('cv:v'); st%V = val; has_v = .true.
                    case ('cv:imin'); st%Imin = val
                    case default; call die(trim(st%kind)//' does not take '''//trim(key)//'''')
                    end select
                end do
                if (st%t /= -1 .and. st%t <= 0) call die('t must be positive')
                if (st%kind == 'cc' .and. .not. has_i) call die('cc needs I=')
                if (st%kind == 'cv' .and. .not. has_v) call die('cv needs V=')
                if (st%kind == 'cv' .and. st%t < 0 .and. st%Imin < 0) call die('cv needs t= or Imin= to end')
                if (st%kind == 'rest' .and. st%t < 0) call die('rest needs t=')
                tmp(n) = st
            end if
            if (len_trim(rest) == 0) exit
        end do
        nsteps = n*max(1, p%cycles)
        allocate(steps(nsteps))
        do sc = 1, max(1, p%cycles)
            steps((sc-1)*n+1:sc*n) = tmp(1:n)
        end do
    contains
        subroutine read_real(s, v)
            character(len=*), intent(in) :: s
            real(dp), intent(out) :: v
            character(len=64) :: buf
            integer :: i, ios
            buf = s
            do i = 1, len_trim(buf)
                if (buf(i:i) == 'd' .or. buf(i:i) == 'D') buf(i:i) = 'e'
            end do
            read(buf, *, iostat=ios) v
            if (ios /= 0) call die('bad number '''//trim(s)//'''')
        end subroutine read_real
    end subroutine parse_protocol

    subroutine lower(s)
        character(len=*), intent(inout) :: s
        integer :: i
        do i = 1, len(s)
            if (s(i:i) >= 'A' .and. s(i:i) <= 'Z') s(i:i) = achar(iachar(s(i:i)) + 32)
        end do
    end subroutine lower

    subroutine die(msg)
        character(len=*), intent(in) :: msg
        write(*, '(A)') 'error: '//msg
        stop 1
    end subroutine die

    ! ================================================================== driver
    real(dp) function margin_of(xn, st, I)
        real(dp), intent(in) :: xn(N,nc), I
        type(step_t), intent(in) :: st
        real(dp) :: v
        ! a discharge ends at Vmin, a charge at Vmax
        v = voltage(xn, I)
        if (I > 0) then
            margin_of = v - st%Vmin
        else if (I < 0) then
            margin_of = st%Vmax - v
        else
            margin_of = min(v - st%Vmin, st%Vmax - v)
        end if
    end function margin_of

    !> Distance to a full host on discharge (theta = 1 - LIMIT) or an empty one on charge (theta = LIMIT), scaled so
    !> that EVENT_DV is THETA_TOL in theta; negative once crossed, huge if none applies (see Stepper.event_margin).
    real(dp) function event_margin(xn, I)
        real(dp), intent(in) :: xn(N,nc), I
        real(dp) :: thmax, thmin, thj, d
        integer :: j
        event_margin = huge(1.0_dp)
        if (.not. p%R3_on .or. I == 0.0_dp) return
        thmax = -huge(1.0_dp); thmin = huge(1.0_dp)
        do j = c0, nc
            thj = sigmoid(xn(TH,j))
            thmax = max(thmax, thj); thmin = min(thmin, thj)
        end do
        if (I > 0) then
            d = (1.0_dp - LIMIT) - thmax
        else
            d = thmin - LIMIT
        end if
        event_margin = d*EVENT_DV/THETA_TOL
    end function event_margin

    !> Advance by dt with backward Euler, halving the sub-step on Newton failure; with a cc step, locate a
    !> cutoff crossing within EVENT_DV. On failure, state and t_done hold the progress made.
    subroutine advance(state, dt, I, use_margin, st, t_done, stopped, fail, host_event)
        real(dp), intent(inout) :: state(N,nc)
        real(dp), intent(in) :: dt, I
        logical, intent(in) :: use_margin
        type(step_t), intent(in) :: st
        real(dp), intent(out) :: t_done
        logical, intent(out) :: stopped, fail
        logical, intent(in), optional :: host_event      ! locate a full (or empty) host instead of a cutoff
        real(dp) :: hh, new(N,nc), m
        integer :: failures
        t_done = 0.0_dp; hh = dt; failures = 0
        stopped = .false.; fail = .false.
        do while (t_done < dt)
            hh = min(hh, dt - t_done)
            call newton_step(state, hh, I, new, fail)
            if (fail) then
                failures = failures + 1
                if (hh/2 < MIN_SUBSTEP .or. failures >= MAX_FAILURES) return
                fail = .false.
                hh = hh/2
                cycle
            end if
            if (use_margin) then
                m = margin_of(new, st, I)
                if (present(host_event)) then
                    if (host_event) m = event_margin(new, I)
                end if
                if (m < 0.0_dp) then
                    if (m < -EVENT_DV .and. hh/2 >= EVENT_MIN_DT) then
                        hh = hh/2
                        cycle
                    end if
                    state = new
                    t_done = t_done + hh
                    stopped = .true.
                    return
                end if
            end if
            state = new
            t_done = t_done + hh
            hh = 2.0_dp*hh
        end do
    end subroutine advance

    !> One time step at constant voltage: find the current I with V(I) = V_set (bracketing, then Illinois).
    subroutine cv_step(state, hh, V_set, I, new, fail)
        real(dp), intent(in) :: state(N,nc), hh, V_set
        real(dp), intent(inout) :: I
        real(dp), intent(out) :: new(N,nc)
        logical, intent(out) :: fail
        real(dp), allocatable :: st_x(:,:,:), st_I(:)
        integer :: nst, it, side
        real(dp) :: fI, grow, a, b, fa, fb
        logical :: have_a, have_b
        allocate(st_x(N,nc,8), st_I(8))
        nst = 0
        fail = .false.
        fI = f(I)
        if (abs(fI) <= CV_TOL) then
            call take(I); return
        end if
        grow = max(abs(I), 1.0e-3_dp*mass)
        have_a = .false.; have_b = .false.
        a = 0; b = 0; fa = 0; fb = 0
        do it = 1, 60
            if (fI > 0) then
                a = I; fa = fI; have_a = .true.
                if (have_b) exit
                if (I < 0) then
                    I = 0.0_dp
                else
                    I = I + grow
                end if
            else
                b = I; fb = fI; have_b = .true.
                if (have_a) exit
                if (I > 0) then
                    I = 0.0_dp
                else
                    I = I - grow
                end if
            end if
            grow = grow*2.0_dp
            fI = f(I)
            if (abs(fI) <= CV_TOL) then
                call take(I); return
            end if
        end do
        if (.not. (have_a .and. have_b)) then
            fail = .true.; return
        end if
        side = 0
        do it = 1, 200
            if (ieee_is_finite(fa) .and. ieee_is_finite(fb)) then
                I = (a*fb - b*fa)/(fb - fa)
                if (.not. (a < I .and. I < b)) I = 0.5_dp*(a + b)
            else
                I = 0.5_dp*(a + b)
            end if
            fI = f(I)
            if (abs(fI) <= CV_TOL .or. (b - a) <= 1.0e-14_dp*mass) then
                ! a collapsed bracket is a solution only at the set voltage (not across a jump in V(I))
                if (abs(fI) <= CV_ACCEPT) then
                    call take(I)
                else
                    fail = .true.
                end if
                return
            end if
            if (fI > 0) then
                a = I; fa = fI
                if (side == 1 .and. ieee_is_finite(fb)) fb = fb*0.5_dp
                side = 1
            else
                b = I; fb = fI
                if (side == -1 .and. ieee_is_finite(fa)) fa = fa*0.5_dp
                side = -1
            end if
        end do
        fail = .true.
    contains
        real(dp) function f(Ic)
            real(dp), intent(in) :: Ic
            real(dp) :: xn(N,nc)
            logical :: fl
            real(dp), allocatable :: tx(:,:,:), tI(:)
            call newton_step(state, hh, Ic, xn, fl)
            if (fl) then
                if (Ic < 0) then
                    f = ieee_value(1.0_dp, ieee_positive_inf)
                else
                    f = ieee_value(1.0_dp, ieee_negative_inf)
                end if
                return
            end if
            if (nst == size(st_I)) then
                allocate(tx(N,nc,2*nst), tI(2*nst))
                tx(:,:,1:nst) = st_x; tI(1:nst) = st_I
                call move_alloc(tx, st_x); call move_alloc(tI, st_I)
            end if
            nst = nst + 1
            st_x(:,:,nst) = xn; st_I(nst) = Ic
            f = voltage(xn, Ic) - V_set
        end function f
        subroutine take(Ic)                       ! the state computed at current Ic (the latest one)
            real(dp), intent(in) :: Ic
            integer :: q
            do q = nst, 1, -1
                if (st_I(q) == Ic) then
                    new = st_x(:,:,q)
                    return
                end if
            end do
            fail = .true.                         ! no feasible current at the set voltage
        end subroutine take
    end subroutine cv_step

    !> A constant-voltage sub-step: cv_step over hh, halved (down to MIN_SUBSTEP) while no current holds V_set
    !> for that long. I is unchanged on failure.
    subroutine cv_advance(state, hh, V_set, I, new, h_done, fail)
        real(dp), intent(in) :: state(N,nc), hh, V_set
        real(dp), intent(inout) :: I
        real(dp), intent(out) :: new(N,nc), h_done
        logical, intent(out) :: fail
        real(dp) :: I_guess
        I_guess = I
        h_done = hh
        do
            I = I_guess
            call cv_step(state, h_done, V_set, I, new, fail)
            if (.not. fail) return
            if (h_done/2 < MIN_SUBSTEP) then
                I = I_guess
                h_done = 0
                return
            end if
            h_done = h_done/2
        end do
    end subroutine cv_advance

    subroutine run_corrected(nml_unit, data_dir, out_file, reason, nsteps_done, mAhg_out)
        integer, intent(in) :: nml_unit
        character(len=*), intent(in) :: data_dir
        character(len=*), intent(out) :: out_file, reason
        integer, intent(out) :: nsteps_done
        real(dp), intent(out) :: mAhg_out
        type(step_t), allocatable :: steps(:)
        real(dp), allocatable :: state(:,:), new(:,:)
        real(dp) :: t, mAhg, I, t_step, hh, h_done
        integer :: nsteps, k, n_done
        logical :: stopped, fail, use_margin, hit, fail2
        real(dp) :: h_event
        real(dp), allocatable :: prog(:,:), redo(:,:)
        character(len=32) :: why
        character(len=512) :: file
        real(dp) :: last_write

        call read_params(nml_unit, p, file)
        out_file = file
        call setup(data_dir)
        call parse_protocol(trim(p%steps), steps, nsteps)
        allocate(state(N,nc), new(N,nc), prog(N,nc), redo(N,nc))
        open(newunit=ou, file=trim(out_file), status='replace', action='write')
        call write_header()
        call initial_state(state, fail)
        if (fail) call die('the initial state could not be solved')
        x0_state = state
        t = 0.0_dp; mAhg = 0.0_dp; n_done = 0
        I = first_current(steps(1))
        call write_row(t, state, mAhg, I, 1)
        last_write = t
        reason = ''
        do k = 1, nsteps
            t_step = 0.0_dp
            if (steps(k)%kind /= 'cv') I = first_current(steps(k))
            do
                hh = p%dt
                if (steps(k)%t >= 0) hh = min(p%dt, steps(k)%t - t_step)
                if (steps(k)%kind == 'cv') then
                    call cv_advance(state, hh, steps(k)%V, I, new, h_done, fail)
                    stopped = steps(k)%Imin >= 0 .and. abs(I) <= steps(k)%Imin*1.0e-3_dp*mass
                    why = 'current_limit'
                else
                    use_margin = steps(k)%kind == 'cc'
                    new = state
                    call advance(new, hh, I, use_margin, steps(k), h_done, stopped, fail)
                    if (.not. fail) then
                        if (stopped .and. voltage(new, I) <= steps(k)%Vmin) then
                            why = 'cutoff_low'
                        else
                            why = 'cutoff_high'
                        end if
                    end if
                end if
                if (fail) then
                    ! keep the sub-steps completed before the failure: the limit is judged where it was reached
                    if (h_done > 0.0_dp) then
                        prog = new
                    else
                        prog = state
                    end if
                    why = limit_reason(prog, I)
                    if ((why == 'insertion_full' .or. why == 'insertion_empty') .and. steps(k)%kind == 'cc') then
                        ! the run ends at a full (or empty) host: redo the time step and stop it where the host
                        ! reaches the limit, a well-determined state (closer to the end the rate no longer depends
                        ! on 1 - theta, and the voltage is not determined by Newton's tolerance)
                        redo = state
                        call advance(redo, hh, I, .true., steps(k), h_event, hit, fail2, host_event=.true.)
                        if (.not. fail2 .and. hit) then
                            new = redo; h_done = h_event
                        end if
                    end if
                    if (h_done > 0.0_dp) then
                        mAhg = mAhg + 1000.0_dp*(I/mass)*h_done/3600.0_dp
                        state = new; t = t + h_done; n_done = n_done + 1
                    end if
                    if (len_trim(why) == 0 .or. p%end_on_cutoff) then
                        if (len_trim(why) == 0) why = 'solver_fail'
                        call finish(why)
                        return
                    end if
                    ! a physical limit ends this step, as a cutoff does; the protocol goes on
                    call write_row(t, state, mAhg, I, k)
                    last_write = t; reason = why
                    exit
                end if
                mAhg = mAhg + 1000.0_dp*(I/mass)*h_done/3600.0_dp
                state = new; t = t + h_done; t_step = t_step + h_done; n_done = n_done + 1
                if (.not. all(ieee_is_finite(state))) then
                    call finish('nan'); return
                end if
                if (stopped .or. (steps(k)%t >= 0 .and. t_step >= steps(k)%t*(1.0_dp - 1.0e-12_dp))) then
                    call write_row(t, state, mAhg, I, k)
                    last_write = t
                    if (stopped) then
                        reason = why
                    else
                        reason = 'duration'
                    end if
                    if (stopped .and. why(1:6) == 'cutoff' .and. p%end_on_cutoff) then
                        call done(why); return
                    end if
                    exit
                end if
                if (t - last_write >= p%write_interval) then
                    call write_row(t, state, mAhg, I, k)
                    last_write = t
                end if
                if (t >= 99.0_dp*3600.0_dp) then
                    call finish('max_time'); return
                end if
            end do
        end do
        if (nsteps /= 1) reason = 'end_of_protocol'
        call done(reason)
    contains
        real(dp) function first_current(st)
            type(step_t), intent(in) :: st
            first_current = 0.0_dp
            if (st%kind == 'cc') first_current = st%I*1.0e-3_dp*mass
        end function first_current
        subroutine finish(r)
            character(len=*), intent(in) :: r
            call write_row(t, state, mAhg, I, k)
            call done(r)
        end subroutine finish
        subroutine done(r)
            character(len=*), intent(in) :: r
            reason = r
            nsteps_done = n_done
            mAhg_out = mAhg
            close(ou)
        end subroutine done
    end subroutine run_corrected

end module zn_corrected
