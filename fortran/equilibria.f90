!> Speciation for the corrected model: free-ion concentrations from the four transported totals
!> (python/znmno2_model/eqchem.py). Ideal activities (mol/L) and the 21 complexation equilibria of
!> Herrmann et al., Adv. Energy Mater. (2023), Table S1, or the SOLUTION_SPECIES of a PHREEQC database.
!> Unknowns: log10 free concentrations of the masters H+, Zn2+, Mn2+, SO4 2-; equations: the balances
!> on the totals H_T, Zn_T, Mn_T, S_T.
!>
!> SPDX-License-Identifier: BSD-3-Clause
module zn_equilibria
    use zn_params, only: dp
    implicit none
    private
    public :: eq_t, eq_init, eq_solve, eq_species, eq_sensitivity, LN10, solve_dense

    real(dp), parameter :: LN10 = log(10.0_dp)
    integer, parameter :: MAXSP = 200

    type :: eq_t
        integer :: ns = 0
        character(len=32) :: names(MAXSP)
        real(dp) :: nu(MAXSP, 4) = 0         ! masters H+, Zn+2, Mn+2, SO4-2
        real(dp) :: logk(MAXSP) = 0
        real(dp) :: z(MAXSP) = 0
        real(dp) :: abs_nu_h(MAXSP) = 0
    end type eq_t

contains

    ! ------------------------------------------------------------------ set-up
    subroutine add_species(e, name, nu, logk)
        type(eq_t), intent(inout) :: e
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: nu(4), logk
        integer :: k
        do k = 1, e%ns
            if (e%names(k) == name) exit
        end do
        if (k > e%ns) then
            if (e%ns == MAXSP) stop 'too many species in the equilibria database'
            e%ns = e%ns + 1
            k = e%ns
        end if
        e%names(k) = name
        e%nu(k,:) = nu
        e%logk(k) = logk
        e%z(k) = charge_of(name)
    end subroutine add_species

    integer function charge_of(name)
        character(len=*), intent(in) :: name
        integer :: n, i, d
        n = len_trim(name)
        charge_of = 0
        i = n
        do while (i >= 1)
            if (index('0123456789', name(i:i)) == 0) exit
            i = i - 1
        end do
        if (i < 1) return
        if (name(i:i) /= '+' .and. name(i:i) /= '-') return
        if (i == n) then
            d = 1
        else
            read(name(i+1:n), *) d
        end if
        charge_of = d
        if (name(i:i) == '-') charge_of = -d
    end function charge_of

    !> Masters, then the complexes of Herrmann et al. (2023) Table S1 or of a PHREEQC database; then the
    !> log K overrides of logk_file ('name value' per line).
    subroutine eq_init(e, database, logk_file)
        type(eq_t), intent(out) :: e
        character(len=*), intent(in) :: database, logk_file
        integer :: k
        e%ns = 0
        call add_species(e, 'H+', [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
        call add_species(e, 'Zn+2', [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
        call add_species(e, 'Mn+2', [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], 0.0_dp)
        call add_species(e, 'SO4-2', [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], 0.0_dp)
        if (len_trim(database) == 0) then
            ! (name, H+, Zn+2, Mn+2, SO4-2, log10 beta); H2O does not enter the ideal balances
            call add_species(e, 'OH-', [-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], -14.0_dp)
            call add_species(e, 'H2SO4', [2.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], 0.0_dp)
            call add_species(e, 'HSO4-', [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], 1.98_dp)
            call add_species(e, 'ZnOH+', [-1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], -7.5_dp)
            call add_species(e, 'Zn(OH)2', [-2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], -16.4_dp)
            call add_species(e, 'Zn(OH)3-', [-3.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], -28.2_dp)
            call add_species(e, 'Zn(OH)4-2', [-4.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], -41.3_dp)
            call add_species(e, 'Zn2OH+3', [-1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp], -9.0_dp)
            call add_species(e, 'Zn2(OH)6-2', [-6.0_dp, 2.0_dp, 0.0_dp, 0.0_dp], -54.3_dp)
            call add_species(e, 'Zn4(OH)4+4', [-4.0_dp, 4.0_dp, 0.0_dp, 0.0_dp], -27.0_dp)
            call add_species(e, 'ZnSO4', [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], 2.37_dp)
            call add_species(e, 'Zn(SO4)2-2', [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp], 3.28_dp)
            call add_species(e, 'Zn(SO4)3-4', [0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp], 1.7_dp)
            call add_species(e, 'Zn(SO4)4-6', [0.0_dp, 1.0_dp, 0.0_dp, 4.0_dp], 1.7_dp)
            call add_species(e, 'MnOH+', [-1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], -10.59_dp)
            call add_species(e, 'Mn(OH)2', [-2.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], -18.54_dp)
            call add_species(e, 'Mn(OH)3-', [-3.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], -34.8_dp)
            call add_species(e, 'Mn(OH)4-2', [-4.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], -48.3_dp)
            call add_species(e, 'Mn2(OH)3+', [-3.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], -23.9_dp)
            call add_species(e, 'Mn2OH+3', [-1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], -10.56_dp)
            call add_species(e, 'MnSO4', [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], 2.25_dp)
        else
            call read_phreeqc(e, database)
        end if
        if (len_trim(logk_file) > 0) call read_overrides(e, logk_file)
        do k = 1, e%ns
            e%abs_nu_h(k) = abs(e%nu(k,1))
        end do
    end subroutine eq_init

    subroutine read_overrides(e, path)
        type(eq_t), intent(inout) :: e
        character(len=*), intent(in) :: path
        character(len=512) :: line, name
        real(dp) :: val
        integer :: u, ios, k, h
        open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'logk_file not found'
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            h = index(line, '#')
            if (h > 0) line = line(1:h-1)
            if (len_trim(line) == 0) cycle
            read(line, *) name, val
            do k = 1, e%ns
                if (e%names(k) == name) exit
            end do
            if (k > e%ns) then
                write(*, '(A)') 'log K override for unknown species '''//trim(name)//''''
                stop 1
            end if
            e%logk(k) = val
        end do
        close(u)
    end subroutine read_overrides

    !> SOLUTION_SPECIES of a PHREEQC database made of H+, H2O, Zn+2, Mn+2, SO4-2 (as speciation.py's
    !> Database.from_phreeqc): log_k, or the analytic expression at 25 C; redox species left out.
    subroutine read_phreeqc(e, path)
        type(eq_t), intent(inout) :: e
        character(len=*), intent(in) :: path
        character(len=1024) :: raw, line, lhs, rhs, opt
        character(len=32) :: cname, sp(20)
        real(dp) :: cf(20), nu(4), logk, a(6), tk
        logical :: in_block, have, ok
        integer :: u, ios, h, eqp, nl, nr, i
        tk = 25.0_dp + 273.15_dp
        open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'equilibria_db not found'
        in_block = .false.; have = .false.; ok = .false.
        logk = 0; nu = 0; cname = ''
        do
            read(u, '(A)', iostat=ios) raw
            if (ios /= 0) exit
            if (.not. in_block) then
                if (trim(raw) == 'SOLUTION_SPECIES' .or. trim(raw) == 'SOLUTION_SPECIES'//achar(13)) in_block = .true.
                cycle
            end if
            if (trim(raw) == 'PHASES') exit
            do i = 1, len_trim(raw)                      ! tabs separate fields as spaces do
                if (raw(i:i) == achar(9)) raw(i:i) = ' '
            end do
            h = index(raw, '#')
            line = raw
            if (h > 0) line = raw(1:h-1)
            if (len_trim(line) == 0) cycle
            opt = adjustl(line)
            eqp = index(line, '=')
            if (eqp > 0 .and. opt(1:1) /= '-' .and. opt(1:5) /= 'log_k') then
                if (have .and. ok) call add_species(e, cname, nu, logk)
                lhs = line(1:eqp-1); rhs = line(eqp+1:)
                call parse_side(lhs, cf, sp, nl)
                call parse_side(rhs, cf(11:), sp(11:), nr)
                cname = sp(11)
                nu = 0
                ok = cf(11) == 1.0_dp .and. index(line, 'e-') == 0 .and. .not. is_master(cname) .and. cname /= 'H2O'
                do i = 1, nl
                    call add_nu(sp(i), cf(i))
                end do
                do i = 2, nr
                    call add_nu(sp(10+i), -cf(10+i))
                end do
                logk = 0
                have = .true.
                cycle
            end if
            if (.not. have) cycle
            opt = adjustl(line)
            if (opt(1:1) == '-') opt = opt(2:)
            call lower(opt)
            if (opt(1:6) == 'log_k ' .or. opt(1:5) == 'logk ') then
                read(opt(index(opt, ' '):), *) logk
            else if (opt(1:8) == 'analytic') then
                a = 0
                read(opt(index(opt, ' '):), *, iostat=ios) a
                logk = a(1) + a(2)*tk + a(3)/tk + a(4)*log10(tk) + a(5)/tk**2 + a(6)*tk**2
            end if
        end do
        if (have .and. ok) call add_species(e, cname, nu, logk)
        close(u)
    contains
        logical function is_master(s)
            character(len=*), intent(in) :: s
            is_master = s == 'H+' .or. s == 'Zn+2' .or. s == 'Mn+2' .or. s == 'SO4-2'
        end function is_master
        subroutine add_nu(s, c)
            character(len=*), intent(in) :: s
            real(dp), intent(in) :: c
            select case (trim(s))
            case ('H+'); nu(1) = nu(1) + c
            case ('Zn+2'); nu(2) = nu(2) + c
            case ('Mn+2'); nu(3) = nu(3) + c
            case ('SO4-2'); nu(4) = nu(4) + c
            case ('H2O')
            case default; ok = .false.
            end select
        end subroutine add_nu
    end subroutine read_phreeqc

    !> 'a A + b B' -> coefficients and species (terms separated by ' + ').
    subroutine parse_side(s, cf, sp, n)
        character(len=*), intent(in) :: s
        real(dp), intent(out) :: cf(:)
        character(len=32), intent(out) :: sp(:)
        integer, intent(out) :: n
        character(len=1024) :: rest, tok
        integer :: p, i
        rest = s
        n = 0
        do
            p = index(rest, ' + ')
            if (p > 0) then
                tok = adjustl(rest(1:p-1)); rest = rest(p+3:)
            else
                tok = adjustl(rest); rest = ''
            end if
            if (len_trim(tok) > 0) then
                n = n + 1
                i = 1
                do while (index('0123456789.', tok(i:i)) > 0)
                    i = i + 1
                end do
                if (i > 1) then
                    read(tok(1:i-1), *) cf(n)
                else
                    cf(n) = 1.0_dp
                end if
                sp(n) = adjustl(tok(i:))
            end if
            if (len_trim(rest) == 0) exit
        end do
    end subroutine parse_side

    subroutine lower(s)
        character(len=*), intent(inout) :: s
        integer :: i
        do i = 1, len(s)
            if (s(i:i) >= 'A' .and. s(i:i) <= 'Z') s(i:i) = achar(iachar(s(i:i)) + 32)
        end do
    end subroutine lower

    ! ------------------------------------------------------------------ dense solve
    !> Solve A X = B (n x n, m right-hand sides) by Gaussian elimination with partial pivoting.
    !> ok = .false. for a zero pivot.
    subroutine solve_dense(n, m, A, B, ok)
        integer, intent(in) :: n, m
        real(dp), intent(inout) :: A(n,n), B(n,m)
        logical, intent(out) :: ok
        integer :: k, p, i, j
        real(dp) :: f, tmp
        ok = .true.
        do k = 1, n
            p = k
            do i = k + 1, n
                if (abs(A(i,k)) > abs(A(p,k))) p = i
            end do
            if (A(p,k) == 0.0_dp) then
                ok = .false.
                return
            end if
            if (p /= k) then
                do j = 1, n
                    tmp = A(k,j); A(k,j) = A(p,j); A(p,j) = tmp
                end do
                do j = 1, m
                    tmp = B(k,j); B(k,j) = B(p,j); B(p,j) = tmp
                end do
            end if
            do i = k + 1, n
                f = A(i,k)/A(k,k)
                do j = k + 1, n
                    A(i,j) = A(i,j) - f*A(k,j)
                end do
                do j = 1, m
                    B(i,j) = B(i,j) - f*B(k,j)
                end do
            end do
        end do
        do k = n, 1, -1
            do j = 1, m
                do i = k + 1, n
                    B(k,j) = B(k,j) - A(k,i)*B(i,j)
                end do
                B(k,j) = B(k,j)/A(k,k)
            end do
        end do
    end subroutine solve_dense

    ! ------------------------------------------------------------------ speciation
    !> Concentrations of every species [mol/L] from x = log10 master concentrations.
    subroutine eq_species(e, x, m)
        type(eq_t), intent(in) :: e
        real(dp), intent(in) :: x(4)
        real(dp), intent(out) :: m(:)
        integer :: s
        real(dp) :: a
        do s = 1, e%ns
            a = e%logk(s) + (((x(1)*e%nu(s,1) + x(2)*e%nu(s,2)) + x(3)*e%nu(s,3)) + x(4)*e%nu(s,4))
            m(s) = 10.0_dp**min(max(a, -300.0_dp), 300.0_dp)
        end do
    end subroutine eq_species

    subroutine residual(e, x, T, F, m, bal, d)
        type(eq_t), intent(in) :: e
        real(dp), intent(in) :: x(4), T(4)
        real(dp), intent(out) :: F(4), m(:), bal(4), d
        integer :: s, k
        call eq_species(e, x, m)
        bal = 0; d = 0
        do s = 1, e%ns
            do k = 1, 4
                bal(k) = bal(k) + m(s)*e%nu(s,k)
            end do
            d = d + m(s)*e%abs_nu_h(s)
        end do
        F(1) = (bal(1) - T(1))/d
        do k = 2, 4
            F(k) = log10(bal(k)) - log10(T(k))
        end do
    end subroutine residual

    subroutine jacobian(e, m, bal, d, T, J)
        type(eq_t), intent(in) :: e
        real(dp), intent(in) :: m(:), bal(4), d, T(4)
        real(dp), intent(out) :: J(4,4)
        real(dp) :: db(4,4), dd(4), dm
        integer :: s, k, q
        db = 0; dd = 0
        do s = 1, e%ns
            do q = 1, 4
                dm = m(s)*e%nu(s,q)*LN10
                do k = 1, 4
                    db(k,q) = db(k,q) + e%nu(s,k)*dm
                end do
                dd(q) = dd(q) + e%abs_nu_h(s)*dm
            end do
        end do
        do q = 1, 4
            J(1,q) = (db(1,q)*d - (bal(1) - T(1))*dd(q))/(d*d)
            do k = 2, 4
                J(k,q) = db(k,q)/(bal(k)*LN10)
            end do
        end do
    end subroutine jacobian

    subroutine initial_guess(T, x)
        real(dp), intent(in) :: T(4)
        real(dp), intent(out) :: x(4)
        integer :: k
        x(1) = log10(max(T(1), 0.0_dp) + 1e-5_dp)
        do k = 2, 4
            x(k) = log10(max(T(k), 1e-300_dp)) - 0.5_dp
        end do
    end subroutine initial_guess

    !> Free-ion log10 concentrations x for totals T = [H_T, Zn_T, Mn_T, S_T] [mol/L]: Newton with the step
    !> scaled to at most one decade and backtracking on |F|, then a bracketed fallback. ok = .false. if it
    !> does not converge. x holds the warm start on entry when warm = .true.
    subroutine eq_solve(e, Tin, x, warm, ok)
        type(eq_t), intent(in) :: e
        real(dp), intent(in) :: Tin(4)
        real(dp), intent(inout) :: x(4)
        logical, intent(in) :: warm
        logical, intent(out) :: ok
        real(dp), parameter :: tol = 1e-13_dp
        integer, parameter :: max_iter = 100
        real(dp) :: T(4), F(4), F1(4), m(e%ns), m1(e%ns), bal(4), bal1(4), d, d1, J(4,4), dx(4,1), err, big, f0, lam
        integer :: it, ib
        logical :: sok
        T = Tin
        T(2:4) = max(T(2:4), 1e-30_dp)
        if (.not. warm) call initial_guess(T, x)
        call residual(e, x, T, F, m, bal, d)
        do it = 1, max_iter
            err = maxval(abs(F))
            if (err < tol) then
                ok = .true.
                return
            end if
            call jacobian(e, m, bal, d, T, J)
            dx(:,1) = -F
            call solve_dense(4, 1, J, dx, sok)
            if (.not. sok) dx = 0
            where (.not. (abs(dx(:,1)) <= huge(1.0_dp))) dx(:,1) = 0.0_dp
            big = maxval(abs(dx(:,1)))
            dx(:,1) = dx(:,1)*min(1.0_dp, 1.0_dp/max(big, 1e-300_dp))
            f0 = norm4(F)
            lam = 1.0_dp
            do ib = 1, 40
                call residual(e, x + lam*dx(:,1), T, F1, m1, bal1, d1)
                if (norm4(F1) < f0) exit
                lam = 0.5_dp*lam
            end do
            x = x + lam*dx(:,1)
            F = F1; m = m1; bal = bal1; d = d1
        end do
        if (maxval(abs(F)) < 1e-10_dp) then
            ok = .true.
            return
        end if
        call bracketed(e, T, x)
        call residual(e, x, T, F, m, bal, d)
        ok = maxval(abs(F)) < 1e-10_dp
    end subroutine eq_solve

    pure real(dp) function norm4(v)
        real(dp), intent(in) :: v(4)
        norm4 = sqrt(((v(1)*v(1) + v(2)*v(2)) + v(3)*v(3)) + v(4)*v(4))
    end function norm4

    !> Zn, Mn and S balances at fixed log10 c_H (xh): Newton on log10 c_Zn, c_Mn, c_SO4.
    subroutine inner(e, xh, T, x)
        type(eq_t), intent(in) :: e
        real(dp), intent(in) :: xh, T(4)
        real(dp), intent(inout) :: x(4)
        real(dp) :: F(4), m(e%ns), bal(4), d, J(4,4), J3(3,3), dx(3,1), big
        integer :: it
        logical :: sok
        x(1) = xh
        do it = 1, 200
            call residual(e, x, T, F, m, bal, d)
            if (all(abs(F(2:4)) < 1e-13_dp)) exit
            call jacobian(e, m, bal, d, T, J)
            J3 = J(2:4, 2:4)
            dx(:,1) = -F(2:4)
            call solve_dense(3, 1, J3, dx, sok)
            if (.not. sok) dx = 0
            where (.not. (abs(dx(:,1)) <= huge(1.0_dp))) dx(:,1) = 0.0_dp
            big = maxval(abs(dx(:,1)))
            x(2:4) = x(2:4) + dx(:,1)*min(1.0_dp, 1.0_dp/max(big, 1e-300_dp))
        end do
    end subroutine inner

    !> Robust fallback: bisection on log10 c_H in [-16, 2], the other balances solved at each trial pH.
    subroutine bracketed(e, T, x)
        type(eq_t), intent(in) :: e
        real(dp), intent(in) :: T(4)
        real(dp), intent(inout) :: x(4)
        real(dp) :: lo, hi, mid, g, m(e%ns)
        integer :: it, s
        lo = -16.0_dp; hi = 2.0_dp
        call initial_guess(T, x)
        do it = 1, 64
            mid = 0.5_dp*(lo + hi)
            call inner(e, mid, T, x)
            call eq_species(e, x, m)
            g = 0
            do s = 1, e%ns
                g = g + m(s)*e%nu(s,1)
            end do
            g = g - T(1)
            if (g < 0) then
                lo = mid
            else
                hi = mid
            end if
            if (hi - lo < 1e-12_dp) exit
        end do
        call inner(e, 0.5_dp*(lo + hi), T, x)
    end subroutine bracketed

    !> d x / d T (4 x 4) at a converged speciation x of totals T (implicit function theorem).
    subroutine eq_sensitivity(e, x, Tin, S)
        type(eq_t), intent(in) :: e
        real(dp), intent(in) :: x(4), Tin(4)
        real(dp), intent(out) :: S(4,4)
        real(dp) :: T(4), F(4), m(e%ns), bal(4), d, J(4,4)
        integer :: k
        logical :: sok
        T = Tin
        T(2:4) = max(T(2:4), 1e-30_dp)
        call residual(e, x, T, F, m, bal, d)
        call jacobian(e, m, bal, d, T, J)
        S = 0
        S(1,1) = 1.0_dp/d
        do k = 2, 4
            S(k,k) = 1.0_dp/(T(k)*LN10)
        end do
        call solve_dense(4, 4, J, S, sok)
    end subroutine eq_sensitivity

end module zn_equilibria
