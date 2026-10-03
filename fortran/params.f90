!> Parameters of the corrected model (units: cm, s, mol, A, V; electrolyte inputs in mol/L), read from the
!> namelist input file shared with the Python package and the C++ program. Defaults and meanings are those of
!> python/znmno2_model/params.py (docs/parameters.md); tests/test_ports.py checks that they agree.
!>
!> SPDX-License-Identifier: BSD-3-Clause
module zn_params
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: dp, params_t, read_params, M_ZMO, M_host, active_mass

    type :: params_t
        real(dp) :: A_cell = 0.178134094_dp                ! cross-section of the cell [cm2], the same in every region (1-D)
        real(dp) :: L_probe = 0.1_dp                       ! probe region length [cm]; 0 = no probe region
        real(dp) :: eps_probe = 1.0_dp                     ! probe region: open electrolyte
        real(dp) :: tau_probe = 1.0_dp                     ! probe region tortuosity
        real(dp) :: L_sep = 0.06_dp                        ! separator thickness [cm]
        real(dp) :: eps_sep = 0.9_dp                       ! separator porosity
        real(dp) :: tau_factor_sep = 2.0_dp                ! tortuosity = tau_factor * eps**bruggeman
        real(dp) :: bruggeman_sep = -0.5_dp                ! separator Bruggeman exponent
        real(dp) :: L_cath = 0.0218_dp                     ! cathode thickness [cm]
        real(dp) :: eps_cath = 0.815_dp                    ! initial cathode porosity
        real(dp) :: tau_factor_cath = 2.0_dp               ! cathode tortuosity = tau_factor_cath * eps**bruggeman_cath
        real(dp) :: bruggeman_cath = -0.5_dp               ! cathode Bruggeman exponent
        integer :: n_probe = 20                            ! finite-volume cells per region
        integer :: n_sep = 30                              ! separator cells
        integer :: n_cath = 40                             ! cathode cells
        real(dp) :: sigma = 0.1_dp                         ! cathode solid conductivity [S/cm]; effective sigma (1 - eps)
        real(dp) :: vf_MnO2 = 0.0001_dp                    ! pristine MnO2 (R1)
        real(dp) :: vf_ZMO = 0.01_dp                       ! Zn_z MnO2, dissolution/deposition (R2)
        real(dp) :: vf_host = 0.03_dp                      ! Zn-insertion host (R3), fixed
        real(dp) :: vf_ZHS = 0.0_dp                        ! zinc hydroxide sulfate
        real(dp) :: M_MnO2 = 86.937_dp                     ! MnO2 molar mass [g/mol]
        real(dp) :: rho_MnO2 = 5.03_dp                     ! MnO2 density [g/cm3]
        real(dp) :: rho_ZMO = 5.0_dp                       ! Zn_z MnO2 density [g/cm3]
        real(dp) :: rho_host = 5.0_dp                      ! insertion host density [g/cm3]
        real(dp) :: M_ZHS = 549.819_dp                     ! Zn4SO4(OH)6 . 5 H2O
        real(dp) :: rho_ZHS = 2.67_dp                      ! ZHS density [g/cm3]
        real(dp) :: r_MnO2 = 0.002_dp                      ! MnO2 particle radius [cm]
        real(dp) :: r_ZMO = 0.002_dp                       ! ZMO particle radius [cm]
        real(dp) :: r_host = 0.002_dp                      ! host particle radius [cm]
        real(dp) :: r_ZHS = 0.002_dp                       ! ZHS (and ZnO, Zn(OH)2) particle radius [cm]
        real(dp) :: z_ZMO = 0.5_dp                         ! Zn per Mn in the dissolving phase
        real(dp) :: zmin = 0.2_dp                          ! insertion range of the host (Zn per Mn)
        real(dp) :: zmax = 0.5_dp                          ! Zn per Mn of the full host
        real(dp) :: theta0 = 0.001_dp                      ! initial insertion fraction (0 = zmin, charged)
        real(dp) :: mass_AM = 0.0_dp                       ! active mass [g] for mAh/g; 0 = MnO2-equivalent mass of the solids
        real(dp) :: c_ZnSO4 = 2.0_dp                       ! initial ZnSO4 [mol/L]
        real(dp) :: c_MnSO4 = 0.05_dp                      ! initial MnSO4 [mol/L]
        real(dp) :: c_H2SO4 = 0.0_dp                       ! initial H2SO4 [mol/L]
        real(dp) :: D_Zn = 7.15e-06_dp                     ! Zn2+ diffusion coefficient [cm2/s]
        real(dp) :: D_Mn = 6.88e-06_dp                     ! Mn2+ diffusion coefficient [cm2/s]
        real(dp) :: D_SO4 = 1.07e-05_dp                    ! SO4 2- diffusion coefficient [cm2/s]
        real(dp) :: D_H = 9e-05_dp                         ! H+ diffusion coefficient [cm2/s]
        real(dp) :: U1 = 1.986_dp                          ! R1 MnO2 + 4H+ + 2e -> Mn2+ + 2H2O
        real(dp) :: k1 = 1e-10_dp                          ! [mol/cm2/s]
        real(dp) :: alpha1 = 0.5_dp                        ! R1 transfer coefficient
        real(dp) :: U2 = 2.49_dp                           ! R2 Zn_z MnO2 + 4H+ + (2-2z)e <-> z Zn2+ + Mn2+ + 2H2O
        real(dp) :: k2 = 1e-09_dp                          ! R2 rate constant [mol/cm2/s]
        real(dp) :: alpha2 = 0.5_dp                        ! R2 transfer coefficient
        real(dp) :: a_seed_R2 = 10.0_dp                    ! deposition area besides ZMO and ZHS [cm2/cm3]
        real(dp) :: k3 = 1e-09_dp                          ! R3 insertion
        real(dp) :: alpha3 = 0.5_dp                        ! R3 transfer coefficient
        real(dp) :: V_at_zmin = 1.75_dp                    ! empirical OCP spline end points [V]
        real(dp) :: V_at_zmax = 1.45_dp                    ! R3 OCP spline value at z_max [V]
        real(dp) :: c_ref3 = 2.0_dp                        ! reference Zn2+ concentration of the R3 OCP [mol/L]
        real(dp) :: k_an = 1e-06_dp                        ! Zn anode, i0 = F k_an sqrt(c_Zn2+)
        real(dp) :: alpha_an = 0.5_dp                      ! anode transfer coefficient
        real(dp) :: logK_ZHS = 28.4_dp                     ! log10([Zn2+]^4 [SO4 2-] / [H+]^6) at saturation (Herrmann et al.)
        real(dp) :: k_ZHS = 1e-07_dp                       ! precipitation/dissolution rate constant [mol/cm2/s]
        real(dp) :: a_seed_ZHS = 10.0_dp                   ! precipitation area besides existing ZHS [cm2/cm3]
        character(len=512) :: ph_mode = 'speciation'       ! pH in the reactions: speciation | zhs_equilibrium | fixed | spline
        real(dp) :: pH_fixed = 4.8_dp                      ! ph_mode = fixed (R2Fixed of Bernard et al. 2025)
        character(len=512) :: species = 'with_H'           ! with_H: H_T transported | no_H: H_T held at its initial value (charge line)
        character(len=512) :: transport = 'ions'           ! ions: each total moves with its ion's D | quasi: species fluxes summed
        character(len=512) :: basis = 'free'               ! Zn, Mn, SO4 in the Nernst and rate terms: free (speciation) | totals
        logical :: R1_on = .true.                          ! R1 (pristine MnO2 dissolution) on
        logical :: R2_on = .true.                          ! R2 (ZMO dissolution/deposition) on
        logical :: R3_on = .true.                          ! R3 (Zn insertion) on
        character(len=512) :: zhs = 'kinetic'              ! kinetic | equilibrium (instantaneous) | lumped (into R1/R2) | off
        real(dp) :: zhs_nucleation = 1.0_dp                ! Zn supersaturation c_Zn/c_sat needed for growth on the seed area (Herrmann: 1.05)
        logical :: ZnO_on = .false.                        ! extra precipitates (kinetic, same law as ZHS)
        logical :: ZnOH2_on = .false.                      ! Zn(OH)2 precipitation on
        real(dp) :: logK_ZnO = 11.17_dp                    ! log10([Zn2+]/[H+]^2) at saturation (Herrmann & Horstmann 2024, Table 1)
        real(dp) :: logK_ZnOH2 = 12.45_dp                  ! log10([Zn2+]/[H+]^2) at Zn(OH)2 saturation
        real(dp) :: k_ZnO = 1e-07_dp                       ! ZnO precipitation rate constant [mol/cm2/s]
        real(dp) :: k_ZnOH2 = 1e-07_dp                     ! Zn(OH)2 precipitation rate constant [mol/cm2/s]
        real(dp) :: M_ZnO = 81.38_dp                       ! ZnO molar mass [g/mol]
        real(dp) :: rho_ZnO = 5.61_dp                      ! ZnO density [g/cm3]
        real(dp) :: M_ZnOH2 = 99.42_dp                     ! Zn(OH)2 molar mass [g/mol]
        real(dp) :: rho_ZnOH2 = 3.05_dp                    ! Zn(OH)2 density [g/cm3]
        character(len=512) :: r3_ocp = 'spline_nernst'     ! spline_nernst: empirical OCP with Nernstian ends | spline | nernst
        real(dp) :: r3_end_width = 0.01_dp                 ! spline_nernst: the end terms act within about this fraction of theta = 0 and 1
        real(dp) :: U3_nernst = 1.55_dp                    ! Herrmann et al. (2024), U_ins,Zn
        character(len=512) :: logk_file = ''               ! log K overrides ('name value' per line); empty: literature values
        character(len=512) :: equilibria_db = ''           ! PHREEQC database for the equilibria; empty: Herrmann et al. (2023) Table S1
        real(dp) :: D_OH = 5.27e-05_dp                     ! species diffusion coefficients for transport = quasi [cm2/s]
        real(dp) :: D_HSO4 = 1.33e-05_dp                   ! HSO4- diffusion coefficient (transport = quasi) [cm2/s]
        real(dp) :: D_complex = 5e-06_dp                   ! every other complex (assumed; plan Q-1)
        real(dp) :: R = 8.314462618_dp                     ! gas constant [J/mol/K]
        real(dp) :: T = 298.15_dp                          ! temperature [K]
        real(dp) :: F = 96485.33212_dp                     ! Faraday constant [C/mol]
        character(len=512) :: steps = 'cc I=100 Vmin=1.0'  ! protocol (docs/protocol.md); I in mA/g, positive = discharge
        integer :: cycles = 1                              ! number of times the protocol is repeated
        logical :: end_on_cutoff = .true.                  ! a voltage cutoff ends the protocol (GITT); False: go to the next step
        real(dp) :: V_min = 1.0_dp                         ! default lower cutoff [V]
        real(dp) :: V_max = 1.9_dp                         ! default upper cutoff [V]
        real(dp) :: dt = 10.0_dp                           ! [s]
        real(dp) :: write_interval = 60.0_dp               ! [s]
        real(dp) :: newton_tol = 1e-09_dp                  ! Newton update tolerance (scaled)
        integer :: newton_max_iter = 30                    ! Newton iterations per (sub-)step
    end type params_t

contains

    !> Read the groups &cell ... &output from an open namelist file (absent groups keep the defaults).
    subroutine read_params(u, p, file)
        integer, intent(in) :: u
        type(params_t), intent(inout) :: p
        character(len=*), intent(out) :: file
        character(len=12), parameter :: groups(9) = [character(len=12) :: 'cell', 'solids', 'electrolyte', &
            'reactions', 'options', 'constants', 'protocol', 'numerics', 'output']
        integer :: ios, g
        character(len=256) :: msg
        real(dp) :: A_cell, L_probe, eps_probe, tau_probe, L_sep, eps_sep, tau_factor_sep, bruggeman_sep
        real(dp) :: L_cath, eps_cath, tau_factor_cath, bruggeman_cath, sigma, vf_MnO2, vf_ZMO, vf_host
        real(dp) :: vf_ZHS, M_MnO2, rho_MnO2, rho_ZMO, rho_host, M_ZHS, rho_ZHS, r_MnO2
        real(dp) :: r_ZMO, r_host, r_ZHS, z_ZMO, zmin, zmax, theta0, mass_AM
        real(dp) :: c_ZnSO4, c_MnSO4, c_H2SO4, D_Zn, D_Mn, D_SO4, D_H, U1
        real(dp) :: k1, alpha1, U2, k2, alpha2, a_seed_R2, k3, alpha3
        real(dp) :: V_at_zmin, V_at_zmax, c_ref3, k_an, alpha_an, logK_ZHS, k_ZHS, a_seed_ZHS
        real(dp) :: pH_fixed, zhs_nucleation, logK_ZnO, logK_ZnOH2, k_ZnO, k_ZnOH2, M_ZnO, rho_ZnO
        real(dp) :: M_ZnOH2, rho_ZnOH2, r3_end_width, U3_nernst, D_OH, D_HSO4, D_complex, R
        real(dp) :: T, F, V_min, V_max, dt, newton_tol, write_interval
        integer :: n_probe, n_sep, n_cath, cycles, newton_max_iter
        character(len=512) :: ph_mode, species, transport, basis, zhs, r3_ocp, logk_file, equilibria_db
        character(len=512) :: steps
        logical :: R1_on, R2_on, R3_on, ZnO_on, ZnOH2_on, end_on_cutoff
        namelist /cell/ A_cell, L_probe, eps_probe, tau_probe, L_sep, eps_sep, tau_factor_sep, bruggeman_sep, &
            L_cath, eps_cath, tau_factor_cath, bruggeman_cath, n_probe, n_sep, n_cath, sigma
        namelist /solids/ vf_MnO2, vf_ZMO, vf_host, vf_ZHS, M_MnO2, rho_MnO2, rho_ZMO, rho_host, M_ZHS, &
            rho_ZHS, r_MnO2, r_ZMO, r_host, r_ZHS, z_ZMO, zmin, zmax, theta0, mass_AM
        namelist /electrolyte/ c_ZnSO4, c_MnSO4, c_H2SO4, D_Zn, D_Mn, D_SO4, D_H
        namelist /reactions/ U1, k1, alpha1, U2, k2, alpha2, a_seed_R2, k3, alpha3, V_at_zmin, V_at_zmax, &
            c_ref3, k_an, alpha_an, logK_ZHS, k_ZHS, a_seed_ZHS
        namelist /options/ ph_mode, pH_fixed, species, transport, basis, R1_on, R2_on, R3_on, zhs, &
            zhs_nucleation, ZnO_on, ZnOH2_on, logK_ZnO, logK_ZnOH2, k_ZnO, k_ZnOH2, M_ZnO, rho_ZnO, M_ZnOH2, &
            rho_ZnOH2, r3_ocp, r3_end_width, U3_nernst, logk_file, equilibria_db, D_OH, D_HSO4, D_complex
        namelist /constants/ R, T, F
        namelist /protocol/ steps, cycles, end_on_cutoff, V_min, V_max
        namelist /numerics/ dt, newton_tol, newton_max_iter
        namelist /output/ file, write_interval

        A_cell = p%A_cell; L_probe = p%L_probe; eps_probe = p%eps_probe; tau_probe = p%tau_probe
        L_sep = p%L_sep; eps_sep = p%eps_sep; tau_factor_sep = p%tau_factor_sep; bruggeman_sep = p%bruggeman_sep
        L_cath = p%L_cath; eps_cath = p%eps_cath; tau_factor_cath = p%tau_factor_cath; bruggeman_cath = p%bruggeman_cath
        n_probe = p%n_probe; n_sep = p%n_sep; n_cath = p%n_cath; sigma = p%sigma
        vf_MnO2 = p%vf_MnO2; vf_ZMO = p%vf_ZMO; vf_host = p%vf_host; vf_ZHS = p%vf_ZHS
        M_MnO2 = p%M_MnO2; rho_MnO2 = p%rho_MnO2; rho_ZMO = p%rho_ZMO; rho_host = p%rho_host
        M_ZHS = p%M_ZHS; rho_ZHS = p%rho_ZHS; r_MnO2 = p%r_MnO2; r_ZMO = p%r_ZMO
        r_host = p%r_host; r_ZHS = p%r_ZHS; z_ZMO = p%z_ZMO; zmin = p%zmin
        zmax = p%zmax; theta0 = p%theta0; mass_AM = p%mass_AM; c_ZnSO4 = p%c_ZnSO4
        c_MnSO4 = p%c_MnSO4; c_H2SO4 = p%c_H2SO4; D_Zn = p%D_Zn; D_Mn = p%D_Mn
        D_SO4 = p%D_SO4; D_H = p%D_H; U1 = p%U1; k1 = p%k1
        alpha1 = p%alpha1; U2 = p%U2; k2 = p%k2; alpha2 = p%alpha2
        a_seed_R2 = p%a_seed_R2; k3 = p%k3; alpha3 = p%alpha3; V_at_zmin = p%V_at_zmin
        V_at_zmax = p%V_at_zmax; c_ref3 = p%c_ref3; k_an = p%k_an; alpha_an = p%alpha_an
        logK_ZHS = p%logK_ZHS; k_ZHS = p%k_ZHS; a_seed_ZHS = p%a_seed_ZHS; ph_mode = p%ph_mode
        pH_fixed = p%pH_fixed; species = p%species; transport = p%transport; basis = p%basis
        R1_on = p%R1_on; R2_on = p%R2_on; R3_on = p%R3_on; zhs = p%zhs
        zhs_nucleation = p%zhs_nucleation; ZnO_on = p%ZnO_on; ZnOH2_on = p%ZnOH2_on; logK_ZnO = p%logK_ZnO
        logK_ZnOH2 = p%logK_ZnOH2; k_ZnO = p%k_ZnO; k_ZnOH2 = p%k_ZnOH2; M_ZnO = p%M_ZnO
        rho_ZnO = p%rho_ZnO; M_ZnOH2 = p%M_ZnOH2; rho_ZnOH2 = p%rho_ZnOH2; r3_ocp = p%r3_ocp
        r3_end_width = p%r3_end_width; U3_nernst = p%U3_nernst; logk_file = p%logk_file; equilibria_db = p%equilibria_db
        D_OH = p%D_OH; D_HSO4 = p%D_HSO4; D_complex = p%D_complex; R = p%R
        T = p%T; F = p%F; steps = p%steps; cycles = p%cycles
        end_on_cutoff = p%end_on_cutoff; V_min = p%V_min; V_max = p%V_max; dt = p%dt
        newton_tol = p%newton_tol; newton_max_iter = p%newton_max_iter; write_interval = p%write_interval
        file = 'znmno2_out.txt'
        do g = 1, size(groups)
            rewind(u)
            select case (g)
            case (1); read(u, nml=cell, iostat=ios, iomsg=msg)
            case (2); read(u, nml=solids, iostat=ios, iomsg=msg)
            case (3); read(u, nml=electrolyte, iostat=ios, iomsg=msg)
            case (4); read(u, nml=reactions, iostat=ios, iomsg=msg)
            case (5); read(u, nml=options, iostat=ios, iomsg=msg)
            case (6); read(u, nml=constants, iostat=ios, iomsg=msg)
            case (7); read(u, nml=protocol, iostat=ios, iomsg=msg)
            case (8); read(u, nml=numerics, iostat=ios, iomsg=msg)
            case (9); read(u, nml=output, iostat=ios, iomsg=msg)
            end select
            if (ios > 0) then
                write(*, '(A)') 'error reading &'//trim(groups(g))//': '//trim(msg)
                stop 1
            end if
        end do
        p%A_cell = A_cell; p%L_probe = L_probe; p%eps_probe = eps_probe; p%tau_probe = tau_probe
        p%L_sep = L_sep; p%eps_sep = eps_sep; p%tau_factor_sep = tau_factor_sep; p%bruggeman_sep = bruggeman_sep
        p%L_cath = L_cath; p%eps_cath = eps_cath; p%tau_factor_cath = tau_factor_cath; p%bruggeman_cath = bruggeman_cath
        p%n_probe = n_probe; p%n_sep = n_sep; p%n_cath = n_cath; p%sigma = sigma
        p%vf_MnO2 = vf_MnO2; p%vf_ZMO = vf_ZMO; p%vf_host = vf_host; p%vf_ZHS = vf_ZHS
        p%M_MnO2 = M_MnO2; p%rho_MnO2 = rho_MnO2; p%rho_ZMO = rho_ZMO; p%rho_host = rho_host
        p%M_ZHS = M_ZHS; p%rho_ZHS = rho_ZHS; p%r_MnO2 = r_MnO2; p%r_ZMO = r_ZMO
        p%r_host = r_host; p%r_ZHS = r_ZHS; p%z_ZMO = z_ZMO; p%zmin = zmin
        p%zmax = zmax; p%theta0 = theta0; p%mass_AM = mass_AM; p%c_ZnSO4 = c_ZnSO4
        p%c_MnSO4 = c_MnSO4; p%c_H2SO4 = c_H2SO4; p%D_Zn = D_Zn; p%D_Mn = D_Mn
        p%D_SO4 = D_SO4; p%D_H = D_H; p%U1 = U1; p%k1 = k1
        p%alpha1 = alpha1; p%U2 = U2; p%k2 = k2; p%alpha2 = alpha2
        p%a_seed_R2 = a_seed_R2; p%k3 = k3; p%alpha3 = alpha3; p%V_at_zmin = V_at_zmin
        p%V_at_zmax = V_at_zmax; p%c_ref3 = c_ref3; p%k_an = k_an; p%alpha_an = alpha_an
        p%logK_ZHS = logK_ZHS; p%k_ZHS = k_ZHS; p%a_seed_ZHS = a_seed_ZHS; p%ph_mode = ph_mode
        p%pH_fixed = pH_fixed; p%species = species; p%transport = transport; p%basis = basis
        p%R1_on = R1_on; p%R2_on = R2_on; p%R3_on = R3_on; p%zhs = zhs
        p%zhs_nucleation = zhs_nucleation; p%ZnO_on = ZnO_on; p%ZnOH2_on = ZnOH2_on; p%logK_ZnO = logK_ZnO
        p%logK_ZnOH2 = logK_ZnOH2; p%k_ZnO = k_ZnO; p%k_ZnOH2 = k_ZnOH2; p%M_ZnO = M_ZnO
        p%rho_ZnO = rho_ZnO; p%M_ZnOH2 = M_ZnOH2; p%rho_ZnOH2 = rho_ZnOH2; p%r3_ocp = r3_ocp
        p%r3_end_width = r3_end_width; p%U3_nernst = U3_nernst; p%logk_file = logk_file; p%equilibria_db = equilibria_db
        p%D_OH = D_OH; p%D_HSO4 = D_HSO4; p%D_complex = D_complex; p%R = R
        p%T = T; p%F = F; p%steps = steps; p%cycles = cycles
        p%end_on_cutoff = end_on_cutoff; p%V_min = V_min; p%V_max = V_max; p%dt = dt
        p%newton_tol = newton_tol; p%newton_max_iter = newton_max_iter; p%write_interval = write_interval
        call check(p%ph_mode, [character(len=16) :: 'speciation', 'zhs_equilibrium', 'fixed', 'spline'], 'ph_mode')
        call check(p%species, [character(len=16) :: 'with_H', 'no_H'], 'species')
        call check(p%transport, [character(len=16) :: 'ions', 'quasi'], 'transport')
        call check(p%basis, [character(len=16) :: 'free', 'totals'], 'basis')
        call check(p%zhs, [character(len=16) :: 'kinetic', 'equilibrium', 'lumped', 'off'], 'zhs')
        call check(p%r3_ocp, [character(len=16) :: 'spline_nernst', 'spline', 'nernst'], 'r3_ocp')
    contains
        subroutine check(val, allowed, name)
            character(len=*), intent(in) :: val, allowed(:), name
            if (.not. any(allowed == val)) then
                write(*, '(A)') name//' must be one of the allowed values, got '''//trim(val)//''''
                stop 1
            end if
        end subroutine check
    end subroutine read_params

    real(dp) function M_ZMO(p)
        type(params_t), intent(in) :: p
        M_ZMO = 65.38_dp*p%z_ZMO + 54.938_dp + 2*15.999_dp
    end function M_ZMO

    real(dp) function M_host(p)
        type(params_t), intent(in) :: p
        M_host = 65.38_dp*p%zmin + 54.938_dp + 2*15.999_dp
    end function M_host

    !> Active mass [g]: mass_AM, or the solids' MnO2-equivalent mass (mol Mn x M_MnO2).
    real(dp) function active_mass(p)
        type(params_t), intent(in) :: p
        real(dp) :: mol_mn
        if (p%mass_AM > 0) then
            active_mass = p%mass_AM
            return
        end if
        mol_mn = p%vf_MnO2*p%rho_MnO2/p%M_MnO2 + p%vf_ZMO*p%rho_ZMO/M_ZMO(p) + p%vf_host*p%rho_host/M_host(p)
        active_mass = mol_mn*p%M_MnO2*p%A_cell*p%L_cath
    end function active_mass

end module zn_params
