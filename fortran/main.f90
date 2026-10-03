!> Zn/MnO2 cell model: usage  znmno2_f [input.nml]   (default input file: znmno2.nml)
!>
!> &run mode = 'corrected' (default), 'faithful_charge' or 'faithful_phcell'; data_dir = folder of the spline tables.
!> The output file is set in &output. See docs/model.md and docs/parameters.md.
!>
!> SPDX-License-Identifier: BSD-3-Clause
program znmno2
    use zn_faithful, only: run_faithful
    use zn_corrected, only: run_corrected
    implicit none
    character(len=32) :: mode = 'corrected'
    character(len=512) :: data_dir = 'python/znmno2_model/data'
    character(len=512) :: file = 'znmno2_out.txt'
    real(8) :: write_interval = 60.0d0
    character(len=512) :: input_file
    character(len=128) :: reason
    character(len=512) :: out_file
    integer :: u, ios, nsteps
    real(8) :: mAhg
    namelist /run/ mode, data_dir
    namelist /output/ file, write_interval

    input_file = 'znmno2.nml'
    if (command_argument_count() >= 1) call get_command_argument(1, input_file)
    open(newunit=u, file=trim(input_file), status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'cannot open the input file'
    read(u, nml=run, iostat=ios)
    if (ios > 0) stop 'error reading &run'
    rewind(u)
    read(u, nml=output, iostat=ios)
    if (ios > 0) stop 'error reading &output'
    select case (trim(mode))
    case ('faithful_charge', 'faithful_phcell')
        call run_faithful(mode(10:15), u, trim(data_dir), trim(file), reason)
        write(*, '(A)') trim(mode)//': '//trim(reason)//'; wrote '//trim(file)
    case ('corrected')
        call run_corrected(u, trim(data_dir), out_file, reason, nsteps, mAhg)
        write(*, '(A,I0,A,F0.1,A)') 'exit '''//trim(reason)//''' after ', nsteps, ' steps, ', mAhg, &
            ' mAh/g; wrote '//trim(out_file)
    case default
        stop 'mode must be corrected, faithful_charge or faithful_phcell'
    end select
    close(u)
end program znmno2
