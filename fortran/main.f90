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
    call check_groups(u)
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
contains

    !> Every &group must be one of the program's, and appear once (otherwise a misspelled group would be skipped
    !> silently, and only the first of two groups with the same name would be read).
    subroutine check_groups(u)
        integer, intent(in) :: u
        character(len=12), parameter :: known(11) = [character(len=12) :: 'run', 'faithful', 'cell', 'solids', &
            'electrolyte', 'reactions', 'options', 'constants', 'protocol', 'numerics', 'output']
        logical :: seen(11)
        character(len=2048) :: line
        character(len=64) :: name
        character(len=1) :: quote
        integer :: ios, i, j, k, g
        seen = .false.
        rewind(u)
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            quote = ' '
            i = 1
            do while (i <= len_trim(line))
                if (quote /= ' ') then
                    if (line(i:i) == quote) quote = ' '
                else if (line(i:i) == '''' .or. line(i:i) == '"') then
                    quote = line(i:i)
                else if (line(i:i) == '!') then
                    exit
                else if (line(i:i) == '&') then
                    j = i + 1
                    do while (j <= len_trim(line))
                        if (index('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_', line(j:j)) == 0) exit
                        j = j + 1
                    end do
                    name = line(i+1:j-1)
                    do k = 1, len_trim(name)
                        if (name(k:k) >= 'A' .and. name(k:k) <= 'Z') name(k:k) = achar(iachar(name(k:k)) + 32)
                    end do
                    g = findloc(known, name, 1)
                    if (g == 0) then
                        write(*, '(A)') 'unknown namelist group &'//trim(name)
                        stop 1
                    end if
                    if (seen(g)) then
                        write(*, '(A)') 'namelist group &'//trim(name)//' appears twice'
                        stop 1
                    end if
                    seen(g) = .true.
                    i = j - 1
                end if
                i = i + 1
            end do
        end do
        rewind(u)
    end subroutine check_groups

end program znmno2
