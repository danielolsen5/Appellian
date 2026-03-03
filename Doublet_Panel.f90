program  Doublet_Panel
    implicit none
    
    ! read in airfoil data from json
    character(len=100) :: naca_file = "naca2412.txt" 
    real, allocatable :: info(:,:)
    integer :: i, ios, n_lines, unit_num
    real :: temp_x, temp_y
    character(len=200) :: line

        unit_num = 10

        ! Count input lines
    open(unit=unit_num, file=trim(naca_file), status='old', action='read')
        n_lines = 0
    do
        read(unit_num, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle ! Skip empty lines
        n_lines = n_lines + 1
    end do
    rewind(unit_num)

        ! read data
    allocate(info(n_lines, 2))
    
    i = 1
    do while (i <= n_lines)
        read(unit_num, *, iostat=ios) temp_x, temp_y
        if (ios == 0) then
                info(i, 1) = temp_x
                info(i, 2) = temp_y
                i = i + 1
        else
            exit
        end if
    end do
    close(unit_num)

    ! Verification
    print *, "Loaded ", n_lines, " points."
    print *, "First point: ", info(1, :)

    
end program  Doublet_Panel  