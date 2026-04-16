program Doublet_Panel
    implicit none

    ! Constants
    real, parameter :: pi = 3.14159265358979323846
    real, parameter :: inv_2pi = 0.159154943
    real, parameter :: inv_pi = 0.318309886

    ! Flow conditions
    real :: alpha_deg = 1.0 
    real :: al

    ! File reading variables
    integer :: unit_num = 10
    integer :: n_lines = 0
    integer :: ios, i, j, M, N_sys
    character(len=100) :: naca_file = "naca2412.txt" 
    real, allocatable :: info(:,:)
    real :: temp_x, temp_y
    character(len=200) :: line

    ! Geometry and Panel arrays 
    real, allocatable, dimension(:,:) :: ep, A_mat
    real, allocatable, dimension(:)   :: pt1_x, pt1_y, pt2_x, pt2_y
    real, allocatable, dimension(:)   :: th, co_x, co_y, b_vec, G, cp, sig, dl, phi

    ! Temporary variables
    real :: xt, zt, x2t, z2t, x_rot, z_rot, x2_rot, r1, r2, th1, th2
    real :: xw, zw, dthw, temp_rhs, r_dist, vel

    ! ------------------------------------------
    ! FILE READING
    ! ------------------------------------------
    open(unit=unit_num, file=trim(naca_file), status='old', action='read')
    do
        read(unit_num, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle 
        n_lines = n_lines + 1
    end do
    rewind(unit_num)
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
    
    ! ------------------------------------------
    ! SETUP & ALLOCATION
    ! ------------------------------------------
    M = n_lines - 1   
    N_sys = M + 1

    allocate(ep(n_lines, 2), pt1_x(M), pt1_y(M), pt2_x(M), pt2_y(M))
    allocate(th(M), co_x(M), co_y(M), cp(M), sig(M), dl(M), phi(M))
    allocate(A_mat(N_sys, N_sys), b_vec(N_sys), G(N_sys))

    al = alpha_deg * (pi / 180.0)

    do i = 1, n_lines
        ep(i, 1) = info(n_lines - i + 1, 1)
        ep(i, 2) = info(n_lines - i + 1, 2)
    end do

    do i = 1, M
        pt1_x(i) = ep(i, 1)
        pt1_y(i) = ep(i, 2)
        pt2_x(i) = ep(i+1, 1)
        pt2_y(i) = ep(i+1, 2)
        
        th(i) = atan2(pt2_y(i) - pt1_y(i), pt2_x(i) - pt1_x(i))
        co_x(i) = (pt1_x(i) + pt2_x(i)) / 2.0
        co_y(i) = (pt1_y(i) + pt2_y(i)) / 2.0
        
        sig(i) = cos(al)*sin(th(i)) - sin(al)*cos(th(i))
    end do

    ! ------------------------------------------
    ! ESTABLISH INFLUENCE COEFFICIENTS
    ! ------------------------------------------
    A_mat = 0.0
    b_vec = 0.0

    do i = 1, M
        temp_rhs = 0.0
        do j = 1, M
            xt = co_x(i) - pt1_x(j)
            zt = co_y(i) - pt1_y(j)
            x2t = pt2_x(j) - pt1_x(j)
            z2t = pt2_y(j) - pt1_y(j)

            x_rot = xt*cos(th(j)) + zt*sin(th(j))
            z_rot = -xt*sin(th(j)) + zt*cos(th(j))
            x2_rot = x2t*cos(th(j)) + z2t*sin(th(j))
            
            if (i == 1) dl(j) = x2_rot

            r1 = sqrt(x_rot**2 + z_rot**2)
            r2 = sqrt((x_rot - x2_rot)**2 + z_rot**2)
            th1 = atan2(z_rot, x_rot)
            th2 = atan2(z_rot, x_rot - x2_rot)

            if (i == j) then
                A_mat(i, j) = 0.5
                temp_rhs = temp_rhs + sig(j) * inv_pi * (x_rot * log(r1))
            else
                A_mat(i, j) = -inv_2pi * (th2 - th1)
                temp_rhs = temp_rhs + (sig(j) / (2.0*pi)) * &
                           (x_rot*log(r1) - (x_rot-x2_rot)*log(r2) + z_rot*(th2-th1))
            end if
        end do

        xw = co_x(i) - pt2_x(M)
        zw = co_y(i) - pt2_y(M)
        dthw = -atan2(zw, xw)
        A_mat(i, N_sys) = -inv_2pi * dthw
        
        b_vec(i) = temp_rhs
    end do

    ! Explicit Kutta Condition
    A_mat(N_sys, :) = 0.0
    A_mat(N_sys, 1) = -1.0
    A_mat(N_sys, M) = 1.0
    A_mat(N_sys, N_sys) = -1.0
    b_vec(N_sys) = 0.0

    call lu_solve(N_sys, A_mat, b_vec, G)

    ! ------------------------------------------
    ! POST-PROCESSING
    ! ------------------------------------------
    do i = 1, M
        phi(i) = co_x(i)*cos(al) + co_y(i)*sin(al) + G(i)
    end do

    print *, "Collocation Point X", "    ", "Pressure Coefficient (Cp)"
    do i = 1, M - 1
        r_dist = (dl(i+1) + dl(i)) / 2.0
        vel = (phi(i) - phi(i+1)) / r_dist
        cp(i) = 1.0 - vel**2
        write(*, '(F15.5, 5X, F15.5)') co_x(i), cp(i)
    end do

    print *, " "
    print *, "LIFT COEFFICIENT = ", G(N_sys)

    deallocate(info, ep, pt1_x, pt1_y, pt2_x, pt2_y, th, co_x, co_y, cp, sig, dl, phi, A_mat, b_vec, G)

contains

    subroutine lu_solve(N, A, b, x)
        integer, intent(in) :: N
        real, dimension(N,N), intent(inout) :: A
        real, dimension(N), intent(in) :: b
        real, dimension(N), intent(out) :: x 
        integer, allocatable, dimension(:) :: indx
        integer :: D, code
        
        allocate(indx(N))
        call lu_decomp(A, N, indx, D, code)
        if (code == 1) then
            write(*,*) 'Matrix is singular.'
            stop
        else
            call lu_back_sub(A, N, indx, b, x)
        end if
        deallocate(indx)
    end subroutine lu_solve

    subroutine lu_decomp(A, N, indx, D, code)
      real,dimension(N,N),intent(inout) :: A
      integer,intent(in) :: N
      integer,dimension(N),intent(out) :: indx
      integer,intent(out) :: code, D
      real,dimension(N) :: vv
      real,parameter :: tiny=1.5e-20
      integer :: i, j, k, imax
      real :: amax, dum, s
    
      D = 1
      code = 0
      do i=1,N
        amax=0.0
        do j=1,N
          if (abs(A(i,j)) > amax) amax = abs(A(i,j))
        end do
        if (amax <= tiny) then
          code = 1 
          return
        end if
        vv(i) = 1.0 / amax
      end do
    
      do j=1,N
        do i=1,j-1
          s = A(i,j)
          do k=1,i-1
            s = s - A(i,k)*A(k,j)
          end do
          A(i,j) = s
        end do
        amax = 0.0
        do i=j,N
          s = A(i,j)
          do k=1,j-1
            s = s - A(i,k)*A(k,j)
          end do
          A(i,j) = s
          dum = vv(i)*abs(s)
          if (dum >= amax) then
            imax = i
            amax = dum
          end if
        end do
        if (j /= imax) then
          do k=1,N
            dum = A(imax,k)
            A(imax,k) = A(j,k)
            A(j,k) = dum
          end do
          D = -D
          vv(imax) = vv(j)
        end if
        indx(j) = imax
        if (j /= N) then
          dum = 1.0 / A(j,j)
          do i=j+1,N
            A(i,j) = A(i,j)*dum
          end do
        end if
      end do
    end subroutine lu_decomp

    subroutine lu_back_sub(A, N, indx, b, x)
      integer,intent(in) :: N
      real,dimension(N,N),intent(in) :: A
      real,dimension(N),intent(in) :: b
      integer,dimension(N),intent(in) :: indx
      real,dimension(N),intent(out) :: x
      real :: sum
      integer :: ii,i,j,ll
    
      x = b
      ii = 0
      do i=1,N
        ll = indx(i)
        sum = x(ll)
        x(ll) = x(i)
        if (ii /= 0) then
          do J=ii,i-1
            sum = sum - A(i,J)*x(J)
          end do
        else if(sum /= 0.0) then
          ii = i
        end if
        x(i) = sum
      end do
      do i=N,1,-1
        sum = x(i)
        do j=i+1,N
          sum = sum - A(i,j)*x(j)
        end do
        x(i) = sum / A(i,i)
      end do
    end subroutine lu_back_sub

end program Doublet_Panel