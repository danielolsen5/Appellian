program Doublet_Panel
    implicit none

    ! Constants
    real, parameter :: pi = 3.14159265358979323846
    
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
    real, allocatable, dimension(:,:) :: ep, A_mat, B_mat
    real, allocatable, dimension(:)   :: pt1_x, pt1_y, pt2_x, pt2_y
    real, allocatable, dimension(:)   :: th, co_x, co_y, b_vec, G, cp

    ! Temporary and Loop variables
    real :: dz, dx, xt, zt, x2t, z2t, x_rot, z_rot, x2_rot
    real :: r1, r2, ul, wl, u_vel, w_vel, r_wake, temp_val, r_dist, vloc, vel
    real :: inv_2pi

    inv_2pi = 1.0 / (2.0 * pi) 

    ! ------------------------------------------
    ! FILE READING 
    ! ------------------------------------------
    open(unit=unit_num, file=trim(naca_file), status='old', action='read')

    ! Count input lines
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

    ! Set matrix sizes dynamically based on the parsed coordinate file
    M = n_lines - 1   ! Number of panels
    N_sys = M + 1     ! System size (M panels + 1 wake panel)

    ! Allocate arrays
    allocate(ep(N_sys, 2))
    allocate(pt1_x(M), pt1_y(M), pt2_x(M), pt2_y(M))
    allocate(th(M), co_x(M), co_y(M), cp(M))
    allocate(A_mat(N_sys, N_sys), B_mat(N_sys, N_sys))
    allocate(b_vec(N_sys), G(N_sys))

    ! Convert angle of attack to radians
    al = alpha_deg * (pi / 180.0) ! Replaces AL=ALPHA/57.2958

    ! Convert paneling to clockwise 
    ! (Following Steven Yon's logic: EP(I,1) = EPT(N-I+1, 1))
    do i = 1, N_sys
        ep(i, 1) = info(n_lines - i + 1, 1)
        ep(i, 2) = info(n_lines - i + 1, 2)
    end do

    ! Establish coordinates of panel end points
    do i = 1, M
        pt1_x(i) = ep(i, 1)
        pt1_y(i) = ep(i, 2)
        pt2_x(i) = ep(i+1, 1)
        pt2_y(i) = ep(i+1, 2)
    end do

    ! Find panel angles TH(j) and collocation points CO(I)
    do i = 1, M
        dz = pt2_y(i) - pt1_y(i)
        dx = pt2_x(i) - pt1_x(i)
        th(i) = atan2(dz, dx)

        co_x(i) = (pt2_x(i) - pt1_x(i)) / 2.0 + pt1_x(i)
        co_y(i) = (pt2_y(i) - pt1_y(i)) / 2.0 + pt1_y(i)
    end do

    ! ------------------------------------------
    ! ESTABLISH INFLUENCE COEFFICIENTS
    ! ------------------------------------------
    A_mat = 0.0
    B_mat = 0.0

    do i = 1, M
        do j = 1, M
            ! Convert collocation point to local panel coordinates
            xt = co_x(i) - pt1_x(j)
            zt = co_y(i) - pt1_y(j)
            x2t = pt2_x(j) - pt1_x(j)
            z2t = pt2_y(j) - pt1_y(j)

            x_rot = xt*cos(th(j)) + zt*sin(th(j))
            z_rot = -xt*sin(th(j)) + zt*cos(th(j))
            x2_rot = x2t*cos(th(j)) + z2t*sin(th(j))
            ! z2_rot evaluates to 0

            r1 = sqrt(x_rot**2 + z_rot**2)
            r2 = sqrt((x_rot - x2_rot)**2 + z_rot**2)

            ! Compute velocity induced at the Ith collocation point by the Jth panel
            if (i == j) then
                ul = 0.0
                wl = -1.0 / (2.0 * pi * x_rot) 
            else
                ul = inv_2pi * (z_rot/(r1**2) - z_rot/(r2**2))
                wl = -inv_2pi * (x_rot/(r1**2) - (x_rot - x2_rot)/(r2**2))
            end if

            ! Rotate back to global frame
            u_vel = ul*cos(-th(j)) + wl*sin(-th(j))
            w_vel = -ul*sin(-th(j)) + wl*cos(-th(j))

            ! Components of velocity normal and tangential to panel I
            A_mat(i, j) = -u_vel*sin(th(i)) + w_vel*cos(th(i))
            B_mat(i, j) = u_vel*cos(th(i)) + w_vel*sin(th(i))
        end do

        ! Include the influence of the wake panel
        r_wake = sqrt((co_x(i) - pt2_x(M))**2 + (co_y(i) - pt2_y(M))**2)

        u_vel = inv_2pi * (co_y(i) / (r_wake**2))
        w_vel = -inv_2pi * (co_x(i) - pt2_x(M)) / (r_wake**2)

        A_mat(i, N_sys) = -u_vel*sin(th(i)) + w_vel*cos(th(i))
        B_mat(i, N_sys) = u_vel*cos(th(i)) + w_vel*sin(th(i))

        ! RHS formulation (replaces A(I, N+1) vector storage from screenshot)
        b_vec(i) = cos(al)*sin(th(i)) - sin(al)*cos(th(i))
    end do

    ! Prepare the matrix for solution by providing Kutta condition on the last row
    do i = 1, N_sys
        A_mat(N_sys, i) = 0.0
    end do
    A_mat(N_sys, 1) = -1.0
    A_mat(N_sys, M) = 1.0
    A_mat(N_sys, N_sys) = -1.0
    b_vec(N_sys) = 0.0

    ! Solve for the solution vector of doublet strengths
    call lu_solve(N_sys, A_mat, b_vec, G)

    ! ------------------------------------------
    ! VELOCITIES AND CP CALCULATIONS
    ! ------------------------------------------
    print *, "Collocation Point X", "    ", "Pressure Coefficient (Cp)"
    print *, "--------------------------------------------------------"
    do i = 1, M
        temp_val = 0.0
        do j = 1, N_sys
            temp_val = temp_val + B_mat(i, j) * G(j)
        end do

        if (i /= 1 .and. i /= M) then
            r_dist = sqrt((co_x(i+1) - co_x(i-1))**2 + (co_y(i+1) - co_y(i-1))**2)
            vloc = (G(i+1) - G(i-1)) / r_dist
        else if (i == 1) then
            r_dist = sqrt((co_x(2) - co_x(1))**2 + (co_y(2) - co_y(1))**2)
            vloc = (G(2) - G(1)) / r_dist
        else if (i == M) then
            r_dist = sqrt((co_x(M) - co_x(M-1))**2 + (co_y(M) - co_y(M-1))**2)
            vloc = (G(M) - G(M-1)) / r_dist
        end if

        vel = cos(al)*cos(th(i)) + sin(al)*sin(th(i)) + temp_val + vloc / 2.0
        cp(i) = 1.0 - vel**2
        
        write(*, '(F15.5, 5X, F15.5)') co_x(i), cp(i)
    end do

    print *, " "
    print *, "LIFT COEFFICIENT = ", G(N_sys)

    ! Cleanup
    deallocate(info, ep, pt1_x, pt1_y, pt2_x, pt2_y, th, co_x, co_y, cp, A_mat, B_mat, b_vec, G)

contains

    ! ------------------------------------------
    ! LU SOLVER
    ! ------------------------------------------
    subroutine lu_solve(N, A, b, x)
        implicit none
        
        integer, intent(in) :: N
        real, dimension(N,N), intent(inout) :: A
        real, dimension(N), intent(in) :: b
        real, dimension(N), intent(out) :: x 
        
        integer, allocatable, dimension(:) :: indx
        integer :: D, info
        
        allocate(indx(N))
        
        ! Compute decomposition
        call lu_decomp(A, N, indx, D, info)
        
        ! if the matrix is nonsingular, then backsolve to find X
        if (info == 1) then
            write(*,*) 'Subroutine lu_decomp() failed. The given matrix is singular. Quitting...'
            stop
        else
            call lu_back_sub(A, N, indx, b, x)
        end if
        
        ! Cleanup
        deallocate(indx)
        
    end subroutine lu_solve

    subroutine lu_decomp(A, N, indx, D, code)
      implicit none
    
      real,dimension(N,N),intent(inout) :: A
      integer,intent(in) :: N
      integer,dimension(N),intent(out) :: indx
      integer,intent(out) :: code, D
    
      real,dimension(N) :: vv
      real,parameter :: tiny=1.5e-20
      integer :: i, j, k, imax
      real :: amax, dum, s
    
      ! Initialize
      D = 1
      code = 0
      imax = 0
    
      ! Loop over rows to get implicit scaling information
      do i=1,N
        ! Get largest element in this row
        amax=0.0
        do j=1,N
          if (abs(A(i,j)) > amax) then
            amax = abs(A(i,j))
          end if
        end do
    
        ! Check the largest element in this row is nonzero
        if (amax <= tiny) then
          code = 1 ! Singular matrix
          return
        end if
    
        ! Store scaling
        vv(i) = 1.0 / amax
      end do
    
      ! Loop over columns of Crout's method
      do j=1,N
        do i=1,j-1
          s = A(i,j)
          do k=1,i-1
            s = s - A(i,k)*A(k,j)
          end do
          A(i,j) = s
        end do
    
        ! Initialize search for largest pivot element
        amax = 0.0
        do i=j,N
          s = A(i,j)
          do k=1,j-1
            s = s - A(i,k)*A(k,j)
          end do
          A(i,j) = s
      
          ! Determine figure of merit for the pivot
          dum = vv(i)*abs(s)
          if (dum >= amax) then
            imax = i
            amax = dum
          end if
        end do
    
        ! Figure out if we need to interchange rows
        if (j /= imax) then
          ! Perform interchange
          do k=1,N
            dum = A(imax,k)
            A(imax,k) = A(j,k)
            A(j,k) = dum
          end do
      
          ! Update the sign of D since a row interchange has occurred
          D = -D
      
          ! Interchange the implicit scaling factor
          vv(imax) = vv(j)
        end if
    
        ! Store pivoting
        indx(j) = imax
    
        ! Divide by pivot element
        if (j /= N) then
          dum = 1.0 / A(j,j)
          do i=j+1,N
            A(i,j) = A(i,j)*dum
          end do
        end if
      end do
    end subroutine lu_decomp

    subroutine lu_back_sub(A, N, indx, b, x)
      implicit none
    
      integer,intent(in) :: N
      real,dimension(N,N),intent(in) :: A
      real,dimension(N),intent(in) :: b
      integer,dimension(N),intent(in) :: indx
      ! Removed the 'allocatable' attribute here
      real,dimension(N),intent(out) :: x
    
      real :: sum
      integer :: ii,i,j,ll
    
      ! Standard array assignment instead of allocate(x, source=b)
      x = b
    
      ! Set tracker to ignore leading zeros in b
      ii = 0
    
      ! Forward substitution
      do i=1,N
        ! Untangle pivoting
        ll = indx(i)
        sum = x(ll)
        x(ll) = x(i)
    
        ! If a nonzero element of b has already been encountered
        if (ii /= 0) then
          do J=ii,i-1
            sum = sum - A(i,J)*x(J)
          end do
        ! Check for first nonzero element of b
        else if(sum /= 0.0) then
          ii = i
        end if
    
        x(i) = sum
      end do
    
      ! Back substitution
      do i=N,1,-1
        sum = x(i)
        do j=i+1,N
          sum = sum - A(i,j)*x(j)
        end do
        x(i) = sum / A(i,i)
      end do
    end subroutine lu_back_sub

end program Doublet_Panel