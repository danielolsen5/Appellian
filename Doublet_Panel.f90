program Doublet_Panel
    implicit none

    ! Constants
    real, parameter :: pi = 3.14159265358979323846
    integer, parameter :: N = 199 ! Number of panels
    
    ! Flow conditions
    real :: alpha, u_inf 
    real :: v_inf = 20.0
    real :: alpha_deg = 2.0
    
    ! Geometry arrays (N+1 nodes define N panels)
    real, dimension(N+1) :: x_node, y_node
    real, dimension(N)   :: x_col, y_col, nx, ny, tx, ty, len_p
    
    ! Singularity and Matrix arrays
    real, dimension(N)   :: sigma, mu, rhs, V_t, Cp
    real, dimension(N,N) :: B, C, A_matrix
    
    ! Loop and temporary variables
    integer :: i, j
    real :: x_loc, y_loc, r1, r2, theta1, theta2, dtheta
    real :: x_te, y_te, theta_w1, theta_w2, dtheta_w, C_wake
    real :: c_lw, mu_w, L_w, c_l
    real :: chord = 1.0
    real :: rho = 1.225

    integer :: unit_num = 10
    integer :: n_lines = 0
    integer :: ios 

    ! Reading Geometry Data
    character(len=100) :: naca_file = "naca2412.txt" 
    real, allocatable :: info(:,:)
    real :: temp_x, temp_y
    character(len=200) :: line

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

    ! Transfer info to nodes (Ensure n_lines matches N+1)
    if (n_lines < N + 1) then
        print *, "Error: Not enough nodes in file for N panels."
        stop
    end if

    do i = 1, N+1
        x_node(i) = info(i, 1)
        y_node(i) = info(i, 2)
    end do

    !------------------------------------------
    ! Calculations
    !------------------------------------------
    alpha = alpha_deg * (pi / 180.0)
    u_inf = V_inf * cos(alpha)
    v_inf = V_inf * sin(alpha)
    
    ! Calculate panel properties
    do i = 1, N
        ! Collocation points (midpoints)
        x_col(i) = 0.5 * (x_node(i) + x_node(i+1))
        y_col(i) = 0.5 * (y_node(i) + y_node(i+1))
        
        ! Tangent vector
        tx(i) = x_node(i+1) - x_node(i)
        ty(i) = y_node(i+1) - y_node(i)
        len_p(i) = sqrt(tx(i)**2 + ty(i)**2)
        
        ! Normalize tangents
        tx(i) = tx(i) / len_p(i)
        ty(i) = ty(i) / len_p(i)
        
        ! Outward Normal vector
        nx(i) = ty(i)
        ny(i) = -tx(i)
        
        ! Source Strengths
        sigma(i) = nx(i)*u_inf + ny(i)*v_inf
    end do
    
    ! Build Influence Matrices
    do i = 1, N
        do j = 1, N
            if (i == j) then
                C(i,j) = 0.5
                B(i,j) = (len_p(j) / (2.0 * pi)) * (log(len_p(j) / 2.0) - 1.0)
            else
                x_loc = (x_col(i) - x_node(j)) * tx(j) + (y_col(i) - y_node(j)) * ty(j)
                y_loc = (x_col(i) - x_node(j)) * nx(j) + (y_col(i) - y_node(j)) * ny(j)
                
                r1 = sqrt(x_loc**2 + y_loc**2)
                r2 = sqrt((x_loc - len_p(j))**2 + y_loc**2)
                
                theta1 = atan2(y_loc, x_loc)
                theta2 = atan2(y_loc, x_loc - len_p(j))
                
                dtheta = theta2 - theta1
                if (dtheta > pi) dtheta = dtheta - 2.0*pi
                if (dtheta < -pi) dtheta = dtheta + 2.0*pi
                
                C(i,j) = -(1.0 / (2.0 * pi)) * dtheta
                B(i,j) = (1.0 / (4.0 * pi)) * &
                         (x_loc * log(r1**2 / r2**2) + len_p(j) * log(r2**2) + &
                          2.0 * y_loc * dtheta - 2.0 * len_p(j))
            end if
        end do
    end do
    
    ! Kutta Condition - Wake influence
    x_te = x_node(N+1)
    y_te = y_node(N+1)
    
    do i = 1, N
        x_loc = (x_col(i) - x_te) * cos(alpha) + (y_col(i) - y_te) * sin(alpha)
        y_loc = -(x_col(i) - x_te) * sin(alpha) + (y_col(i) - y_te) * cos(alpha)
        
        theta_w1 = atan2(y_loc, x_loc)
        theta_w2 = atan2(y_loc, -10000.0) ! Representative "infinity"
        dtheta_w = theta_w2 - theta_w1
        
        C_wake = -(1.0 / (2.0 * pi)) * dtheta_w
        
        C(i,1) = C(i,1) - C_wake
        C(i,N) = C(i,N) + C_wake
    end do
    
    ! Solve System
    do i = 1, N
        rhs(i) = 0.0
        do j = 1, N
            rhs(i) = rhs(i) - B(i,j) * sigma(j)
        end do
    end do
    
    A_matrix = C
    call gauss_elimination(N, A_matrix, rhs, mu)
    
    ! ------------------------------------------
    ! Velocities and Cp (Fixed to include panels 1 and N)
    ! ------------------------------------------
    do i = 1, N
        if (i == 1) then
            ! Forward difference for panel 1
            r1 = sqrt((x_col(2) - x_col(1))**2 + (y_col(2) - y_col(1))**2)
            V_t(i) = (mu(2) - mu(1)) / r1 + (u_inf * tx(i) + v_inf * ty(i))
        else if (i == N) then
            ! Backward difference for panel N
            r1 = sqrt((x_col(N) - x_col(N-1))**2 + (y_col(N) - y_col(N-1))**2)
            V_t(i) = (mu(N) - mu(N-1)) / r1 + (u_inf * tx(i) + v_inf * ty(i))
        else
            ! Central difference for interior panels
            r1 = sqrt((x_col(i+1) - x_col(i-1))**2 + (y_col(i+1) - y_col(i-1))**2)
            V_t(i) = (mu(i+1) - mu(i-1)) / r1 + (u_inf * tx(i) + v_inf * ty(i))
        end if
        
        Cp(i) = 1.0 - (V_t(i) / V_inf)**2
    end do
    
    print *, "Doublet solution completed successfully."
    
    ! ------------------------------------------
    ! Lift Calc 1: Kutta-Joukowski (Wake Doublet)
    ! ------------------------------------------
    ! By Katz & Plotkin convention: Gamma = mu_upper - mu_lower. 
    ! Assuming counter-clockwise node ordering starting at TE-lower:
    mu_w = mu(1) - mu(N) 
    L_w = rho * V_inf * mu_w
    c_lw = (2.0 * L_w) / (rho * V_inf**2 * chord)
    print *, "Section lift coefficient from wake (Kutta-Joukowski): ", c_lw
    
    ! ------------------------------------------
    ! Lift Calc 2: Pressure Integration
    ! ------------------------------------------
    c_l = 0.0
    do i = 1, N
        ! Project pressure force normal vector onto the lift axis
        ! F_x = -Cp * nx * len_p,  F_y = -Cp * ny * len_p
        ! Lift = -F_x * sin(alpha) + F_y * cos(alpha)
        c_l = c_l + Cp(i) * len_p(i) * (nx(i) * sin(alpha) - ny(i) * cos(alpha)) / chord
    end do
    
    print *, "Section lift coefficient from Cp integration: ", c_l

contains

    subroutine lu_solve(N, A, b, x)
        ! Solves a general Ax=b on an nxn matrix
        ! This replaces A (in place) with its LU decomposition (permuted row-wise)
        implicit none
        
        integer, intent(in) :: N
        real, dimension(N,N), intent(inout) :: A
        real, dimension(N), intent(in) :: b
        ! Removed 'allocatable' below so it matches the fixed-size 'mu' array
        real, dimension(N), intent(out) :: x 
        
        integer, allocatable, dimension(:) :: indx
        integer :: D, info
        
        allocate(indx(N))
        
        ! Compute decomposition
        call lu_decomp(A, N, indx, D, info)
        
        ! if the matrix is nonsingular, then backsolve to find X
        if (info == 1) then
            write(*,*) 'Subroutine lu_decomp() failed. The given matrix is singular (i.e. no unique solution). Quitting...'
            stop
        else
            call lu_back_sub(A, N, indx, b, x)
        end if
        
        ! Cleanup
        deallocate(indx)
        
    end subroutine lu_solve

    subroutine lu_decomp(A, N, indx, D, code)
  ! Given an N x N matrix A, this routine replaces it by the LU
  ! decomposition of a rowwise permutation of itself. A and N  
  ! are input. indx is an output vector which records the row  
  ! permutation effected by the partial pivoting; D is output  
  ! as -1 or 1, depending on whether the number of row inter-  
  ! changes was even or odd, respectively. This routine is used
  ! in combination with lu_back_sub to solve linear equations or to 
  ! invert a matrix. Return code is 1 if matrix is singular.  

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
  ! Solves the set of N linear equations Ax = b.  Here A is     
  ! input, not as the matrix A but rather as its LU decomposition, 
  ! determined by the routine LUDCMP. indx is input as the permuta-
  ! tion vector returned by LUDCMP. b is input as the right-hand   
  ! side vector b. The solution vector is x. A, N, b and
  ! indx are not modified by this routine and can be used for suc- 
  ! cessive calls with different right-hand sides. This routine is 
  ! also efficient for plain matrix inversion.                     

  implicit none

  integer,intent(in) :: N
  real,dimension(N,N),intent(in) :: A
  real,dimension(N),intent(in) :: b
  integer,dimension(N),intent(in) :: indx
  real,dimension(:),allocatable,intent(out) :: x

  real :: sum
  integer :: ii,i,j,ll

  ! Initialize solution
  allocate(x, source=b)

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