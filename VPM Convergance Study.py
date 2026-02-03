import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

# --- CONFIGURATION ---
Debug = False       
Plot = True         
rho = 1.225         
R = 1.0             
xi0 = -0.10         
eta0 = 0.0         
V_inf = 10          
alpha = np.radians(5)
eps = R - np.sqrt(R**2 - eta0**2) - xi0

# Convergence Settings
n_iterations = 7  # How many times to double the points

# --- HELPER FUNCTIONS ---

def transform_joukowski_cylinder(n_up, n_lo, R, xi0, eta0, eps, 
                                 xl, xt, theta_le, theta_te,
                                 method=2, smoothing=0.0):
    # Cylinder Center
    zeta_0 = xi0 + 1j * eta0
    
    # --- DISTRIBUTION GENERATOR ---
    def get_distribution(num_points, method_type, w_smooth):
        beta = np.linspace(0, np.pi, num_points + 1)
        
        if method_type == 1: # UNIFORM
            d_pure = beta / np.pi
        elif method_type == 2: # COSINE
            d_pure = (1 - np.cos(beta)) / 2
        else:
            d_pure = beta / np.pi
            
        d_uniform = beta / np.pi
        d_final = (1 - w_smooth) * d_pure + w_smooth * d_uniform
        return d_final

    # Calculate Upper and Lower distributions
    dist_upper = get_distribution(n_up, method, smoothing)
    dist_lower = get_distribution(n_lo, method, smoothing)

    # Mapping
    theta_upper = theta_te + (theta_le - theta_te) * dist_upper
    theta_lower = theta_le + (theta_te + 2*np.pi - theta_le) * dist_lower
    
    # Combine (Drop duplicate LE point from upper)
    theta = np.concatenate((theta_upper[:-1], theta_lower))
    
    # Generate Cylinder Surface
    zeta_surface = R * np.exp(1j * theta) + zeta_0
    
    # Joukowski Transform
    z_surface = zeta_surface + (R - eps)**2 / zeta_surface
    
    return zeta_surface, z_surface

def build_A(x, y, xc, yc, dx, dy, l, n, DEBUG=False):
    A = np.zeros((n, n), dtype=float)
    
    for i in range(n-1):
        # outward unit normal at control point i 
        nx = dx[i] / l[i]
        ny = -dy[i] / l[i]

        for j in range(n-1):
            dxj = dx[j]; dyj = dy[j]; lj = l[j]

            # 1.6.20 xi, eta in panel-j coordinates
            M = np.array([[dxj,  dyj], [-dyj, dxj]]) 
            N = np.array([xc[i] - x[j], yc[i] - y[j]])

            xi_eta = (1.0 / lj) * M.dot(N) 
            xi = float(xi_eta[0])
            eta = float(xi_eta[1])

            # phi 1.6.21
            phi = np.arctan2(eta * lj, (eta**2 + xi**2 - xi * lj))
            
            # psi 1.6.22
            num = xi**2 + eta**2
            den = (xi - lj)**2 + eta**2
            psi = 0.5 * np.log(num / (den if den > 1e-12 else 1e-12)) # avoid div 0

            # P 1.6.23
            W = np.array([
                [(lj - xi) * phi + (eta * psi),       xi * phi - eta * psi],
                [eta * phi - (lj - xi) * psi - lj,  -eta * phi - xi * psi + lj]
            ])
            
            V = np.array([[dxj,  -dyj], [dyj, dxj]])

            P = (1.0 / (2.0 * np.pi * lj**2)) * (np.matmul(V, W))

            # build A
            A[i, j] += nx * P[1, 0] + ny * P[0, 0]
            A[i, (j + 1) ] += nx * P[1, 1] + ny * P[0, 1]
            
    return A

def calculate_lift_coefficient(l, gamma, c, V_inf):
    C_L = 0 
    n = len(l)
    for i in range(n):
        segment_contribution = (l[i] / c) * ((gamma[i] + gamma[i+1]) / V_inf)
        C_L += segment_contribution 
    return C_L

# --- MAIN CONVERGENCE LOOP ---

# Storage for plotting
n_values_list = []
error_values_list = []

# Geometric Angles (Constant for all N)
xl = xi0 - np.sqrt(R**2 - eta0**2) 
xt = xi0 + np.sqrt(R**2 - eta0**2) 
theta_le = np.pi - np.arcsin(eta0/(xl - xi0)) 
theta_te = -np.arcsin(-eta0/(xi0 -  xt)) 

print(f"{'Points (N)':<15} | {'CL (VPM)':<15} | {'Error (%)':<15}")
print("-" * 50)

for i in range(n_iterations):
    
    # 1. Update Number of Points
    n = 20 * (2**i) # Starts at 20, then 40, 80...
    
    # 2. Recalculate Upper/Lower split based on N
    l_u = R * np.abs((theta_le-theta_te)) 
    l_t = np.pi * 2 * R 
    p_u = l_u/l_t
    n_upper = round(p_u * n) 
    n_lower = n - n_upper 
    
    # 3. Generate Geometry
    zeta_surf, z_surf = transform_joukowski_cylinder(
        n_upper, n_lower, R, xi0, eta0, eps,
        xl, xt, theta_le, theta_te,
        method=2, smoothing=0.6 
    )
    
    # Normalize Geometry
    z_real_raw = z_surf.real
    z_imag_raw = z_surf.imag
    x_min = np.min(z_real_raw)
    z_real_shifted = z_real_raw - x_min # Shift to 0
    chord_scale = 1.0 / np.max(z_real_shifted)
    
    # Coordinates for VPM
    x = z_real_shifted[::-1] * chord_scale # Flip to clockwise
    y = z_imag_raw[::-1] * chord_scale     # Flip to clockwise
    
    # 4. VPM Setup
    n_nodes = len(x)
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    l = np.sqrt(dx**2 + dy**2)
    
    # Control points
    xc = 0.5 * (x[:-1] + x[1:])
    yc = 0.5 * (y[:-1] + y[1:])
    
    # RHS Vector B
    B = V_inf * ((dy * np.cos(alpha)) - (dx * np.sin(alpha))) / l
    B = np.append(B, 0) # Kutta condition RHS
    
    # LHS Matrix A
    A = build_A(x, y, xc, yc, dx, dy, l, n_nodes, DEBUG=Debug)
    
    # Apply Kutta Condition at the end of A
    A[-1,:] = 0.0
    A[-1,0] = 1.0
    A[-1,-1] = 1.0
    
    # 5. Solve
    try:
        gamma = np.linalg.solve(A, B)
        
        # 6. Post-Process (Lift)
        c = 1.0
        C_L = calculate_lift_coefficient(l, gamma, c, V_inf)
        
        # Analytical Solution (Kutta-Joukowski)
        # Note: Recalculating exact analytical CL for reference
        num = (2 * np.pi) * (np.sin(alpha) + ((eta0 * np.cos(alpha)) / np.sqrt(R**2 - eta0**2)))
        den = 1 + (xi0 / (np.sqrt(R**2 - eta0**2) - xi0))
        C_Lkj = num / den
        
        # Calculate Error
        error_pct = np.abs(((C_L - C_Lkj) / C_Lkj) * 100)
        
        # Store
        n_values_list.append(n)
        error_values_list.append(error_pct)
        
        print(f"{n:<15} | {C_L:<15.5f} | {error_pct:<15.5f}")

    except np.linalg.LinAlgError:
        print(f"{n:<15} | Singular Matrix Error")

# --- PLOTTING ---
if Plot:
    plt.figure(figsize=(10, 6))
    plt.loglog(n_values_list, error_values_list, 'o-', linewidth=2, color='b')
    plt.title(f"VPM Convergence (Alpha={np.degrees(alpha):.1f}°)")
    plt.xlabel("Number of Panels (N)")
    plt.ylabel("Absolute Error in CL (%)")
    plt.show()