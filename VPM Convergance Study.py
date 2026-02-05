import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.integrate import quad

# --- CONFIGURATION ---
Debug = False        
Plot = True          
rho = 1.225          
R = 1.0              
xi0 = -0.10          
V_inf = 10           
alpha = np.radians(5)

# Iteration Settings
eta0_list = [0.0, 0.05, 0.10, 0.15] # Iterating through camber
n_iterations = 8   # How many times to double the points

# Storage for all plotting data
all_run_data = {}

# --- HELPER FUNCTIONS (UPDATED FROM VISUALIZATION SOURCE) ---

def get_point_at_theta(t, R, xi0, eta0, eps):
    """Helper to get a single (x, y) point for a given theta on the physical plane."""
    zeta = R * np.exp(1j * t) + (xi0 + 1j * eta0)
    z = zeta + (R - eps)**2 / zeta
    return z.real, z.imag

def arc_length_integrand(t, R, xi0, eta0, eps):
    """Calculates ds/dt for arc length integration."""
    dt = 1e-6
    z_m1_r, z_m1_i = get_point_at_theta(t - dt, R, xi0, eta0, eps)
    z_p1_r, z_p1_i = get_point_at_theta(t + dt, R, xi0, eta0, eps)
    
    dz_dt_r = (z_p1_r - z_m1_r) / (2 * dt)
    dz_dt_i = (z_p1_i - z_m1_i) / (2 * dt)
    
    return np.sqrt(dz_dt_r**2 + dz_dt_i**2)

def length_clustering(len_up, len_lo, n_total):
    """
    Calculates arc-length distributions (s) for upper and lower surfaces.
    """
    total_len = len_lo + len_up
    # Proportional point count
    n_up = int(round((len_up / total_len) * n_total))
    n_lo = n_total - n_up
    
    # Cosine spacing in Arc Length Domain
    i = np.arange(n_up + 1)
    s_up = (len_up / 2) * (1 - np.cos((i * np.pi) / n_up))
    
    j = np.arange(n_lo + 1)
    s_lo = (len_lo / 2) * (1 - np.cos((j * np.pi) / n_lo))
    
    return s_up, s_lo, n_up, n_lo

def generate_joukowski_from_s(s_up_arr, s_lo_arr, len_up, len_lo, 
                              R, xi0, eta0, eps, theta_le, theta_te):
    """
    Maps arc length (s) distributions back to physical coordinates (z).
    """
    # Linear approximation of s vs theta
    t_up = theta_te + (s_up_arr / len_up) * (theta_le - theta_te)
    t_lo = theta_le + (s_lo_arr / len_lo) * (theta_te + 2*np.pi - theta_le)
    
    all_t = np.concatenate([t_up, t_lo[1:]]) # Drop duplicate LE
    
    # Generate Z
    mu = xi0 + 1j * eta0
    zeta = mu + R * np.exp(1j * all_t)
    z = zeta + (R - eps)**2 / zeta 
    
    return z.real, z.imag

def rearrange_points(x, y):
    """
    Rearranges coordinates to ensure CW winding and start at Trailing Edge.
    """
    # 1. Check Winding (Shoelace formula)
    # Positive area = CCW (Top First), Negative = CW (Bottom First)
    area = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
    
    # We want CW (Negative Area) for VPM usually, but let's just ensure consistent ordering
    # If the input generates CCW, we flip to CW.
    if area > 0:
        x = x[::-1]
        y = y[::-1]

    # 2. Find Trailing Edge (Max X)
    idx_te = np.argmax(x[:-1])
    
    # 3. Roll arrays to start at TE
    if idx_te != 0:
        x_unique = np.roll(x[:-1], -idx_te)
        y_unique = np.roll(y[:-1], -idx_te)
        
        # Close the loop
        x = np.append(x_unique, x_unique[0])
        y = np.append(y_unique, y_unique[0])
        
    return x, y

# --- VPM SOLVER FUNCTIONS (UNCHANGED) ---

def build_A(x, y, xc, yc, dx, dy, l, n, DEBUG=False):
    A = np.zeros((n, n), dtype=float)
    
    for i in range(n-1):
        nx = dx[i] / l[i]
        ny = -dy[i] / l[i]

        for j in range(n-1):
            dxj = dx[j]; dyj = dy[j]; lj = l[j]
            M = np.array([[dxj,  dyj], [-dyj, dxj]]) 
            N = np.array([xc[i] - x[j], yc[i] - y[j]])
            xi_eta = (1.0 / lj) * M.dot(N) 
            xi = float(xi_eta[0])
            eta = float(xi_eta[1])
            phi = np.arctan2(eta * lj, (eta**2 + xi**2 - xi * lj))
            num = xi**2 + eta**2
            den = (xi - lj)**2 + eta**2
            psi = 0.5 * np.log(num / (den if den > 1e-12 else 1e-12)) 

            W = np.array([
                [(lj - xi) * phi + (eta * psi),       xi * phi - eta * psi],
                [eta * phi - (lj - xi) * psi - lj,  -eta * phi - xi * psi + lj]
            ])
            V = np.array([[dxj,  -dyj], [dyj, dxj]])
            P = (1.0 / (2.0 * np.pi * lj**2)) * (np.matmul(V, W))

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

# --- MAIN CONVERGENCE LOOPS ---

print(f"{'Camber (eta0)':<15} | {'Points (N)':<15} | {'CL (VPM)':<15} | {'Error (%)':<15}")
print("-" * 65)

# OUTER LOOP: Iterate through Camber values
for eta0 in eta0_list:
    
    # Storage for this specific eta0 run
    n_values_list = []
    error_values_list = []
    
    # 1. Recalculate Geometry Constants
    # Calculate epsilon automatically based on geometry (Updated from Source 1)
    eps = R - np.sqrt(R**2 - eta0**2) - xi0
    
    xl = xi0 - np.sqrt(R**2 - eta0**2) 
    xt = xi0 + np.sqrt(R**2 - eta0**2) 
    theta_le = np.pi - np.arcsin(eta0/(xl - xi0)) 
    theta_te = -np.arcsin(-eta0/(xi0 - xt)) 

    # 2. Calculate Exact Arc Lengths (Z-Plane) for Clustering
    true_len_upper, _ = quad(arc_length_integrand, theta_te, theta_le, args=(R, xi0, eta0, eps))
    true_len_lower, _ = quad(arc_length_integrand, theta_le, theta_te + 2*np.pi, args=(R, xi0, eta0, eps))

    # INNER LOOP: Iterate through N points
    for i in range(n_iterations):
        
        # 1. Update Number of Points
        n = 20 * (2**i) 
        
        # 2. Generate Clustering Distribution (Updated)
        s_upper, s_lower, _, _ = length_clustering(true_len_upper, true_len_lower, n)

        # 3. Generate Coordinates from s (Updated)
        z_real_final, z_imag_final = generate_joukowski_from_s(
            s_upper, s_lower, true_len_upper, true_len_lower,
            R, xi0, eta0, eps, theta_le, theta_te
        )
        
        # 4. Normalize Geometry
        # Shift so min x is at 0, then scale by chord
        z_real_shifted = z_real_final - np.min(z_real_final)
        chord_scale = 1.0 / np.max(z_real_shifted)
        
        x_raw = z_real_shifted * chord_scale
        y_raw = z_imag_final * chord_scale

        # 5. Rearrange Points
        x, y = rearrange_points(x_raw, y_raw)
        
        # --- VPM SOLVER START ---
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
        
        # Solve
        try:
            gamma = np.linalg.solve(A, B)
            
            # Post-Process (Lift)
            c = 1.0
            C_L = calculate_lift_coefficient(l, gamma, c, V_inf)
            
            # Analytical Solution (Kutta-Joukowski)
            num = (2 * np.pi) * (np.sin(alpha) + ((eta0 * np.cos(alpha)) / np.sqrt(R**2 - eta0**2)))
            den = 1 + (xi0 / (np.sqrt(R**2 - eta0**2) - xi0))
            C_Lkj = num / den
            
            # Calculate Error
            error_pct = np.abs(((C_L - C_Lkj) / C_Lkj) * 100)
            
            # Store
            n_values_list.append(n)
            error_values_list.append(error_pct)
            
            print(f"{eta0:<15.2f} | {n:<15} | {C_L:<15.5f} | {error_pct:<15.5f}")

        except np.linalg.LinAlgError:
            print(f"{eta0:<15.2f} | {n:<15} | Singular Matrix Error")
            
    # Save data for this eta0
    all_run_data[eta0] = (n_values_list, error_values_list)

# --- PLOTTING ---
if Plot:
    plt.figure(figsize=(10, 6))
    
    # Iterate through stored results and plot
    for eta0_val, (n_vals, err_vals) in all_run_data.items():
        plt.loglog(n_vals, err_vals, 'o-', linewidth=2, label=f'eta0 = {eta0_val}')
        
    plt.title(f"VPM Convergence")
    plt.xlabel("Number of Panels")
    plt.ylabel("Absolute Error in CL (%)")
    plt.legend()
    plt.show()