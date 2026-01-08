import numpy as np
import matplotlib.pyplot as plt

"""
Joukowski Airfoil Generation
"""

# --- CONFIGURATION ---
Debug = True       # Toggle debug outputs
Plot = True         # Toggle plotting
Export = False      # Toggle export to file
Number = False      # Toggle printing point indices on plot
Curvature = False   # Toggle curvature calculation/plotting

# Clustering Settings
# 1: No clustering (Uniform)
# 2: Cosine clustering 
clustering_method = 2 

# Smoothing Factor (0.0 to 1.0)
smoothing_factor = .6

# Geometry Parameters
n = 200 # number of points
rho = 1.225 # density
R = 1.0       # Radius of cylinder
xi0 = -0.10  # Real center "x" (thickness)
eta0 = 0.01   # Imaginary center "y" (camber)
eps = R - np.sqrt(R**2 - eta0**2) - xi0 # Eccentricity (TE sharpness)
V_inf = 10    # Freestream velocity 
alpha = np.radians(0) # Angle of attack


# Geometric Angles for Leading Edge (LE) and Trailing Edge (TE)
xl = xi0 - np.sqrt(R**2 - eta0**2) 
xt = xi0 + np.sqrt(R**2 - eta0**2) 
theta_le = np.pi - np.arcsin(eta0/(xl - xi0)) 
theta_te = -np.arcsin(-eta0/(xi0 -  xt)) 

# Point Distribution Settings
l_u = np.abs(R * (xl-xt)) # length of upper surface in seta plane
l_t = np.pi * 2 * R # circumfrance in zeta
p_u = l_u/l_t
n_upper = p_u * n # scale number of points on the upper surface with its length
n_upper = round(n_upper) # force n_upper to be a whole number
n_lower = n - n_upper # Points on lower surface

if Debug:
        print(f"Target LE Angle: {np.degrees(theta_le):.2f}°")
        print(f"Target TE Angle: {np.degrees(theta_te):.2f}°")

# --- AERODYNAMIC CALCULATIONS ---
gamma_k = 4 * np.pi * V_inf * (np.sqrt(R**2 - eta0**2)*np.sin(alpha) + eta0*np.cos(alpha))

if abs(gamma_k) <= (4 * np.pi * V_inf * R):
    theta_stag_aft = alpha - np.arcsin(gamma_k / (4 * np.pi * V_inf * R))
    theta_stag_fwd = np.pi - theta_stag_aft + 2 * alpha
else:
    theta_stag_fwd = alpha - np.arcsin(gamma_k / (4 * np.pi * V_inf * R))
    theta_stag_aft = np.pi - theta_stag_fwd + 2 * alpha

# --- MAIN TRANSFORM FUNCTION ---
def transform_joukowski_cylinder(n_up, n_lo, R, xi0, eta0, eps, 
                                 method=2, smoothing=0.0):
    # Cylinder Center
    zeta_0 = xi0 + 1j * eta0
    
    # Geometric Angles for Leading Edge (LE) and Trailing Edge (TE)
    theta_le = np.pi - np.arcsin(eta0/(xl - xi0)) 
    theta_te = -np.arcsin(-eta0/(xi0 -  xt)) 
    
    
    # --- DISTRIBUTION GENERATOR ---
    def get_distribution(num_points, method_type, w_smooth):
        # 1. Base Uniform Vector (0 to pi)
        beta = np.linspace(0, np.pi, num_points + 1)
        
        # 2. Calculate Pure Distribution
        if method_type == 1: # UNIFORM
            d_pure = beta / np.pi
            
        elif method_type == 2: # COSINE
            d_pure = (1 - np.cos(beta)) / 2
            
        else: # Default to uniform if unknown
            print("Unknown method. Defaulting to Uniform.")
            d_pure = beta / np.pi
            
        # 3. Calculate Uniform Distribution (for blending)
        d_uniform = beta / np.pi
        
        # 4. Apply Smoothing (Blending)
        # w_smooth = 0.0 -> Pure Cosine
        # w_smooth = 1.0 -> Pure Uniform
        d_final = (1 - w_smooth) * d_pure + w_smooth * d_uniform
        
        return d_final

    # Calculate Upper and Lower distributions independently
    dist_upper = get_distribution(n_up, method, smoothing)
    dist_lower = get_distribution(n_lo, method, smoothing)

    # --- MAPPING ---
    # Segment 1: Upper Surface (TE -> LE)
    theta_upper = theta_te + (theta_le - theta_te) * dist_upper
    
    # Segment 2: Lower Surface (LE -> TE) 
    theta_lower = theta_le + (theta_te + 2*np.pi - theta_le) * dist_lower
    
    # Combine (Drop duplicate LE point from upper)
    theta = np.concatenate((theta_upper[:-1], theta_lower))
    
    # Generate Cylinder Surface
    zeta_surface = R * np.exp(1j * theta) + zeta_0
    
    # Joukowski Transform
    z_surface = zeta_surface + (R - eps)**2 / zeta_surface
    
    return zeta_surface, z_surface


# --- EXECUTION ---

# 1. Generate Raw Coordinates
zeta_surf, z_surf = transform_joukowski_cylinder(
    n_upper, n_lower, R, xi0, eta0, eps, 
    method=clustering_method, 
    smoothing=smoothing_factor
)

# 2. Extract and Normalize
z_real_raw = z_surf.real
z_imag_raw = z_surf.imag

# Shift LE to (0,0)
x_min = np.min(z_real_raw)
if x_min < 0:
    z_real_shifted = z_real_raw + abs(x_min)
else:
    z_real_shifted = z_real_raw - abs(x_min)

# Scale Chord to 1.0
chord_scale = 1.0 / np.max(z_real_shifted)
z_real = z_real_shifted[::-1] * chord_scale # Flip to clockwise
z_imag = z_imag_raw[::-1] * chord_scale     # Flip to clockwise

"""
VPM Solver
"""
x = z_real
y = z_imag
n = len(x)
n_pan = n - 1

#panel lengths 
dx = x[1:] - x[:-1]
dy = y[1:] - y[:-1]
l  = np.sqrt(dx**2 + dy**2)

# B
B = V_inf * ((dy * np.cos(a)) - (dx * np.sin(a))) / l
B = np.append(B, 0)

# control points
xc = 0.5 * (x[:-1] + x[1:])
yc = 0.5 * (y[:-1] + y[1:])
if Debug:
    print("xc:", xc)
    print("yc:", yc)

# build A 1.6.25
def build_A(x, y, xc, yc, dx, dy, l, DEBUG=False):
    n_pan = len(dx)
    A = np.zeros((n, n), dtype=float)
    eps = 1e-15

    for i in range(n-1):
        
        # outward unit normal at control point i 
        nx = dx[i] / l[i]
        ny = -dy[i] / l[i]

        for j in range(n-1):
            dxj = dx[j]; dyj = dy[j]; lj = l[j]

            #1.6.20  xi, eta in panel-j coordinates
            # use panel j to build the local frame 
            M = np.array([[dxj,  dyj],
                           [-dyj, dxj]])   # transform for panel j

            #1.6.20
            # vector from panel-j start to control point i
            N = np.array([xc[i] - x[j], yc[i] - y[j]])

            xi_eta = (1.0 / lj) * M.dot(N)   # [xi, eta]
            xi = float(xi_eta[0])
            eta = float(xi_eta[1])

            # phi 1.6.21
            phi = np.arctan2(eta * lj, (eta**2 + xi**2 - xi * lj))
            
            # psi 1.6.22
            num = xi**2 + eta**2
            den = (xi - lj)**2 + eta**2
            psi = 0.5 * np.log(num / (den ))

            # P 1.6.23
            W = np.array([
                [(lj - xi) * phi + (eta * psi),       xi * phi - eta * psi],
                [eta * phi - (lj - xi) * psi - lj,  -eta * phi - xi * psi + lj]
            ])
            
            V = np.array([[dxj,  -dyj],
                           [dyj, dxj]])

            # P 
            P = (1.0 / (2.0 * np.pi * lj**2)) * (np.matmul(V, W))

            # build A
            A[i, j] += nx * P[1, 0] + ny * P[0, 0]
            A[i, (j + 1) ] += nx * P[1, 1] + ny * P[0, 1]

            if DEBUG and (i, j) in ((0,0), (0,1), (0,2)):
                print("\nDEBUG: influence of panel j={} on cp i={}".format(j, i))
                print(" panel j start (x,y) = ({:.12e}, {:.12e})".format(x[j], y[j]))
                print(" panel j end   (x,y) = ({:.12e}, {:.12e})".format(x[j+1], y[j+1]))
                print(" cp (x,y) = ({:.12e}, {:.12e})".format(xc[i], yc[i]))
                print(" panel length l_j = {:.12e}".format(lj))
                print("    xi = {:.15e}".format(xi))
                print("   eta = {:.15e}".format(eta))
                print("   psi = {:.15e}".format(psi))
                print("   phi = {:.15e}".format(phi))
                print("   M =\n", M)
                print("   N =\n", N)
                print("   P =\n", P)
    return A

A = build_A(x, y, xc, yc, dx, dy, l, DEBUG=False)
A[-1,:] = 0.0
A[-1,0] = 1.0
A[-1,-1] = 1.0


gamma = np.linalg.solve(A, B)


"""
Lift calculation

"""

# Section lift 1.6.32
def calculate_lift_coefficient(l, gamma, c, V_inf):

    C_L = 0  # Initialize the lift coefficient 
    n = len(l)

    # Loop through each segment from i=0 to n-1
    for i in range(n):
        segment_contribution = (l[i] / c) * ((gamma[i] + gamma[i+1]) / V_inf)
        C_L += segment_contribution  # Add the segment's contribution to the total
    return C_L
c = 1
C_L = calculate_lift_coefficient(l, gamma, c, V_inf)
print("C_L,VPM: ", C_L)

#analytical solution

del_z = 4 * ((R**2 - eta0**2) / (np.sqrt(R**2 - eta0**2) - xi0))
L_kj = rho * V_inf * gamma_k # Kutta - Joukowski theorem Eq 8
#C_Lkj = L_kj / (rho * V_inf**2 * x_max * del_z) # non-dimentionalize section lift to get lift coefficient 
num = (2 * np.pi) * (np.sin(a) + ((eta0 * np.cos(a)) / np.sqrt(R**2 - eta0**2)))
den = 1 + (xi0 / (np.sqrt(R**2 - eta0**2) - xi0))
C_Lkj = num / den
c_vpm = np.max(x) + abs(np.min(x))

#print("del z: ", del_z)
#print("vpm chord: ", c_vpm)
print("Kutta Joukowski Lift: ", C_Lkj)
#print("gamma_K: ", gamma_k)

e = ((C_L - C_Lkj) / C_L) * 100
#print("percent error: ", e)


# Section lift 1.6.32
def calculate_lift_coefficient(l, gamma, c, V_inf):

    C_L = 0  # Initialize the lift coefficient 
    n = len(l)

    # Loop through each segment from i=0 to n-1
    for i in range(n):
        segment_contribution = (l[i] / c) * ((gamma[i] + gamma[i+1]) / V_inf)
        C_L += segment_contribution  # Add the segment's contribution to the total
    return C_L

# leading edge moment 
def calculate_moment(x, y, gamma, c, l, V_inf, a):
    C_mle = 0 # initialize value
    for i in range(n - 1):  # stop one short to avoid index overflow
        sx = 2*x[i]*gamma[i] + x[i]*gamma[i+1] + x[i+1]*gamma[i] + 2*x[i+1]*gamma[i+1]
        sy = 2*y[i]*gamma[i] + y[i]*gamma[i+1] + y[i+1]*gamma[i] + 2*y[i+1]*gamma[i+1]
        
        panel_contribution = (l[i] / c) * (
            (sx * np.cos(a) / (V_inf * c)) +
            (sy * np.sin(a) / (V_inf * c))
        )
        C_mle += (-1/3) * panel_contribution

    return C_mle

"""
Calculation printouts
"""
C_mle = calculate_moment(x, y, gamma, c, l, V_inf, a)
print(f"Leading-edge Moment Coefficient: {C_mle:.5f}")

C_qc = C_mle + (.25 * C_L)
print(f"Quarter Chord Moment Coefficient: {C_qc:.5f}")