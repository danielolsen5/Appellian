import numpy as np
import json
import matplotlib.pyplot as plt


#--------------------------------
# Global variables
#--------------------------------

v_inf = 1          # freestream velocity
alpha = 0           # AOA in degrees
a = np.radians(alpha)
rho = 1.225         # air density
c = 1               # chord length
T = 2               # simulaion time
n_steps = 10       # number of time steps
dt = T/n_steps      # time step 

#-------------------------------
#Read JSON geometry
#-------------------------------

filename = "json_input.json"
with open(filename, "r") as fh:
    input_dict = json.load(fh)

airfoil = input_dict["geometry"]
info = []
with open(airfoil, "r") as f:
    for line in f:
        if line.strip() == "":
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        info.append([float(parts[0]), float(parts[1])])

info = np.array(info)
x = info [:, 0]
y = info [:, 1]

n = len(x) # number of geometry points
n_pan = n-1 # number of panels

#-------------------------------------
# VPM formulation
#-------------------------------------

#calculate panel lengths
dx = x[1:] - x[:-1]
dy = y[1:] - y[:-1]
l = np.sqrt(dx**2 +dy**2)

# find control point coordinates
xc = 0.5 * (x[1:] - x[:-1])
yc = 0.5 * (y[1:] - y[:-1])

#-------------------------------
# A matrix
#-------------------------------

def build_A(x, y, xc, yx, dx, dy, l, n_pan):
    A = np.zeros ((n, n), dtype=float)
    for i in range(n-1):

        #outward unit normal vector at each control point
        nx = dx[i]/l[i]
        ny = dy[i]/l[i]
        for j in range(n-1):
            dxj = dx[j]; dyj = dy[j]; lj=l[j] # define lengths for j index

            #define local frame from panel j
            M = np.array([[dxj, dyj],
                          [-dyj, dxj]])
            
            # vector from panel j to start of control point i
            N = np.array([xc[i]-x[j], yc[i] - y[j]])

            xi_eta = (1 / lj) * M.dot(N)
            xi = float(xi_eta[0])
            eta = float(xi_eta[1])

            phi = np.arctan2(eta * lj, (eta**2 + xi**2 -xi * lj))

            psi_num = xi**2 + eta**2
            psi_den = (xi - lj)**2 + eta**2
            psi = 0.5 * np.log(psi_num / psi_den)

            W = np.array([
                         [(lj - xi) * phi + (eta * psi),      xi * phi - eta * psi],
                         [eta * phi - (lj - xi) * psi - lj,   -eta * phi - xi * psi + lj]
            ])

            V = np.array([
                         [dxj, -dyj],
                         [dyj, dxj]
            ])

            P = (1.0 / (2.0 * np.pi * lj**2)) * np.matmul(V, W)

            # Build A 
            A[i, j] += nx * P[1, 0] + ny * P[0, 0]
            A[i, (j+1)] += nx * P[1, 1] + ny * P[0, 1]


    return A

A = build_A(x, y, xc, yc, dx, dy, l, n_pan)

# Apply Kutta Condition
A[-1,:] = 0
A[-1,0] = 1
A[-1, -1] = 1

#------------------------
# B vector
#------------------------

def build_B(v_inf, dy, dx, a):
    B = v_inf * ((dy*np.cos(a)) - (dx * np.sin(a)))
    B = np.append(B, 0)
    return B

B = build_B(v_inf, dy, dx, a)


#----------------------
# circulation vector
#----------------------

gamma = np.linalg.solve(A, B)


#-------------------------
# Lift coefficient
#-------------------------

def calc_lift(l, gamma, c, v_inf):
    C_L = 0
    n = len(l)
    for i in range(n):
        contribution = (l[i] / c) * ((gamma[i] + gamma[i+1]) / v_inf)
        C_L += contribution
    return C_L


#-------------------------
# Leading edge moment 
#-------------------------
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

#-------------------------
# Time Stepping
#-------------------------

# wake properties
wake_x = []       # X-positions of shed vortices
wake_y = []       # Y-positions of shed vortices
wake_gamma = []  # Strengths of shed vortices

# Pre-build A matrix 
A_matrix = build_A(x, y, xc, yc, dx, dy, l, n_pan)

# Apply Kutta Condition to the static matrix
A_matrix[-1, :] = 0
A_matrix[-1, 0] = 1
A_matrix[-1, -1] = 1

# Record each C_l
cl_history = []

#-------------------------
# Time Stepping Loop
#-------------------------
wake_x, wake_y, wake_gamma = [], [], []
cl_history, time_history = [], []

def get_total_cl(l, gamma, v_inf, c):
    Gamma = np.sum(0.5 * (gamma[:-1] + gamma[1:]) * l)
    return 2 * Gamma / (v_inf * c)

def build_B_fs(v_inf, dy, dx, a, l):
    nx, ny = -dy/l, dx/l
    u_inf, v_inf_y = v_inf * np.cos(a), v_inf * np.sin(a)
    return -(u_inf * nx + v_inf_y * ny)

# Initial steady solve
B_fs = build_B_fs(v_inf, dy, dx, a, l)
gamma = np.linalg.solve(A_matrix, np.append(B_fs, 0))
gamma_airfoil_prev = np.sum(0.5 * (gamma[:-1] + gamma[1:]) * l)

def get_freestream_conditions(t, v_base, alpha_base):
    # Time-varying parameters
    # Gradually increase velocity and AOA over the first half of simulation
    ramp_end = T / 2 
    
    if t < ramp_end:
        # Linear ramp for velocity and AOA
        v_t = v_base + (2 * (t / ramp_end))  # Increase by 10 m/s
        alpha_deg_t = alpha_base + (10 * (t / ramp_end))  # Increase by 5 degrees
    else:
        # Hold constant after ramp
        v_t = v_base + 10
        alpha_deg_t = alpha_base + 5
        
    return v_t, np.radians(alpha_deg_t)

for step in range(n_steps):
    t = step * dt
    
    # 1. Update time-varying freestream conditions
    v_inf_t, a_t = get_freestream_conditions(t, v_inf, alpha)
    
    # 2. Re-calculate the freestream RHS (Quasi-steady contribution)
    # This represents the lift produced if conditions were maintained [cite: 242]
    B_fs = v_inf_t * ((dy * np.cos(a_t)) - (dx * np.sin(a_t)))
    
    # 3. Convect existing wake at current v_inf_t
    for i in range(len(wake_x)):
        wake_x[i] += v_inf_t * np.cos(a_t) * dt
        wake_y[i] += v_inf_t * np.sin(a_t) * dt
    
    # 4. Shed new vortex (Kelvin's Theorem)
    gamma_airfoil_current = np.sum(0.5 * (gamma[:-1] + gamma[1:]) * l)
    dg_shed = -(gamma_airfoil_current - gamma_airfoil_prev)
    
    wake_gamma.append(dg_shed)
    wake_x.append(x[0] + 0.5 * v_inf * np.cos(a) * dt) # Shed at Trailing Edge
    wake_y.append(y[0] + 0.5 * v_inf * np.sin(a) * dt)
    
    gamma_airfoil_prev = gamma_airfoil_current
    cl_history.append(2 * gamma_airfoil_current / (v_inf * c))
    time_history.append(step * dt)

#-------------------------
# Visualization
#-------------------------
plt.figure(figsize=(10, 4))
plt.plot(x, y, 'k-'); plt.scatter(wake_x, wake_y, c=wake_gamma, s=5, cmap='coolwarm')
plt.title("Wake Path"); plt.axis('equal'); plt.xlabel("x"); plt.ylabel("y")
plt.show()