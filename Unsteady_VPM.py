import numpy as np
import json
import matplotlib.pyplot as plt


#--------------------------------
# Global variables
#--------------------------------

v_inf = 20          # freestream velocity
alpha = 0           # AOA in degrees
a = np.radians(alpha)
rho = 1.225         # air density
c = 1               # chord length
T = 5               # simulaion time
n_steps = 100       # number of time steps
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

A = build_A(x, y, xc, yc, dx, dy, l)

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
# Time Stepping Setup
#-------------------------

# 1. Initialize Wake Containers
wake_x = []       # X-positions of wake vortices
wake_y = []       # Y-positions of wake vortices
wake_gamma = []  # Strengths of wake vortices

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

for step in range(n_steps):
    # 1. Convect existing wake (vortices stay stationary relative to fluid)
    for i in range(len(wake_x)):
        wake_x[i] += v_inf * np.cos(a) * dt
        wake_y[i] += v_inf * np.sin(a) * dt
        
    # 2. Wake influence on airfoil (Section II)
    B_wake = np.zeros(n_pan)
    for i in range(n_pan):
        v_ind_n = 0
        nx, ny = -dy[i]/l[i], dx[i]/l[i]
        for wx, wy, wg in zip(wake_x, wake_y, wake_gamma):
            rx, ry = xc[i] - wx, yc[i] - wy
            r2 = max(rx**2 + ry**2, 1e-6)
            u_w, v_w = (wg / (2 * np.pi)) * (ry / r2), -(wg / (2 * np.pi)) * (rx / r2)
            v_ind_n += (u_w * nx + v_w * ny)
        B_wake[i] = -v_ind_n 
        
    # 3. Solve for new surface circulation
    B_total = np.append(B_fs + B_wake, 0)
    gamma = np.linalg.solve(A_matrix, B_total)
    
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
plt.subplot(1, 2, 1)
plt.plot(x, y, 'k-'); plt.scatter(wake_x, wake_y, c=wake_gamma, s=5, cmap='coolwarm')
plt.title("Wake Path"); plt.axis('equal'); plt.xlabel("x"); plt.ylabel("y")

plt.subplot(1, 2, 2)
plt.plot(time_history, cl_history)
plt.title("Lift Coefficient over Time"); plt.xlabel("Time (s)"); plt.ylabel("$C_L$")
plt.tight_layout(); plt.savefig('wake_analysis.png')