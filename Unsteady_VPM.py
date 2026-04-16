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

for i in range(n_steps):

    gamma_surface = np.zerros(n_steps) # circulation on the airfoil surface at each time step
    gamma_wake = np.zerros(n_steps-1) # circulation of each wake vortex, there is no wake at T=0\
