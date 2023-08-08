import math
import numpy as np
from Control_MYO.Helpers import clamp, pi2pi

class SMC(object):
    # Constructor
    def __init__(self, lam, k, u_bounds):
        self.lam   = np.array([lam])
        self.k     = np.array([[k]]) 
        self.u_min = u_bounds[0]
        self.u_max = u_bounds[1]
        print("SMC Object Created!")
    def __call__(self, x, x_d):

        K = 0.1705
        T = 7.1167
        epsilon = 3e-2*1 # Chattering depends on epsilon and dt

        A = np.array([[0, 1],
                      [0, -1/T]])
        B = np.array([[0],
                      [K/T]])
        # Trim x and x_d, only 1st two needed
        x   = np.array([x[:2,]]).T
        x_d = np.array([x_d[:2,]]).T

        u_equiva = np.linalg.inv(self.lam@B)@self.lam@(x_d[1] - A@x)

        error = x - x_d
        sigma = self.lam@error
        # u_switch = -np.linalg.inv(self.lam@B)@self.k@np.sign(sigma)
        u_switch = -np.linalg.inv(self.lam@B)@self.k@(sigma/(abs(sigma)+epsilon))
        u = u_equiva + u_switch
        u = clamp(u, self.u_min, self.u_max)
        return u, sigma

class SMC_NL(object):
    # Constructor
    def __init__(self, lam, k, u_bounds):
        self.lam   = lam
        self.k     = k
        self.u_min = u_bounds[0]
        self.u_max = u_bounds[1]
        print("SMC with Nonlinear Nomoto Model Object Created!")
    def __call__(self, x, x_d):

        K  = 0.1705
        T  = 7.1167
        n1 = -1.25
        n3 = 3

        epsilon = 3e-2*1 # Chattering depends on epsilon and dt

        f = (-1/T)*(n3*x[1]**3 + n1*x[1])
        b = K/T

        # # Trim x and x_d, only 1st two needed
        # x   = np.array([x[:2,]]).T
        # x_d = np.array([x_d[:2,]]).T

        u_equiva = (b)**(-1)*(-f + x_d[1] - self.lam*(x[1] - x_d[1]))

        sigma = self.lam*(x[0] - x_d[0]) + (x[1] - x_d[1])
        # u_switch = (b)**(-1)*(-self.k*np.sign(sigma))
        u_switch = (b)**(-1)*(-self.k*(sigma/(abs(sigma)+epsilon)))
        u = u_equiva + u_switch
        u = clamp(u, self.u_min, self.u_max)
        return u, sigma

# class SMC_Asym(object):
#     # Constructor
#     def __init__(self, lam, k, u_bounds, u_rate_bounds):
#         self.lam   = lam
#         self.lam_bar = lam_bar
#         self.k     = k
#         self.u_min = u_bounds[0]
#         self.u_max = u_bounds[1]
#         self.u_rate_max = u_rate_bounds[0]
#         self.u_rate_min = u_rate_bounds[1]
#         print("SMC with Nonlinear Nomoto Model Object Created!")
#     def __call__(self, x, x_d, u_vec):

#         K  = 0.1705
#         T  = 7.1167
#         n1 = -1.25
#         n3 = 3

#         epsilon = 3e-2*1 # Chattering depends on epsilon and dt

#         x_ddot = # r_dot, psi_ddot
#         f     = (-1/T)*(n3*x[1]**3 + n1*x[1])
#         f_dot = (-1/T)*(3*n3*x[1]**2*x_ddot + n1*x_ddot)
#         b = K/T

#         u     = u_vec[0]
#         u_dot = u_vec[1]

#         # # Trim x and x_d, only 1st two needed
#         # x   = np.array([x[:2,]]).T
#         # x_d = np.array([x_d[:2,]]).T

#         u_equiva = (b)**(-1)*(-f + x_d[1] - self.lam*(x[1] - x_d[1]))

#         e      = x[0] - x_d[0]
#         e_dot  = x[1] - x_d[1]
#         e_ddot = f_dot + b*u_dot - x_d[2]
#         sigma     = self.lam*e + e_dot
#         sigma_dot = self.lam*e_dot + e_ddot
#         s = self.lam_bar*sigma + sigma_dot
        
#         u_switch = (b)**(-1)*(-self.k*np.sign(s))

#         u_dot = u_equiva + u_switch
#         # u_dot = clamp(u_dot, self.u_rate_min, self.u_rate_max)
#         u = u + u_dot*dt    # Euler integration
#         u = clamp(u, self.u_min, self.u_max)
#         return u, sigma