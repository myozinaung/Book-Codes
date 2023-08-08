import numpy as np
import math
from Control_MYO.Helpers import clamp

# 1st Order Nomoto Model
# T*r_dot + r = K*delta
def Nomoto1st(t,x,u,p):
    ### INPUTS ###
    # t     = time [for time-dependent terms]
    # x     = state vector
    # u     = control vector
    # p     = parameters [for cosnstant and time-varying]
    
    ### OUTPUT ###
    # x_dot = derivative of x (same length as x)
    x_dot = np.zeros(x.size)
    
    # ------Start of own code------
    K = 0.1705
    T = 7.1167
    
    r     = x[1]
    delta = u[0]
 
    x_dot[0] = r
    x_dot[1] = (K*delta-r)/T + p
    # ------End of own code--------
    
    return x_dot


# Nonlinear Nomoto Model
# T*r_dot + (n3*r^3 + n1*r) = K*delta
# n1 = 1 (course stable), n1 = -1 (course unstable), n3 = nonlinear term
def Nomoto1stNL(t,x,u,p):
    x_dot = np.zeros(x.size)
    
    # ------Start of own code------
    K = 0.1705
    T = 7.1167

    n1 = -1.25 # (1:course stable, -1:course unstable)
    n3 = 3
    
    r     = x[1]
    delta = u[0]
 
    x_dot[0] = r
    x_dot[1] = (K*delta - n3*r**3 - n1*r)/T + p
    # ------End of own code--------
    
    return x_dot

def Nomoto1stNL_RudDyn2(t,x,u,p):
    x_dot = np.zeros(x.size)
    
    # ------Start of own code------
    K = 0.1705
    T = 7.1167
    TR = 5

    n1 = -1.25 # (1:course stable, -1:course unstable)
    n3 = 3
    
    r     = x[1]
    r_dot = x[2]
    delta = u[0]
 
    x_dot[0] = r
    x_dot[1] = r_dot + p
    x_dot[2] = (K*delta - (T+TR)*r_dot - n3*r**3 - n1*r)/(T*TR)
    # ------End of own code--------
    
    return x_dot

def Nomoto1stNL_RudDyn1(t,x,u,p): # state 3 is delta, not psi_dddot
    x_dot = np.zeros(x.size)
    
    # ------Start of own code------
    K = 0.1705
    T = 7.1167
    KR = 1
    TR = 5

    n1 = -1.25 # (1:course stable, -1:course unstable)
    n3 = 3
    
    delta = u[0]
 
    x_dot[0] = x[1]
    x_dot[1] = (-1/T)*(n3*x[1]**3 + n1*x[1]) + (K/T)*x[2] + p
    x_dot[2] = (-1/TR)*x[2] + (KR/TR)*delta
    # ------End of own code--------
    
    return x_dot

def Nomoto1st_FF(x_d, K_add):
    K = 0.1705
    T = 7.1167

    r_d     = x_d[1]
    r_dot_d = x_d[2]
    delta_FF = (1/K + K_add)*(T*r_dot_d + r_d)
    return delta_FF

def Nomoto1stNL_FF(x_d, K_add):
    K = 0.1705
    T = 7.1167
    n1 = -1.25
    n3 = 3

    r_d     = x_d[1]
    r_dot_d = x_d[2]
    delta_FF = (1/K + K_add)*(T*r_dot_d + n3*r_d**3 + n1*r_d)
    return delta_FF

# Nonlinear Mapping of delta_linear to delta_nonlinear
# used with Feedbacl Linearization
def Nomoto1st_FL(x, delta_linear):
    K = 0.1705
    T = 7.1167
    n3 = 3

    r = x[1]

    delta_NL = (n3/K)*r**3 + (T/K)*delta_linear
    return delta_NL

# 3rd Order ODE Model, obtained from 2nd Order Nomoto Model
def BechWenger3rdNL(t,x,u,p): # Not working yet
    # A Tanker
    L = 350 # [m]
    V = 5 # [m/s]
    
    # Tanker (Ballast Condition)
    K0 = 5.88
    T10 = -16.91*-1
    T20 = 0.45
    T30 = 1.43

    # # Tanker (Fully Loaded)
    # K0 = 0.83
    # T10 = -2.88
    # T20 = 0.38
    # T30 = 1.07

    K  = K0*(V/L)
    T1 = T10*(L/V)
    T2 = T20*(L/V)
    T3 = T30*(L/V)
     
    n3 = 1# alpha
    n1 = 1# beta

    a = (1/T1) + (1/T2)
    b = 1/(T1*T2)
    c = (K*T3)/(T1*T2)
    d = K/(T1*T2)

    H = n3*x[1]**3 + n1*x[1]
    delta = u[0]
    x_dot = np.zeros(x.size)
    delta_dot = p # can be obtained by Steering Gear Dynamics
    x_dot[0] = x[1] # psi_dot
    x_dot[1] = x[2] # psi_ddot
    x_dot[2] = - a*x[2] - b*H + c*delta_dot + d*delta
    
    # # Alternative, with delta_dot
    # x_dot[0] = x[1] # psi_dot
    # x_dot[1] = x[2] + c*delta# psi_ddot
    # x_dot[2] = - a*(x[2]+c*delta) - b*H + d*delta # psi_dddot

    return x_dot

def SteeringGear(delta, delta_d, KR, TR, delta_max, delta_dot_max, dt): # Tanker

    delta_d = clamp(delta_d, -delta_max, delta_max)
    # delta_e = delta_d - delta
    delta_dot = (-1/TR)*delta + delta_d
    delta_dot = clamp(delta_dot, -delta_dot_max, delta_dot_max)

    delta = delta + delta_dot*dt # Euler Integration
    delta = delta*(KR/TR)
    delta = np.array([delta])
    return delta, delta_dot