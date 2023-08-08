# Control Methods
def P(error, PID_Gains, u_bound):
    KP = PID_Gains[0]
    u_min = u_bound[0]
    u_max = u_bound[1]
    
    u = KP*error
    u = clamp(u,u_min,u_max)
    return u

def PD(error, error_last, dt, PID_Gains, u_bound):
    KP = PID_Gains[0]
    KD = PID_Gains[2]
    u_min = u_bound[0]
    u_max = u_bound[1]

    Proportional = KP*error
    Differential = KD*(error-error_last)/dt
    u = Proportional + Differential
    u = clamp(u,u_min,u_max)
    return u

def PID(error, error_last, Integral, dt, PID_Gains, u_bound):
    KP = PID_Gains[0]
    KI = PID_Gains[1]
    KD = PID_Gains[2]
    
    u_min = u_bound[0]
    u_max = u_bound[1]

    Proportional = KP*error
    Integral    += KI*error*dt                 # Forward Euler
    Integral     = clamp(Integral,u_min,u_max) # Intergral Anti-Windup (Method 1)
    Differential = KD*(error-error_last)/dt
    u = Proportional + Integral + Differential
    u = clamp(u,u_min,u_max)
    return u, Integral

    # Integral = integral[KI*error - KW*(u_unclamped - u_clamped)]dt # Anti-windup (method 2)
    
def PIDAutopilot(error, u, x_rate, Integral, dt, PID_Gains, u_bound, u_rate_bound):
    KP = PID_Gains[0]
    KI = PID_Gains[1]
    TD = PID_Gains[2]
    TE = PID_Gains[3]
    
    u_min = u_bound[0]
    u_max = u_bound[1]
    
    u_rate_min = u_rate_bound[0]
    u_rate_max = u_rate_bound[1]

    error = -error
    Proportional = KP*error
    Integral    += KI*error*dt                 # Forward Euler
    Integral     = clamp(Integral,u_min,u_max) # Intergral Anti-Windup (Method 1)
    ## General Integrater Implementation ##
    def f(t,x,u,p):
        u_dot = -(Proportional + Integral + TD*x_rate + x)/TE
        u_dot = clamp(u_dot, u_rate_min, u_rate_max) # Rate limit
        return u_dot
    u = Integrators.RKGill(f,[],u,[],[],dt)

    ## Euler Integrator Implementation ##
    # u_dot = -(Proportional + Integral + TD*x_rate + u)/TE
    # u_dot = clamp(u_dot, u_rate_min, u_rate_max) # Rate limit
    # u = u + u_dot*dt # Euler Integration (could be changed to RKGill)
    
    u = clamp(u,u_min,u_max)
    return u, Integral

import Integrators
import numpy as np    
# State-Space Version of PID with derivative filter
# if "s" is used for derivative, the order of numerator will be larger(not implementable) 
def FilteredPID(x, u, dt, PID_Gains, N, y_bound): # x: internal state, u:error
    KP = PID_Gains[0]
    KI = PID_Gains[1]
    KD = PID_Gains[2]
    # N =  Derivative Filter Coefficient
    
    y_min = y_bound[0]
    y_max = y_bound[1]

    
    def f(t,x,u,p):
        u = np.array([[u]])
        # Denominator Coefficients
        s2 = 1
        s1 = N
        s0 = 0
        A = np.array([[-s1, -s0],
                      [1, 0]])
        B = np.array([[1],[0]])
        x = A.dot(x) + B.dot(u)
        return x
    x = Integrators.RKGill(f,[],x,u,[],dt)
    
    # Numerator Coefficients
    n2 = (KP + KD*N)
    n1 = (KP*N + KI)
    n0 = KI*N 
    C = np.array([[n1-N*n2, n0]]) # s1=N and s0=0, C = np.array([[n1-s1*n2, n0-s0*n2]])
    D = np.array([[n2]])
    y = C.dot(x) + D.dot(u)

    y = clamp(y[0,0],y_min,y_max)
    return y, x

# Helper Functions
# Clapm x between x_min and x_max
def clamp(x,x_min,x_max):
    if x <= x_min:
        return x_min
    elif x >= x_max:
        return x_max
    else:
        return x