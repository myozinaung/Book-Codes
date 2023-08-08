from Control_MYO import Integrators
import math
import numpy as np

def Lowpass1st(x, u, omega_f, dt): # x: filter ode state, u:signal to be filtered, T_filter: cut-off frequency time
    def f(t,x,u,p):
        return -(omega_f)*x + u
    x = Integrators.RKGill(f,[],x,u,[],dt)
    y = (omega_f)*x
    return y, x

def Lowpass2nd(x, u, omega_f, dt):
    def f(t,x,u,p):
        u = np.array([[u]])
        zeta = math.sin(45*math.pi/180)
        # Denominator Coefficients
        s2 = 1
        s1 = 2*zeta*omega_f
        s0 = omega_f**2
        A = np.array([[-s1, -s0],
                      [1, 0]])
        B = np.array([[1],[0]])
        x = A.dot(x) + B.dot(u)
        return x
    x = Integrators.RKGill(f,[],x,u,[],dt)
    
    n0 = omega_f**2 # Numerator Coefficient
    C = np.array([[0, n0]])
    y = C.dot(x)
    return y[0,0], x    

def Lowpass3rd(x, u, omega_f, dt):
    def f(t,x,u,p):
        u = np.array([[u]])
        zeta = math.sin(30*math.pi/180)
        # Denominator Coefficients
        s3 = 1
        s2 = (1+2*zeta)*omega_f
        s1 = (1+2*zeta)*omega_f**2
        s0 = omega_f**3
        A = np.array([[-s2, -s1, -s0],
                      [1, 0, 0],
                      [0, 1, 0]])
        B = np.array([[1],[0],[0]])
        x = A.dot(x) + B.dot(u)
        return x
    x = Integrators.RKGill(f,[],x,u,[],dt)
    
    n0 = omega_f**3 # Numerator Coefficient
    C = np.array([[0, 0, n0]])
    y = C.dot(x)
    return y[0,0], x    

def Lowpass4th(x, u, omega_f, dt):
    def f(t,x,u,p):
        u = np.array([[u]])
        zeta1 = math.sin(22.5*math.pi/180)
        zeta2 = math.sin(67.5*math.pi/180)
        # Denominator Coefficients
        s4 = 1
        s3 = 2*(zeta1+zeta2)*omega_f
        s2 = 2*(2*zeta1*zeta2+1)*omega_f**2
        s1 = 2*(zeta1+zeta2)*omega_f**3
        s0 = omega_f**4
        A = np.array([[-s3, -s2, -s1, -s0],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0]])
        B = np.array([[1],[0],[0],[0]])
        x = A.dot(x) + B.dot(u)
        return x
    x = Integrators.RKGill(f,[],x,u,[],dt)
    
    n0 = omega_f**4 # Numerator Coefficient
    C = np.array([[0, 0, 0, n0]])
    y = C.dot(x)
    return y[0,0], x  
    
def NotchFilter(x, u, omega_f, zeta, dt):
    zeta = -zeta + 1 # to reverse (from 0-1 to 1-0)
    s1 = 2*omega_f
    s0 = omega_f**2
    def f(t,x,u,p):
        u = np.array([[u]])

        # Denominator Coefficients
        s2 = 1 
        s1 = 2*omega_f
        s0 = omega_f**2
        A = np.array([[-s1, -s0],
                      [1, 0]])
                      
        B = np.array([[1],[0]])
        x = A.dot(x) + B.dot(u)
        return x
    x = Integrators.RKGill(f,[],x,u,[],dt)
    
    # Numerator Coefficients
    n2 = 1 
    n1 = 2*zeta*omega_f
    n0 = omega_f**2 
    C = np.array([[n1-s1*n2, n0-s0*n2]])
    D = np.array([[n2]])
    y = C.dot(x) + D.dot(u)
    return y[0,0], x   

def Highpass1st(x, u, omega_f, zeta, dt): # x: filter ode state, u:signal to be filtered, T_filter: cut-off frequency time
    def f(t,x,u,p):
        A = np.array([[-omega_f]])
        B = np.array([[1]])
        x = A.dot(x) + B.dot(u)
        return x
    x = Integrators.RKGill(f,[],x,u,[],dt)
    
    C = np.array([[-zeta*omega_f]])
    D = np.array([[zeta]])
    y = C.dot(x) + D.dot(u)

    return y[0,0], x    

# s alone not OK, so cascade with lowpass filter = s*omega_f/(s+omega_f)
# if omega_f is large enough --> LP filter becomes "1"
def Derivative(x, u, omega_f, dt): 
    def f(t,x,u,p):
        return -omega_f*x + u
    # x = Integrators.euler(f,[],x,u,[],dt)
    x = Integrators.RKGill(f,[],x,u,[],dt)
    y = -omega_f**2*x + omega_f*u
    return y, x

def DerivativeBackward(x_last, x, dt): 
    return (x - x_last)/dt

# Reference Generator Full State Feedback control and Feedforward Term
# Generate velocity and acceleration reference (v_ref and a_ref) from the position ref (x_ref)
# It also smooth(filter) sudden change in x_ref
def RefGen3rd(x, x_ref, omega_n, zeta, dt, v_max):
# Input: x_ref (position only)
# Parameters: omega_n(natural frequency), zeta(damping coefficient)
# Output: x (desired position, velocity and acceleration) filtered with 3rd order LP
    x_dot = np.zeros(3)
    x_dot[0] = x[1]
    x_dot[1] = x[2]
    if abs(x[1]) > v_max:
        x[1] = np.sign(x[1])*v_max
    x_dot[2] = -(2*zeta+1)*omega_n*x[2] - (2*zeta+1)*omega_n**2*x[1] + omega_n**3*(x_ref - x[0])
    x = x + x_dot*dt # Euler Integration
    # Saturation can be imposed on "x" if needed (e.g. Velocity saturation)
    
    return x

def RefGen2nd(x, x_ref, omega_n, zeta, dt, v_max):
    x_dot = np.zeros(2)
    x_dot[0] = x[1]
    x_dot[1] = - 2*zeta*omega_n*x[1] + omega_n**2*(x_ref - x[0])
    x = x + x_dot*dt
    if abs(x[1]) > v_max:
        x[1] = np.sign(x[1])*v_max
    return x

def RefGen4th(x, x_ref, omega_n, zeta1, zeta2, dt, v_max):
    x_dot = np.zeros(4)
    x_dot[0] = x[1]
    x_dot[1] = x[2]
    if abs(x[1]) > v_max:
        x[1] = np.sign(x[1])*v_max
    x_dot[2] = x[3]
    x_dot[3] = -2*(zeta1+zeta2)*omega_n*x[3] - 2*(2*zeta1*zeta2+1)*omega_n**2*x[2] - 2*(zeta1+zeta2)*omega_n**3*x[1] + omega_n**4*(x_ref - x[0])
    x = x + x_dot*dt # Euler Integration
    # Saturation can be imposed on "x" if needed (e.g. Velocity saturation)
    
    return x

