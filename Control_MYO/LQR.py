import scipy.linalg as LA # Linear Algebra
import math
import numpy as np
from Control_MYO.Helpers import clamp, pi2pi

# Continuous-Time Linear Quadratic Regulator (LQR)
def CLQR_Gain(A, B, Q, R): # Gain need to be solved only once if system is static
    X = np.matrix(LA.solve_continuous_are(A, B, Q, R))
    K = np.matrix(LA.inv(R)*(B.T*X))
    eigVals, eigVecs = LA.eig(A-B*K)
    print("LQR Gains:-K",-K)
    return K

class LQR(object):
    # Constructor
    def __init__(self, A, B, Q, R, u_bounds):
        self.K = CLQR_Gain(A, B, Q, R)
        self.u_min = u_bounds[0]
        self.u_max = u_bounds[1]
    def __call__(self, x):
        u = -self.K @ x.T
        u = pi2pi(u)
        u = clamp(u, self.u_min, self.u_max)
        return u, self.K
