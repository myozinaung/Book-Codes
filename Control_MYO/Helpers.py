# Helper Functions
import math
import numpy as np
# Clamp x between x_min and x_max
def clamp(x,x_min,x_max):
    if x <= x_min:
        return x_min
    elif x >= x_max:
        return x_max
    else:
        return x
        
# Rewrite an angle(rad) within -pi to pi
def pi2pi(x):
    x = math.fmod(x+np.sign(x)*math.pi,2*math.pi)- np.sign(x)*math.pi
    # Note: Python x % y doesn't work same as math.fmod
    return x    