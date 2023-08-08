from Control_MYO.Helpers import clamp
class PID(object):
    # Constructor
    def __init__(self,
                 KP = 1.0, KI = 0.0, KD = 0.0,
                 u_bounds = (None, None)):
        self.KP = KP
        self.KI = KI
        self.KD = KD
        self.u_min, self.u_max = u_bounds

        self.Integral = 0   # To keep the Intergral State
        self.error_last = 0 # Since there is error_last at intial point
    
    def __call__(self, error, dt):
        
        self.Proportional = self.KP*error
        self.Integral    += self.KI*error*dt
        self.Integral     = clamp(self.Integral,self.u_min,self.u_max)
        self.Differential = self.KD*(error-self.error_last)/dt
        u = self.Proportional + self.Integral + self.Differential
        u = clamp(u,self.u_min,self.u_max)

        self.error_last = error
        return u

class PID_OnRate(object):
    # Constructor
    def __init__(self,
                 KP = 1.0, KI = 0.0, KD = 0.0,
                 u_bounds = (None, None)):
        self.KP = KP
        self.KI = KI
        self.KD = KD
        self.u_min, self.u_max = u_bounds

        self.Integral = 0   # To keep the Intergral State
    
    def __call__(self, error, x_rate, dt):
        
        self.Proportional = self.KP*error
        self.Integral    += self.KI*error*dt
        self.Integral     = clamp(self.Integral,self.u_min,self.u_max)
        self.Differential = self.KD*x_rate
        u = self.Proportional + self.Integral + self.Differential
        u = clamp(u,self.u_min,self.u_max)

        return u