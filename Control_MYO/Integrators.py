import math
# Integrators (Fixe time step) # NO CHANGE NECESSARY
# Euler Integration Method (1st Order)
def euler(f,t,x,u,p,dt):

    x_next = x + f(t,x,u,p)*dt
    
    return x_next
    
# Range-Kutta Integration Method (4th Order)
def RKGill(f,t,x,u,p,dt):
    k1  = f(t,x,u,p)
    x_temp = x + 0.5*dt*k1
    
    k2  = f(t,x_temp,u,p)
    x_temp = x + 0.5*dt*(math.sqrt(2)-1)*k1 + dt*(1-(1/math.sqrt(2)))*k2
    
    k3  = f(t,x_temp,u,p)
    x_temp = x - 0.5*dt*math.sqrt(2)*k2 + dt*(1+(1/math.sqrt(2)))*k3
    
    k4  = f(t,x_temp,u,p)
    x_next = x + dt*(k1+(2-math.sqrt(2))*k2+(2+math.sqrt(2))*k3+k4)/6
    
    return x_next