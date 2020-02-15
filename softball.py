from vpython import * 
import numpy

# Roslyn Parker
# Phsyical Simulations
# 13 February 2020

def fGravity():
    """Returns the force of gravity on the object"""
    return m * numpy.array([0, -g, 0])

def fDrag(v: numpy.array):
    """Returns the rayleigh drag force on the object"""
    return -0.5 * rho * (pi * r**2) * cD * v * numpy.linalg.norm(v)

def force(x: numpy.array, v: numpy.array):
    """Returns the sum of all forces on the object"""
    return fGravity() + fDrag(v)

def euler(xn: numpy.array, vn: numpy.array, dt: float, m: float):
    """Predicts the positon and velocity a time dt in the future using
       Euler's method"""
    xnp1 = xn + dt * vn # x(t + dt)
    vnp1 = vn + dt * force(xn, m) / m
    return xnp1, vnp1

def implicitEuler(xn: numpy.array, vn: numpy.array, dt: float, m: float):
    """Predicts the positon and velocity a time dt in the future using
       semi-implicit Euler integration"""
    vnp1 = vn + dt * force(xn, m) / m
    xnp1 = xn + dt * vnp1 # x(t + dt)
    return xnp1, vnp1

def f(phi: numpy.array):
    x = phi[0:3]
    v = phi[3:6]

    dphidt = numpy.zeros(6) # [0, 0, 0, 0, 0, 0]
    dphidt[0:3] = v #dx/dt = v
    dphidt[3:6] = force(x, m) / m #dv/dt = F / m
    return dphidt

def rk2(xn: numpy.array, vn: numpy.array, dt: float, m: float):
    """Predicts the position and velocity a time dt in the future using the
       second-order Runge-Kutta method"""
    phin = numpy.concatenate((xn, vn))

    k1 = f(phin) * dt
    k2 = f(phin + k1 / 2) * dt
    phinp1 = phin + k2

    xnp1 = phinp1[0:3]
    vnp1 = phinp1[3:6]

    return xnp1, vnp1

def rk4(xn: numpy.array, vn: numpy.array, dt: float, m: float):
    """Predicts the position and velocity a time dt in the future using the
       fourth-order Runge-Kutta method"""
    phin = numpy.concatenate((xn, vn))

    k1 = f(phin) * dt
    k2 = f(phin + k1 / 2) * dt
    k3 = f(phin + k2 / 2) * dt
    k4 = f(phin + k3) * dt
    phinp1 = phin + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    xnp1 = phinp1[0:3]
    vnp1 = phinp1[3:6]

    return xnp1, vnp1

# Throw Constants
tSpeed = 25 # m/s
tAngle = radians(45) # 45 degrees to radians
x0 = numpy.array([0, 2, 0]) # starting position in m
v0 = numpy.array([25, 0, 0])

# Constants
m = 0.185 # kg
r = 0.037 # m 
cD = 0.4 # drag coefficient set to 0 for no drag
g = 9.8 # gravity
rho = 1.225

tf = 10 # simulation endtime
dt = 0.033 # s

integrators = {'Euler': euler, 'Semi-Implicit Euler': implicitEuler, 'RK2': rk2, 'RK4': rk4}

# Loop over integrators
for method in integrators:
    xgraph = graph(title=method, xtitle='X', ytitle='Y')
    xSim = gcurve(graph=xgraph, color=color.magenta)
    
    # Run simulation using this integrator
    t = 0
    x = x0
    v = v0
    xSim.plot(t, x[1])

    while t < tf:
        t = t + dt
        x, v = integrators[method](x, v, dt, m)
        xSim.plot(t, x[1])