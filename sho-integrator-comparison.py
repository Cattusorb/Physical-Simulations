from vpython import *
import numpy

# Roslyn Parker
# Physical Simulations
# 7 February 2020

# "gravity is real" - Jeremy Ouellette 2020

def force(pos: numpy.array, m: float):
    """Force for this object as a function of position and its mass"""
    return -(k * pos)

def euler(xn: numpy.array, vn: numpy.array, dt: float, m: float):
    """Predicts the positon and velocity a time dt in the future using
       Euler's method"""
    xnp1 = xn + dt * vn # x(t + dt)
    vnp1 = vn + dt * force(xn, k) / m
    return xnp1, vnp1

def implicitEuler(xn: numpy.array, vn: numpy.array, dt: float, m: float):
    """Predicts the positon and velocity a time dt in the future using
       semi-implicit Euler integration"""
    vnp1 = vn + dt * force(xn, k) / m
    xnp1 = xn + dt * vnp1 # x(t + dt)
    return xnp1, vnp1

def leapfrog(xn: numpy.array, vn: numpy.array, dt: float, m: float):
    """Predicts the positon and velocity a time dt in the future using
       the leapfrog method"""
    xnp1 = xn + vn * dt + force(xn, k) / m * dt**2 / 2
    vnp1 = vn + (force(xn, k) + force(xnp1, k)) / m * dt / 2
    return xnp1, vnp1

def f(phi: numpy.array):
    x = phi[0:3]
    v = phi[3:6]

    dphidt = numpy.zeros(6) # [0, 0, 0, 0, 0, 0]
    dphidt[0:3] = v #dx/dt = v
    dphidt[3:6] = force(x, k) / m #dv/dt = F / m
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

def energy(x: numpy.array, v: numpy.array, m: float):
    """Returns the energy of an object in free fall along the x axis"""
    return k * x[0]**2 / 2 + m * v[0]**2 / 2

# Constants
force
g = 9.8 # gravitational acceleration
m = 2 # mass in kg
k = 8 * pi**2
x0 = numpy.array([0.25, 0, 0]) # Initial position in m
v0 = numpy.array([0, 0, 0]) # Initial velocity in m/s
tf = 30 # Simulation end time
dt = 0.05 # timestep

# Calculate total energy
eTot = k * x0[0]**2 / 2

# Define integrators
integrators = {'Euler': euler, 'Semi-Implicit Euler': implicitEuler, 'Leapfrog': leapfrog, 'RK2': rk2, 'RK4': rk4}

# Loop over integrators
for method in integrators:
    xgraph = graph(title=method, ytitle='y (m)')
    xSim = gcurve(graph=xgraph, color=color.blue)
    xExact = gcurve(graph=xgraph, color=color.black)

    egraph = graph(xtitle='t (s)', ytitle='E (J)')
    eSim = gcurve(graph=egraph, color=color.blue)
    eExact = gcurve(graph=egraph, color=color.black)
    
    # Run simulation using this integrator
    t = 0
    x = x0
    v = v0
    xSim.plot(t, x[0])
    xExact.plot(t, x[0])
    eSim.plot(t, energy(x, v, m))
    eExact.plot(t, eTot)

    while t < tf:
        t = t + dt
        x, v = integrators[method](x, v, dt, m)
        xSim.plot(t, x[0])
        xExact.plot(t, x0[0] * cos(2 * pi * t))
        eSim.plot(t, energy(x, v, m))
        eExact.plot(t, eTot)