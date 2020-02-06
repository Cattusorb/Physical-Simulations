from vpython import *
import numpy

# Roslyn Parker
# Physical Simulations
# 4 February 2020

# "gravity is real" - Jeremy Ouellette 2020

def force(pos: numpy.array, m: float):
    """Force for this object as a function of position and its mass"""
    return m * numpy.array([0, -g, 0])

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

def leapfrog(xn: numpy.array, vn: numpy.array, dt: float, m: float):
    """Predicts the positon and velocity a time dt in the future using
       the leapfrog method"""
    xnp1 = xn + vn * dt + force(xn, m) / m * dt**2 / 2
    vnp1 = vn + (force(xn, m) + force(xnp1, m)) / m * dt / 2
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

def energy(x: numpy.array, v: numpy.array, m: float):
    """Returns the energy of an object in free fall along the y axis"""
    return m * v[1]**2 / 2 + m * g * x[1]

# Constants
force
g = 9.8 # gravitational acceleration
m = 2 # mass in kg
x0 = numpy.array([0, 20, 0]) # Initial position in m
v0 = numpy.array([0, 0, 0]) # Initial velocity in m/s
tf = 10 # Simulation end time
dt = 1 # timestep


# Calculate total energy
eTot = energy(x0, v0, m)

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
    xSim.plot(t, x[1])
    xExact.plot(t, x[1])
    eSim.plot(t, energy(x, v, m))
    eExact.plot(t, eTot)

    while t < tf:
        t = t + dt
        x, v = integrators[method](x, v, dt, m)
        xSim.plot(t, x[1])
        xExact.plot(t, x0[1] - g * t**2 / 2)
        eSim.plot(t, energy(x, v, m))
        eExact.plot(t, eTot)
