from vpython import *
import numpy

# Roslyn Parker
# Physical Simulations
# 13 February 2020

def fGravity():
    """Returns the force of gravity on the object"""
    return m * numpy.array([0, -g, 0])

def fDrag(v: numpy.array):
    """Returns the rayleigh drag force on the object"""
    return -0.5 * rho * (pi * r**2) * cD * v * numpy.linalg.norm(v)

def fSpring(x: numpy.array):
    """Returns the spring force on the object"""
    xhat = x / numpy.linalg.norm(x)
    springRestPos = l0 * xhat
    return -k * (x - springRestPos)

def force(x: numpy.array, v: numpy.array):
    """Returns the sum of all forces on the object"""
    return fGravity() + fSpring(x) + fDrag(v)

def f(phi: numpy.array):
    """Returns the time derivative f(phi) for a set of coupled differential
       equations of the form dphi/dt = f(phi)."""
    x = phi[0:3]
    v = phi[3:6]
    dphidt = numpy.zeros(6)
    dphidt[0:3] = v
    dphidt[3:6] = force(x, v) / m
    return dphidt

def rk4(xn: numpy.array, vn: numpy.array, dt: float):
    """Predicts the position and velocity a time dt in the future using the
       fourth-order Runge-Kutta method"""
    # Form state vector from xn and vn
    phin = numpy.zeros(6)
    phin[0:3] = xn
    phin[3:6] = vn
    # Do RK4 Step
    k1 = dt * f(phin)
    k2 = dt * f(phin + k1 / 2)
    k3 = dt * f(phin + k2 / 2)
    k4 = dt * f(phin + k3)
    phinp1 = phin + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    # Return the relevant portions of the phinp1
    return phinp1[0:3], phinp1[3:6]

# Environmental parameters
g = 9.8
rho = 1.225
zeta = 2 # run with 1, and 2 as well

l0 = 0.5 # Natural length of spring in m
k = pi**2 / 10 # spring constant in N/m
m = 0.100 # mass in kg
r = 0.2 # sphere radius in m
cD = zeta * 2 * sqrt(m * k) # drag coefficient | to turn off, set to 0

# Define starting position
x = numpy.array([l0 / 4, -l0 - m * g / k - 0.1, 0])
v = numpy.array([0, 0, 0])

# Graph
graph(title='damped sho', xtitle='time', ytitle='position')
crv = gcurve(color=color.magenta)

t = 0 
dt = 0.067

while t <= 10:
    x, v = rk4(x, v, dt)
    crv.plot(t, x[1])
    t = t + dt