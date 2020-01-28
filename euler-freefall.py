from vpython import * 
import numpy

# Roslyn Parker
# Phsyical Simulations
# 28 January 2020

g = numpy.array([0, -9.8, 0]) # force of gravity

def force(x: numpy.array, m: float):
    return m * g

def integrate(x: numpy.array, v: numpy.array, dt: float, m: float):
    xNext = x + dt * v # x(t + dt)
    vNext = v + dt * force(x, m) / m
    return xNext, vNext

# Simulation Parameters
x = numpy.array([0, 20, 0]) # position
v = numpy.array([0, 0, 0]) # velocity
m = 2 # mass
dt = 0.1

print(x, v)

t = 0
while t < 5:
    t = t + dt
    x, v = integrate(x, v, dt, m) 
    print(x, v)
    if x[1] <= 0:
        break