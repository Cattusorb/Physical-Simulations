from vpython import * 
import numpy

# Roslyn Parker
# Phsyical Simulations
# 30 January 2020

g = 9.8

def force(x: numpy.array, k: float):
    return -(k * x)

def integrate(x: numpy.array, v: numpy.array, dt: float, m: float, k: float):
    xNext = x + dt * v # x(t + dt)
    vNext = v + dt * force(x, k) / m
    return xNext, vNext

# Simulation Parameters
x0 = numpy.array([.25, 0, 0]) # initial position
v0 = numpy.array([0, 0, 0]) # initial velocity
m = 2 # mass
k = 2 # spring constant
dt = 0.01
x = x0 # current position
v = v0 # current velocity
w = sqrt(k/m)
E = m * v**2 / 2 + k * x**2 / 2

# Euler mass Spring vs. Exact Mass Spring
grph = graph(width=600, height=400, title='Euler Mass-Spring Simulation', xtitle='Time', ytitle='Position')
curve = gcurve(graph=grph, color=color.magenta, markers=True, label='curve')
ecurve = gcurve(graph=grph, color=color.green, label='exact curve')
mechcurve = gcurve(graph=grph, color=color.yellow, markers=True, label='mechanical energy curve')
emechcurve = gcurve(graph=grph, color=color.orange, label='exact mechanical energy curve')

t = 0
while t < 20:
    t = t + dt
    x, v = integrate(x, v, dt, m, k) 
    E = m * v**2 / 2 + k * x**2 / 2

    # Exact Values
    exactX = x0[0] * cos(w * t)
    exactV = -w * x0 * sin(w * t)
    exactE = m * exactV**2 / 2 + k * exactX**2 / 2
    
    #plot
    curve.plot(t, x[0])
    ecurve.plot(t, exactX)
    mechcurve.plot(t, E[0])
    emechcurve.plot(t, exactE[0])
