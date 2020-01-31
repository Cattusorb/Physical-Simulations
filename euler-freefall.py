from vpython import * 
import numpy

# Roslyn Parker
# Phsyical Simulations
# 30 January 2020

g = 9.8

def force(x: numpy.array, m: float):
    return m * numpy.array([0, -g, 0]) # force of gravity

def integrate(x: numpy.array, v: numpy.array, dt: float, m: float):
    xNext = x + dt * v # x(t + dt)
    vNext = v + dt * force(x, m) / m
    return xNext, vNext

# Simulation Parameters
x0 = numpy.array([0, 20, 0]) # initial position
v0 = numpy.array([0, 0, 0]) # initial velocity
m = 2 # mass
dt = .25
x = x0 # current position
v = v0 # current velocity

exactE = m * g * x0[1]

grph = graph(width=600, height=400, title='Euler Freefall Simulation', xtitle='Time', ytitle='Position')
curve = gcurve(graph=grph, color=color.magenta, markers=True, label='curve')
ecurve = gcurve(graph=grph, color=color.green, label='exact curve')
mechcurve = gcurve(graph=grph, color=color.yellow, label='mechanical energy curve')
emechcurve = gcurve(graph=grph, color=color.orange, label='exact mechanical energy curve')

t = 0
while t < 5:
    t = t + dt
    x, v = integrate(x, v, dt, m) 
    E = m * g * x + m * v**2 / 2
   
    # Exact Values
    y = x0[1] - g * t**2 / 2
    exactV = -(g * t)

    #plot
    curve.plot(t, x[1])
    ecurve.plot(t, y)
    mechcurve.plot(t, E[1])
    emechcurve.plot(t, exactE)
    
