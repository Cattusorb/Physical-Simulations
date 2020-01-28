from vpython import * 
import numpy

def double(x: float):
    return 2*x

def firstElement(x: list):
    return x[0]

def f(x: float):
    return x**2, x**3 - 2*x


grph = graph(width=600, height=400, title='Graph 0')
crv = gcurve(graph=grph, color=color.magenta, markers=True, label='f(x)')
crv2 = gcurve(graph=grph, color=color.purple, label='g(x)')

grph2 = graph(width=600, height=400, title='Graph 1', xtitle='x', ytitle='y')
sinCurve = gcurve(graph=grph2, color=color.orange, markers=True, label='sin(x)')
cosCurve = gcurve(graph=grph2, color=color.red, label='cos(x)')

# Plot x vs. f(x) for 0 to 5 step by 0.1
# Add plots of sin(x) and cos(x) on grph2

for x in numpy.arange(-5, 5.1, 0.1):
    y1, y2 = f(x)
    crv.plot(x, y1)
    crv2.plot(x, y2)
    sinCurve.plot(x, sin(x))
    cosCurve.plot(x, cos(x))


