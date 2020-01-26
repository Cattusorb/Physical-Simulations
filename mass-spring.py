from vpython import * 
import numpy as numpy
import math 

# Roslyn Parker
# Physical Simulations
# 26 January 2020


# Animation parameters
framerate = 60
dt = 1.0 / framerate

l0 = 100 # equilibrium length of spring in cm

# Set up camera
scene.forward = vector(-1, -0.8, -1)

# Set up ceiling, mass and spring
table = box(pos=vector(0, 0, 0), length=400, width=200, height=0.5, color=color.blue)
wall = box(pos=vector(0, 25, 0), length=0.5, width=150, height=50, color=color.orange)
spring = helix(pos=vector(0, 25,0), axis=vector(1, 0, 0), length=l0, radius=15, thickness=1, color=color.red)
mass = box(pos=vector(l0, 25, 0), length=40, width=40, height=50)

m = 1 #mass
k = 20 #spring constant 
t = 0 #time
x = l0 #position
A = 40 #initial compression

while t < 100:
    rate(framerate)
    t = t + dt

    w = sqrt(k/m)
    x =  l0 + (cos(w * t) * (l0 - A))
    # x is equal to this because it has to be occulating around the equilibrium 
    # point, which is why it is l0 + and not l0 *. The cos(w * t) is multipled by (l0 - A)
    # because we have to know how much the spring is compressed at any given point.
    mass.pos.x =  x
    spring.length = x
