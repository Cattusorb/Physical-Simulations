from vpython import *
import numpy
import random

# Roslyn Parker
# Physical Simulations
# 18 February 2020

# Firework Simulation

# Launch firework from the ground
# Fly through the air
# Detonate into fragments

# You can make one firework, or a show!

# Set the make_trail attribute for the fragments

# Make the initial launch invisible visible=False

def force(x: numpy.array, v: numpy.array):
    """Returns the force of gravity on the object"""
    return m * numpy.array([0, -g, 0])

def f(phi: numpy.array):
    x = phi[0:3]
    v = phi[3:6]

    dphidt = numpy.zeros(6) # [0, 0, 0, 0, 0, 0]
    dphidt[0:3] = v #dx/dt = v
    dphidt[3:6] = force(x, v) / m #dv/dt = F / m
    return dphidt

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

def launchFirework(start: int):
    firework = sphere(pos=vector(start,0,start), visible=False)
    x0 = numpy.array([start, 0, start]) # starting position in m
    v0 = fSpeed * numpy.array([cos(fAngle), sin(fAngle), 0])
    
    t = 0
    x = x0
    v = v0

    # Launch rocket to explode height, average 15-22 m
    while x[1] <= random.randrange(15, 30):
        rate(framerate)
        x, v = rk4(x, v, dt, m)
        firework.pos = vector(x[0], x[1], x[2])
        t = t + dt

    frags = []
    fragM = 1
    fragColor = random.choice(colors)
    for angle in range(0, 361, 20):
        speed = random.randrange(200, 500)
        fx0 = numpy.array([0, x[1], 0]) # starting position in m
        fv0 = speed * numpy.array([cos(angle), sin(angle), 0])
        frag = sphere(pos=vector(0,fx0[1],0), radius=0.2, make_trail=True, color=fragColor)
        frags.append([fx0, fv0, frag])

    t = 0
    while t < 1:
        rate(framerate)
        for frag in frags: 
            fx = frag[0]
            fv = frag[1]
            fragment = frag[2]
            fx, fv = rk4(fx, fv, dt, fragM)
            fragment.pos = vector(fx[0], fx[1], fx[2])
        t = t + dt

# Set-up Scene
scene.forward = vector(-1, -0.5, -1)
lw = 50
floor = box(pos=vector(0, 0, 0), length=lw, width=lw, height=0.5, color=color.magenta)

# Firework Constants
fSpeed = 20 # m/s ~80mph
fAngle = radians(90) # 90 degrees to radians, straight up

# Constants
g = 9.8 # gravity
m = 4.535 # kg 
framerate = 60
dt = 1.0 / framerate
colors = [color.red, color.orange, color.yellow, color.green, color.blue, color.purple]

start = [0, -20, 20, -40, 40]
for firework in range(0, 5):
    launchFirework(start[0])
    start.pop(0)