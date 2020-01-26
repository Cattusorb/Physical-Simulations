from vpython import * 
import numpy as numpy

# Roslyn Parker
# Physical Simulations
# 24 January 2020

# Animation parameters
framerate = 60
dt = 1.0 / framerate

# Gravity
g = numpy.array([0, -9.8, 0])

# Projectile parameters
nProjectiles = 5
launchSpeed = 30

# Floor
floorWidth = 50
floorLength = 100
floorThickness = 0.5
floor = box(center=vector(0,0,0), length=floorLength, width=floorWidth, 
height=floorThickness)

# Set up camera
scene.forward = vector(1, -1, -1)
scene.range = 60

# Set up initial conditions for projectiles
x0 = numpy.array([-floorLength/2, 0, 0])

# Specify velocity for each projectile  
launchAngleR = radians(80) # converts degrees to radians
launchAngleO = radians(60)
launchAngleY = radians(40)
launchAngleG = radians(30)
launchAngleB = radians(20)
v0R = launchSpeed * numpy.array([cos(launchAngleR), sin(launchAngleR), 0])
v0O = launchSpeed * numpy.array([cos(launchAngleO), sin(launchAngleO), 0])
v0Y = launchSpeed * numpy.array([cos(launchAngleY), sin(launchAngleY), 0])
v0G = launchSpeed * numpy.array([cos(launchAngleG), sin(launchAngleG), 0])
v0B = launchSpeed * numpy.array([cos(launchAngleB), sin(launchAngleB), 0])

def getPosition(x0, v0, g, t):
    '''
        Gets the postion of the ball, returns a numpy array
    '''
    position = x0 + v0 * t + g * t**2 / 2
    return position

def launch(position, ball, t, dt, framerate, x0, v0, g):
    '''
        launches a projectile til it hits the floor
    '''
    while position[1] >=0:
        rate(framerate)
        t = t + dt
        position = getPosition(x0, v0, g, t)
        ball.pos = vector(position[0], position[1], position[2])

t = 0
positionR = getPosition(x0, v0R, g, t)
positionO = getPosition(x0, v0O, g, t)
positionY = getPosition(x0, v0Y, g, t)
positionG = getPosition(x0, v0G, g, t)
positionB = getPosition(x0, v0B, g, t)

ballR = sphere(pos=vector(positionR[0], positionR[1], positionR[2]), color=color.red, make_trail=True)
ballO = sphere(pos=vector(positionO[0], positionO[1], positionO[2]), color=color.orange, make_trail=True)
ballY = sphere(pos=vector(positionY[0], positionY[1], positionY[2]), color=color.yellow, make_trail=True)
ballG = sphere(pos=vector(positionG[0], positionG[1], positionG[2]), color=color.green, make_trail=True)
ballB = sphere(pos=vector(positionB[0], positionB[1], positionB[2]), color=color.blue, make_trail=True)

launch(positionR, ballR, t, dt, framerate, x0, v0R, g)
launch(positionO, ballO, t, dt, framerate, x0, v0O, g)
launch(positionY, ballY, t, dt, framerate, x0, v0Y, g)
launch(positionG, ballG, t, dt, framerate, x0, v0G, g)
launch(positionB, ballB, t, dt, framerate, x0, v0B, g)