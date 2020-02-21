from vpython import *
import numpy

# Roslyn Parker
# Physical Simulations
# 21 February 2020

# Launch firework from the ground
# Fly through the air
# Detonate into fragment

def fGravity():
    """Returns the force of gravity on the object"""
    return m * numpy.array([0, -g, 0])

def fDrag(v: numpy.array):
    """Returns the rayleigh drag force on the object"""
    return -0.5 * rho * (pi * r**2) * cD * v * numpy.linalg.norm(v)

def force(x: numpy.array, v: numpy.array):
    """Returns the force of gravity on the object"""
    return fGravity() + fDrag(v)

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

def initialLaunch(startX: float, startZ: float):
    '''
        Launches the firework up into the air to a certain height
    '''
    firework = sphere(pos=vector(startX,0.5,startZ), radius=1, make_trail=True, trail_type="points")
    x0 = numpy.array([startX, 0.5, startZ]) # starting position in m
    v0 = fSpeed * numpy.array([cos(fAngle), sin(fAngle), 0])
    
    t = 0
    x = x0
    v = v0

    # Launch rocket to explode height, average 15-22 m
    while x[1] <= numpy.random.randint(15, 30):
        rate(framerate)
        x, v = rk4(x, v, dt, m)
        firework.pos = vector(x[0], x[1], x[2])
        t = t + dt
    firework.visible = False

    return x, v

def createFrags(x: numpy.array, v: numpy.array):
    '''
        Creates fragments to explode!
    '''
    phi = tan(x[0] / x[1])

    frags = [] # create list of frags
    fragM = 0.005 # frag mass
    step = numpy.random.randint(5, 20) # random step for the angle for variation
    fragColor = numpy.random.choice(colors)
    for angle in range(0, 361, step):
        speed = fSpeed / 4
        fx0 = numpy.array([x[0], x[1], x[2]])
        vx = sin(angle) * cos(phi)
        vy = sin(angle) * sin(phi)
        vz = cos(angle)
        fv0 = speed * numpy.array([vx, vy, vz])
        frag = sphere(pos=vector(fx0[0],fx0[1],fx0[2]), radius=0.05,  make_trail=True, color=fragColor)
        frags.append([fx0, fv0, frag]) # add fragment into to the list

    return frags, fragM

def explode(fragments, m): 
    '''
        Explodes the fragments out in all directions
    '''
    t = 0
    while t < 1:
        rate(framerate)
        for frag in fragments:
            frag[0], frag[1] = rk4(frag[0], frag[1], dt, m) # get x and v
            frag[2].pos = vector(frag[0][0], frag[0][1], frag[0][2]) # set x and v on object
        t = t + dt

def launchFirework(x: float, z: float):
    x, v = initialLaunch(x, z)
    frags, m = createFrags(x, v)
    explode(frags, m)

# Constants
framerate = 60
dt = 1.0 / framerate
colors = [color.red, color.orange, color.yellow, color.green, color.magenta, color.purple]

# Environmental Constants
g = 9.8 # gravitational acceleration in m/s**2
rho = 1.225 # air density in kg/m**3

# Firework Constants
m = 2.267 # firework mass in kg ~5lb
r = 0.075 # firework diameter in m
fSpeed = 35 # m/s ~80mph
fAngle = radians(80) # slight angle
cD = 1 # drag coefficient 

# Set-up Scene
scene.forward = vector(-10, -5, -10)
lw = 50
floor = box(pos=vector(0, 0, 0), length=lw, width=lw, height=0.5, color=color.blue)

startingPts = [[0, 0], [20, -20], [-20, 20], [40, -40], [-40, 40]]
for fireworks in range(0, 5):
    launchFirework(startingPts[0][0], startingPts[0][1])
    startingPts.pop(0)
    sleep(2) # wait 2 seconds