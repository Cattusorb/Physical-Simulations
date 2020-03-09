from types import FunctionType
from vpython import *
import numpy

'''
Roslyn Parker
Physical Simulations
9 March 2020
'''

# Standalone functions to handle numerics
def rk4(phin: numpy.array, dt: float, f: FunctionType):
    """Finds the state vector dt in the future using an RK4 integration step"""
    k1 = f(phin) * dt
    k2 = f(phin + k1 / 2) * dt
    k3 = f(phin + k2 / 2) * dt
    k4 = f(phin + k3) * dt
    phinp1 = phin + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return phinp1

def assemblePhi(x: numpy.array, p: numpy.array, theta: float, L: float):
    """Assembles the state vector from its component parts"""
    phi = numpy.zeros(8)
    phi[0:3] = x
    phi[3:6] = p
    phi[6] = theta
    phi[7] = L
    return phi

def parsePhi(phi: numpy.array):
    """Disassembles the state vector into its component parts"""
    x = phi[0:3]
    p = phi[3:6]
    theta = phi[6]
    L = phi[7]
    return x, p, theta, L

class RigidBody():
    '''Defines a basic interface of a rigid body which can retate about a single axis'''

    def __init__(self, m: float, I: float, x: numpy.array, v: numpy.array):
        '''Initialize the rigid body's state'''
        self.m = m
        self.I = I # moment of inertia
        self.x = x
        self.p = m * v # p is momentum
        self.theta = 0
        self.L = 0
        self.force = m * numpy.array([0, -9.8, 0])
        self.torque = 0

    def addForce(self, F: numpy.array, r: numpy.array = numpy.array([0,0,0])):
        '''Adds a constant force to the object'''
        self.force = self.force + F
        self.torque = r[0] * F[1] - r[1] * F[0]

    def addImpluse(self, J: numpy.array, r: numpy.array):
        '''Adds an impulse to the object'''
        # will be used for football problem
        self.p = self.p + J 
        self.L = self.L + r[0] * J[1] - r[1] * J[0]

    def update(self, dt: float):
       '''Updates the object's state to a time dt in the future'''
       phin = assemblePhi(self.x, self.p, self.theta, self.L)
       phinp1 = rk4(phin, dt, self.ddt)
       self.x, self.p, self.theta, self.L = parsePhi(phinp1)

    def ddt(self, phi: numpy.array):
        x, p, theta, L = parsePhi(phi)

        dxdt = p / self.m # dx/dt = p/m
        dpdt = self.force # dp/dt = F
        dthetadt = L / self.I # dtheta/dt = L / I
        dLdt = self.torque # dL/dt = torque

        return assemblePhi(dxdt, dpdt, dthetadt, dLdt)

class Ball(RigidBody):
    '''A sphere-like rigid body'''
    def __init__(self, m: float, r: float, x: numpy.array, v: numpy.array):
        I = 2 * m * r**2 / 5 
        super(Ball, self).__init__(m, I, x, v)
        self.graphic = sphere(pos=vector(x[0], x[1], x[2]), radius=r, axis=vector(0, 0, 1), texture=textures.stones)

    def render(self):
        '''Updates the object's graphical representation to reflect 
            it's current position'''
        self.graphic.pos = vector(self.x[0], self.x[1], self.x[2])
        self.graphic.up = vector(cos(self.theta), sin(self.theta), 0)

class Cylinder(RigidBody):
    '''A cylindrical solid rigid body'''
    def __init__(self, m: float, r: float, x: numpy.array, v: numpy.array):
        I = m * r**2 / 2
        super(Cylinder, self).__init__(m, I, x, v)
        self.graphic = cylinder(pos=vector(x[0], x[1], x[2]), radius=r, axis=vector(0, 0, 1), texture=textures.stones)

    def render(self):
        '''Updates the object's graphical representation to reflect 
            it's current position'''
        self.graphic.pos = vector(self.x[0], self.x[1], self.x[2])
        self.graphic.up = vector(cos(self.theta), sin(self.theta), 0)

class Hoop(RigidBody):
    '''A cylindrical solid rigid body'''
    def __init__(self, m: float, r: float, x: numpy.array, v: numpy.array):
        a = r/10
        I = m * (4 * (2 * r)**2 + 3 * (a/2))
        super(Hoop, self).__init__(m, I, x, v)
        self.graphic = ring(pos=vector(x[0], x[1], x[2]), radius=r, thickness=a, axis=vector(0, 0, 1), texture=textures.stones)

    def render(self):
        '''Updates the object's graphical representation to reflect 
            it's current position'''
        self.graphic.pos = vector(self.x[0], self.x[1], self.x[2])
        self.graphic.up = vector(cos(self.theta), sin(self.theta), 0)

class Football(RigidBody):
    '''A cylindrical solid rigid body'''
    def __init__(self, m: float, r: float, x: numpy.array, v: numpy.array, l: float):
        b = r * 2
        a = r / 10
        I = (m *  b**2 / 5) * (1 + a**2 / b**2)
        super(Football, self).__init__(m, I, x, v)
        self.graphic = ellipsoid(pos=vector(x[0], x[1], x[2]), length=l, radius=r, texture=textures.metal, make_trail=True)

    def render(self):
        '''Updates the object's graphical representation to reflect 
            it's current position'''
        self.graphic.pos = vector(self.x[0], self.x[1], self.x[2])
        self.graphic.up = vector(cos(self.theta), sin(self.theta), 0)

framerate = 30
dt = 1.0 / framerate

t = 0

# Set up camera
scene.forward = vector(1, -1, -5)
scene.range = 20
# Floor
floorWidth = 50
floorLength = 100
floorThickness = 0.5
floor = box(center=vector(0,0,0), length=floorLength, width=floorWidth, 
height=floorThickness, texture=textures.gravel)

#ball = Ball(1, 3.0, numpy.array([0, 0, 0]), numpy.array([0, 0, 0]))
#ball.addForce(numpy.array([0, 9.8, 0]), numpy.array([3, 0, 0]))

#cyl = Cylinder(1, 3.0, numpy.array([0, 0, 5]), numpy.array([0, 0, 0]))
#cyl.addForce(numpy.array([0, 9.8, 0]), numpy.array([3, 0, 0]))

#hoop = Hoop(1, 3.0, numpy.array([0, 0, 10]), numpy.array([0, 0, 0]))
#hoop.addForce(numpy.array([0, 9.8, 0]), numpy.array([3, 0, 0]))

'''
When running the simulation of the three rotating objects,
the ring looks like it is rotating the slowest and the cylinder
looks like it is rotating the fastest, but the ball also looks like
it is rotating at about the same speed; maybe slightly slower. 

This makes sense because the ball and cylinder are solid objects and the 
ring is somewhat like a bicycle wheel. I don't know for 
sure why the ring looks slower, but when a bike wheel 
is spinning the inside part of the wheel is rotating 
slower than the tire part. This is because there is more 
surface to be covered on the outer edge of the wheel. 

I made the textures all the same for the three objects because
when I first ran it with different textures for each, it was
hard to tell which one was rotating faster. The textures are still
not uniform for all objects, which may make a difference in 
judgement, but overall I feel like the cylinder is rotating
the fastest. 
'''

football = Football(0.415, 0.17, numpy.array([0, 0, 0]), numpy.array([0, 0, 0]), 0.283)
football.addImpluse(numpy.array([8, 8, 0]), numpy.array([0, radians(40), 0]))

'''
Do you get a realistic range for your kick and a realistic tumbling action?
I feel like the result I got from running this simulation 
was a little bit weird. At first, it goes up at an angle, then
that angle decreases slightly, but it is still a straight line. Then 
gravity takes over and it curves off back down to the ground. 
This simulation looks very weird because it's 2 lines with a curve
on the end. It doesn't feel as smooth as when watching an 
actual football be kicked into the air. 

I cannot really see the tumbling action that goes on 
in this simulation because it's pretty small and there is
quite the distance for the football to go. 
'''

'''while t <= 5:
    rate(framerate)
    ball.update(dt)
    ball.render()
    cyl.update(dt)
    cyl.render()
    hoop.update(dt)
    hoop.render()
'''
while football.x[1] >= 0: 
    rate(framerate)
    football.update(dt)
    football.render()
    t = t + dt
