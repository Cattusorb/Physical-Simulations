from types import FunctionType
from vpython import *
import numpy

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
    def __init__(self, m: float, r: float, x: numpy.array, v: numpy.array, axis: numpy.array):
        I = m * r**2 / 2
        super(Cylinder, self).__init__(m, I, x, v)
        self.graphic = cylinder(pos=vector(x[0], x[1], x[2]), radius=r, axis=vector(axis[0], axis[1], axis[2]))

    def render(self):
        '''Updates the object's graphical representation to reflect 
            it's current position'''
        self.graphic.pos = vector(self.x[0], self.x[1], self.x[2])

class Hoop(RigidBody):
    '''A cylindrical solid rigid body'''
    def __init__(self, m: float, r: float, x: numpy.array, v: numpy.array):
        Rthinkness = r/10
        I = m * (4 * (2 * r)**2 + 3 * (Rthinkness/2))
        super(Hoop, self).__init__(m, I, x, v)
        self.graphic = ring(pos=vector(x[0], x[1], x[2]), radius=r, thickness=Rthinkness)

    def render(self):
        '''Updates the object's graphical representation to reflect 
            it's current position'''
        self.graphic.pos = vector(self.x[0], self.x[1], self.x[2])

class Football(RigidBody):
    '''A cylindrical solid rigid body'''
    def __init__(self, m: float, r: float, x: numpy.array, v: numpy.array, l: float, w: float, h: float):
        I = 0
        super(Football, self).__init__(m, I, x, v)
        self.graphic = ellipsoid(pos=vector(x[0], x[1], x[2]), size=(l, h, w), radius=r)

    def render(self):
        '''Updates the object's graphical representation to reflect 
            it's current position'''
        self.graphic.pos = vector(self.x[0], self.x[1], self.x[2])

framerate = 30
dt = 1.0 / framerate

t = 0
body = Ball(1, 3.0, numpy.array([0, 0, 0]), numpy.array([0, 0, 0]))
body.addForce(numpy.array([0, 9.8, 0]), numpy.array([1, 0, 0]))

# moments of inertia are different for each shape, google it
# J (impluse) = integral of F * dt
# change in p = J 

while t < 5:
    rate(framerate)
    body.update(dt)
    body.render()
    t = t + dt
