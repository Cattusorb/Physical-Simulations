from __future__ import annotations
from typing import Tuple, Union
from types import FunctionType
from vpython import *
import numpy

g = 9.8 # Gravitational acceleration in m/s^2

# Standalone vector function
def star(a: numpy.array) -> numpy.array:
    return numpy.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])

# Standalone functions to handle numerics
def rk4(phin: numpy.array, dt: float, f: FunctionType) -> numpy.array:
    """Finds the state vector dt in the future using an RK4 integration step"""
    k1 = f(phin) * dt
    k2 = f(phin + k1 / 2) * dt
    k3 = f(phin + k2 / 2) * dt
    k4 = f(phin + k3) * dt

    phinp1 = phin + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return phinp1

def assemblePhi(x: numpy.array, p: numpy.array, R: numpy.array, L: numpy.array) -> numpy.array:
    """Assembles the state vector from its component parts"""
    phi = numpy.zeros(18)
    phi[0:3] = x
    phi[3:6] = p
    phi[6:15] = R.flatten()
    phi[15:18] = L

    return phi

def parsePhi(phi: numpy.array) -> Tuple[numpy.array, numpy.array, numpy.array, numpy.array]:
    """Disassembles the state vector into its component parts"""
    x = phi[0:3]
    p = phi[3:6]
    R = phi[6:15].reshape(3, 3)
    L = phi[15:18]
    return x, p, R, L
  
class RigidBody:
    """Base class for a 3D rigid body, defines the basic methods for specifying
       and updating position and orientation in 3D space"""

    def __init__(self,
    m: float,
    Ibody: numpy.array,
    x: numpy.array = numpy.array([0., 0., 0.]),
    v: numpy.array = numpy.array([0., 0., 0.]),
    axis: numpy.array = numpy.array([1, 0, 0]),
    up: numpy.array = numpy.array([0, 1, 0]),
    omega: numpy.array = numpy.array([0., 0., 0.])):
        # Initialize object properties
        # Mass
        self.m = m

        # Moment of inertia
        self.Ibody = Ibody
        self.Ibodyinv = numpy.linalg.inv(Ibody)

        # Initialize state information
        # Position
        self.x = x

        # Momentum
        self.p = m * v

        # Rotation
        # Normalize axis and up
        axis = axis / numpy.linalg.norm(axis)
        up = up / numpy.linalg.norm(up)

        self.R = numpy.zeros([3, 3])
        self.R[0:3][0] = axis
        self.R[0:3][1] = up
        self.R[0:3][2] = numpy.cross(axis, up)

        # Angular momentum
        self.L = self.R.dot(Ibody.dot(self.R.transpose())).dot(omega)

        # Initialize world-centered force and torque accumulator
        self.force = self.m * numpy.array([0, -g, 0])

        # Gravity induces no torque about the center of mass so initial torque is 0
        self.torque = numpy.zeros(3)

    def v(self) -> numpy.array:
        """Returns the object's velocity in world coordinates"""
        return self.p / self.m

    def applyImpulse(self, J: numpy.array, r: numpy.array):
        """Applies an impulse in the world frame to the body at optional
           position in the world frame"""
        self.p += J
        self.L += numpy.cross(r - self.x, J)

    def ddt(self, phi: numpy.array) -> numpy.array:
        """Returns a vector with time derivatives of all state variables"""
        # Get individual components from phi
        x, p, R, L = parsePhi(phi)

        # dx/dt = p / m
        dxdt = p / self.m

        # dp/dt = force
        dpdt = self.force

        # dR/dt = omega* * R
        omega = R.dot(self.Ibodyinv.dot(R.transpose())).dot(L)
        dRdt = star(omega).dot(R)

        #dL/dt = torque
        dLdt = self.torque

        # Assemble return vector
        return assemblePhi(dxdt, dpdt, dRdt, dLdt)
    
    def update(self, dt: float):
        """Updates the object's state to a time dt in the future"""
        
        phin = assemblePhi(self.x, self.p, self.R, self.L)
        phinp1 = rk4(phin, dt, self.ddt)
        self.x, self.p, self.R, self.L = parsePhi(phinp1)

    def collidesWith(self, other: RigidBody) -> Tuple[bool, Contact]:
        """Determines if this body collides with another and returns a contact
           point if it does"""
        raise NotImplementedError("Need to implement for this RigidBody subtype")

    def collidesWithBoundary(self, boundary: Boundary) -> Tuple[bool, Contact]:
        """Determines if this body collides with a boundary and returns a
           contact point if it does"""
        raise NotImplementedError("Need to implement for this RigidBody subtype")

class Sphere(RigidBody):
    """A solid sphere"""
    def __init__(self,
    m: float, r: float,
    x: numpy.array = numpy.zeros(3),
    v: numpy.array = numpy.zeros(3)):
        self.r = r
        I = (2 * m * r**2 / 5) * numpy.identity(3)
        super(Sphere, self).__init__(m, I, x, v)

        self.graphic = sphere(pos=vector(x[0], x[1], x[2]), radius=r)
        return
    
    def render(self):
        # Set axis from rotation matrix
        xhat = numpy.array([self.graphic.radius, 0, 0])
        axis = self.R.dot(xhat)

        self.graphic.pos = vector(self.x[0], self.x[1], self.x[2])
        self.graphic.axis = vector(axis[0], axis[1], axis[2])

    def collidesWith(self, other: RigidBody) -> Tuple[bool, Contact]:
        # Find closest point on ther object
        dx = self.x - other.x
        magdx = numpy.linalg.norm(dx)

        # If this point lies within the sphere's radius we have a collision
        if magdx < self.r + other.r:
            # Generate the contact
            # Contact normal is just the vector pointing from this sphere's
            # center to the closest point on the other object, divided by its 
            # magnitude
            nhat = dx / magdx

            # Contact position will be defined as the midpoint of the two
            # object's overlap
            xContact = ((self.x - self.r * nhat) + (other.x + other.r * nhat)) / 2

            return True, Contact(self, other, xContact, nhat)
        else:
            return False, None

    def collidesWithBoundary(self, boundary: Boundary) -> Tuple[bool, Contact]:
        # Get object's displacement from plane's center point in plane-centered
        # coordinates
        dx = self.x - boundary.x
        dl = abs(dx.dot(boundary.axis)) - self.r
        dh = dx.dot(boundary.normal)
        dw = abs(dx.dot(numpy.cross(boundary.axis, boundary.normal))) - self.r
        
        if dh < self.r and dl < boundary.length / 2 and dw < boundary.width / 2:
            contact = Contact(self, boundary, self.x - self.r * boundary.normal, boundary.normal)
            return True, contact
        else:
            return False, None


class Boundary():
    """Represents an immovable plane"""
    def __init__(self,
    x: numpy.array,
    length: float,
    width: float,
    axis: numpy.array = numpy.array([1, 0, 0]),
    normal: numpy.array = numpy.array([0, 1, 0])):
        self.x = x # Center position

        self.axis = axis / numpy.linalg.norm(axis) # Axis vector
        self.normal = normal / numpy.linalg.norm(normal) # Normal vector

        self.length = length # Extent along axis
        self.width = width # Extent along axis cross normal

        self.m = numpy.inf

        # Default thickness
        depth = 0.1

        # Set the graphical representation
        self.graphic = box(
            pos=vector(x[0], x[1] - depth / 2, x[2]),
            size=vector(length, depth, width),
            axis=length*vector(self.axis[0], self.axis[1], self.axis[2]),
            up=vector(self.normal[0], self.normal[1], self.normal[2]),
            color=color.blue)

    def v(self) -> numpy.array:
        return numpy.zeros(3)

    def applyImpulse(self, J: numpy.array, r: numpy.array):
        return
        
class Contact:
    """Represents a contact point for a collision betweeen two rigid bodies"""
    def __init__(self, a: RigidBody, b: Union[RigidBody, Boundary], x: numpy.array, normal: numpy.array):
        self.a = a # First body
        self.b = b # Second body or boundary
        self.x = x # Position in world coordinates where contact occurs
        self.normal = normal # Contact normal vector, points outwards from B

    def isColliding(self) -> bool:
        """Determines if the two bodies involved in the contact are actively
           colliding."""
        COLLISION_THRESHOLD_SPEED = 0.01
        return (self.a.v() - self.b.v()).dot(self.normal) < -COLLISION_THRESHOLD_SPEED
        
    def resolveCollision(self, epsilon: float):
        """Performs collision resolution for a contact given a coefficient of
           restitution"""
        
        # Calculate the impulse that should be applied
        # First, get relative velocity of the two centers of mass
        vclose = (self.a.v() - self.b.v()).dot(self.normal)

        # Next, get terms in the denominator
        invMassA = 1 / self.a.m
        invMassB = 1 / self.b.m

        # Put it all together
        magj = -(1 + epsilon) * vclose / (invMassA + invMassB)
        j = magj * self.normal

        self.a.applyImpulse(j, self.x)
        self.b.applyImpulse(-j, self.x)