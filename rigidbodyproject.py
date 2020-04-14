from __future__ import annotations
from typing import Tuple, Union
from types import FunctionType
from vpython import *
import numpy

g = 0 # Gravitational acceleration in m/s^2

# Define unit vectors for later reference
xhat = numpy.array([1, 0, 0])
yhat = numpy.array([0, 1, 0])
zhat = numpy.array([0, 0, 1])

def arr2vec(a: numpy.array) -> vector:
    """Converts a numpy array into a VPython vector"""
    return vector(a[0], a[1], a[2])

def normalize(v: numpy.array) -> numpy.array:
    """Returns a unit vector in the direction of v"""
    return v / numpy.linalg.norm(v)

def star(a: numpy.array) -> numpy.array:
    """Returns the 'star' matrix for a vector"""
    return numpy.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])

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
        axis = normalize(axis)
        up = normalize(up)
        self.R = numpy.array([axis, up, numpy.cross(axis, up)]).transpose()

        # Angular momentum
        self.L = self.R.dot(Ibody.dot(self.R.transpose())).dot(omega)

        # Initialize world-centered force and torque accumulator
        self.force = self.m * numpy.array([0, -g, 0])

        # Gravity induces no torque about the center of mass so initial torque is 0
        self.torque = numpy.zeros(3)

    def v(self) -> numpy.array:
        """Returns the object's velocity in world coordinates"""
        return self.p / self.m

    def Iinv(self) -> numpy.array:
        """Returns the current moment of inertia tensor in world coordinates"""
        return self.R.dot(self.Ibodyinv.dot(self.R.transpose()))

    def omega(self) -> numpy.array:
        """Returns the current angular velocity in world coordinates"""
        return self.Iinv().dot(self.L)

    def worldPositionOf(self, x: numpy.array) -> numpy.array:
        """Returns the world coordinates of a postion in body-centered
           coordinates"""
        return self.R.dot(x) + self.x

    def bodyPositionOf(self, x: numpy.array) -> numpy.array:
        """Returns the body-centered coordinates of a position in world
           coordinates"""
        return self.R.transpose().dot(x - self.x)

    def velocityOf(self, x: numpy.array) -> numpy.array:
        """Returns the velocity of a point on the object in given in
           world-centered coordinates"""
        return self.v() + numpy.cross(self.omega(), x - self.x)

    def transform(self, v: numpy.array) -> numpy.array:
        """Transforms a vector from body coordinates to world coordinates"""
        return self.R.dot(v)

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

    def distanceFrom(self, other: RigidBody) -> float:
        """Returns the straight line distance between the centers of mass of
           this object and another"""
        return numpy.linalg.norm(self.x - other.x)

    def projectionOn(self, axis: numpy.array):
        """Returns the extent of the object along a given axis"""
        raise NotImplementedError("Need to implement for this RigidBody subtype")

    def closestPointTo(self, other: RigidBody) -> numpy.array:
        """Returns the point on this object closest to another object's center
           of mass"""
        raise NotImplementedError("Need to implement for this RigidBody subtype")

    def boundingRadius(self) -> float:
        """Returns the radius of the bounding sphere for this object"""
        raise NotImplementedError("Need to implement for this RigidBody subtype")

    def collidesWith(self, other: RigidBody) -> Tuple[bool, Contact]:
        """Determines if this body collides with another and returns a contact
           point if it does"""
        raise NotImplementedError("Need to implement for this RigidBody subtype")

    def collidesWithBlock(self, other: Block) -> Tuple[bool, Contact]:
        """Determines if this body collides with a Block and returns a contact
           point if it does"""
        raise NotImplementedError("Need to implement for this RigidBody subtype")

    def collidesWithBoundary(self, other: Boundary) -> Tuple[bool, Contact]:
        """Determines if this body coolides with a Boundary and returns all
           contact points if it does"""
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

        self.graphic = sphere(pos=arr2vec(x), radius=r)
        return
    
    def render(self):
        # Set axis from rotation matrix
        axis = self.R.dot(numpy.array([self.graphic.radius, 0, 0]))

        self.graphic.pos = arr2vec(self.x)
        self.graphic.axis = arr2vec(axis)

    def projectionOn(self, axis: numpy.array):
        return self.r

    def closestPointTo(self, other: RigidBody) -> numpy.array:
        # Get vector pointing from this object to the other
        dx = other.x - self.x

        # Closest point is this sphere's center of mass position plus a distance
        # r along this vector
        return self.x + self.r * dx / numpy.linalg.norm(dx)        

    def boundingRadius(self) -> float:
        return self.r

    def collidesWith(self, other: RigidBody) -> Tuple[bool, Contact]:
        # Find closest point on ther object
        closestPoint = other.closestPointTo(self)
        dx = self.x - closestPoint
        magdx = numpy.linalg.norm(dx)
        overlap = self.r - magdx
        
        # If this point lies within the sphere's radius we have a collision
        if overlap > 0:
            # Generate the contact
            # Contact normal is just the vector pointing from this sphere's
            # center to the closest point on the other object, divided by its 
            # magnitude
            nhat = normalize(dx)

            # Contact position will be defined as the midpoint of the two
            # object's overlap
            xContact = closestPoint + overlap * nhat / 2

            return True, Contact(self, other, xContact, nhat)
        else:
            return False, None

    def collidesWithBlock(self, other: Block) -> Tuple[bool, Contact]:
        # For a sphere, a collision with a block is the same as a collision with
        # any other object
        return self.collidesWith(other)

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

############ HELPER METHODS FOR BLOCK #####################

def sat(a: Block, b: Block) -> Tuple[bool, Contact]:
    """Performs the separating axis test to determine if two blocks overlap,
       returning the Contact"""
    # We need to check for overlapping projections along the following axes:
    #   xhat, yhat, zhat for body a
    #   xhat, yhat, zhat for body b
    #   the 9 cross product combinations of these

    dx = a.x - b.x

    minOverlap = inf

    edgeCollision = False

    # Do all three faces of a
    for aAxis in [xhat, yhat, zhat]:
        axis = a.transform(aAxis)
        distanceBetween = abs(dx.dot(axis))
        aProjection = a.projectionOn(axis)
        bProjection = b.projectionOn(axis)
        overlap = aProjection / 2 + bProjection / 2 - distanceBetween

        if overlap < 0:
            return False, None
        elif overlap < minOverlap:
            minOverlap = overlap
            faceNormal = axis
            faceBlock = a
            vertexBlock = b

    # Now do three faces of b
    for bAxis in [xhat, yhat, zhat]:
        axis = b.transform(bAxis)
        distanceBetween = abs(dx.dot(axis))
        aProjection = a.projectionOn(axis)
        bProjection = b.projectionOn(axis)
        overlap = aProjection / 2 + bProjection / 2 - distanceBetween

        if overlap < 0:
            return False, None
        elif overlap < minOverlap:
            minOverlap = overlap
            faceNormal = axis
            faceBlock = b
            vertexBlock = a

    # Now do nine edge cross product combinations
    for aEdge in [xhat, yhat, zhat]:
        for bEdge in [xhat, yhat, zhat]:
            axis = numpy.cross(a.transform(aEdge), b.transform(bEdge))

            if numpy.linalg.norm(axis) < 0.001:
                continue

            axis = normalize(axis)
            distanceBetween = abs(dx.dot(axis))
            aProjection = a.projectionOn(axis)
            bProjection = b.projectionOn(axis)
            overlap = aProjection / 2 + bProjection / 2 - distanceBetween

            if overlap < 0:
                return False, None
            elif overlap < minOverlap:
                minOverlap = overlap
                aEdgeBody = aEdge
                bEdgeBody = bEdge
                edgeCollision = True

    # If we get here, we have a collision
    if edgeCollision:
        return True, edgeEdgeContact(a, b, aEdgeBody, bEdgeBody, minOverlap)
    else:
        return True, vertexFaceContact(vertexBlock, faceBlock, minOverlap, faceNormal)

def vertexFaceContact(
    vertexBlock: Block,
    faceBlock: Block,
    minOverlap: float,
    faceNormal: numpy.array) -> Contact:
    """Returns the contact for an edge-edge collision"""

    # See if we need to flip the direction of the contact axis
    toCenter = vertexBlock.x - faceBlock.x
    if toCenter.dot(faceNormal) < 0:
        faceNormal = -faceNormal

    # Find which vertex on a is in contact with b
    vertex = vertexBlock.size / 2

    if vertexBlock.transform(xhat).dot(faceNormal) > 0:
        vertex[0] = -vertex[0]

    if vertexBlock.transform(yhat).dot(faceNormal) > 0:
        vertex[1] = -vertex[1]

    if vertexBlock.transform(zhat).dot(faceNormal) > 0:
        vertex[2] = -vertex[2]

    xContact = vertexBlock.worldPositionOf(vertex)
    return Contact(vertexBlock, faceBlock, xContact, faceNormal)

def edgeEdgeContact(
    a: Block,
    b: Block,
    aEdgeBody: numpy.array,
    bEdgeBody: numpy.array,
    minOverlap: float) -> Contact:
    """Returns the contact for an edge-edge collision"""

    # Get contact direction
    aEdgeWorld = normalize(a.transform(aEdgeBody))
    bEdgeWorld = normalize(b.transform(bEdgeBody))
    contactDirection = normalize(numpy.cross(aEdgeWorld, bEdgeWorld))

    # See if we need to flip the direction of the contact axis
    toCenter = a.x - b.x
    if toCenter.dot(contactDirection) < 0:
        contactDirection = -contactDirection

    # Pick points midway on colliding edges
    aEdgePoint = a.size / 2
    bEdgePoint = b.size / 2

    axes = [xhat, yhat, zhat]
    for i in range(len(axes)):
        axis = axes[i]

        if numpy.array_equal(axis, aEdgeBody):
            aEdgePoint[i] = 0
        elif a.transform(axis).dot(contactDirection) > 0:
            aEdgePoint[i] = -aEdgePoint[i]


        if numpy.array_equal(axis, bEdgeBody):
            bEdgePoint[i] = 0
        elif b.transform(axis).dot(contactDirection) < 0:
            bEdgePoint[i] = -bEdgePoint[i]

    # Now do math to get points on edge closest to collision
    aPoint = a.worldPositionOf(aEdgePoint)
    bPoint = b.worldPositionOf(bEdgePoint)

    toSt = aPoint - bPoint

    dpStaA = aEdgeWorld.dot(toSt)
    dpStaB = bEdgeWorld.dot(toSt)

    smA = numpy.sum(aEdgeWorld**2)
    smB = numpy.sum(bEdgeWorld**2)

    dpEdges = aEdgeWorld.dot(bEdgeWorld)

    denom = smA * smB - dpEdges**2
    aDist = (dpEdges * dpStaB - smB * dpStaA) / denom
    bDist = (smA * dpStaB - dpEdges * dpStaA) / denom
    
    nearestPointA = aPoint + aEdgeWorld * aDist
    nearestPointB = bPoint + bEdgeWorld * bDist
    xContact = (nearestPointA + nearestPointB) / 2

    return Contact(a, b, xContact, contactDirection)

class Block(RigidBody):
    """A solid block (three-dimensional rectangular prism)"""
    def __init__(self,
    m: float,
    size: numpy.array,
    x: numpy.array = numpy.zeros(3),
    v: numpy.array = numpy.zeros(3),
    axis: numpy.array = numpy.array([1, 0, 0]),
    up: numpy.array = numpy.array([0, 1, 0]),
    omega: numpy.array = numpy.zeros(3)):
        # Define moment of inertia tensor
        dx = size[0]
        dy = size[1]
        dz = size[2]
        Ixx = (m / 12) * (dy**2 + dz**2)
        Iyy = (m / 12) * (dx**2 + dz**2)
        Izz = (m / 12) * (dx**2 + dy**2)
        I = numpy.array([[Ixx, 0, 0], [0, Iyy, 00], [0, 0, Izz]])

        # Normalize axis and up
        axis = normalize(axis)
        up = normalize(up)

        super(Block, self).__init__(m, I, x, v, axis, up, omega)
        self.size = size
        self.graphic = box(pos=arr2vec(x), length=dx, height=dy, width=dz, axis=arr2vec(axis), up=arr2vec(up))
        return
    
    def render(self):
        # Set axis from rotation matrix
        axis = self.R.dot(numpy.array([self.size[0], 0, 0]))
        up = self.R.dot(numpy.array([0, 1, 0]))

        self.graphic.pos = vector(self.x[0], self.x[1], self.x[2])
        self.graphic.axis = vector(axis[0], axis[1], axis[2])
        self.graphic.up = vector(up[0], up[1], up[2])

    def projectionOn(self, axis: numpy.array) -> float:
        """Returns the projection (extent) of the block along an axis"""
        # Smallest coordinate value projected onto the axis
        minProjection = inf
        maxProjection = -inf

        # Get half width of object in body coordinates
        dx = self.size[0] / 2
        dy = self.size[1] / 2
        dz = self.size[2] / 2

        for x in [-dx, dx]:
            for y in [-dy, dy]:
                for z in [-dz, dz]:
                    xvBody = numpy.array([x, y, z])
                    xvWorld = self.worldPositionOf(xvBody)
                    projection = xvWorld.dot(axis)

                    if projection < minProjection:
                        minProjection = projection

                    if projection > maxProjection:
                        maxProjection = projection

        return maxProjection - minProjection

    def closestPointTo(self, other: RigidBody) -> numpy.array:
        otherBody = self.bodyPositionOf(other.x)
        xclose = max(-self.size[0] / 2, min(self.size[0] / 2, otherBody[0]))
        yclose = max(-self.size[1] / 2, min(self.size[1] / 2, otherBody[1]))
        zclose = max(-self.size[2] / 2, min(self.size[2] / 2, otherBody[2]))

        close = numpy.array([xclose, yclose, zclose])
        return self.worldPositionOf(close)

    def boundingRadius(self) -> float:
        # Code for this function added by Roslyn Parker
        x = self.x[0] + self.size[0]
        y = self.x[1] + self.size[1]
        z = self.x[2] + self.size[2]

        r = (self.x[0] - x)**2 + (self.x[1] - y)**2 + (self.x[2] - z)**2
        return r

    def collidesWith(self, other: RigidBody) -> Tuple[bool, Contact]:
        return other.collidesWithBlock(self)
    
    def collidesWithBlock(self, other: Block) -> Tuple[bool, Contact]:
        return sat(self, other)

    def collidesWithBoundary(self, boundary: Boundary) -> Tuple[bool, [Contact]]:
        # Loop over all vertices (corners) of the block. When we find that one
        # of the vertices has made contact with the boundary, we create the
        # Contact for that point and return it.

        # Find the formula that tells you when a corner on the block has hit
        # the boundary

        # Loop over all vertices of the block
        # When we find a collision, we return contact
        # Get half width of object in body coordinates
        dx = self.size[0] / 2
        dy = self.size[1] / 2
        dz = self.size[2] / 2

        contacts = [] # Added this list

        for x in [-dx, dx]:
            for y in [-dy, dy]:
                for z in [-dz, dz]:
                    xvBody = numpy.array([x, y, z])
                    xvWorld = self.worldPositionOf(xvBody)

                    if (xvWorld - boundary.x).dot(boundary.normal) < 0:
                        contacts.append(Contact(self, boundary, xvWorld, boundary.normal))
        # Added this code
        if len(contacts) > 0:
            return True, contacts
        
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

    def Iinv(self) -> numpy.array:
        """Returns the current moment of inertia tensor in world coordinates"""
        return numpy.zeros([3, 3])

    def omega(self) -> numpy.array:
        """Returns the current angular velocity in world coordinates"""
        return numpy.zeros(3)

    def velocityOf(self, x: numpy.array) -> numpy.array:
        """Returns the velocity of a point on the object in given in
           world-centered coordinates"""
        return numpy.zeros(3)

    def applyImpulse(self, J: numpy.array, r: numpy.array):
        return
        
class Contact:
    """Represents a contact point for a collision betweeen two rigid bodies"""

    # Threshold velocity 
    COLLIDING_THRESHOLD_SPEED = 0.001

    def __init__(self, a: RigidBody, b: Union[RigidBody, Boundary], x: numpy.array, normal: numpy.array):
        self.a = a # First body
        self.b = b # Second body or boundary
        self.x = x # Position in world coordinates where contact occurs
        self.normal = normal # Contact normal vector, points outwards from B

    def isColliding(self) -> bool:
        """Determines if the two bodies involved in the contact are actively
           colliding."""
        xadot = self.a.velocityOf(self.x)
        xbdot = self.b.velocityOf(self.x)
        vrel = xadot - xbdot
        vclose = vrel.dot(self.normal)

        return vclose < 0
        
    def resolveCollision(self, epsilon: float):
        """Performs collision resolution for a contact given a coefficient of
           restitution"""
        
        # Calculate the impulse that should be applied
        # First, get relative velocity of the two centers of mass
        xadot = self.a.velocityOf(self.x)
        xbdot = self.b.velocityOf(self.x)
        vrel = self.normal.dot(xadot - xbdot)

        # Next, get terms in the denominator
        # Inverse masses
        invMassA = 1 / self.a.m
        invMassB = 1 / self.b.m

        # Angular terms
        n = self.normal
        ra = self.x - self.a.x
        rb = self.x - self.b.x
        angTermA = n.dot(numpy.cross(self.a.Iinv().dot(numpy.cross(ra, n)), ra))
        angTermB = n.dot(numpy.cross(self.b.Iinv().dot(numpy.cross(rb, n)), rb))

        # Put it all together
        magj = -(1 + epsilon) * vrel / (invMassA + invMassB + angTermA + angTermB)
        j = magj * self.normal

        self.a.applyImpulse(j, self.x)
        self.b.applyImpulse(-j, self.x)

class World:
    """Represents the entire simulation"""
    def __init__(self):
        self.bodies = [] # List of bodies in the simulation
        self.boundaries = [] # List of boundaries
        self.contacts = [] # List of contacts between bodies

    def generateContacts(self):
        """Generates contacts between objects in the simulation"""
        # Check for boundary collisions
        for body in self.bodies:
            for boundary in self.boundaries:
                collisionDetected, ncontacts = body.collidesWithBoundary(boundary)

                if collisionDetected:
                    for contact in ncontacts: 
                        self.contacts.append(contact)


        # Collision detection
        for ia in range(len(self.bodies)):
            for ib in range(ia+1, len(self.bodies)):
                a = self.bodies[ia]
                b = self.bodies[ib]
                
                # Code added here by Roslyn Parker
                # If the bounding radius of the block a 
                # is overlapping with the bounding radius of 
                # block b then collisionDetected should be aware

                collisionDetected, contact = a.collidesWith(b)

                if collisionDetected:
                    self.contacts.append(contact)

    def resolveCollisions(self):
        """Resolves any collisions identified in the simulation"""
        hadCollision = True
        EPSILON = 0.4

        while hadCollision:
            hadCollision = False

            for contact in self.contacts:
                if contact.isColliding():
                    # Try to resolve the collision
                    contact.resolveCollision(EPSILON)
                    hadCollision = True
                else:
                    # Collision resolved, remove from contacts list
                    self.contacts.remove(contact)

    
    def update(self, dt: float):
        """Propagates all objects to their new positions, then identifies and
           resolves any collisions"""
        # Update positions of bodies
        for body in self.bodies:
            body.update(dt)

        # Generate contacts between colliding bodies
        self.generateContacts()

        # Resolve collisions
        self.resolveCollisions()

    def render(self):
        """Updates the graphical representation of all objects in the
         simulation"""
        for body in self.bodies:
            body.render()

    def removeObjectsBelow(self, y):
        """Removes all objects below a certain height to keep the simulation
           from getting bogged down"""
        for i in range(len(self.bodies) - 1, -1, -1):
            if self.bodies[i].x[1] < y:
                self.bodies[i].graphic.visible = False
                del self.bodies[i]