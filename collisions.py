from rigidbody import *

class World:
    """Represents our entire simulation"""
    def __init__(self):
        self.bodies = [] # List of RigidBodies in the simulation
        self.boundaries = [] # List of Boundaries in the simulation
    
    def update(self, dt):
        """Advances all objects dt into the future using our integrator,
           performs collision detection and resolution if necessary"""
        for body in self.bodies:
            body.update(dt)

        contacts = []
        for body in self.bodies: 
            for boundary in self.boundaries: 
                collisionDetected, contact = body.collidesWithBoundary(boundary)

                if collisionDetected: 
                    contacts.append(contact)

        for ia in range(len(self.bodies) - 1):
            for ib in range(ia + 1, len(self.bodies)):
                a = self.bodies[ia]
                b = self.bodies[ib]

                collisionDetected, contact = a.collidesWith(b)
                if collisionDetected:
                    contacts.append(contact)
        
        for contact in contacts:
            if contact.isColliding():
                contact.resolveCollision(0.85)
            else:
                contacts.remove(contact)
    
    def render(self):
        """Updates the graphical representation of all objects to reflect their
           new positions"""
        for body in self.bodies:
            body.render()

    def removeObjectsBelow(self, y):
        """Removes all objects below a certain height to keep the simulation
           from getting bogged down"""
        for i in range(len(self.bodies) - 1, -1, -1):
            if self.bodies[i].x[1] < y:
                self.bodies[i].graphic.visible = False
                del self.bodies[i]


framerate = 60
dt = 1 / framerate

world = World()
#world.bodies.append(Sphere(0.1, 0.2, numpy.array([3, 10, 0]))) # ball 1
#world.bodies.append(Sphere(0.1, 0.2, numpy.array([3, 12, 0]))) # ball 2
#world.boundaries.append(Boundary(numpy.array([0, 0, 0]), 10, 10)) # ball surface
world.boundaries.append(Boundary(numpy.array([-4, 0, 0]), 7, 7, numpy.array([0, 0, 1]), numpy.array([1, 0.8, 0])))
world.boundaries.append(Boundary(numpy.array([4, 0, 0]), 7, 7, numpy.array([0, 0, 1]), numpy.array([1, -0.8, 0])))


t = 0
spawn = 0
while t < 5:
    rate(framerate)
    world.update(dt)
    world.render()
    world.removeObjectsBelow(-2)
    if spawn % 5 == 0: 
        world.bodies.append(Sphere(0.1, 0.2, numpy.array([numpy.random.randint(-5, 5), 5, numpy.random.randint(-5, 5)])))
    spawn += 1
    t += dt