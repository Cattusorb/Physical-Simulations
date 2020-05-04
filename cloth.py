from vpython import *
import numpy

g = numpy.array([0, -9.8, 0])

def arr2vec(a: numpy.array) -> vector:
    return vector(a[0], a[1], a[2])

class Particle:
    """Represents a single point mass on the cloth"""
    def __init__(
        self,
        m: float,
        x: numpy.array,
        v: numpy.array = numpy.zeros(3),
        damping: float = 0.9,
        isMovable: bool = True):

        self.m = m
        self.x = x
        self.xold = x
        self.isMovable = isMovable
        self.damping = damping
        self.windForce = numpy.zeros(3)

    def update(self, dt: float):
        """Uses a force model to update the position of the particle"""
        if self.isMovable:
            xold = self.x
            self.x = self.x + (self.x - self.xold) * self.damping + (g + self.windForce / self.m)* dt**2
            self.xold = xold
            self.windForce = numpy.zeros(3)

class Face():
    """Defines a triangular face of our fabric. The entire sheet of fabric can
       can be described as a collection of these faces"""
    def __init__(self, a: Particle, b: Particle, c: Particle):
        """Constructs a new face. The area vector is defined by half the cross
           product (b.x - a.x) cross (c.x - a.x)"""
        self.a = a
        self.b = b
        self.c = c

    def applyWindForce(self, wind: numpy.array):
        """Applies the force of the wind to the particles on this face"""
        area = numpy.cross(self.b.x - self.a.x, self.c.x - self.a.x)
        fWind = wind.dot(area) * area / numpy.linalg.norm(area)
        self.a.windForce += fWind
        self.b.windForce += fWind
        self.c.windForce += fWind


class Constraint():
    """Represents a pair of particles whose distance must be maintained"""
    def __init__(self, a: Particle, b: Particle):
        self.a = a
        self.b = b
        self.restLength = numpy.linalg.norm(self.a.x - self.b.x)

    def enforce(self):
        """Updates the positions of the particles so they satisfy the constraint"""
        dx = self.a.x - self.b.x
        length = numpy.linalg.norm(dx)
        delta = length - self.restLength

        if self.a.isMovable:
            self.a.x = self.a.x - (delta / 2) * (dx / length)

        if self.b.isMovable:
            self.b.x = self.b.x + (delta / 2) * (dx / length)

class Cloth():
    """Represents a sheet of fabric"""
    def __init__(
        self,
        m: float,
        length: float,
        width: float,
        nLength: int,
        nWidth: int,
        x: numpy.array = numpy.zeros(3),
        lengthDirection: numpy.array = numpy.array([1, 0, 0]),
        widthDirection: numpy.array = numpy.array([0, 0, 1]),
        damping: float = 0.9,
        constraintCycles: int = 5):

        self.damping = damping
        self.constraintCycles = constraintCycles

        lhat = lengthDirection / numpy.linalg.norm(lengthDirection)
        what = widthDirection / numpy.linalg.norm(widthDirection)
        dl = lhat * length / (nLength - 1)
        dw = what * width / (nWidth - 1)

        mParticle = m / (nLength * nWidth)
        self.nLength = nLength
        self.nWidth = nWidth
        self.particles = [[None for i in range(nWidth)] for j in range(nLength)]
        x0 = x - length * lhat / 2 - width * what / 2

        for i in range(nLength):
            for j in range(nWidth):
                xParticle = x0 + i * dl + j * dw
                self.particles[i][j] = Particle(mParticle, xParticle, numpy.zeros(3), damping)

        # Create all constraints
        self.constraints = []

        # Nearest neighbors along width
        for i in range(nLength):
            for j in range(nWidth - 1):
                pa = self.particles[i][j]
                pb = self.particles[i][j + 1]
                self.constraints.append(Constraint(pa, pb))

        # Nearest neighbors along length
        for i in range(nLength - 1):
            for j in range(nWidth):
                pa = self.particles[i][j]
                pb = self.particles[i + 1][j]
                self.constraints.append(Constraint(pa, pb))

        # Diagonals
        for i in range(nLength - 1):
            for j in range(nWidth - 1):
                pa = self.particles[i][j]
                pb = self.particles[i + 1][j + 1]
                self.constraints.append(Constraint(pa, pb))

        for i in range(1, nLength):
            for j in range(nWidth - 1):
                pa = self.particles[i][j]
                pb = self.particles[i - 1][j + 1]
                self.constraints.append(Constraint(pa, pb))

        # Secondary neighbors along width
        for i in range(nLength):
            for j in range(nWidth - 2):
                pa = self.particles[i][j]
                pb = self.particles[i][j + 2]
                self.constraints.append(Constraint(pa, pb))

        # Secondary neighbors along length
        for i in range(nLength - 2):
            for j in range(nWidth):
                pa = self.particles[i][j]
                pb = self.particles[i + 2][j]
                self.constraints.append(Constraint(pa, pb))

        # Diagonals
        for i in range(nLength - 2):
            for j in range(nWidth - 2):
                pa = self.particles[i][j]
                pb = self.particles[i + 2][j + 2]
                self.constraints.append(Constraint(pa, pb))

        for i in range(2, nLength):
            for j in range(nWidth - 2):
                pa = self.particles[i][j]
                pb = self.particles[i - 2][j + 2]
                self.constraints.append(Constraint(pa, pb))

        # Create triangular faces for use in calculating wind force
        self.faces = []
        for i in range(nLength - 1):
            for j in range(nWidth - 1):
                pij = self.particles[i][j]
                pip1j = self.particles[i + 1][j]
                pijp1 = self.particles[i][j + 1]
                pip1jp1 = self.particles[i + 1][j + 1]

                self.faces.append(Face(pij, pip1j, pijp1))
                self.faces.append(Face(pip1jp1, pijp1, pip1j))


    def update(self, dt: float, wind: numpy.array = numpy.zeros(3)):
        for face in self.faces:
            face.applyWindForce(wind)

        for i in range(self.nLength):
            for j in range(self.nWidth):
                self.particles[i][j].update(dt)

        for n in range(self.constraintCycles):
            for constraint in self.constraints:
                constraint.enforce()

    def getFrame(self) -> [[vector]]:
        """Returns a 2D array of vpython vectors reflecting the positions of all
           particles on the mesh"""
        particlePositions = [[None for i in range(self.nWidth)] for j in range(self.nLength)]

        for i in range(self.nLength):
            for j in range(self.nWidth):
                particlePositions[i][j] = arr2vec(self.particles[i][j].x)

        return particlePositions

class MeshRendering():
    """Graphical rendering of the fabric as a mesh"""
    def __init__(self, particlePositions: [[vector]], particleSize: float = 0.1):
        self.ni = len(particlePositions)
        self.nj = len(particlePositions[0])

        self.particles = \
            [[sphere(radius=particleSize) for i in range(self.nj)] for j in range(self.ni)]

        self.update(particlePositions)

    def update(self, particlePositions: [[vector]]):
        """Updates the mesh's rendering to match a new set of particle positions"""
        # Update vertex positions
        for i in range(self.ni):
            for j in range(self.nj):
                self.particles[i][j].pos = particlePositions[i][j]