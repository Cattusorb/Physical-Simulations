from rigidbodycollisions import * 

# Roslyn Parker
# Physical Simulations
# 13 April 2020

framerate = 60
dt = 1.0 / framerate
world = World()

m = 0.5
size = numpy.array([0.6, 0.2, 0.4])
x = numpy.array([0, 5, 0])
v = numpy.zeros(3)
axis = normalize(numpy.array([1, 1, -1]))
up = normalize(numpy.cross(axis, numpy.array([0, 0, 1])))
block = Block(m, size, x, v, axis, up)
world.bodies.append(block)

floor = Boundary(numpy.zeros(3), 10, 10)
world.boundaries.append(floor)

t = 0
while (t < 3):
    rate(framerate)
    world.update(dt)
    world.render()
    t += dt