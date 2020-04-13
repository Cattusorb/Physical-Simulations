from rigidbodycollisions import * 

# Roslyn Parker
# Physical Simulations
# 13 April 2020

'''
The simulation does look fairly acctruate, 
even when gravity is turned off. The tumbling 
action that occurs it mostly accurate and looks cool! 
'''

framerate = 60
dt = 1.0 / framerate
world = World()

m = 0.5
size = numpy.array([0.6, 0.2, 0.4])
axis = normalize(numpy.array([0, 1, 0]))
up = normalize(numpy.cross(axis, numpy.array([0, 1, 1])))
block = Block(m, size, numpy.array([0, 0, 0]), numpy.zeros(3), axis, up)
world.bodies.append(block)

sphere = Sphere(1.0, 0.2, numpy.array([-2, 0.2, 0.3]), numpy.array([2, 0, 0]))
world.bodies.append(sphere)

scene.autoscale = False

t = 0
while (t < 3):
    rate(framerate)
    world.update(dt)
    world.render()
    t += dt