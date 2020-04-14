 # Roslyn Parker
 # Physical Simulations
 # Project 2
 # 14 April 2020

from rigidbodyproject import * 

'''
● Performing the SAT is relatively expensive 
● Don’t want to do it for each pair of objects 
● Do initial broad-phase collision detection 
● Draw a bounding sphere around each box 
● Perform SAT on objects whose bounding spheres overlap 
● Axis-aligned bounding boxes (AABBs) are better
'''
# then set-up broad-phase collision (listed above in block comment)

def createGrid(n: int): 
    '''Created a grid of size n x n'''
    size = numpy.array([0.6, 0.6, 0.6]) # square
    axis = normalize(numpy.array([0, 1, 0]))
    up = normalize(numpy.cross(axis, numpy.array([0, 1, 1])))
    m = 0.5
    
    blocks = []
    for x in range (0, n):
        for y in range (0, n):
            block = Block(m, size, numpy.array([x, y, 0]), numpy.array([0, 0, 0]), axis, up)
            blocks.append(block)
    
    return blocks

framerate = 60
dt = 1.0 / framerate
world = World()

gridSize = 5 # set grid size

# x and y sphere starting point 
# (for best sim)
if gridSize % 2 == 0: 
    xy = gridSize / 3
else: 
    xy = gridSize / 2

sphere = Sphere(1.0, 0.4, numpy.array([xy, xy, 5]), numpy.array([0, 0, -5]))
world.bodies.append(sphere)

# create grid and add blocks to world
gridBlocks = createGrid(gridSize)
for block in gridBlocks:
    world.bodies.append(block)

t = 0
while (t < 3):
    rate(framerate)
    world.update(dt)
    world.render()
    t += dt