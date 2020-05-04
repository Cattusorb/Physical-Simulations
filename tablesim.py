from cloth import *

# Roslyn Parker
# Physical Simulations
# 4 May 2020

def createTable(cloth: Cloth, particalr: float, centerl: int, centerw: int):
    ''' creates a "table" with the radius r '''
    # Currently is a square table.. 
    # set all particles within the radius isMovable = False
    # if the particle is within the radius the circle
    # then make it imovable, else, no!
    for i in range(0, cloth.nLength - 1):
        for j in range(0, cloth.nWidth - 1):
            if numpy.linalg.norm(cloth.particles[i][j].x) < particalr:
                cloth.particles[i][j].isMovable = False

    '''
    for i in range(centerl, centerl + particalr):
        for j in range(centerw, centerw + particalr):
            cloth.particles[i][j].isMovable = False

    for i in range(centerl - particalr, centerl):
        for j in range(centerw, centerw + particalr):
            cloth.particles[i][j].isMovable = False

    for i in range(centerl, centerl + particalr):
        for j in range(centerw - particalr, centerw):
            cloth.particles[i][j].isMovable = False

    for i in range(centerl - particalr, centerl):
        for j in range(centerw - particalr, centerw):
            cloth.particles[i][j].isMovable = False'''

plw = 30
lw = 5
cloth = Cloth(1, lw, lw, plw, plw, numpy.zeros(3), numpy.array([1, 0, 0]), numpy.array([0, 0, 1]))
createTable(cloth, 1, 14, 14)

framerate = 30
dt = 1 / framerate

# Pre compute simulation
t = 0
frames = [cloth.getFrame()]
while t < 5:
    cloth.update(dt, numpy.array([0, 0, 0]))
    frames.append(cloth.getFrame())
    t += dt
    print(t)

input("Press Enter to render...")

# Render
#scene.autoscale = False
mesh = MeshRendering(frames[0], 0.02)
for frame in frames:
    rate(framerate)
    mesh.update(frame)

print("Done!")