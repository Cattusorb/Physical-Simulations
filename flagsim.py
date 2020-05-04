from cloth import *

# Roslyn Parker
# Physical Simulations
# 4 May 2020

def flagpole(cloth: Cloth, width: float):
    ''' binds the left edge of a cloth to a "flagpole" '''
    for i in range(0, width):
        cloth.particles[0][i].isMovable = False
        print(i)

width = 15
cloth = Cloth(1, 5, 3, 25, width, numpy.zeros(3), numpy.array([1, 0, 0]), numpy.array([0, 1, 0]))
flagpole(cloth, width)

framerate = 30
dt = 1 / framerate

# Pre compute simulation
t = 0
frames = [cloth.getFrame()]
while t < 5:
    cloth.update(dt, numpy.array([0.3, 0, 0.2]))
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