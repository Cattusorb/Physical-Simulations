from cloth import *

# Roslyn Parker
# Physical Simulations
# 4 May 2020

width = 30
cloth = Cloth(1, 5, 5, width, width, numpy.zeros(3), numpy.array([1, 0, 0]), numpy.array([0, 0, 1]))
cloth.particles[0][0].isMovable = False
cloth.particles[0][-1].isMovable = False
cloth.particles[-1][0].isMovable = False
cloth.particles[-1][-1].isMovable = False

framerate = 30
dt = 1 / framerate

# Pre compute simulation
t = 0
frames = [cloth.getFrame()]
while t < 5:
    if t > 2: 
        cloth.particles[0][-1].isMovable = True
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