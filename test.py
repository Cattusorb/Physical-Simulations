from vpython import *
#sphere()
ball = sphere(pos=vector(-5, 0, 0), radius=0.5, color=color.yellow, make_trail=True)
wallR = box(pos=vector(6, 0, 0), size=vector(0.2, 12, 12), color=color.purple)
wallL = box(pos=vector(-6, 0, 0), size=vector(0.2, 12, 12), color=color.purple)
wallT = box(pos=vector(0, 6, 0), size=vector(12, 0.2, 12), color=color.purple)
wallB = box(pos=vector(0, -6, 0), size=vector(12, 0.2, 12), color=color.purple)
wallZ = box(pos=vector(0, 0, -6), size=vector(12, 12, 0.2), color=color.purple)
ball.velocity = vector(25,10,10)

vscale = 0.1
varr = arrow(pos=ball.pos, axis=ball.velocity, color=color.yellow)

deltat = 0.005
t = 0

scene.autoscale = False

while True:
    rate(50)
    if ball.pos.x > wallR.pos.x:
        ball.velocity.x = -ball.velocity.x
    if ball.pos.x < wallL.pos.x: 
        ball.velocity.x = -ball.velocity.x
    if ball.pos.y > wallT.pos.y:
        ball.velocity.y = -ball.velocity.y
    if ball.pos.y < wallB.pos.y:
        ball.velocity.y = -ball.velocity.y
    if ball.pos.z < wallZ.pos.z:
        ball.velocity.z = -ball.velocity.z
    if ball.pos.z > -wallZ.pos.z:
        ball.velocity.z = -ball.velocity.z
    ball.pos = ball.pos + ball.velocity*deltat
    t = t + deltat
    varr.pos = ball.pos
    varr.axis = ball.velocity*vscale


