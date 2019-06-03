import numpy as np
import random

SOFTENING = 1e-9
class Body:
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]

    def randomizeBodies(self,n):
        for i in range(n):
           self.x.append(random.uniform(-1.0,1.0)) # generates a random number between -1.0 to 1.0
           self.y.append(random.uniform(-1.0,1.0))
           self.z.append(random.uniform(-1.0,1.0))
           self.vx.append(random.uniform(-1.0,1.0))
           self.vy.append(random.uniform(-1.0,1.0))
           self.vz.append(random.uniform(-1.0,1.0))


    def bodyForce(p, dt, n):
        for i in range(n):
            Fx=0.0
            Fy=0.0
            Fz=0.0
            for j in range(n):
                dx=p[j].x - p[i].x
                dy=p[j].y - p[i].y
                dz=p[j].z - p[i].z
                distSqr = dx*dx + dy*dy + dz*dz + SOFTENING
                invDist = np.sqrt(distSqr)
                invDist3 = invDist**3

                Fx += dx*invDist3
                Fy += dy*invDist3
                Fz += dz*invDist3
            p[i].vx += dt*Fx
            p[i].vy += dt*Fy
            p[i].vz += dt*Fz


if __name__ == "__main__":
    nBodies = 2<<11;
    dt = 0.01
    nIters  = 10

    B = Body();
    B.randomizeBodies(nBodies) 