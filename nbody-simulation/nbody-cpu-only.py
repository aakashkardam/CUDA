import numpy as np
import random
import time

SOFTENING = 1e-9


# Class Definition for the system of particles as Body
class Body:
  x=[]
  y=[]
  z=[]
  vx=[]
  vy=[]
  vz=[]


# Generates random data for the position (x,y,z) and the velocity (vx,vy,vz) for the nbodies
  def randomizeBodies(self,n):
    for i in range(n):
      self.x.append(random.uniform(-1.0,1.0)) # generates a random number between -1.0 to 1.0
      self.y.append(random.uniform(-1.0,1.0))
      self.z.append(random.uniform(-1.0,1.0))
      self.vx.append(random.uniform(-1.0,1.0))
      self.vy.append(random.uniform(-1.0,1.0))
      self.vz.append(random.uniform(-1.0,1.0))


# Calculates the bodyForces on a body due to all other bodies
  def bodyForce(self, dt, n):
    for i in range(n):
      Fx=0.0
      Fy=0.0
      Fz=0.0
      for j in range(n):
        dx=self.x[j] - self.x[i]
        dy=self.y[j] - self.y[i]
        dz=self.z[j] - self.z[i]
        distSqr = dx*dx + dy*dy + dz*dz + SOFTENING
        invDist = np.sqrt(distSqr)
        invDist3 = invDist**3

        Fx += dx*invDist3
        Fy += dy*invDist3
        Fz += dz*invDist3
      self.vx[i] += dt*Fx
      self.vy[i] += dt*Fy
      self.vz[i] += dt*Fz


# Prints the calculated forces on the bodies after nIters
  def printInteractions(self):
    for i in range(np.shape(B.x)[0]):
      print("BDY ",i,"\tB.x = ","{0:.2}".format(self.x[i]), "\tB.y = ", "{0:.2}".format(self.y[i]),\
              "\tB.z = ", "{0:.2}".format(self.z[i]), "\tB.vx = ","{0:.2}".format(self.vx[i]), "\tB.vy = ",\
              "{0:.2}".format(self.vy[i]), "\tB.vz = ", "{0:.2}".format(self.vz[i]))


#if __name__ == "__main__":
nBodies = 2<<11;
dt = 0.01
nIters  = 10

B = Body();
B.randomizeBodies(nBodies)
start=time.time()
B.bodyForce(dt,nBodies)
end=time.time()
print("time elapsed in calculating the body forces = ",end-start)
B.printInteractions()
