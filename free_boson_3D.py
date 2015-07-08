import comparisons as cmp
import numpy as np
from scipy import linalg
import time
import dbm
import math
import glob
start_time = time.time()
import os

def getCornerEnt(Lx,Ly,Lz):
  perx = pery = perz = False  #PBC or not along the x, y and z directions
  massterm = 1

  Ns = Lx * Ly * Lz

  K = np.zeros((Ns, Ns))
  #Loop to calculate the matrix K:
  for x in range(Lx):
    for y in range(Ly):
      for z in range(Lz):
        K[site(Lx,Ly,Lz,x,y,z),site(Lx,Ly,Lz,x,y,z)]= 6.0 + (massterm ** (2))
        xp = (x+1)%Lx
        yp = (y+1)%Ly
        zp = (z+1)%Lz
        if (xp > x) or perx:
          K[site(Lx,Ly,Lz,xp,y,z), site(Lx,Ly,Lz,x, y,z)] = -1.0 
          K[site(Lx,Ly,Lz,x, y,z), site(Lx,Ly,Lz,xp,y,z)] = -1.0
        if (yp > y) or pery:
          K[site(Lx,Ly,Lz,x,yp,z), site(Lx,Ly,Lz,x,y, z)] = -1.0 
          K[site(Lx,Ly,Lz,x,y, z), site(Lx,Ly,Lz,x,yp,z)] = -1.0
        if (zp > z) or perz:
          K[site(Lx,Ly,Lz,x,y,zp), site(Lx,Ly,Lz,x,y,z )] = -1.0
          K[site(Lx,Ly,Lz,x,y,z ), site(Lx,Ly,Lz,x,y,zp)] = -1.0
  
  Eval,Evec = np.linalg.eigh(K) #, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False, check_finite=True)
  #print Evec
  P = 1./2. * np.matrix(Evec) * np.matrix(np.diag(np.sqrt(Eval))) * np.matrix(Evec.T)
  X = 1./2. * np.matrix(Evec) * np.matrix(np.diag(1. / np.sqrt(Eval))) * np.matrix(Evec.T)
  
  result=0
  for x0 in range(1,Lx):
    for y0 in range(1,Ly):
      for z0 in range(1,Lz):
        #corners:
        result += 1.0/4.0*(getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.geq,cmp.geq,cmp.geq),X,P) 
                         + getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.lt, cmp.lt, cmp.geq),X,P) 
                         + getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.geq,cmp.lt, cmp.lt ),X,P)
                         + getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.lt, cmp.geq,cmp.lt ),X,P))
        #edges:
        result -= 1.0/4.0*(getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.geq,cmp.geq,cmp.any),X,P) 
                         + getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.lt, cmp.lt ,cmp.any),X,P) 
                         + getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.geq,cmp.any,cmp.geq),X,P)
                         + getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.lt ,cmp.any,cmp.lt ),X,P)
                         + getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.any,cmp.geq,cmp.geq),X,P)
                         + getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.any,cmp.lt ,cmp.lt ),X,P))
        #planes:
        result += 1.0/4.0*(getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.geq,cmp.any,cmp.any),X,P) 
                         + getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.any,cmp.geq,cmp.any),X,P) 
                         + getEntropy(Lx,Ly,Lz,getRegionA(Lx,Ly,Lz,x0,y0,z0,cmp.any,cmp.any,cmp.geq),X,P))
      #loop over z0
    #loop over y0
  #loop over x0
  return result
#.......................................END getCornerEnt.......................................

def getEntropy(Lx,Ly,Lz,regA,X,P):
  sitesA = []
  for x in range(Lx):
    for y in range(Ly):
      for z in range(Lz):
        if regA[x,y,z] == True: sitesA.append(site(Lx,Ly,Lz,x,y,z))
  NsA=len(sitesA)
  Pred = np.zeros((NsA, NsA))
  Xred = np.zeros((NsA, NsA))

  # ............ Following step differs from Mathematica ...........
  #Get X and P for sites inside A only:
  for a in range(0, NsA):
    for b in range(0, NsA):
        Pred[a,b] = P[sitesA[a] , sitesA[b]]
        Xred[a,b] = X[sitesA[a] , sitesA[b]]
  #....end of the step .......

  Csquared = (Xred.T).dot(Pred)
  Ev = np.sqrt(np.linalg.eigvals(Csquared))
  #print np.linalg.eigvals(Csquared)
  #print
  Sn = 0. 
  for j in range(0, NsA):
    if Ev[j] > 0.5:
      Sn += (Ev[j]+1./2) * np.log(abs(Ev[j]+1./2.))-(Ev[j]-1./2.)*np.log(abs(Ev[j]-1./2)) 
  return Sn
#........................................END getEntropy........................................

def getRegionA(Lx,Ly,Lz,x0,y0,z0,fx,fy,fz):
  regA = np.zeros( (Lx,Ly,Lz), dtype='bool' )
  
  for x in range(Lx):
    for y in range(Ly):
      for z in range(Lz):
        if (fx(x,x0) and fy(y,y0) and fz(z,z0)):
          regA[x,y,z] = True
        #endif
      #end for loop over z
    #end for loop over y
  # end for loop over x
  
  #if the size of region A is more than half of the lattice, swap regions A and B:
  if( (regA==True).sum() > (Lx*Ly*Lz/2) ): regA = np.logical_not(regA)
  
  return regA
  
def site(Lx,Ly,Lz,x,y,z):
    return x + (y*Lx) + (z*Lx*Ly) #convert (x,y,z) pair to a single site number
