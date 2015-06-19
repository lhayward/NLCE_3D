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

  def site(x, y, z):
    return x + (y*Lx) + (z*Lx*Ly) #convert (x,y,z) pair to a single site number

  K = np.zeros((Ns, Ns))
  #Loop to calculate the matrix K:
  for x in range(Lx):
    for y in range(Ly):
      for z in range(Lz):
        K[site(x,y,z),site(x,y,z)]= 6.0 + (massterm ** (2))
        xp = (x+1)%Lx
        yp = (y+1)%Ly
        zp = (z+1)%Lz
        if (xp > x) or perx:
          K[site(xp,y,z), site(x, y,z)] = -1.0 
          K[site(x, y,z), site(xp,y,z)] = -1.0
        if (yp > y) or pery:
          K[site(x,yp,z), site(x,y, z)] = -1.0 
          K[site(x,y, z), site(x,yp,z)] = -1.0
        if (zp > z) or perz:
          K[site(x,y,zp), site(x,y,z )] = -1.0
          K[site(x,y,z ), site(x,y,zp)] = -1.0
  
  Eval,Evec = np.linalg.eigh(K) #, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False, check_finite=True)
  print Evec
  P = 1./2. * np.matrix(Evec) * np.matrix(np.diag(np.sqrt(Eval))) * np.matrix(Evec.T)
  X = 1./2. * np.matrix(Evec) * np.matrix(np.diag(1. / np.sqrt(Eval))) * np.matrix(Evec.T)
  
#   filename=file.split(".")[0]
#   db = dbm.open(filename,'w')
#   for key in db.keys():
#       if db[key] == str(-99):
#           clust = np.array(kg.decompress(key))
#           #print clust
#           sitesA = []
#           for y in range(0, Ly):
#               for x in range(0, Lx):
#                   if clust[y, x] == 1:
#                       sitesA.append(site(x, y))
#           NsA=len(sitesA)
#           Pred = np.zeros((NsA, NsA))
#           Xred = np.zeros((NsA, NsA))
#           # ............ Following step differs from Mathematica ...........
#           #Get X and P for sites inside A only:
#           for a in range(0, NsA):
#               for b in range(0, NsA):
#                   Pred[a,b] = P[sitesA[a] , sitesA[b]]
#                   Xred[a,b] = X[sitesA[a] , sitesA[b]]
#                   #....end of the step .......
# 
#           Csquared = (Xred.T).dot(Pred)
#           Ev = np.sqrt(np.linalg.eigvals(Csquared))
#           Sn = 0. 
#           for j in range(0, NsA):
#               if Ev[j] > 0.5:
#                   Sn += (Ev[j]+1./2) * np.log(abs(Ev[j]+1./2.))-(Ev[j]-1./2.)*np.log(abs(Ev[j]-1./2)) 
#           db[key] = str("%.15f"%Sn)
#   db.close()
# 
#   newfilename = filename.split("/")[0]+"/"+filename.split("/")[1][2:]
#   os.rename(filename+".db", newfilename+".db")
#   print file , "done"
#     
#   print "Free boson calculation complete"

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
