import comparisons as cmp
import numpy as np
#from scipy import linalg
#import time
#start_time = time.time()

#Ev_tol = 10e-16

###############################################################################################
#######################################  getCornerEnt  ########################################
###############################################################################################
def getCornerEnt(Lx,Ly,Lz,alpha,massterm):
  perx = pery = perz = False  #PBC or not along the x, y and z directions
  
  Ns = Lx * Ly * Lz
  L  = (Lx,Ly,Lz)
  
  K = np.zeros((Ns, Ns))
  #Loop to calculate the matrix K:
  for x in range(Lx):
    for y in range(Ly):
      for z in range(Lz):
        K[site(L,x,y,z),site(L,x,y,z)]= 6.0 + (massterm ** (2))
        xp = (x+1)%Lx
        yp = (y+1)%Ly
        zp = (z+1)%Lz
        if (xp > x) or perx:
          K[site(L,xp,y,z), site(L,x, y,z)] = -1.0 
          K[site(L,x, y,z), site(L,xp,y,z)] = -1.0
        if (yp > y) or pery:
          K[site(L,x,yp,z), site(L,x,y, z)] = -1.0 
          K[site(L,x,y, z), site(L,x,yp,z)] = -1.0
        if (zp > z) or perz:
          K[site(L,x,y,zp), site(L,x,y,z )] = -1.0
          K[site(L,x,y,z ), site(L,x,y,zp)] = -1.0
  
  #use eigh because we know the matrix is symmetric
  Eval,Evec = np.linalg.eigh(K) #, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False, check_finite=True)
  #print Evec
  P = 1./2. * np.matrix(Evec) * np.matrix(np.diag(np.sqrt(Eval))) * np.matrix(Evec.T)
  X = 1./2. * np.matrix(Evec) * np.matrix(np.diag(1. / np.sqrt(Eval))) * np.matrix(Evec.T)
  
  #result=0
  result = np.zeros(len(alpha))
#   #loop over all possible locations of the corner:
#   for x0 in range(1,Lx):
#     for y0 in range(1,Ly):
#       for z0 in range(1,Lz):
#         r0 = (x0,y0,z0)
#         
#         ################################# CORNER CONTRIBUTION #################################
#         #corners:
#         result += 1.0/4.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.geq,cmp.geq),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.lt, cmp.lt, cmp.geq),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.lt, cmp.lt ),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.lt, cmp.geq,cmp.lt ),X,P))
#         #result += getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.geq,cmp.geq),X,P)
#         
#         #edges:
#         result -= 1.0/4.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.geq,cmp.any),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.lt, cmp.lt ,cmp.any),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.any,cmp.geq),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.lt ,cmp.any,cmp.lt ),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.geq,cmp.geq),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.lt ,cmp.lt ),X,P))
#         #planes:
#         result += 1.0/4.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.any,cmp.any),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.geq,cmp.any),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.any,cmp.geq),X,P))
#                          
# #         ################################## EDGE CONTRIBUTION ##################################                         
# #         result += 1.0/6.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.geq,cmp.any),X,P)
# #                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.lt, cmp.lt ,cmp.any),X,P)
# #                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.any,cmp.geq),X,P)
# #                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.lt ,cmp.any,cmp.lt ),X,P)
# #                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.geq,cmp.geq),X,P)
# #                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.lt ,cmp.lt ),X,P))
# #         #planes:
# #         result -= 1.0/3.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.any,cmp.any),X,P)
# #                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.geq,cmp.any),X,P)
# #                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.any,cmp.geq),X,P))
#       #loop over z0
#     #loop over y0
#   #loop over x0
  
  ################################# CORNER CONTRIBUTION #################################
  #loop over all possible locations of the corner:
  #corners:
  for x0 in range(1,Lx):
    for y0 in range(1,Ly):
      for z0 in range(1,Lz):
        r0 = (x0,y0,z0)
        result += getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.geq,cmp.geq),X,P)
  
  #edges:
  for x0 in range(1,Lx):
    for y0 in range(1,Ly):
      r0 = (x0,y0,1)
      result -= (Lz-1)/2.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.geq,cmp.any),X,P))
  for x0 in range(1,Lx):
    for z0 in range(1,Lz):
      r0 = (x0,1,z0)
      result -= (Ly-1)/2.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.any,cmp.geq),X,P))
  for y0 in range(1,Ly):
    for z0 in range(1,Lz):
      r0 = (1,y0,z0)
      result -= (Lx-1)/2.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.geq,cmp.geq),X,P))
      
  #planes:
  for x0 in range(1,Lx):
    r0 = (x0,1,1)
    result += (Ly-1)*(Lz-1)/4.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.any,cmp.any),X,P))
  for y0 in range(1,Ly):
    r0 = (1,y0,1)
    result += (Lx-1)*(Lz-1)/4.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.geq,cmp.any),X,P))
  for z0 in range(1,Lz):
    r0 = (1,1,z0)
    result += (Lx-1)*(Ly-1)/4.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.any,cmp.geq),X,P))
  
  return result
#.......................................END getCornerEnt.......................................

def getEntPart1(regA,P,X):
  sitesA = [i for i in range(len(regA)) if regA[i]]   
  Pred = P[sitesA][:,sitesA]
  Xred = X[sitesA][:,sitesA]
  return Pred,Xred
  

###############################################################################################
########################################  getEntropy  #########################################
###############################################################################################
def getEntropy((Lx,Ly,Lz),alpha,regA,X,P):
#   sitesA = []
#   for x in range(Lx):
#     for y in range(Ly):
#       for z in range(Lz):
#         if regA[x,y,z] == True: sitesA.append(site((Lx,Ly,Lz),x,y,z))
#   NsA=len(sitesA)
#   Pred = np.zeros((NsA, NsA))
#   Xred = np.zeros((NsA, NsA))
#   
#   # ............ Following step differs from Mathematica ...........
#   #Get X and P for sites inside A only:
#   for a in range(0, NsA):
#     for b in range(0, NsA):
#       Pred[a,b] = P[sitesA[a] , sitesA[b]]
#       Xred[a,b] = X[sitesA[a] , sitesA[b]]
#   #....end of the step .......

  sitesA = [i for i in range(len(regA)) if regA[i]]   
  Pred = P[sitesA][:,sitesA]
  Xred = X[sitesA][:,sitesA]
  #Pred,Xred = getEntPart1(regA,P,X)
  
  #Csquared = Xred.dot(Pred)
  #Csquared = (Xred.T).dot(Pred)
  Csquared = np.matrix(Xred)*np.matrix(Pred)
  
  Ev = np.sqrt(np.linalg.eigvals(Csquared))
  
#   Sn = np.zeros(len(alpha))  
#   for j in range(0, NsA):
#     #if Ev[j] > 0.5:
#     if Ev[j].real > 0.5:
#       for i,n in enumerate(alpha):
#         if n == 1:
#           Sn[i] += (Ev[j]+1./2)*np.log(abs(Ev[j]+1./2.)) - (Ev[j]-1./2.)*np.log(abs(Ev[j]-1./2))
#         else:
#           Sn[i] += 1.0/(n-1.0) * np.log( (Ev[j]+1./2)**n - (Ev[j]-1./2.)**n )

  #Sn = 0. 
  Sn = np.zeros(len(alpha))
  Ev_new = np.array([e.real for e in Ev if (e.real - 1./2.)>0])
  for i,n in enumerate(alpha):
    if n == 1:
      Sn[i] = np.sum( (Ev_new+1./2)*np.log(Ev_new+1./2.) - (Ev_new-1./2.)*np.log(Ev_new-1./2) )
    else:
      Sn[i] = 1.0/(n-1.0)*np.sum( np.log( (Ev_new+1./2)**n - (Ev_new-1./2.)**n ) )
  
  return Sn
#........................................END getEntropy........................................

###############################################################################################
########################################  getRegionA  #########################################
###############################################################################################
def getRegionA((Lx,Ly,Lz),(x0,y0,z0),fx,fy,fz):
#   regA = np.zeros( (Lx,Ly,Lz), dtype='bool' )
#   
#   for x in range(Lx):
#     for y in range(Ly):
#       for z in range(Lz):
#         if (fx(x,x0) and fy(y,y0) and fz(z,z0)):
#           regA[x,y,z] = True
#         #endif
#       #end for loop over z
#     #end for loop over y
#   # end for loop over x
#   
#   #if the size of region A is more than half of the lattice, swap regions A and B:
#   if( (regA==True).sum() > (Lx*Ly*Lz/2.0) ): regA = np.logical_not(regA)
  NTot=Lx*Ly*Lz
  regA = np.zeros( NTot, dtype='bool' )
  
  for x in range(Lx):
    for y in range(Ly):
      for z in range(Lz):
        if (fx(x,x0) and fy(y,y0) and fz(z,z0)):
          regA[site((Lx,Ly,Lz),x,y,z)] = True
        #endif
      #end for loop over z
    #end for loop over y
  # end for loop over x
  
  #if the size of region A is more than half of the lattice, swap regions A and B:
  if( (regA==True).sum() > (NTot/2.0) ): regA = np.logical_not(regA)
  
  return regA

###############################################################################################
###########################################  site  ############################################
###############################################################################################
def site((Lx,Ly,Lz),x,y,z):
  return x + (y*Lx) + (z*Lx*Ly) #convert (x,y,z) pair to a single site number
