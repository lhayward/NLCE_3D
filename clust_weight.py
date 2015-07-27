import numpy
import free_boson_3D
#import scipy as sp
#import property_gen as pg

def weight(Lx,Ly,Lz,w,alpha,massterm):
  w_clust_name = clust_name(Lx,Ly,Lz)
  
  if w_clust_name not in w:
    #First term in weight of this cluster is the property of the cluster:
    w[w_clust_name] = free_boson_3D.getCornerEnt(Lx,Ly,Lz,alpha,massterm)
  
    wformula = "W(" + clust_name(Lx,Ly,Lz) + ") = P(" + clust_name(Lx,Ly,Lz) + ")"
    #print wformula
  
    # Subtract weights of all subclusters:
    for x in range(1, Lx+1):
      for y in range(1, Ly+1):
        for z in range(1, Lz+1):
          # Check that we are not at the Lx x Ly x Lz cluster, and also check that we are not at
          # the 1 x 1 x 1 cluster (both of these cases don't contribute:
          if (x!=Lx or y!=Ly or z!=Lz) and max(x,y,z)>1:
        
            coeff = (Lx-x+1)*(Ly-y+1)*(Lz-z+1)
          
            #weight is stored such that x<=y<=z, so sort the current x,y,z
            [xs,ys,zs] = sorted([x,y,z])
          
            wformula += "+ %d*W("%(-coeff) + clust_name(xs,ys,zs) + ")"
            w[w_clust_name] -= coeff * w[clust_name(xs,ys,zs)]
          
    #print wformula
    #print
  #end if
  
  return w
  
def clust_name(Lx,Ly,Lz):
  return "%02dx%02dx%02d"%(Lx,Ly,Lz)
