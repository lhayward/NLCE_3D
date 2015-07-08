import numpy
#import scipy as sp
#import property_gen as pg

def weight(m,n,w):
  w_mxn_name = '%02d%02d'%(m,n)
  
  # First term in weight of mxn is property of mxn
  ### !!! WRITE CODE TO CALCULATE PROPERTY:
  ### w[w_mxn_name] = pg.property(m,n,Lx,Ly,lattice)
  w[w_mxn_name] = 0
  
  wformula = "W%02d%02d=P%02d%02d"%(m,n,m,n)
  print wformula
  print
  
  for y in range (1,n+1):
    for x in range (y,m+1):
      if (y < n or x < m) and x > 1:
        print wformula
        
        if x > n or x==y: coeff = (m-x+1)*(n-y+1)
        else: coeff = (m-x+1)*(n-y+1)+(m-y+1)*(n-x+1)
        
        wformula += "%+d*W%02d%02d"%(-coeff,x,y)
        
        ### !!! Write NLCE base case# 
        ### w[w_mxn_name] -= coeff * w['%02d%02d'%(x,y)]
        
  ###print w[w_mxn_name]
  print wformula
  
  return w
