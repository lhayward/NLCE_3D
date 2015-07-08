#import numpy as np
import mxn_order
import mxn_weight

#############################
# User settings

order_min = 1
order_max = 3 
order = mxn_order.Max()
#############################

# filename = "Results_%02d"%(R,)
# print "Writing file ",filename
# f = open(filename, 'a')
total = 0
w = {} # weights

clusters = []

for ord in range(order_min,order_max+1):
  print
  print "Order %02d" %ord
  for Lx,Ly,Lz in order.clusters(ord):
    print "%02d x %02d x %02d" %(Lx,Ly,Lz)
    # w = mxn_weight.weight(m,n,R,w) # performs cluster weight calculations
    
    # #Embedding factor (1 for squares, 2 for rectangles):
    # Lc = 1
    # if m != n: Lc = 2   #see line 24, mxn_weight.py
    # # cannot use total += w['%02d%02d'%(m,n)] or else W0202 somehow gets changed every iteration
    # total = total + Lc*w['%02d%02d'%(m,n)]
    
  # # Save result to file
  # f.write("%d %.15f"%(ord,total)+'\n')

# f.close()

# print "Order done: ",str(order_max)
