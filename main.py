#import numpy as np
import clust_order
import clust_weight

#############################
# User settings

order_min = 2
order_max = 7
order = clust_order.Max()
#############################

filename = "results.txt"
print "Writing file ",filename
f = open(filename, 'w')

total = 0
w = {} # weights
clusters = []

for ord in range(order_min,order_max+1):
  print
  print "Order %02d:" %ord
  for Lx,Ly,Lz in order.clusters(ord):
    curr_clust_name = clust_weight.clust_name(Lx,Ly,Lz)
    print "  " + curr_clust_name
    
    w = clust_weight.weight(Lx,Ly,Lz,w) # performs cluster weight calculations
    
    #Embedding factor:
    ef = 6
    if Lx==Ly and Ly==Lz: ef = 1
    elif Lx==Ly or Lx==Lz or Ly==Lz: ef = 3
    
    total = total + ef*w[curr_clust_name]
    
  # Save result to file
  f.write("%d %.15f"%(ord,total)+'\n')

f.close()
print "\nOrder done: ",str(order_max)
print
