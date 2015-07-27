import clust_order
import clust_weight
import numpy as np
import sys  #for sys.stdout.flush()
import time

def decimalStr(num):
  res = str(num)
  length = len(res)
  index = res.find('.')
  if index >= 0:
    res = res[0:index] + '-' + res[(index+1):length]
  return res

#############################
# User settings

order_min = 2
order_max = 4
order = clust_order.Max()
massterm = 0
#############################

t1 = time.clock()

w = {} # weights
clusters = []
#alpha=np.array( np.linspace(0.4,10,49).tolist() + [20,50,100,200,500,1000] )
alpha = np.array( [0.5, 1, 2])
total = np.zeros(len(alpha))

f=[0 for i in alpha]
for i,n in enumerate(alpha):
  filename = "results_mass0-0_alpha" + decimalStr(n) + ".txt"
  f[i] = open(filename, 'w')

for ord in range(order_min,order_max+1):
  print
  print "Order %02d:" %ord
  for Lx,Ly,Lz in order.clusters(ord):
    curr_clust_name = clust_weight.clust_name(Lx,Ly,Lz)
    print "  " + curr_clust_name
    sys.stdout.flush()
    
    w = clust_weight.weight(Lx,Ly,Lz,w,alpha,massterm) # performs cluster weight calculations
    #Embedding factor:
    ef = 6
    if Lx==Ly and Ly==Lz: ef = 1
    elif Lx==Ly or Lx==Lz or Ly==Lz: ef = 3
    
    total = total + ef*w[curr_clust_name]
    
  # Save result to file
  for i in range(len(alpha)):
    f[i].write("%d %.15f"%(ord,total[i])+'\n')
    f[i].flush()

for ff in f:
  ff.close()
print "\nOrder done: ",str(order_max)
print

t2 = time.clock()
print "Total time: " + str(t2-t1) + " sec."
