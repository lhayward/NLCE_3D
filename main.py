import clust_order
import clust_weight
import numpy as np
import os.path  #to check if file exists
import sys  #for sys.stdout.flush()
import time

alphaEqTol=1e-10

def arrayEq(a1, a2):
  eq = True
  
  if len(a1) != len(a2):
    eq = False
  else:
    for elem1, elem2 in zip(a1, a2):
      #if elem1 != elem2:
      #Check if the alpha values are equal (within tolerance):
      if abs(elem1-elem2) > alphaEqTol:
        eq = False
        print "%.20f neq %.20f" %(elem1,elem2)
  #end if-else
  
  return eq
####### end arrayEq(a1, a2) function #######

def decimalStr(num):
  res = str(num)
  length = len(res)
  index = res.find('.')
  if index >= 0:
    res = res[0:index] + '-' + res[(index+1):length]
  return res
####### end decimalStr(num) function #######

def readArray(line):
  start = max( 0, line.find('[') )
  end = min( len(line), line.find(']') )
  return line[start+1:end] 
####### end readArray(line) function #######

def readWeights(alpha, massterm):
  w={}
  filename = "weights_mass" + decimalStr(massterm) + ".txt"
  if os.path.isfile(filename):
    fin = open(filename,'r')
    line = fin.readline()
    alpha_read = readArray(line).split()
    alpha_read = [float(a) for a in alpha_read]
    
    #if the alpha list is the same, then can use these previous stored weights:
    if arrayEq(alpha, alpha_read):
      line = fin.readline() #read blank line
      
      line = fin.readline()
      while line != "":
        #Read the name of the cluster:
        clust_name = line[ 0 : min(line.find(':'),len(line)) ]
        
        #Read the weights of that cluster for each alpha:
        weightArr = readArray(line).split()
        weightArr = np.array([float(elem) for elem in weightArr])
        
        #Store these weights:
        w[clust_name] = weightArr
        
        line = fin.readline()
      #end loop over lines
    else:
      print "WARNING: alpha list for the weights stored in " + filename + " does not match " + \
      "the current alpha list"
    
    fin.close()
    
  return w
####### end readWeights(alpha, massterm) #######

#############################
# User settings

order_min = 2
order_max = 4
order = clust_order.Max()
massterm = 0.0
#############################

t1 = time.clock()

clusters = []
alpha=np.array( np.linspace(0.4,10,49).tolist() + [20,50,100,200,500,1000] )
#alpha = np.array( [0.5, 1, 2])
total = np.zeros(len(alpha))
w = readWeights(alpha,massterm) #try to read in weights (if there are any stored)
print "\nInitial weights:"
print [key for key in w]
#print w

#Save the weights to file:
filename = "weights_mass" + decimalStr(massterm) + ".txt"
fout_w = open(filename, 'w')
#Write the alpha array:
fout_w.write("alpha = [ ")
for n in alpha:
  fout_w.write("%.3f\t" %n)
fout_w.write(" ]\n\n")
#Write the existing weights:
for key in w:
  fout_w.write(key + ": [ ")
  for i in range(len(alpha)):
    fout_w.write("%.20e\t" %w[key][i])
  fout_w.write(" ]\n")
fout_w.flush()

fout_res=[0 for i in alpha]
for i,n in enumerate(alpha):
  filename = "results_mass" + decimalStr(massterm) + "_alpha" + decimalStr(n) + ".txt"
  fout_res[i] = open(filename, 'w')

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
    
    #Write the new weight to file:
    fout_w.write(curr_clust_name + ": [ ")
    for i in range(len(alpha)):
      fout_w.write("%.20e\t" %w[curr_clust_name][i])
    fout_w.write(" ]\n")
    fout_w.flush()
    
  #end loop over loop over clusters
    
  # Save result to file
  for i in range(len(alpha)):
    fout_res[i].write("%d %.15f"%(ord,total[i])+'\n')
    fout_res[i].flush()
#end loop over orders

for fout in fout_res:
  fout.close()

fout_w.close()

print "\nFinal weights:"
print [key for key in w]
#print w

#Write the weights:
#for key in w:
#  fout_w.write(key + ": [ ")
#  for i in range(len(alpha)):
#    fout_w.write("%.20e\t" %w[key][i])
#  fout_w.write(" ]\n")

print "\nOrder done: ",str(order_max)
print

t2 = time.clock()
print "Total time: " + str(t2-t1) + " sec."
