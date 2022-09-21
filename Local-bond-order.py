import math
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.special import sph_harm as sh
#from numba import jit

#@jit(nopython=True)
def bond_order(x,y,z,l):
              Y = np.zeros((natom,2*l+1) , dtype =np.complex) #Y_lm array for ith particle for nearest neighbour
              neighbour = np.zeros(natom)
              for i in range(natom):
                       #mod = complex(0,0)
                       for j in range(i+1,natom):
                                    xij = x[j] - x[i]
                                    yij = y[j] - y[i]
                                    zij = z[j] - z[i]
                                    xij = xij - L*round(xij/L)
                                    yij = yij - L*round(yij/L)
                                    zij = zij - L*round(zij/L)
                                    rij = math.sqrt(xij*xij+yij*yij+zij*zij)
                                    if (rij < rcut):
                                        thetaij = math.acos(zij/rij)*180./math.pi
                                        if (xij > 0. and yij > 0.):
                                                phiij = math.atan(yij/xij)*180./math.pi
                                        if (xij < 0. and (yij > 0. or yij < 0.)):
                                                phiij = 180. + math.atan(yij/xij)*180./math.pi
                                        if (xij > 0. and yij < 0.):
                                                phiij = 360. + math.atan(yij/xij)*180./math.pi
                                        m = - l
                                        while m < l :
                                                    Y[i][m+l] = Y[i][m+l] + sh(m,l,phiij*math.pi/180.,thetaij*math.pi/180.)
                                                    Y[j][m+l] = Y[j][m+l] + sh(m,l,phiij*math.pi/180.,thetaij*math.pi/180.) 
                                                    m = m + 1 
                                        neighbour[i] = neighbour[i] + 1
                                        neighbour[j] = neighbour[j] + 1
                       if neighbour[i]!=0 :
                                       Y[i] = Y[i]/neighbour[i]
                                       #print (Y[i])     
                                       
      
              
              q6 = np.zeros((natom,2*l+1) , dtype =np.complex) #Y_lm array for ith particle for nearest neighbour 
              Q = []        
              for i in range(natom):
                       mod = 0.0
                       q6[i] = q6[i] + Y[i]
                       for j in range(i+1 , natom):
                                xij = x[j] - x[i]
                                yij = y[j] - y[i]
                                zij = z[j] - z[i]
                                xij = xij - L*round(xij/L)
                                yij = yij - L*round(yij/L)
                                zij = zij - L*round(zij/L)
                                rij = math.sqrt(xij*xij+yij*yij+zij*zij)          
                                if rij < rcut :       
                                        q6[i] = q6[i] + Y[j]                                       
                                        q6[j] = q6[j] + Y[i]
                       q6[i] = q6[i]/(neighbour[i]+1)
                       for m in range(2*l+1):
                              mod = mod + q6[i][m]*np.conjugate(q6[i][m])
                       Q.append(mod.real)                       
              Q = [math.sqrt(number*4*math.pi/(2*l+1))  for number in Q]
                                      
              return Q                      
#-----------------------------------------------------------------------------------------------------------------
#frequency calculation
def frequency(Q1,atom):
    cwd = os.getcwd()
    file1 = open(cwd + '/freq_eq_q6.dat','w')
    #file1 = open('freq_eq_q6-32.dat','w')
    step = 0.01
    minimum = min(Q1)
    maximum = max(Q1)
    print(maximum , minimum,atom)
    array = np.arange(minimum,maximum,step)
    binc = int((maximum-minimum)/step)+1
    hist = np.zeros(binc)
    count = 0
    for i in range(atom):
        bin1 = int((Q1[i]-minimum)/step)
        if (bin1 > 0):
           hist[bin1] = hist[bin1] + 1
           count = count +1
        #print(x[i])

    for i in range(binc):
        prob = hist[i]/count
        q = array[i]
        print(prob ,'\t' ,q)
        file1.write("%f %f\n" % (prob, q))
    plt.plot(array,hist/count,'r')
    plt.title('Local-bond-order-parameter-frequency-equi')
    plt.xlabel('avg q6')
    plt.ylabel('probability')
    plt.legend()
    plt.savefig(cwd + '/freq_q6(i)_equi32.png')
    file1.close()                                                               
#-------------------------------------------------------------------------------------
natom = 4000
L = 20.0
ln = 8
rcut = 1.5
  
                  
  
        
time =  [800000]

folder = 10
nf= 1
        
X = np.zeros((len(time),folder,natom*nf))
Y = np.zeros((len(time),folder,natom*nf))
Z = np.zeros((len(time),folder,natom*nf))

cwd = os.getcwd()
for fi in range(folder): #storing data for particular time fames of interest
	f = open(cwd +'/../' +str(fi+11) +'/Steadystate-Coordinates1.dat' , 'r')
	for t in range(len(time)):
		tr = time[t]  + 1
		nm = str(int(t*0.0001))
		t_stop = np.arange(time[t],tr,500)
		#print(t_stop)
		k = len(t_stop) - 1
		i = 0 # initiale point of particle array
		while True :
			s = f.readline().split()
			if not s: break 
			if int(s[4]) in t_stop :
				X[t,fi,i] = float(s[0]) - L*round(float(s[0])/L)
				Y[t,fi,i] = float(s[1]) - L*round(float(s[1])/L)
				Z[t,fi,i] = float(s[2]) - L*round(float(s[2])/L)
				d = int(s[3])
				i = i + 1
			if int(s[4]) > t_stop[k] : break
	f.close()

x = []
y = []
z = []
for t in range(len(time)):
	#Q6cor = np.zeros(max_bin)
	#RDF = np.zeros(max_bin)
	#corr_avg = []
	print(time[t])
	f2 = open(cwd + '/qs' +'_'+ str(time[t])+'.dat','w')
	for fi in range(folder):
		i = 0
		while True:
			if i == natom*nf : break
			x.append(X[t,fi,i])
			y.append(Y[t,fi,i])
			z.append(Z[t,fi,i])	
			i = i +1
			if i % natom == 0 :
				bndOrddatset = np.zeros((natom,ln))
				for li in range(1,ln+1):
				      cor = []
				      cor = bond_order(x,y,z,li)
				      for atom in range(natom):
					      bndOrddatset[atom,li-1] = cor[atom]
				for atom in range(natom):
					f2.write("%d %f" % (atom+1,x[atom]))
					for li in range(1,ln+1):
					     f2.write("\t %f" % (bndOrddatset[atom,li-1]))
					f2.write('\n')
				x.clear()
				y.clear()
				z.clear()  
	f2.close()        
