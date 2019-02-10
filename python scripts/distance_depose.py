# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 09:52:18 2017

@author: pablojensen
"""

#import matplotlib.pyplot as plt
import xmltodict
import os
import matplotlib.pyplot as plt
import numpy
import numpy.matlib
from scipy import stats
import networkx as nx
import scipy.sparse
    
def buildDiagonalMatrix(adjacencyMatrix):
    DiagonalMatrix=numpy.zeros(adjacencyMatrix.shape)
    for i in range(adjacencyMatrix.shape[1]):
        DiagonalMatrix[i,i]=adjacencyMatrix[i,:].sum()
    return numpy.matrix(DiagonalMatrix)

def buildDistance(Lplus):
    return numpy.matrix(n)
                 
####### DEMARRE, CALCUL COMMUTING TIME A PARTIR ADJ MATRIX

adjacencyMatrix=numpy.genfromtxt('Karate_FA_LinLog_grav0.csv',delimiter=';')

#print(adjacencyMatrix)
dim=adjacencyMatrix.shape[1]
label=numpy.zeros(dim-1)
#montre relation indice -> label
for i in range(0,dim-1):
   label[i]=adjacencyMatrix[0,i+1]

adjacencyMatrix=adjacencyMatrix[1:dim,1:dim] #enleve 1ere colonne/ligne (labels)
dim=adjacencyMatrix.shape[1]
print(dim)

#boucle sur noeud

###### calcule matrice distances (avec sqrt) 
vg=numpy.sum(adjacencyMatrix)
diagonalMatrix=buildDiagonalMatrix(adjacencyMatrix)
LaplacianMatrix=diagonalMatrix-adjacencyMatrix
minus0=numpy.matlib.ones((dim,dim))
minus=minus0/dim
invl=numpy.linalg.inv(LaplacianMatrix-minus)
Lplus=invl+minus
n=numpy.matlib.zeros((Lplus.shape[1],Lplus.shape[1]))
ndir=numpy.matlib.zeros((Lplus.shape[1],Lplus.shape[1]))

for i in range(Lplus.shape[1]):
    print 'i', i
    for j in range(Lplus.shape[1]):
        n[i,j]=numpy.sqrt(Lplus[i,i]+Lplus[j,j]-2*Lplus[i,j])


#        for k in range(Lplus.shape[1]):
#           ndir[j,i]+=(Lplus[j,k]-Lplus[j,i]-Lplus[i,k]+Lplus[i,i])*diagonalMatrix[k,k]

n=n*numpy.sqrt(vg)
#print n

#file=open("distn.csv","w")
#file.write(n)
#file.close()

########compute distance = shortest path length
print 'calcule plus courts chemins'
G = nx.DiGraph(adjacencyMatrix) 
#distsp=numpy.matlib.zeros((Lplus.shape[1],Lplus.shape[1]))

#il faudrait calculer j>i et symétriser 
distsp=numpy.zeros((dim,dim))
for i in range(Lplus.shape[1]):
    print i
    for j in range(i+1,Lplus.shape[1]):
        distsp[i,j] = nx.shortest_path_length(G, i, j)
#        print i,j,distsp[i,j]

print 'complete plus courts chemins'
for i in range(1,Lplus.shape[1]):
    print i
    for j in range(0,i):
        distsp[i,j] = distsp[j,i]


#file=open("distsp.csv","w")
#file.write(distsp)
#file.close()

src="Karate_FA_LinLog_grav0.gexf"
#G=nx.read_gexf(os.path.join(src))

with open(src) as fd:
   doc = xmltodict.parse(fd.read())

nodes=doc['gexf']['graph']['nodes']['node']


coord=numpy.zeros((dim,2))
i=0
for _n in nodes:
   coord[i,0]=float(_n['viz:position']['@x'])
   coord[i,1]=float(_n['viz:position']['@y'])
   i=i+1


dist=numpy.zeros((dim,dim))

#il faudrait calculer j>i et symétriser 
print 'calcule distances gephi'
for i in range(0,dim):
    print i  
    for j in range(i+1,dim):
          dist[i,j]=numpy.sqrt((coord[i,0]-coord[j,0])**2+(coord[i,1]-coord[j,1])**2)
           

print 'complete dist'
for i in range(1,dim):
    print i
    for j in range(0,i):
        dist[i,j] = dist[j,i]


# normalise distances pour pouvoir les comparer           
#distymean = numpy.mean(disty)
#disty /= distymean
#file=open("disteuclidien.csv","w")
#file.write(dist)
#file.close()


#for i in range(0,dim):
#    print i, 'dist11min',distrw[11,i],'dist11sym',n[11,i],'distYH',disty[11,i]

#agrege distances
#ICI CHOISIR Y ET Y1 DISTANCES A COMPARER (SP, RW, EUCLI)

length=dim*(dim-1)
dmax=numpy.amax(dist)
x=numpy.zeros(length)
y=numpy.zeros(length)
y1=numpy.zeros(length)
count=0
count1=0
count2=0
x1=numpy.zeros(length)
x2=numpy.zeros(length)

for i in range(0,dim):
   print i 
   for j in range(0,dim):
        if j!=i:
#            y[count]=distsp[i,j]
#            print 'dist', i,j,disty[i,j],ndir[i,j],distrw[i,j],n[i,j]
            x[count]=int(dist[i,j]/dmax*10)
            y[count]=n[i,j]
            y1[count]=distsp[i,j]
#	    if label[i]==33:
#            	print label[i],label[j],dist[i,j],n[i,j],distsp[i,j]
            if distsp[i,j]==1: 
                x1[count1]=dist[i,j]
                count1+=1
            if distsp[i,j]==2: 
                x2[count2]=dist[i,j]
                count2+=1
            count=count+1
           

csp=numpy.zeros(100)
isp=numpy.zeros(100)
isp2=numpy.zeros(100)
 
for count in range(0,length):
	isp[int(x[count])]=isp[int(x[count])]+y1[count]
	csp[int(x[count])]=csp[int(x[count])]+1
	isp2[int(x[count])]=isp2[int(x[count])]+y1[count]*y1[count]

for i in range(100):
	if csp[i] > 0:
		isp[i]=isp[i]/csp[i]
		isp2[i]=isp2[i]/csp[i]

for i in range(100):
	if csp[i] > 0:
	   print i, isp[i],numpy.sqrt(isp2[i]),numpy.sqrt(isp2[i]-isp[i]*isp[i]),csp[i]

slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
print ("R2 xy",r_value**2)
plt.plot(x, intercept + slope*x, 'b', label='fitted line')
plt.plot(x,y,'+r')
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y1)
print ("R2 xy1",r_value**2)
plt.plot(x,y1,'ob')
#plt.xlabel('euclidian distances')
#plt.ylabel('mean commuting time')
#plt.show()
slope, intercept, r_value, p_value, std_err = stats.linregress(y,y1)
print ("R2 yy1",r_value**2)
plt.plot(y, intercept + slope*y, 'r', label='fitted line')
plt.plot(y,y1,'or')
plt.show()
