import numpy as np
import matplotlib.pyplot as plt
from AnalysisLib import *
from mpl_toolkits import mplot3d
from sklearn.decomposition import PCA
from math import *

pdbfile = "last-r3.pdb"
group = "WAT"
plot = "False"
binx = 10
biny = 15

Xfull, Yfull, Zfull = findcentre(pdbfile, group, visible=False)

n, bins, _ = plt.hist(Zfull, 20)

#select the smallest diameter of the cylinder, find the boundaries 
delta = [ abs(i - j) for i,j in zip(n[:-1], n[1:])]
minBound = 10000
maxBound = -10000
for i, j in enumerate(delta):
	if j < 7.0 and i in range(4,16):
		minBin = bins[i]
		maxBin = bins[i+2]
		if minBound > minBin:
			minBound = minBin
		if maxBound < maxBin:
			maxBound = maxBin
plt.vlines(minBound, 0, 200)
plt.vlines(maxBound, 0, 200)
#print(n)
#print(delta)
#print(minBound)
#print(maxBound)

# select the water molecules X,Y,Z in X/Y/Zfilter in the pore -> within the range
Xfilter = []
Yfilter = []
Zfilter = []
detZ = maxBound - minBound 
for i, Z in enumerate(Zfull):
	if Z > minBound and Z < maxBound:
		Zfilter.append(Z)
		Xfilter.append(Xfull[i])
		Yfilter.append(Yfull[i])

#centre the pore
Xfilter = Xfilter - np.mean(Xfilter)
Yfilter = Yfilter - np.mean(Yfilter)
Zfilter = Zfilter - np.mean(Zfilter)

#stack them into one matrix, x,y,z coordinates in three different columns
X = np.vstack((Xfilter, Yfilter, Zfilter))
X = X.transpose()


#plot selected WAT in 3D - centred
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(X[:,0], X[:,1], X[:,2])


#PCA analysis to get PCs to describe the pore - limits: assuming the 3rd PC is the direction of the pore, as the pore is very round, flat and fat
pca = PCA(n_components=3)
pca.fit(X)
PCvector  = pca.components_
print("PCvector = %s " % PCvector)

for i in range(0,3):
	ax.plot([0, 10*PCvector[i][0]], [0, 10*PCvector[i][1]], [0, 10*PCvector[i][2]])
	if i == 2:
		COSangle = np.sqrt((PCvector[i][0])**2 + (PCvector[i][1])**2  ) / np.sqrt( PCvector[2].dot(PCvector[2]))

#get the angle from the direction vector
degreeAng = degrees(acos(COSangle))
print("Cylinder degree = %s" % degreeAng)

hCylinder = detZ / sin(radians(degreeAng))
print("Cylinder height = %s" % hCylinder)

#the equation is the number of WAT * Volume of WAT = pi * r^2 * height of cylinder
VolCynWAT = 30.52111 * len(Zfilter)  #Volume of WAT molecule * number of WAT molecule
radiusCyn = np.sqrt( VolCynWAT / (pi * hCylinder) )
diaCyn = radiusCyn * 2
print("Cylinder diameter = %s" % diaCyn)

if plot == "True":
	plt.show()
