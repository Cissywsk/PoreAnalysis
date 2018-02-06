import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
n = 0
Xtot = []
Ytot = []
Ztot = []
Xfull = []
Yfull = []
Zfull = []

#open original pdb, select all water molecules
with open("last-r3.pdb") as inp:
	inp.next()
	for line in inp:
		if "WAT" in line and "TER" not in line:
			data = line.split()
			Xcor = float(data[5])
			Ycor = float(data[6])
			Zcor = float(data[7])
			Xtot.append(Xcor)
			Ytot.append(Ycor)
			Ztot.append(Zcor)
			n = n + 1
			if n == 3:
				Xavg = 16 * Xtot[0] + 1 * Xtot[1] + 1 * Xtot[2]
				Yavg = 16 * Ytot[0] + 1 * Ytot[1] + 1 * Ytot[2]
				Zavg = 16 * Ztot[0] + 1 * Ztot[1] + 1 * Ztot[2]
				data[5] = str(Xavg / 18)
				data[6] = str(Yavg / 18)
				data[7] = str(Zavg / 18)
				Xfull.append(Xavg/18)
				Yfull.append(Yavg/18)
				Zfull.append(Zavg/18)
				Xtot = []
				Ytot = []
				Ztot = []
				n = 0

#draw histogram, understand the water distribution
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
print(n)
print(delta)
print(minBound)
print(maxBound)

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
from mpl_toolkits import mplot3d
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(X[:,0], X[:,1], X[:,2])


#PCA analysis to get PCs to describe the pore - limits: assuming the 3rd PC is the direction of the pore, as the pore is very round, flat and fat
from sklearn.decomposition import PCA
pca = PCA(n_components=3)
pca.fit(X)
PCvector  = pca.components_
print("PCvector = %s " % PCvector)

for i in range(0,3):
	ax.plot([0, 10*PCvector[i][0]], [0, 10*PCvector[i][1]], [0, 10*PCvector[i][2]])
	if i == 2:
		COSangle = np.sqrt((PCvector[i][0])**2 + (PCvector[i][1])**2  ) / np.sqrt( PCvector[2].dot(PCvector[2]))

#get the angle from the direction vector
from math import *
degreeAng = degrees(acos(COSangle))
print("Cylinder degree = %s" % degreeAng)

hCylinder = detZ / sin(radians(degreeAng))
print("Cylinder height = %s" % hCylinder)

#the equation is the number of WAT * Volume of WAT = pi * r^2 * height of cylinder
VolCynWAT = 30.52111 * len(Zfilter)  #Volume of WAT molecule * number of WAT molecule
radiusCyn = np.sqrt( VolCynWAT / (pi * hCylinder) )
diaCyn = radiusCyn * 2
print("Cylinder diameter = %s" % diaCyn)
plt.show()
