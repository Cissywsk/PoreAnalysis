from AnalysisLib import *
from mpl_toolkits import mplot3d
import numpy as np

def APLwithHole(pdbfile):
	group = "PC"
	lowBound = 20
	upBound = 40
	plot = "False"
	binx = 10
	biny = 15

	centres = findcentre(pdbfile, group, visible=False)

	Xfull = centres[0]
	Yfull = centres[1]
	Zfull = centres[2]

	lowZ = [z for z in Zfull if z < lowBound]
	lowX = [x for x, z in zip(Xfull, Zfull) if z < lowBound]
	lowY = [y for y, z in zip(Yfull, Zfull) if z < lowBound]

	upZ = [z for z in Zfull if z > upBound]
	upX = [x for x, z in zip(Xfull, Zfull) if z > upBound]
	upY = [y for y, z in zip(Yfull, Zfull) if z > upBound]

	cellNum = binx * biny

	hLow, xedgesLow, yedgesLow, image = plt.hist2d(lowX, lowY, bins=(binx, biny))
	xEdLow = xedgesLow[-1] - xedgesLow[0]
	yEdLow = yedgesLow[-1] - yedgesLow[0]
	lowArea = xEdLow * yEdLow
	print("Lower Leaflet - HeadGroup Distribution")
	print("x-width: %s" % xEdLow)
	print("y-width: %s" % yEdLow)
	print("area: %s" % lowArea)
	print("HeadGroup Distribution")
	print(hLow)

	hUp, xedgesUp, yedgesUp, image = plt.hist2d(upX, upY, bins=(binx, biny))
	xEdUp = xedgesUp[-1] - xedgesUp[0]
	yEdUp = yedgesUp[-1] - yedgesUp[0]
	upArea = xEdUp * yEdUp
	print("Upper Leaflet - HeadGroup Distribution")
	print("x-width: %s" % xEdUp)
	print("y-width: %s" % yEdUp)
	print("area: %s" % upArea)
	print("HeadGroup Distribution")
	print(hUp)

	if plot == "True":
		ax = plt.axes(projection='3d')

		#overall PC scatter
		#ax.scatter3D(Xfull, Yfull, Zfull)

		#lower leaflet PC scatter - 3D
		ax.scatter3D(lowX, lowY, lowZ, color='r')
		ax.set_zlim(0,50)
		ax.set_title("Lower leaflet")

		#lower leaflet PC heapmap -- checking PC density
		fig = plt.figure()
		plt.hist2d(lowX, lowY, bins=(binx, biny))
		plt.title("Lower leaflet HeadGroup density heatmap")

		#upper leaflet PC scatter
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		ax.scatter3D(upX, upY, upZ, color='g')
		ax.set_zlim(0,50)
		ax.set_title("Upper leaflet")

		#upper leaflet PC heapmap -- checking PC density
		fig = plt.figure()
		plt.hist2d(upX, upY, bins=(binx, biny))
		plt.title("Upper leaflet HeadGroup density heatmap")

		plt.show()

	#Area per lipid of Lower Leaflet
	result = []
	for layer, matrix in zip(["Lower", "Upper"], [hLow, hUp]):
		if layer == "Lower":
			MatrxArea = lowArea
		elif layer == "Upper":
			MatrxArea = upArea
		MaxZero, MaxBorder, lipNum = findHole(matrix)
		TotLip = sum(map(sum, matrix))

		if len(MaxZero) == 1:
			print("No significant hole found")
			AreaPerLip = MatrxArea / TotLip
			result.append(AreaPerLip)
			CloseHoleArea = MatrxArea / TotLip
			result.append(CloseHoleArea)
			print(layer + " Leaflet - Area Per Lipid - Away from Hole: %s" % AreaPerLip)
			print(layer + "Leaflet - Area Per Lipid - Close to Hole: %s" % CloseHoleArea)
		else:
			cellArea = MatrxArea / cellNum
			AreaPerLip = (MatrxArea - cellArea * (len(MaxZero)+len(MaxBorder)))/ (TotLip - lipNum) 
			result.append(AreaPerLip)
			CloseHoleArea = (cellArea * len(MaxBorder)) / lipNum
			result.append(CloseHoleArea)
			print(layer + " Leaflet - Area Per Lipid - Away from Hole: %s" % AreaPerLip)
			print(layer + "Leaflet - Area Per Lipid - Close to Hole: %s" % CloseHoleArea)
	return result

def PoreDiameter(filename):
	inppdb = filename

	Xfull, Yfull, Zfull = findcentre(inppdb, "WAT")

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
		elif j == min(delta[4:16]) and i in range(4,16):
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
#	plt.show()

	return diaCyn
