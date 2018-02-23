from AnalysisLib import *
from mpl_toolkits import mplot3d
import numpy as np


pdbfile = "last-r3.pdb"
group = "PC"
lowBound = 20
upBound = 40
plot = "True"
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
		CloseHoleArea = MatrxArea / TotLip
		print(layer + " Leaflet - Area Per Lipid - Away from Hole: %s" % AreaPerLip)
		print(layer + "Leaflet - Area Per Lipid - Close to Hole: %s" % CloseHoleArea)
	else:
		cellArea = MatrxArea / cellNum
		AreaPerLip = (MatrxArea - cellArea * (len(MaxZero)+len(MaxBorder)))/ (TotLip - lipNum) 
		CloseHoleArea = (cellArea * len(MaxBorder)) / lipNum
		print(layer + " Leaflet - Area Per Lipid - Away from Hole: %s" % AreaPerLip)
		print(layer + "Leaflet - Area Per Lipid - Close to Hole: %s" % CloseHoleArea)
