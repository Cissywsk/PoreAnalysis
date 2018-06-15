import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def findcentre(pdbfile, group, visible=False):
	"""
	This function is to find the centre of mass of a particular residue e.g. WAT or PC and give a overall 3D visualisation if visible=True. 
	
	Usage: data = findcentre("last-r3.pdb", "PC", visible=True), 
	pdbfile: pdbfile name you need to analyse
	group: the particular functional group
	visible: either True or False, if you need to visualise the distribution of your centred residues

	Output:
	The output is a list [Xfull, Yfull, Zfull]
	"""

	n = 0
	sumX = 0
	sumY = 0
	sumZ = 0
	TotMass = 0
	Xfull = []
	Yfull = []
	Zfull = []

	if group == "PC":
		groupAtomNum = 38
	elif group == "WAT":
		groupAtomNum = 3
	else:
		print("Unknown Group, please update the total number of atoms here.")

	with open(pdbfile) as inp:
		inp.next()
		for line in inp:
			if group in line and "TER" not in line:
				data = line.split()
				Xcor = float(data[5])
				Ycor = float(data[6])
				Zcor = float(data[7])
				AtomType = data[10]
				if AtomType == 'C':
					mass = 12
				elif AtomType == 'H':
					mass = 1
				elif AtomType == 'O':
					mass = 16
				elif AtomType == "N":
					mass = 14
				elif AtomType == "P":
					mass = 31
				else:
					print("Cannot find AtomType")
				sumX = sumX + Xcor * mass
				sumY = sumY + Ycor * mass
				sumZ = sumZ + Zcor * mass
				TotMass = TotMass + mass
	#			print(mass)
				n += 1
				if n == groupAtomNum: # the total number of atoms in PC head group
					TotMass = TotMass - mass
					sumX = sumX - Xcor * mass
					sumY = sumY - Ycor * mass
					sumZ = sumZ - Zcor * mass
					centreX = sumX / TotMass
					centreY = sumY / TotMass
					centreZ = sumZ / TotMass

					n = 0
					TotMass = 0
					sumX = 0
					sumY = 0
					sumZ = 0
					Xfull.append(centreX)
	#				print(Xfull)
					Yfull.append(centreY)
					Zfull.append(centreZ)
	result = [Xfull, Yfull, Zfull]

	if visible == True:
		from mpl_toolkits import mplot3d
		fig = plt.figure()
		ax = plt.axes(projection='3d')

		#overall PC scatter
		ax.scatter3D(Xfull, Yfull, Zfull)
		plt.show()

#	print("Find x,y,z of the centre of mass of the %s" % group)
#	print("Output shape is [Xfull, Yfull, Zfull]")
	return result

def findNeighbour(p, q, matrix, full=False):
	"""
	This function is to find the 4 (default: up, down, left, right) or 8 (if full = True) neighbours of a particular cell(p,q) in a 2D matrix
	
	This function is on a periodic boundary condition, which means if the neighbour is on an edge, it uses the values from the otherside. 

	Usage: data = findNeighbour(p, q, full=True)
	p: the row number of the particular cell
	q: the column number of the particular cell
	full: optional. True: if you need all 8 neighbours. Default: 4 neighbours

	Output: 
	[neighbrVal, neighbrCor]
	neighbrVal: the values of the neighbours. 
				default: [Up,Left, Right, Down]
				full = True: [Up, Left, Right, Down, UpLeft, UpRight, DownLeft, DownRight]
	neighbrCor: the row/column number of the neighbours. 

	"""
	r1 = p - 1
	r2 = p
	r3 = p + 1
	c1 = q - 1
	c2 = q
	c3 = q + 1

	if r3 == len(matrix):
		r3 = 0
	if c3 == len(matrix[p, :]):			
		c3 = 0

	UpLeft = matrix[r1, c1]
	Up = matrix[r1, c2]
	UpRight = matrix[r1, c3]
	Left = matrix[r2, c1]
	Right = matrix[r2, c3]
	DownLeft = matrix[r3, c1]
	Down = matrix[r3, c2]
	DownRight = matrix[r3, c3]

	if full == True:
		neighbrVal = [Up, Left, Right, Down, UpLeft, UpRight,DownLeft, DownRight]
		neighbrCor = [[r1, c2], [r2, c1], [r2, c3],[r3, c2], [r1, c1], [r1, c3], [r3, c1], [r3, c3]]
	else:
		neighbrVal = [Up, Left, Right, Down]
		neighbrCor = [[r1, c2], [r2, c1], [r2, c3],[r3, c2]]
			
	return [neighbrVal, neighbrCor]

def findHole(matrix):
	"""
	Use the First Depth Search to find the biggest continuous zero cells in a 2D matrix


	"""

	def findNeighbour(p, q, matrix2D, full=False):
		r1 = p - 1
		r2 = p
		r3 = p + 1
		c1 = q - 1
		c2 = q
		c3 = q + 1

		if r3 == len(matrix2D):
			r3 = 0
		if c3 == len(matrix2D[p, :]):			
			c3 = 0

		UpLeft = matrix2D[r1, c1]
		Up = matrix2D[r1, c2]
		UpRight = matrix2D[r1, c3]
		Left = matrix2D[r2, c1]
		Right = matrix2D[r2, c3]
		DownLeft = matrix2D[r3, c1]
		Down = matrix2D[r3, c2]
		DownRight = matrix2D[r3, c3]

		if full == True:
			neighbrVal = [Up, Left, Right, Down, UpLeft, UpRight,DownLeft, DownRight]
			neighbrCor = [[r1, c2], [r2, c1], [r2, c3],[r3, c2], [r1, c1], [r1, c3], [r3, c1], [r3, c3]]
		else:
			neighbrVal = [Up, Left, Right, Down]
			neighbrCor = [[r1, c2], [r2, c1], [r2, c3],[r3, c2]]
				
		return [neighbrVal, neighbrCor]

	MaxSize = 0

	RefMatrix = np.zeros((len(matrix), len(matrix[0])))

	for i in range(0, len(matrix)):
		for j in range(0, len(matrix[i,:])):
			if RefMatrix[i, j] == 0:
				print RefMatrix
				RefMatrix[i, j] = 1

				TotSize = 0
				ZeroCell = []
				if matrix[i, j] == 0:
					neighbour = findNeighbour(i, j, matrix)
					neighbrVal = neighbour[0]
					neighbrCor = neighbour[1]
					CheckCell = [[i, j]]
					border = findNeighbour(i, j, matrix, full=True)[1][4:]


					for val, cor in zip(neighbrVal, neighbrCor):
						if val == 0:
							if cor not in ZeroCell and cor not in CheckCell:
								ZeroCell.append(cor)
						else:
							if cor not in border:
								border.append(cor)


					while len(ZeroCell) != 0:

						currentCell = ZeroCell[0]
						CheckCell.append(currentCell)

						RefMatrix[currentCell] = 1
						print(RefMatrix)
						ZeroCell = ZeroCell[1:]
						neighbour = findNeighbour(currentCell[0], currentCell[1], matrix)
						bord = findNeighbour(currentCell[0], currentCell[1], matrix, full=True)[1][4:]
						border.extend(bord)
						neighbrVal = neighbour[0]
						neighbrCor = neighbour[1]

						for val, cor in zip(neighbrVal, neighbrCor):
							if val == 0:
								if cor not in ZeroCell and cor not in CheckCell:
									ZeroCell.append(cor)
							else:
								if cor not in border:
									border.append(cor)
					
					UniqueBorder = []
					lipNum = 0

					for k in border:
						if k not in CheckCell and k not in UniqueBorder:
							UniqueBorder.append(k)
							lipNum += matrix[k[0], k[1]]
					TotSize = len(CheckCell) + len(UniqueBorder)

				if TotSize > MaxSize:
					MaxSize = TotSize
					MaxZero = CheckCell
					MaxBorder = UniqueBorder
					MaxLip = lipNum
	return MaxZero, MaxBorder, MaxLip

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def findcentre(pdbfile, group, visible=False):
	"""
	This function is to find the centre of mass of a particular residue e.g. WAT or PC and give a overall 3D visualisation if visible=True. 
	
	Usage: data = findcentre("last-r3.pdb", "PC", visible=True), 
	pdbfile: pdbfile name you need to analyse
	group: the particular functional group
	visible: either True or False, if you need to visualise the distribution of your centred residues

	Output:
	The output is a list [Xfull, Yfull, Zfull]
	"""

	n = 0
	sumX = 0
	sumY = 0
	sumZ = 0
	TotMass = 0
	Xfull = []
	Yfull = []
	Zfull = []

	if group == "PC":
		groupAtomNum = 38
	elif group == "WAT":
		groupAtomNum = 3
	else:
		print("Unknown Group, please update the total number of atoms here.")

	with open(pdbfile) as inp:
		inp.next()
		for line in inp:
			if group in line and "TER" not in line:
				data = line.split()
				Xcor = float(data[5])
				Ycor = float(data[6])
				Zcor = float(data[7])
				AtomType = data[10]
				if AtomType == 'C':
					mass = 12
				elif AtomType == 'H':
					mass = 1
				elif AtomType == 'O':
					mass = 16
				elif AtomType == "N":
					mass = 14
				elif AtomType == "P":
					mass = 31
				else:
					print("Cannot find AtomType")
				sumX = sumX + Xcor * mass
				sumY = sumY + Ycor * mass
				sumZ = sumZ + Zcor * mass
				TotMass = TotMass + mass
	#			print(mass)
				n += 1
				if n == groupAtomNum: # the total number of atoms in PC head group
					TotMass = TotMass - mass
					sumX = sumX - Xcor * mass
					sumY = sumY - Ycor * mass
					sumZ = sumZ - Zcor * mass
					centreX = sumX / TotMass
					centreY = sumY / TotMass
					centreZ = sumZ / TotMass

					n = 0
					TotMass = 0
					sumX = 0
					sumY = 0
					sumZ = 0
					Xfull.append(centreX)
	#				print(Xfull)
					Yfull.append(centreY)
					Zfull.append(centreZ)
	result = [Xfull, Yfull, Zfull]

	if visible == True:
		from mpl_toolkits import mplot3d
		fig = plt.figure()
		ax = plt.axes(projection='3d')

		#overall PC scatter
		ax.scatter3D(Xfull, Yfull, Zfull)
		plt.show()

#	print("Find x,y,z of the centre of mass of the %s" % group)
#	print("Output shape is [Xfull, Yfull, Zfull]")
	return result

def findNeighbour(p, q, matrix, full=False):
	"""
	This function is to find the 4 (default: up, down, left, right) or 8 (if full = True) neighbours of a particular cell(p,q) in a 2D matrix
	
	This function is on a periodic boundary condition, which means if the neighbour is on an edge, it uses the values from the otherside. 

	Usage: data = findNeighbour(p, q, full=True)
	p: the row number of the particular cell
	q: the column number of the particular cell
	full: optional. True: if you need all 8 neighbours. Default: 4 neighbours

	Output: 
	[neighbrVal, neighbrCor]
	neighbrVal: the values of the neighbours. 
				default: [Up,Left, Right, Down]
				full = True: [Up, Left, Right, Down, UpLeft, UpRight, DownLeft, DownRight]
	neighbrCor: the row/column number of the neighbours. 

	"""
	r1 = p - 1
	r2 = p
	r3 = p + 1
	c1 = q - 1
	c2 = q
	c3 = q + 1

	if r3 == len(matrix):
		r3 = 0
	if c3 == len(matrix[p, :]):			
		c3 = 0

	UpLeft = matrix[r1, c1]
	Up = matrix[r1, c2]
	UpRight = matrix[r1, c3]
	Left = matrix[r2, c1]
	Right = matrix[r2, c3]
	DownLeft = matrix[r3, c1]
	Down = matrix[r3, c2]
	DownRight = matrix[r3, c3]

	if full == True:
		neighbrVal = [Up, Left, Right, Down, UpLeft, UpRight,DownLeft, DownRight]
		neighbrCor = [[r1, c2], [r2, c1], [r2, c3],[r3, c2], [r1, c1], [r1, c3], [r3, c1], [r3, c3]]
	else:
		neighbrVal = [Up, Left, Right, Down]
		neighbrCor = [[r1, c2], [r2, c1], [r2, c3],[r3, c2]]
			
	return [neighbrVal, neighbrCor]

def findHole(matrix):
	"""
	Use the First Depth Search to find the biggest continuous zero cells in a 2D matrix


	"""

	def findNeighbour(p, q, matrix2D, full=False):
		r1 = p - 1
		r2 = p
		r3 = p + 1
		c1 = q - 1
		c2 = q
		c3 = q + 1

		if r3 == len(matrix2D):
			r3 = 0
		if c3 == len(matrix2D[p, :]):			
			c3 = 0

		UpLeft = matrix2D[r1, c1]
		Up = matrix2D[r1, c2]
		UpRight = matrix2D[r1, c3]
		Left = matrix2D[r2, c1]
		Right = matrix2D[r2, c3]
		DownLeft = matrix2D[r3, c1]
		Down = matrix2D[r3, c2]
		DownRight = matrix2D[r3, c3]

		if full == True:
			neighbrVal = [Up, Left, Right, Down, UpLeft, UpRight,DownLeft, DownRight]
			neighbrCor = [[r1, c2], [r2, c1], [r2, c3],[r3, c2], [r1, c1], [r1, c3], [r3, c1], [r3, c3]]
		else:
			neighbrVal = [Up, Left, Right, Down]
			neighbrCor = [[r1, c2], [r2, c1], [r2, c3],[r3, c2]]
				
		return [neighbrVal, neighbrCor]

	MaxSize = 0

	for i in range(0, len(matrix)):
		for j in range(0, len(matrix[i,:])):
			TotSize = 0
			ZeroCell = []
			if matrix[i, j] == 0:
				neighbour = findNeighbour(i, j, matrix)
				neighbrVal = neighbour[0]
				neighbrCor = neighbour[1]
				CheckCell = [[i, j]]
				border = findNeighbour(i, j, matrix, full=True)[1][4:]


				for val, cor in zip(neighbrVal, neighbrCor):
					if val == 0:
						if cor not in ZeroCell and cor not in CheckCell:
							ZeroCell.append(cor)
					else:
						if cor not in border:
							border.append(cor)


				while len(ZeroCell) != 0:

					currentCell = ZeroCell[0]
					CheckCell.append(currentCell)
					ZeroCell = ZeroCell[1:]
					neighbour = findNeighbour(currentCell[0], currentCell[1], matrix)
					bord = findNeighbour(currentCell[0], currentCell[1], matrix, full=True)[1][4:]
					border.extend(bord)
					neighbrVal = neighbour[0]
					neighbrCor = neighbour[1]

					for val, cor in zip(neighbrVal, neighbrCor):
						if val == 0:
							if cor not in ZeroCell and cor not in CheckCell:
								ZeroCell.append(cor)
						else:
							if cor not in border:
								border.append(cor)
				
				UniqueBorder = []
				lipNum = 0

				for k in border:
					if k not in CheckCell and k not in UniqueBorder:
						UniqueBorder.append(k)
						lipNum += matrix[k[0], k[1]]
				TotSize = len(CheckCell) + len(UniqueBorder)

			if TotSize > MaxSize:
				MaxSize = TotSize
				MaxZero = CheckCell
				MaxBorder = UniqueBorder
				MaxLip = lipNum
	return MaxZero, MaxBorder, MaxLip

