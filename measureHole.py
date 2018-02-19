def getHole(p, q):
	r1 = p - 1
	r2 = p
	r3 = p + 1
	c1 = q - 1
	c2 = q
	c3 = q + 1

	if r3 == len(data):
		r3 = 0
	if c3 == len(data[i, :]):			
		c3 = 0

	UpLeft = data[r1, c1]
	Up = data[r1, c2]
	UpRight = data[r1, c3]
	Left = data[r2, c1]
	Right = data[r2, c3]
	DownLeft = data[r3, c1]
	Down = data[r3, c2]
	DownRight = data[r3, c3]
	


	neighbrVal = [Up,Left, Right, Down]
	neighbrCor = [[r1, c2], [r2, c1], [r2, c3],[r3, c2]]

	return [neighbrVal, neighbrCor]

def getBorder(p, q):
	r1 = p - 1
	r2 = p
	r3 = p + 1
	c1 = q - 1
	c2 = q
	c3 = q + 1
	
	if r3 >= len(data):
		r3 = 0
	if c3 >= len(data[i, :]):			
		c3 = 0

	UpLeft = data[r1, c1]
	Up = data[r1, c2]
	UpRight = data[r1, c3]
	Left = data[r2, c1]
	Right = data[r2, c3]
	DownLeft = data[r3, c1]
	Down = data[r3, c2]
	DownRight = data[r3, c3]


	neighbrVal = [UpLeft, UpRight,DownLeft, DownRight]
	neighbrCor = [[r1, c1], [r1, c3], [r3, c1],  [r3, c3]]

	return [neighbrVal, neighbrCor]



import numpy as np 

data = np.array([[ 1,  0,  2,  1,  2,  0,  2,  1,  1,  0,  1,  1,  1,  1,  1],
[0, 1, 0, 1, 1, 2, 1, 1, 1, 1, 1, 1, 0, 0, 1],
[1, 0, 1, 1, 2, 1, 0, 2, 1, 1, 0, 1, 0, 1, 1],
[2, 2, 1, 1, 1, 0, 0, 0, 2, 0, 1, 1, 0, 2, 1],
[2, 0, 1, 2, 1, 1, 0, 0, 1, 1, 2, 2, 1, 0, 0],
[1, 1, 2, 1, 0, 0, 1, 1, 0, 2, 1, 1, 1, 0, 2],
[2, 0, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1],
[1, 1, 1, 2, 1, 1, 1, 2, 0, 1, 3, 0, 1, 1, 1],
[3, 1, 2, 1, 3, 0, 2, 1, 1, 0, 0, 1, 1, 2, 0],
[0, 2, 1, 0, 0, 2, 0, 2, 2, 1, 1, 1, 1, 2, 1]])

MaxSize = 0


for i in range(0, len(data)):
	for j in range(0, len(data[i,:])):
		TotSize = 0
		ZeroCell = []
		if data[i, j] == 0:
			neighbour = getHole(i, j)
			neighbrVal = neighbour[0]
			neighbrCor = neighbour[1]
			CheckCell = [[i, j]]
			border = getBorder(i, j)[1]

			for val, cor in zip(neighbrVal, neighbrCor):
				if val == 0:
					if cor not in ZeroCell and cor not in CheckCell:
						ZeroCell.append(cor)
				else:
					if cor not in border:
						border.append(cor)

			while len(ZeroCell) != 0:
#				print(ZeroCell)
				currentCell = ZeroCell[0]
#				print(currentCell)
#				print(ZeroCell)
				CheckCell.append(currentCell)
				ZeroCell = ZeroCell[1:]
				neighbour = getHole(currentCell[0], currentCell[1])
				bord = getBorder(currentCell[0], currentCell[1])
				border.extend(bord[1])
#				print(neighbour)
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
#				print(k)
				if k not in CheckCell and k not in UniqueBorder:
					UniqueBorder.append(k)
					lipNum += data[k[0], k[1]]
			TotSize = len(CheckCell) + len(UniqueBorder)

		if TotSize > MaxSize:
			MaxSize = TotSize
			MaxZero = CheckCell
			MaxBorder = UniqueBorder
			print(MaxZero)
			print(MaxBorder)
			print(lipNum)