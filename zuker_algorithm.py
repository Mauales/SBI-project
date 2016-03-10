import numpy as np
import sys

seq = input("Introduce sequence: ")

pathmatrix=np.zeros((len(seq),len(seq),2))

stacking_energy= [	[-0.9, -1.8, -2.3, -1.1, -1.1, -0.8],
					[-1.7, -2.9, -3.4, -2.3, -2.1, -1.4],
					[-2.1, -2.0, -2.9, -1.8, -1.9, -1.2],
					[-0.9, -1.7, -2.1, -0.9, -1.0, -0.5],
					[-0.5, -1.2, -1.4, -0.8, -0.4, -0.2],
					[-1.0, -1.9, -2.1, -1.1, -1.5, -0.4]]
# A/U = 0, C/G = 1, G/C = 2, U/A = 3, G/U = 4, U/G = 5.

inf = 999
not_calculated = 0.01

destabilizing_energy=[	[inf, 5.3, 6.6, 7.0, 7.4],
						[3.9, 4.8, 5.5, 6.3, 6.7],
						[inf, 4.4, 5.3, 6.1, 6.5]]
			            # 1 ,   5,  10,  20,  30				
# interior[0], 
# bulge[1], 
# hairpin[2].

V = np.nan*np.ones((len(seq),len(seq)))
W = np.nan *np.ones((len(seq),len(seq)))
loop_type = np.zeros((len(seq),len(seq)))

def initialize_everything():
	for i in range(0,len(seq)):
		for j in range(0,len(seq)):
			V[i][j]= not_calculated
			W[i][j]= not_calculated
			loop_type[i][j]= -5 # -5 don't pair
			pathmatrix[i][j][0]= -1
			pathmatrix[i][j][1]= -1
	for i in range(0,len(seq)):
		V[i][i] = inf
		W[i][i] = inf
	for i in range(0,len(seq)-1):
		V[i][i+1] = inf
		W[i][i+1] = inf

def do_basepair(i,j):
	if i == "A" and j == "U":
		return 1
	elif i == "U" and j == "A":
		return 1
	elif i == "G" and j == "U":
		return 1
	elif i == "U" and j == "G":
		return 1
	elif i == "G" and j == "C":
		return 1
	elif i == "C" and j == "G":
		return 1	
	else:
		return 0

def hairpin_loop(i, j):
	hairpin_nucleotides = abs(i-j) + 1
	if hairpin_nucleotides <= 4:
		return inf
	if hairpin_nucleotides == 5:
		hairpin_energy = destabilizing_energy[2][1]

	elif hairpin_nucleotides <= 10:
		interpolation = ((10 - hairpin_nucleotides)*((destabilizing_energy[2][2]- destabilizing_energy[2][1])/5))
		hairpin_energy = destabilizing_energy[2][2] - interpolation

	elif hairpin_nucleotides <= 20:
		interpolation = ((20 - hairpin_nucleotides)*((destabilizing_energy[2][3]- destabilizing_energy[2][2])/10))
		hairpin_energy = destabilizing_energy[2][3] - interpolation

	elif hairpin_nucleotides <= 30:
		interpolation = ((30 - hairpin_nucleotides)*((destabilizing_energy[2][4]- destabilizing_energy[2][3])/10))
		hairpin_energy = destabilizing_energy[2][4] - interpolation
	return hairpin_energy	

def getindex_stacking(row,column):
	if seq[row] == "A" and seq[column] == "U":
		return 0
	elif seq[row] == "C" and seq[column] == "G":
		return 1
	elif seq[row] == "G" and seq[column] == "C":
		return 2
	elif seq[row] == "U" and seq[column] == "A":
		return 3
	elif seq[row] == "G" and seq[column] == "U":
		return 4
	elif seq[row] == "U" and seq[column] == "G":
		return 5
	else:
		return -1

def stacking_loop (row, column):
	column_pointer = getindex_stacking(row,column)
	row_pointer = getindex_stacking(row+1, column-1)

	stacking_energy_value = stacking_energy[row_pointer][column_pointer]

	return stacking_energy_value

def bulge_loop(row1, column1, row2, column2):
	if abs(row2-row1)>abs(column1-column2) or abs(row2-row1)<abs(column1-column2):
		bulge_nucleotide = abs(abs(row2-row1)-abs(column1-column2))

		if bulge_nucleotide <= 5:
			interpolation = ((5-bulge_nucleotide)*(destabilizing_energy[1][1]-destabilizing_energy[1][0])/5)
			bulge_energy = destabilizing_energy[1][1]-interpolation
		elif bulge_nucleotide <= 10:
			interpolation = ((10-bulge_nucleotide)*(destabilizing_energy[1][2]-destabilizing_energy[1][1])/5)
			bulge_energy = destabilizing_energy[1][2]-interpolation
		elif bulge_nucleotide <= 20:
			interpolation = ((20-bulge_nucleotide)*(destabilizing_energy[1][3]-destabilizing_energy[1][2])/10)
			bulge_energy = destabilizing_energy[1][3]-interpolation
		elif bulge_nucleotide <= 30:
			interpolation = ((30-bulge_nucleotide)*(destabilizing_energy[1][4]-destabilizing_energy[1][3])/10)
			bulge_energy = destabilizing_energy[1][4]-interpolation
	return bulge_energy

def interior_loop(row1, column1, row2, column2):
	interior_nucleotides = (row2-row1)+(column1-column2)+2

	if interior_nucleotides == 5:
		interior_energy = destabilizing_energy[0][1]
	elif interior_nucleotides <= 10:
		interpolation = ((10 - interior_nucleotides)*(destabilizing_energy[0][2]-destabilizing_energy[0][1])/5)
		interior_energy=destabilizing_energy[0][2] - interpolation
	elif interior_nucleotides <= 20:
		interpolation = ((20 - interior_nucleotides)*(destabilizing_energy[0][3]-destabilizing_energy[0][2])/10)
		interior_energy=destabilizing_energy[0][3] - interpolation	
	elif interior_nucleotides <= 30:
		interpolation = ((30 - interior_nucleotides)*(destabilizing_energy[0][4]-destabilizing_energy[0][3])/10)
		interior_energy=destabilizing_energy[0][4] - interpolation	
	return interior_energy

def calculate_V(i, j):
	if do_basepair(seq[i],seq[j]) == 0:
		V[i][j] = inf
		return inf
	elif V[i][j] != not_calculated:
		return V[i][j]
	#case 1 FH(i,j)
	FH = minimum_energy = hairpin_loop(i,j)

	#case 2 min[FL(i,j,h,k) + V(h,k)]

	for h in range (i+1,j):
		for k in range (j-1, h, -1):
			if do_basepair(seq[h], seq[k]):
				#print (h,seq[h],"--",seq[k],k)
				if h == i+1 and k == j-1:
					temp_energy = stacking_loop(i,j)+ calculate_V(h,k)
					if (temp_energy < minimum_energy):
						minimum_energy = temp_energy
						pathmatrix[i][j][0] = h
						pathmatrix[i][j][1] = k
						loop_type[i][j] = -1 #stacking
				elif h == i+1 or k == j-1:
					temp_energy = bulge_loop(i,j,h,k) + calculate_V(h,k)
					if (temp_energy < minimum_energy):
						minimum_energy = temp_energy
						pathmatrix[i][j][0] = h
						pathmatrix[i][j][1] = k
						loop_type[i][j] = -2 #for bulge
				else:
					temp_energy = interior_loop(i,j,h,k) + calculate_V(h,k)
					if (temp_energy < minimum_energy):
						minimum_energy = temp_energy
						pathmatrix[i][j][0] = h
						pathmatrix[i][j][1] = k
						loop_type[i][j] = -3 #for interior	
	
	#case 3: min[W(i+1,k) + W(k+1, j-1)] condition i+1 < k < j-1

	for k in range (i+2, j-1):
		temp_energy = W[i+1][k] + W[k+1][j-1]
		if (temp_energy < minimum_energy):
			minimum_energy = temp_energy
			pathmatrix[i][j][0] = k
			pathmatrix[i][j][1] = k
			loop_type[i][j] = -4 #bifurcation

	V[i][j] = minimum_energy

	return minimum_energy

def calculate_W(i, j):
	if W[i][j] != not_calculated:
		return W[i][j]
	minimum_energy = calculate_V(i,j)

	if minimum_energy > calculate_W(i+1, j):
		minimum_energy = W[i+1][j]
		loop_type[i][j] = -5 # don't basepair
	elif minimum_energy > calculate_W(i,j-1):
		minimum_energy = W[i][j-1]
		loop_type[i][j] = -5 # don't basepair
	else:
		for k in range(i+1, j):
			if minimum_energy > (calculate_W(i,k)+calculate_W(k+1,j)):
				minimum_energy = W[i][k] + W[k+1][j]
				loop_type[i][j]

	W[i][j] = minimum_energy
	return minimum_energy			

initialize_everything()
for n in range(1,len(seq)):
		for j in range(n,len(seq)):
			i=j-n
			calculate_V(i,j)
			calculate_W(i,j)
np.set_printoptions(suppress=True)
sys.stdout = open("zuker.txt","w")
print (V)
print (W)