#!/usr/bin/python3
# -*- coding: utf-8 -*

import sys
from numpy import *
from matplotlib import *
from tkinter import *
from check_RNA import *


##PROCESSING INPUT AND OUTPUT

def process():
	
	global text_area
	seq=[]
	iden=[]

	#Processing the user input

	lines = text_area.get("1.0", END).splitlines()
	for things in FASTA_iterator(lines):
		seq.append(things.get_sequence())
		iden.append(things.get_identifier())

	#Calling the Nussinov algotihm and printing the results		

	for q in range(0,len(seq)):
		pair=traceback(fill_matrix(seq[q]),seq[q],0,len(seq[q])-1,[])
		print ("max number of folding pairs: ",len(pair))
		dot_bracket = list("."*len(seq[q]))
		for x in range(0,len(pair)):
			print ('%d %d %s==%s' % (pair[x][0],pair[x][1],pair[x][2],pair[x][3]))
			dot_bracket[pair[x][0]]="("
			dot_bracket[pair[x][1]]=")"

		#Printing dot-bracket
		print ("\nDot-bracket representation\n",str(seq[q])+"\n",''.join(dot_bracket))
	print ("\n---")
	
	window.destroy()

def FASTA_iterator(lines):
 	sequence = ""
 	identifier = "query"
 	for line in lines:
 		if line.startswith(">"):
 			if len(sequence)>0:
 				try:
 					yield RNASequence(identifier, sequence)
 				except IncorrectSequenceLetter as e:
 					sys.stderr.write(str(e)+"\n")
 			identifier = line[1:].strip()
 			sequence = ""
 		else:
 			sequence+=line.strip()		
 	if len(sequence)>0:
 		try:
 			yield RNASequence(identifier, sequence)
 		except IncorrectSequenceLetter as e:
 			sys.stderr.write(str(e)+"\n")

##NUSSINOV ALGORITHM

#Delta scores for pairs

def delta(l,m):
	delta=0
	if l=='A' and m=='U':
		return 2
	elif l=='U' and m=='A':
		return 2
	elif l=='G' and m=='C':
		return 3
	elif l=='C' and m=='G':
		return 3
	elif l=='G' and m=='U':
		return 1
	elif l=='C' and m=='U':
		return 1
	else:
		return 0

def fill_matrix(seq):

#Initialization
	L=len(seq)
	s=zeros((L,L))

#Recursion
	for n in range(1,L):
		for j in range(n,L):
			i=j-n
			paired=s[i+1,j-1]+delta(seq[i],seq[j])
			unpair_i=s[i+1,j]
			unpair_j=s[i,j-1]
			if i+3<=j:
				tmp=[]
				for k in range(i+1,j):
					tmp.append(s[i,k]+s[k+1,j])
				bifu=max(tmp)
				s[i,j]=max(paired,unpair_i,unpair_j,bifu)
			else:
				s[i,j]=max(paired,unpair_i,unpair_j)
	return s



#Traceback
def traceback(s,seq,i,j,pair):
	if i<j:
		if s[i,j]==s[i+1,j]:
			traceback(s,seq,i+1,j,pair)
		elif s[i,j]==s[i,j-1]:
			traceback(s,seq,i,j-1,pair)
		elif s[i,j]==s[i+1,j-1]+delta(seq[i],seq[j]):
			if (j-i >= 4):
				pair.append([i,j,str(seq[i]),str(seq[j])])
			traceback(s,seq,i+1,j-1,pair)
		else:
			for k in range(i+1,j):
				if s[i,j]==s[i,k]+s[k+1,j]:
					traceback(s,seq,i,k,pair)
					traceback(s,seq,k+1,j,pair)
					break
	return pair

#GRAPHICAL INTERFACE INPUT

window = Tk()
window.title("Title of the program")

label = Label(window, text="Enter a RNA sequence:")
label.pack( side = TOP)

frame=Frame(window)
frame.pack()

text_area = Text(frame)
text_area.pack()

b = Button(window, text="Submit", command=process)
b.pack( side=LEFT )

window.mainloop()
