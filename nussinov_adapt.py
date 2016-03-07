import sys
from numpy import *
from matplotlib import *
import tkinter 

#read sequences
#sequences are stored in many lines
#f=open(sys.argv[1], 'r')


def delta(l,m):
 delta=0
 if l=='A' and m=='U':
  return 1
 elif l=='U' and m=='A':
  return 1
 elif l=='G' and m=='C':
  return 1
 elif l=='C' and m=='G':
  return 1
 else:
  return 0

def process():
 global entry_var

 seq=[]
 seq.append(entry_var.get())
 
 for q in range(0,len(seq)):
 	pair=traceback(fill_matrix(seq[q]),seq[q],0,len(seq[q])-1,[])
 	print ("max # of folding pairs: ",len(pair))
 	for x in range(0,len(pair)):
 		print ('%d %d %s==%s' % (pair[x][0],pair[x][1],pair[x][2],pair[x][3]))
 print ("---")

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
   pair.append([i,j,str(seq[i]),str(seq[j])])
   traceback(s,seq,i+1,j-1,pair)
  else:
   for k in range(i+1,j):
    if s[i,j]==s[i,k]+s[k+1,j]:
     traceback(s,seq,i,k,pair)
     traceback(s,seq,k+1,j,pair)
     break
 return pair


top = tkinter.Tk()

label = tkinter.Label(top, text="Enter a RNA sequence:")
label.pack( side=tkinter.LEFT )
entry_var = tkinter.StringVar()
entry = tkinter.Entry(top, textvariable=entry_var)
entry.pack( side=tkinter.LEFT )
b = tkinter.Button(top, text="Submit", command=process)
b.pack( side=tkinter.LEFT )

top.mainloop()