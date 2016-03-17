#!/usr/bin/python3
# -*- coding: utf-8 -*

import argparse
from tkinter import *
import subprocess as sub
from PIL import ImageTk
from PIL import Image

parser = argparse.ArgumentParser(description="Given raw RNA sequences or Fasta RNA sequences can execute Nussinov algorithm or Zuker algorithm when specified and returns output")

parser.add_argument('-i', '--input', 
                    dest = "algorithm", 
                    action="store", 
                    default =None, 
                    required=True, 
                    help="Input must be a string writed in command line specifying the algorithm like: -i Nussinov or -i Zuker")
parser.add_argument('-o', '--output', 
                    dest = "output_file", 
                    action="store", 
                    default =None, 
                    required=True, 
                    help="Output must be a string writed in command line specifying the name where the 2D structure image will be stored, with a '.jpg' extension")

args = parser.parse_args()
alg = args.algorithm
out_file = args.output_file

if alg == 'Nussinov':
	p = sub.Popen('./nussinov_algorithm.py', stdout=sub.PIPE, stderr=sub.PIPE)
	output, errors = p.communicate()
if alg == 'Zuker':
	p = sub.Popen('./zuker_algorithm.py', stdout=sub.PIPE, stderr=sub.PIPE)
	output, errors = p.communicate()

out = output.decode('utf-8')
for element in out.split("\n"):
	if element.startswith((" A"," C"," G", " U")):
		sequence = element[1:]
	elif element.startswith((" .", " (")):
		dot_bracket = element[1:]


command = "\'"+"java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN \""+sequence+"\" -structureDBN \""+dot_bracket+"\" -o \""+out_file+"\"\'"
#print(command)
sub.call("'"+command+"'", shell=True)
root = Tk()
root.title('Results')

im = Image.open(out_file)
image_width = im.size[0]
image_height = im.size[1]

canvas = Canvas(root, width=image_width, height=image_height, bg= "white")
canvas.pack()


canvas.image = ImageTk.PhotoImage(im)
canvas.create_image(0, 0, image=canvas.image, anchor='nw')

frame=Frame(root)
frame.pack()

text_area = Text(frame)
text_area.pack()

text_area.insert(END, output, errors)


root.mainloop()

#sub.call('javaws ./VARNA-WebStart.jnlp', shell=True)

 