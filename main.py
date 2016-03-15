import argparse
from tkinter import *
import subprocess as sub

parser = argparse.ArgumentParser(description="Given raw RNA sequences or Fasta RNA sequences can execute Nussinov algorithm or Zuker algorithm when specified and returns output")

parser.add_argument('-i', '--input', 
                    dest = "algorithm", 
                    action="store", 
                    default =None, 
                    required=True, 
                    help="Input must be a string writed in command line specifing the algorithm like: -i Nussinov or -i Zuker")

args = parser.parse_args()
alg=args.algorithm


if alg == 'Nussinov':
	p = sub.Popen('./nussinov_algorithm.py', stdout=sub.PIPE, stderr=sub.PIPE)
	output, errors = p.communicate()
if alg == 'Zuker':
	p = sub.Popen('./zuker_algorithm.py', stdout=sub.PIPE, stderr=sub.PIPE)
	output, errors = p.communicate()


root = Tk()
root.title('Results')
text = Text(root)
text.pack()
text.insert(END, output)
root.mainloop()

sub.call('javaws ./VARNA-WebStart.jnlp', shell=True)

