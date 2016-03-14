from tkinter import *
import subprocess as sub

p = sub.Popen('./nussinov_algorithm.py',stdout=sub.PIPE,stderr=sub.PIPE)
output, errors = p.communicate()

#sub.call('./VARNA-WebStart.jnlp')

print(errors)

root = Tk()
text = Text(root)
text.pack()
text.insert(END, output)
root.mainloop()