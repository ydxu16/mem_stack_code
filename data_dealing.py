#!/usr/bin/python
import sys
from subprocess import call
from os import path

num_mem = 5
Nx = 64
Ny = 64

if not path.exists("./matlab_data_align"):
  call("mkdir matlab_data_align", shell = True)
if not path.exists("./tecplot_data_movie"):
  call("mkdir tecplot_data_movie", shell = True)

for i in range(num_mem):
  name_i = "file_membrane_" + str(i) + "_psi"
  file = open(name_i, 'r')
  mat_name ="./matlab_data_align/"+name_i
  tec_name ="./tecplot_data_movie/"+name_i
  matlab_file = open(mat_name, "w+")
  tecplot_file = open(tec_name, "w+")

  for line in file:
    words = line.split(',')
    if words[0] == "t":
      matlab_file.writelines(",\n")
      tecplot_file.writelines("ZONE "+ "I = "+ str(Nx) + ", J = "+ str(Ny)+"\n")
    else:
      matlab_file.writelines(line)
      tecplot_file.writelines(line)
  matlab_file.close()
  tecplot_file.close()


