import numpy as np
import subprocess
from subprocess import call
from analysis import *

#inpf = open("kinks",'r')
#lengthf = open("cluster-size",'w')


startCA = 8
endCA = 25
count = 'Hmol'
start_frame = 4
stop_frame = 4000
step_frame = 4

result = open("results.dat", "w")
result.write("Frame")
result.write("	")
result.write("Diameter in A")
result.write("\n")

for f in range(start_frame, stop_frame+1, step_frame):
	with open('trajpdb', 'w') as outf:
		outf.write("parm ../%(count)s.prmtop\n"%locals())
		outf.write("trajin ../analysis-wrap_Prod20-isotropic_0025.nc %(f)s %(f)s" % locals())
		outf.write('\n')
		outf.write("trajout output.pdb pdb\n")
		outf.write("go\n")

	call(["cpptraj","-i","trajpdb"])

	dia = PoreDiameter("output.pdb")

	result.write("%s" % f)
	result.write("	")
	result.write("%s" % dia)
	result.write("\n")

result.close()