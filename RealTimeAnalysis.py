import numpy as np
import subprocess
from subprocess import call
from BilayerLib import *

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
result.write("		")
result.write("Diameter (A)")
result.write("		")
result.write("LowAPL away Hole A^2")
result.write("		")
result.write("LowAPL close Hole A^2")
result.write("		")
result.write("UpAPL away Hole A^2")
result.write("		")
result.write("UpAPL close Hole A^2")
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
	APLsummary = APLwithHole("output.pdb")


	result.write("%s" % f)
	result.write("		")
	result.write("%s" % dia)
	result.write("		")
	result.write("%s" % APLsummary[0])
	result.write("		")
	result.write("%s" % APLsummary[1])
	result.write("		")
	result.write("%s" % APLsummary[2])
	result.write("		")
	result.write("%s" % APLsummary[3])
	result.write("\n")

result.close()