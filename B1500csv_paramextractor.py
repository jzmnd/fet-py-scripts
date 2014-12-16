# Reads B1500 csv VthMeasure files and summarises parameters
# Jeremy Smith
# Northwestern University
# Version 1.4


from numpy import *
import os
import sys
from myfunctions import *


data_path = os.path.dirname(__file__)    # Path name for location of script

print "\n"
print data_path
print "\n"

files = os.listdir(data_path)   # All files in directory
data_summary = []

skipinit = 5   # Initial data points to skip

# Loops through each Vth device file
for d in files:
	print d
	if "VthMeasure" not in d:
		continue
	datad, headd = csvImport(d, data_path, 12)
	l = float(headd[11][1])
	w = float(headd[11][2])
	ci = float(headd[11][5])
	# Data smoothing
	modIdrain_smoothed = adjAvSmooth(datad["modIdrain"], N=2)
	# Finding max transconductance
	diff_sqrt_modIdrain_smoothed = numDiff(sqrt(modIdrain_smoothed), datad["Vgate"])
	tmaxarg = argmax(diff_sqrt_modIdrain_smoothed[skipinit:-1]) + skipinit
	# Saturation mobility (max transconductance)
	satmob_tmax = (2*l/(w*ci))*(diff_sqrt_modIdrain_smoothed[tmaxarg])**2
	# Threshold Voltage (max transconductance)
	vth_tmax = datad["Vgate"][tmaxarg] - sqrt(modIdrain_smoothed)[tmaxarg]/diff_sqrt_modIdrain_smoothed[tmaxarg]
	# On-off ratio
	onoffratio = log10(max(modIdrain_smoothed[skipinit:-1])/min(modIdrain_smoothed[skipinit:-1]))
	# Subtheshold slope
	sts = min(abs(1/array(numDiff([log10(abs(x)) for x in modIdrain_smoothed[skipinit:-1]], datad["Vgate"][skipinit:-1]))))

	data_summary.append([headd[4][1], datad["SATMOB"][0], datad["Vth"][0], onoffratio, sts, satmob_tmax, vth_tmax])

outfile = open(os.path.join(data_path, "summarylist.txt"), "w")
outfile.write("Device\tSatMobility\tVthreshold\tLogOnOff\tSubSlope\tSatMobility(Tmax)\tVthreshold(Tmax)\n")
for a in data_summary:
	outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(a[0], a[1], a[2], a[3], a[4], a[5], a[6]))
outfile.close()
