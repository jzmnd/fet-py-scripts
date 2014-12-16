# Reads B1500 csv Calculates the log of the effective mobility in the linear regime
# Jeremy Smith
# Northwestern University
# Version 1.3


from numpy import *
from scipy.interpolate import interp1d
import os
import sys
from myfunctions import *


data_path = os.path.dirname(__file__)    # Path name for location of script

print "\n"
print data_path
print "\n"

files = os.listdir(data_path)   # All files in directory
data_summary_VG = []
data_summary_VGVT = []

vonset = paramImport("DEVICEPROP_temp.txt", data_path, param_no=1)[0]   # Onset taken voltage from DEVICEPROP file
vglist = range(40, 82, 2)                                          # Choice of gate voltages to extract mobility
vgvtlist = array([1, 5, 10, 15, 20, 25, 30, 40, 50])               # Choice of Vg-Vt values to extract mobility
sweep = "fwd"                                                      # Extract from forward (fwd) or reverse (rev) sweep
vt = paramImport("summarylist.txt", data_path, param_no=4)[1]      # Threshold voltages from satmob vt data

# Loops through each transfer device file
for d in files:
	print d
	if "Transfer" not in d:
		continue
	datad, headd = csvImport(d, data_path, 12)
	
	s = len(datad["Vgate"])/2                # Length of data set (linear only)
	vg = array(datad["Vgate"][:s])           # Gate voltages
	idl = array(datad["modIdrain"][:s])      # Linear Id
	ids = array(datad["modIdrain"][s:])      # Saturation Id

	vd = float(datad["Vdrain"][0])           # Drain voltage in linear
	chW = float(headd[11][2])                # Channel W
	chL = float(headd[11][1])                # Channel L
	Ci = float(headd[11][5])                 # Geometric capacitance of dielectric

	mobEFF = log(abs(chL*idl/(chW*Ci*vd*(vg - vonset[headd[4][1]]))))    # Effective linear mobility calculation
	print "\tmobility at 70V =", mobEFF[where(vg == 70)[0]]
	print "\tvonset =", vonset[headd[4][1]]
	print "\tvt =", vt[headd[4][1]]

	vg_vt = vg - vt[headd[4][1]]             # Vg minus Vt array
	
	mobEFFlist_VG = [headd[4][1]]            # Adds device name to first item of list
	mobEFFlist_VGVT = [headd[4][1]]

	if sweep == "fwd":
		vg_vt_interp = vg_vt[:s/2]           # Forward scan data
		mobEFF_interp = mobEFF[:s/2]
		for v in vglist:
			mobEFFlist_VG.append(mobEFF[where(vg == v)[0][0]])
	elif sweep == "rev":
		vg_vt_interp = vg_vt[s/2:][::-1]     # Reverse scan data (and reversed)
		mobEFF_interp = mobEFF[s/2:][::-1]
		for v in vglist:
			mobEFFlist_VG.append(mobEFF[where(vg == v)[0][1]])
	else:
		sys.exit("ERROR: fwd or rev not chosen")

	vg_vt_interp = delete(vg_vt_interp, where(mobEFF_interp == inf)[0])        # Removes inf values for interpolation
	mobEFF_interp = delete(mobEFF_interp, where(mobEFF_interp == inf)[0])

	mobEFFinterpfunc = interp1d(vg_vt_interp, mobEFF_interp, kind="cubic", bounds_error=False)         # Interpolates mobEFF wrt vg_vt
	for m in mobEFFinterpfunc(vgvtlist):                                                               # Calculates interpolated mobilities
		mobEFFlist_VGVT.append(m)
	print "\tmobility at Vg-Vt=5V =", mobEFFlist_VGVT[1]

	data_summary_VG.append(mobEFFlist_VG)
	data_summary_VGVT.append(mobEFFlist_VGVT)

	# Plots transfer and mobility curves against VG
	quickPlot("mobilityplot_%s"%headd[4][1], data_path, [vg, mobEFF], xlabel="Vg (V)", ylabel="ln(muEFF (cm2/Vs))", yrange=[-12,1], yscale="linear")
	quickPlot("transfer_%s"%headd[4][1], data_path, [vg, idl, ids], xlabel="Vg (V)", ylabel="Id (A)", yrange=[1e-12,1e-2], yscale="log")

# Writes out effective mobilities to mobEFFlist files
outfile = open(os.path.join(data_path, "mobEFFlist_VG.txt"), "w")
outfile.write("Device\tEffective mobilities\n")
for a in data_summary_VG:
	l = ""
	for i in a:
		l += "%s\t"%i
	outfile.write(l[:-1] + "\n")
outfile.close()

outfile = open(os.path.join(data_path, "mobEFFlist_VGVT.txt"), "w")
outfile.write("Device\tEffective mobilities\n")
for a in data_summary_VGVT:
	l = ""
	for i in a:
		l += "%s\t"%i
	outfile.write(l[:-1] + "\n")
outfile.close()
