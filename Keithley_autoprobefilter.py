# Reads Keithley .xls and filters data from autoprober
# Jeremy Smith
# Northwestern University
# Version 1.1


from numpy import *
import xlrd
import os
import sys
from myfunctions import *
from scipy import stats


data_path = os.path.dirname(__file__)    # Path name for location of script

print "\n"
print data_path
print "\n"

files = os.listdir(data_path)   # All files in directory
data_summary = []

skipinit = 5   # Initial data points to skip
l = 50.
w = 150.
ci = 497e-9

# Loops through each device file
for d in files:
	print d
	if "Vgs-Id" not in d:
		continue

	workbook = xlrd.open_workbook(d, logfile=open(os.devnull, 'w'))
	datasheet = workbook.sheet_by_index(0)
	data = {}
	colhead = []
	for h in range(datasheet.ncols):
		data[datasheet.cell_value(0,h)] = []
		colhead.append(datasheet.cell_value(0,h))
	for r in range(datasheet.nrows-1):
		for c in range(datasheet.ncols):
			data[colhead[c]].append(float(datasheet.cell_value(r+1,c)))
	settingsheet = workbook.sheet_by_index(2)

	idrain = array(data["DrainI"])
	igate = array(data["GateI"])
	vgate = array(data["GateV"])
	vdrain = float(settingsheet.cell_value(19,2))

	onoff_at_zero = abs(idrain[-1]/idrain[where(vgate==0)[0][0]])
	leakage_ratio = abs(idrain/igate)

	# Initial filtering steps
	filtered = False
	if onoff_at_zero < 1e2:
		filtered = True
	if leakage_ratio[where(vgate==2)[0][0]] < 10:
		filtered = True
	if abs(idrain[-1]) < 1e-10:
		filtered = True
	if filtered == True:
		print "  FILTERED"
		continue

	# Smoothing the drain current data
	idrain_smoothed = adjAvSmooth(abs(idrain), N=1)

	# Finding max transconductance
	sqrtidrain = sqrt(idrain_smoothed)
	diff_sqrt_idrain_smoothed = array(numDiff(sqrtidrain, vgate))
	tmaxarg = argmax(diff_sqrt_idrain_smoothed[skipinit:-1]) + skipinit
	# Saturation mobility (max transconductance)
	satmob_tmax = (2*l/(w*ci))*(diff_sqrt_idrain_smoothed[tmaxarg])**2
	# Threshold Voltage (max transconductance)
	vth_tmax = vgate[tmaxarg] - sqrt(idrain_smoothed)[tmaxarg]/diff_sqrt_idrain_smoothed[tmaxarg]
	# On-off ratio
	onoffratio = log10(max(idrain_smoothed[skipinit:-1])/min(idrain_smoothed[skipinit:-1]))

	# Finds range of data that lies within the minimum+15% and the maximum-15% and also has a positive transconductance
	fitrange_id = [0.85*min(sqrtidrain[skipinit:-1]) + 0.15*max(sqrtidrain[skipinit:-1]), 0.85*max(sqrtidrain[skipinit:-1]) + 0.15*min(sqrtidrain[skipinit:-1])]
	fitrange_bool = bitwise_and(sqrtidrain > fitrange_id[0], sqrtidrain < fitrange_id[1], diff_sqrt_idrain_smoothed > 0)

	# Checks that there are at least 3 data points to fit
	if sum(fitrange_bool) < 3:
		filtered = True
		print "  FILTERED"
		continue
	
	# Linear Fitting to sqrt(Idrain)
	slope, intercept, r_value, p_value, std_err = stats.linregress(vgate[fitrange_bool], sqrtidrain[fitrange_bool])
	fitline = slope*vgate + intercept

	# Saturation mobility (from slope of sqrt(Idrain) fit)
	satmob_FITTED = (2*l/(w*ci))*slope**2
	# Threshold Voltage (from slope of sqrt(Idrain) fit)
	vth_FITTED = -intercept/slope

	# Second filtering steps
	if abs(vth_FITTED) > 3.0:
		filtered = True
	if satmob_FITTED < 0.1:
		filtered = True
	if satmob_FITTED > 250.:
		filtered = True
	if r_value**2 < 0.9:
		filtered = True
	if filtered == True:
		print "  FILTERED"
		continue

	data_summary.append([d[:-4], satmob_tmax, vth_tmax, log10(onoff_at_zero), onoffratio, log10(leakage_ratio[where(vgate==2)[0][0]]), satmob_FITTED, vth_FITTED, r_value**2])
	quickPlot(d[:-4]+"_SQRTplot", data_path, [vgate, sqrtidrain, fitline], xlabel="VG [V]", ylabel="sqrt(Id) [A^0.5]", yrange=[0, 'auto'])
	quickPlot(d[:-4]+"_TRANSFERplot", data_path, [vgate, idrain_smoothed, abs(igate)], xlabel="VG [V]", ylabel="Id,g [A]", yscale="log", yrange=[1e-12, 1e-3])

outfile = open(os.path.join(data_path, "summarylist.txt"), "w")
outfile.write("Device\tSatMobility(maxtrans)\tVthreshold(maxtrans)\tLogOnOffAtZero\tLogOnOffRatio\tLogLeakageRatioAt2V\tSatMobility(fitted)\tVthreshold(fitted)\tFITr2\n")
for a in data_summary:
	outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8]))
outfile.close()
