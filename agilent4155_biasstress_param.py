#! /usr/bin/env python
"""
agilent4155_biasstress_param.py
Parameter extractor for EasyExpert bias stress files

Created by Jeremy Smith on 2015-08-13
University of California, Berkeley
j-smith@ecs.berkeley.edu
"""

import os
import sys
import re
import numpy as np
import myfunctions as mf
from scipy import stats

__author__ = "Jeremy Smith"
__version__ = "1.2"

data_path = os.path.dirname(__file__)    # Path name for location of script

files = os.listdir(data_path)   # All files in directory
data_summary = []
skipinit = 5

chl = 80
chw = 350
tox = 60
kox = 51
ci = 8.85418782e-7*kox/tox
plot = True
doublesweep = True


def main():
	"""Main function"""

	print "\nBatch importing .csv files..."
	print data_path, '\n'

	for f in files:
		print f
		# Loops through all transfer files
		if "Bias Stress" not in f:
			continue
		datatransfer, datastress = mf.csvBiasStressImport(f, data_path)
		no_scans = len(datatransfer)
		no_stress = len(datastress)
		stress_time = datastress[0]['Time'][-1]

		for s in xrange(no_scans):
			print "    Scan number:", s+1
			outname = "{:s}_{:03d}".format(re.sub("[^\[]*\[|_;|]", "", f[:-4].replace(' ', '_')), s+1)
			vd = np.array(datatransfer[s]['Vdrain'][0])
			vg = np.array(datatransfer[s]['Vgate'])
			ids = np.array(datatransfer[s]['absIdrain'])
			igs = np.array(datatransfer[s]['absIgate'])

			if doublesweep:
				ids_forward = ids[:len(ids)/2]
				vg_forward = vg[:len(vg)/2]
			else:
				ids_forward = ids
				vg_forward = vg

			ids_smoothed = mf.adjAvSmooth(abs(ids_forward), N=1)
			diff_ids_smoothed = np.array(mf.numDiff(ids_smoothed, vg_forward))
			tlmaxarg = np.argmax(diff_ids_smoothed[skipinit:-1]) + skipinit
			# Linear mobility (max transconductance)
			linmob_tmax = (chl/(chw*ci*vd))*(diff_ids_smoothed[tlmaxarg])
			# Linear threshold voltage (max transconductance)
			vthlin_tmax = vg[tlmaxarg] - ids_smoothed[tlmaxarg]/diff_ids_smoothed[tlmaxarg]
			# Linear mobility at vth + 1 V
			linmob_constn = (chl/(chw*ci*vd))*np.interp(vthlin_tmax + 1.0, vg_forward, diff_ids_smoothed)
			# Subthreshold slope
			lin_sts = min(abs(1/np.array(mf.numDiff([np.log10(abs(x)) for x in ids_smoothed[skipinit:-1]], vg_forward[skipinit:-1]))))

			data_summary.append([outname, linmob_tmax, vthlin_tmax, linmob_constn, lin_sts])

			mf.dataOutput(outname + "_transfer.txt",  data_path, [vg, ids, igs], format="%.2f\t %.5e\t %.5e\n")
			
			# Plot transfer
			if plot:
				mf.quickPlot(outname + "_TRANSFERplot", data_path, [vg, ids, igs], xlabel="Vg [V]", ylabel="Id,g [A]", yscale="log", yrange=[1e-12, 1e-3])
				mf.quickPlot(outname + "_TCplot", data_path, [vg_forward, ids_smoothed], xlabel="Vg [V]", ylabel="Id [A]")

			if s < no_stress:
				time = np.array(datastress[s]['Time'])
				idsst = np.array(datastress[s]['absIdStress'])
				igsst = np.array(datastress[s]['absIgStress'])
				mf.dataOutput(outname + "_stress.txt",  data_path, [time, idsst, igsst], format="%.2f\t %.5e\t %.5e\n")

	mf.dataOutput("BIAS_STRESS_SUMMARY.txt", data_path, map(list, zip(*data_summary)), format="%s\t %.5e\t %.5f\t %.5e\t %.5e\n")

	return


if __name__ == "__main__":
	sys.exit(main())
