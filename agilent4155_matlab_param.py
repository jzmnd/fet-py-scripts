#! /usr/bin/env python
"""
agilent4155_matlab_param.py
Prarmeter extractor for matlab generated .xlsx files

Created by Jeremy Smith on 2015-07-07
University of California, Berkeley
j-smith@eecs.berkeley.edu
"""

import os
import sys
import xlrd
import numpy as np
import myfunctions as mf
from scipy import stats

__author__ = "Jeremy Smith"
__version__ = "1.5"

data_path = os.path.dirname(__file__)    # Path name for location of script

files = os.listdir(data_path)   # All files in directory
data_summary = []
colheads = ['VG', 'ID1', 'ID2', 'IG1', 'IG2']
skipinit = 3          # skip initial data points 
rangefittop = 20.0    # % from max for fit range
rangefitbot = 20.0    # % from min for fit range
mobilitycorrection = 1.000     # K_actual/K_file
summary_list_header = [["filename", "satmob_tmax", "vthsat_tmax", "satmob_rev_tmax", "vthsat_rev_tmax",
						"linmob_tmax", "vthlin_tmax", "linmob_rev_tmax", "vthlin_rev_tmax",
						"hysteresis_lin", "hysteresis_sat",
						"onoffratio_lin", "onoffratio_sat",
						"leakage_ratio_lin", "leakage_ratio_sat",
						"Slin", "Slin_rev", "Ssat", "Ssat_rev",
						"satmob_FITTED", "vthsat_FITTED", "r_value"]]
sweepfwddirection = True       # True if sweeping from depletion to accumulation

def main():
	"""Main function"""

	print "\nBatch importing .xlsx files..."
	print data_path, '\n'

	for f in files:
		print f
		# Loops through all transfer files
		if "IDVG.xlsx" in f:

			workbook = xlrd.open_workbook(f, logfile=open(os.devnull, 'w'))

			for dev in workbook.sheet_names():
				if "Sheet" in dev:
					continue
				print "  - device {:s}".format(dev)
				datasheet = workbook.sheet_by_name(dev)
				run_numbers = [str(int(x)) for x in datasheet.row_values(2) if x]

				for i, run in enumerate(run_numbers):
					print "    - run {:s}".format(run)
					data = {}
					# File name for outputs
					outname = f[:-5] + '_' + dev + '_' + run
					# Constant parameters taken from header of .xlsx file
					vd1 = float(datasheet.cell_value(3, 6*i + 3))
					vd2 = float(datasheet.cell_value(4, 6*i + 3))
					chl = float(datasheet.cell_value(1, 1))
					chw = float(datasheet.cell_value(0, 1))
					tox = float(datasheet.cell_value(1, 3))
					kox = float(datasheet.cell_value(0, 3))
					ldr = float(datasheet.cell_value(1, 5))
					lso = float(datasheet.cell_value(0, 5))
					# Calculation of geometric capacitance
					ci = 8.85418782e-7*kox/tox
					# Extract data
					for h in colheads:
						data[h] = []
					for row in range(datasheet.nrows - 9):
						for col, h in enumerate(colheads):
							if datasheet.cell_type(9 + row, 6*i + col) is 0:
								continue
							data[h].append(float(datasheet.cell_value(9 + row, 6*i + col)))
					p = len(data['VG'])/2
					# Forward scan
					if sweepfwddirection:
						vg = np.array(data['VG'][:p])
						id1 = np.array(data['ID1'][:p])
						id2 = np.array(data['ID2'][:p])
						ig1 = np.array(data['IG1'][:p])
						ig2 = np.array(data['IG2'][:p])
					else:
						vg = np.array(data['VG'][:p][::-1])
						id1 = np.array(data['ID1'][:p][::-1])
						id2 = np.array(data['ID2'][:p][::-1])
						ig1 = np.array(data['IG1'][:p][::-1])
						ig2 = np.array(data['IG2'][:p][::-1])
					# Reverse scan
					if sweepfwddirection:
						vg_r = np.array(data['VG'][p:][::-1])
						id1_r = np.array(data['ID1'][p:][::-1])
						id2_r = np.array(data['ID2'][p:][::-1])
						ig1_r = np.array(data['IG1'][p:][::-1])
						ig2_r = np.array(data['IG2'][p:][::-1])
					else:
						vg_r = np.array(data['VG'][p:])
						id1_r = np.array(data['ID1'][p:])
						id2_r = np.array(data['ID2'][p:])
						ig1_r = np.array(data['IG1'][p:])
						ig2_r = np.array(data['IG2'][p:])

					# Smoothing Id for fitting
					id1_smoothed = mf.adjAvSmooth(abs(id1), N=1)
					id2_smoothed = mf.adjAvSmooth(abs(id2), N=1)
					id1_r_smoothed = mf.adjAvSmooth(abs(id1_r), N=1)
					id2_r_smoothed = mf.adjAvSmooth(abs(id2_r), N=1)
					# On-off ratio
					onoffratio1 = np.log10(max(id1[skipinit:-1])/min(abs(id1[skipinit:-1])))
					onoffratio2 = np.log10(max(id2[skipinit:-1])/min(abs(id2[skipinit:-1])))
					# Leakage ratio
					leakage_ratio1 = np.log10(abs(id1/ig1))
					leakage_ratio2 = np.log10(abs(id2/ig2))

					# Finding max saturation transconductance
					sqrtid2 = np.sqrt(id2_smoothed)
					sqrtid2_r = np.sqrt(id2_r_smoothed)
					diff_sqrt_id2_smoothed = np.array(mf.numDiff(sqrtid2, vg))
					diff_sqrt_id2_r_smoothed = np.array(mf.numDiff(sqrtid2_r, vg_r))
					tsmaxarg = np.argmax(diff_sqrt_id2_smoothed[skipinit:-1]) + skipinit
					tsmaxarg_r = np.argmax(diff_sqrt_id2_r_smoothed[skipinit:-1]) + skipinit
					# Saturation mobility (max transconductance)
					satmob_t = mobilitycorrection * (2*chl/(chw*ci))*(diff_sqrt_id2_smoothed)**2
					satmob_r_t = mobilitycorrection * (2*chl/(chw*ci))*(diff_sqrt_id2_r_smoothed)**2
					satmob_tmax = satmob_t[tsmaxarg]
					satmob_r_tmax = satmob_r_t[tsmaxarg_r]
					# Saturation threshold voltage (max transconductance)
					vthsat_tmax = vg[tsmaxarg] - sqrtid2[tsmaxarg]/diff_sqrt_id2_smoothed[tsmaxarg]
					vthsat_r_tmax = vg_r[tsmaxarg_r] - sqrtid2_r[tsmaxarg_r]/diff_sqrt_id2_r_smoothed[tsmaxarg_r]
					# Hysteresis
					hysteresissat = vthsat_tmax - vthsat_r_tmax
					# Calculate subthreshold slopes
					sts_sat = min(abs(1/np.array(mf.numDiff([np.log10(abs(x)) for x in id2_smoothed[skipinit:-1]], vg[skipinit:-1]))))
					sts_r_sat = min(abs(1/np.array(mf.numDiff([np.log10(abs(x)) for x in id2_r_smoothed[skipinit:-1]], vg_r[skipinit:-1]))))

					# Finding max linear transconductance
					diff_id1_smoothed = np.array(mf.numDiff(id1_smoothed, vg))
					diff_id1_r_smoothed = np.array(mf.numDiff(id1_r_smoothed, vg_r))
					tlmaxarg = np.argmax(diff_id1_smoothed[skipinit:-1]) + skipinit
					tlmaxarg_r = np.argmax(diff_id1_r_smoothed[skipinit:-1]) + skipinit
					# Linear mobility (max transconductance)
					linmob_t = mobilitycorrection * (chl/(chw*ci*vd1))*(diff_id1_smoothed)
					linmob_r_t = mobilitycorrection * (chl/(chw*ci*vd1))*(diff_id1_r_smoothed)
					linmob_tmax = linmob_t[tlmaxarg]
					linmob_r_tmax = linmob_r_t[tlmaxarg_r]
					# Linear threshold voltage (max transconductance)
					vthlin_tmax = vg[tlmaxarg] - id1_smoothed[tlmaxarg]/diff_id1_smoothed[tlmaxarg]
					vthlin_r_tmax = vg_r[tlmaxarg_r] - id1_r_smoothed[tlmaxarg_r]/diff_id1_r_smoothed[tlmaxarg_r]
					# Hysteresis
					hysteresislin = vthlin_tmax - vthlin_r_tmax
					# Calculate subthreshold slopes
					sts_lin = min(abs(1/np.array(mf.numDiff([np.log10(abs(x)) for x in id1_smoothed[skipinit:-1]], vg[skipinit:-1]))))
					sts_r_lin = min(abs(1/np.array(mf.numDiff([np.log10(abs(x)) for x in id1_r_smoothed[skipinit:-1]], vg_r[skipinit:-1]))))

					# Finds range of data that lies within the minimum+x% and the maximum-x% and also has a positive transconductance
					fitrange_id_lo = (1-rangefitbot/100.0)*min(sqrtid2[skipinit:-1]) + (rangefitbot/100.0)*max(sqrtid2[skipinit:-1])
					fitrange_id_hi = (1-rangefittop/100.0)*max(sqrtid2[skipinit:-1]) + (rangefittop/100.0)*min(sqrtid2[skipinit:-1])
					fitrange_bool = np.bitwise_and(np.bitwise_and(sqrtid2 > fitrange_id_lo, sqrtid2 < fitrange_id_hi), diff_sqrt_id2_smoothed > 0)
					# Checks that there are at least 3 data points to fit
					if sum(fitrange_bool) < 3:
						print "      NOT ENOUGH DATA TO FIT"
						satmob_FITTED = np.nan
						vthsat_FITTED = np.nan
						r_value = np.nan
					else:
						# Linear Fitting to sqrt(Idrain)
						slope, intercept, r_value, p_value, std_err = stats.linregress(vg[fitrange_bool][skipinit:-1], sqrtid2[fitrange_bool][skipinit:-1])
						fitline = slope*vg + intercept
						# Saturation mobility (from slope of sqrt(Idrain) fit)
						satmob_FITTED = mobilitycorrection * (2*chl/(chw*ci))*slope**2
						# Threshold Voltage (from slope of sqrt(Idrain) fit)
						vthsat_FITTED = -intercept/slope
						# Plot sqrt(Isd)
						mf.quickPlot(outname+"_SQRTplot", data_path, [vg, sqrtid2, fitline],
						xlabel="VG [V]", ylabel="sqrt(Id) [A^0.5]", yrange=[0, 'auto'])

					# Output data
					data_summary.append([outname,
						satmob_tmax, vthsat_tmax,
						satmob_r_tmax, vthsat_r_tmax,
						linmob_tmax, vthlin_tmax,
						linmob_r_tmax, vthlin_r_tmax,
						hysteresislin, hysteresissat,
						onoffratio1, onoffratio2,
						leakage_ratio1[-skipinit],
						leakage_ratio2[-skipinit],
						sts_lin, sts_r_lin, sts_sat, sts_r_sat,
						satmob_FITTED, vthsat_FITTED, r_value**2])

					# Ouput files
					mf.dataOutputHead(outname+"_transfer.txt", data_path, [np.array(data['VG']), abs(np.array(data['ID1'])), abs(np.array(data['ID2'])), abs(np.array(data['IG1'])), abs(np.array(data['IG2'])),
						np.concatenate((linmob_t, linmob_r_t[::-1])), np.concatenate((satmob_t, satmob_r_t[::-1]))],
						[["vg", "idlin", "idsat", "iglin", "igsat", "LINMOB", "SATMOB"]], 
						format_d="%.3f\t %.5e\t %.5e\t %.5e\t %.5e\t %.5e\t %.5e\n", 
						format_h="%s\t")
					
					# Plot transfer
					mf.quickPlot(outname+"_TRANSFERplot", data_path, [vg, id1_smoothed, id1_r_smoothed, abs(ig1), id2_smoothed, id2_r_smoothed, abs(ig2)],
						xlabel="VG [V]", ylabel="Id,g [A]", yscale="log", yrange=[1e-12, 1e-2], col=["r", "r", "r", "b", "b", "b"])

	mf.dataOutputHead("SUMMARY.txt", data_path, map(list, zip(*data_summary)), summary_list_header,
		format_d="%s\t %.5e\t %.5f\t %.5e\t %.5f\t %.5e\t %.5f\t %.5e\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5f\t %.5e\t %.5f\t %.6f\n", 
		format_h="%s\t")

	return


if __name__ == "__main__":
	sys.exit(main())
