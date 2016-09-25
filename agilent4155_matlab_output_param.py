#! /usr/bin/env python
"""
agilent4155_matlab_output_param.py
Prarmeter extractor for matlab generated .xlsx idvd files

Created by Jeremy Smith on 2015-10-29
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
__version__ = "1.0"

data_path = os.path.dirname(__file__)    # Path name for location of script

files = os.listdir(data_path)   # All files in directory
data_summary = []
summary_list_header = [["filename", "channelL", "channelW"]]
fitrange = 3


def main():
	"""Main function"""

	print "\nBatch importing .xlsx files..."
	print data_path, '\n'

	for f in files:
		print f
		# Loops through all transfer files
		if "IDVD.xlsx" in f:

			workbook = xlrd.open_workbook(f, logfile=open(os.devnull, 'w'))

			for dev in workbook.sheet_names():
				if "Sheet" in dev:
					continue
				print "  - device {:s}".format(dev)
				datasheet = workbook.sheet_by_name(dev)
				run_numbers = [str(int(x)) for x in datasheet.row_values(2) if x]
				stepvg = 0

				for i, run in enumerate(run_numbers):
					print "    - run {:s}".format(run)
					data = {}
					gdlin = []
					gdsat = []
					vg_list = []
					# File name for outputs
					outname = f[:-5] + '_' + dev + '_' + run
					# Constant parameters taken from header
					vgmin = float(datasheet.cell_value(3, (stepvg + 2)*i + 1))
					vgmax = float(datasheet.cell_value(4, (stepvg + 2)*i + 1))
					stepvg_prev = stepvg
					stepvg = int(datasheet.cell_value(5, (stepvg + 2)*i + 1))
					chl = float(datasheet.cell_value(1, 1))
					chw = float(datasheet.cell_value(0, 1))
					tox = float(datasheet.cell_value(1, 3))
					kox = float(datasheet.cell_value(0, 3))
					ldr = float(datasheet.cell_value(1, 5))
					lso = float(datasheet.cell_value(0, 5))
					ci = 8.85418782e-7*kox/tox

					colheads = ['VDS'] + ["ID{:d}".format(x + 1) for x in range(stepvg)]

					for h in colheads:
						data[h] = []
					for row in range(datasheet.nrows - 11 - stepvg):
						for col, h in enumerate(colheads):
							if datasheet.cell_type(9 + row, (stepvg_prev + 2)*i + col) is 0:
								continue
							data[h].append(float(datasheet.cell_value(9 + row, (stepvg_prev + 2)*i + col)))

					vds = np.array(data['VDS'])
					output_list = [vds]
					
					for j in range(stepvg):
						ids = np.array(data["ID{:d}".format(j + 1)])
						# Fits to first data points (given by fitrange) i.e. linear
						slope, intercept, r_value, p_value, std_err = stats.linregress(vds[:fitrange], ids[:fitrange])
						gdlin.append(slope)
						# Fits to last data points (given by fitrange) i.e. saturation
						slope, intercept, r_value, p_value, std_err = stats.linregress(vds[-fitrange:], ids[-fitrange:])
						gdsat.append(slope)
						# Update lists
						output_list.append(ids)
						vg_list.append(vgmin + j*(vgmax-vgmin)/(stepvg-1))

					# Output data
					data_summary.append([outname, chl, chw])

					# Ouput files
					mf.dataOutputGen(outname+"_output.txt", data_path, map(list, zip(*output_list)))
					mf.dataOutputHead(outname+"_gm.txt", data_path, [vg_list, gdlin, gdsat], [['VG', 'GDlin', 'GDsat']], 
						format_d="%.2f\t %.5e\t %.5e\n",
						format_h="%s\t")


	mf.dataOutputHead("SUMMARY_OUT.txt", data_path, map(list, zip(*data_summary)), summary_list_header,
		format_d="%s\t %.1f\t %.1f\n", 
		format_h="%s\t")

	return


if __name__ == "__main__":
	sys.exit(main())
