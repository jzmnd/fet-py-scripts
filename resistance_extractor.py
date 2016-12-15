#! /usr/bin/env python
"""
resistance_extractor.py

Created by Jeremy Smith on 2015-07-30
University of California, Berkeley
j-smith@ecs.berkeley.edu
"""

import os
import sys
import numpy as np
import myfunctions as mf
from scipy import stats

__author__ = "Jeremy Smith"
__version__ = "1.0"

data_path = os.path.dirname(__file__)    # Path name for location of script

files = os.listdir(data_path)   # All files in directory
data_summary = []


def main():
	"""Main function"""

	print data_path, '\n'

	for f in files:
		print f
		# Loops through all files
		if "I_V Sweep" not in f:
			continue
		datad, headd = mf.csvImport(f, data_path, 104)
		deviceID = headd[4][1]
		print "   Device:", deviceID
		voltage = np.array(datad["V1"])
		current = np.array(datad["I1"])
		slope, intercept, r_value, p_value, std_err = stats.linregress(voltage, current)
		resistance = 1.0/slope

		data_summary.append([deviceID, resistance])

	mf.dataOutput("SUMMARY.txt", data_path, map(list, zip(*data_summary)), format="%s\t %.5e\n")

	return


if __name__ == "__main__":
	sys.exit(main())