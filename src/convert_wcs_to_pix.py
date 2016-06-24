#! /Users/Austin/progs/Ureka/variants/common/bin/python
#----------------------------
#   NAME
#----------------------------
# convert_wcs_to_pix.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# The purpose of this program is to convert a strong lensing catalog
# given in RA,DEC in wcs degrees to image coordinates
# so that imremap does not have to do this clunkily in C.
# 
# This program will output a file that can then be read in to remap
#
#----------------------------
#   INPUTS:
#----------------------------
# input_strongcat      strong lensing catalog of the lenstool format, i.e.:
#				       ID RA (degrees) DEC (degrees) major_axis minor_axis angles z magnitude
#
# reference_fitsfile   The fits file whose wcs information will be used for the conversion.
#					   This needs to have the same wcs information as the file that will
#					   be used in imremap      
#                    
#                    
#----------------------------
#   OUTPUTS: 
#----------------------------
# output_strongcat       : File (named same as input file, but with _image prepended before the extension)
#					     This file will have the exact same format as input_strongcat, 
#						 but RA and DEC in wcs degrees will be saved in image coordinates
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# bash> chmod +x convert_wcs_to_pix (only required once)
#
# bash> convert_wcs_to_pix stronglensing_example.cat hlsp_frontier_model_macs0416_cats_v3_x-pixels-deflect.fits
#
#----------------------------
from astropy.wcs import WCS
import numpy as np
import sys


def convert_file(input_strongcat,reference_fitsfile):
	''' Save a new file that is identical to the input_filename 
	except ra,dec in degrees are replaced with pixel x and y
	in the image coordinate system of reference_fitsfile'''
	data = np.genfromtxt(input_strongcat,dtype='S20')
	W = WCS(reference_fitsfile) 
	ras_deg = map(float,data[:,1])
	decs_deg = map(float,data[:,2])
	xs,ys = W.wcs_world2pix(ras_deg,decs_deg,1)
	data[:,1] = xs
	data[:,2] = ys
	prefix,extension = input_strongcat.split('.')
	output_strongcat = prefix + '_image' + '.' + extension
	# print output_strongcat
	np.savetxt(output_strongcat,data,fmt='%s')
	# savearray = np.hstack((ids,xs,ys,))
	# w.wcs_world2pix(image_ra,image_dec,1) # position of object in x,y

if __name__ == '__main__': # Execute the stuff below when this program is called from the command line
	input_strongcat = sys.argv[1] 
	reference_fitsfile = sys.argv[2]
	convert_file(input_strongcat,reference_fitsfile)