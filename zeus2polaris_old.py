#!/usr/bin/env python3

import os
import sys
import time
import struct
import argparse
import numpy as np
from tqdm import tqdm


class physical_quantity:
	"  "

	n_quantities = 0
	n_gas_species = 0

	descriptions = {
		0:	'gas number density',
		1:	'dust number density',
		2:	'dust temperature',
		3:	'gas temperature',
		4:	'mag. field x-direction',
		5:	'mag. field y-direction',
		6:	'mag. field z-direction',
		7:	'vel. field x-direction',
		8:	'vel. field y-direction',
		9:	'vel. field z-direction',
		10:	'rad. pressure x-direction',
		11:	'rad. pressure y-direction',
		12:	'rad. pressure z-direction',
		13:	'dust grain alignment radius',
		14:	'min. dust grain radius',
		15:	'max. dust grain radius',
		16:	'exponent of the grain size distribution',
		17:	'number fraction of a gas species',
		18:	'turbulent gas velocity',
		19:	'PDA photon count',
		20:	'ID for the opiate database',
		21:	'ID for dust composition',
		22:	'number density of thermal electrons',
		23:	'Tgaserature of thermal electrons',
		24:	'number density of cosmic ray electrons',
		25:	'min. Lorentz factor',
		26:	'max. Lorentz factor',
		27:	'power-law index of cosmic ray electrons',
		28:	'gas mass density',
		29:	'dust mass density',
		30:	'rad. field density in x-direction',
		31:	'rad. field density in y-direction',
		32:	'rad. field density in z-direction',
		33:	'total rad. field density',
		34:	'average angle between rad. and mag. field',
		35:	'anisotropy factor of the rad. field'
	}	

	def __init__(self, ID, enabled=True):
		physical_quantity.n_quantities += 1
		self.ID = ID
		self.enabled = enabled
		self.description = descriptions[ID]

def print_(string, verbose):
	if verbose:
		print(string)

def read_data(filename, shape=None, unshape=True, trim_ghosts=None, entries=-1):

	func_name = "["+sys._getframe().f_code.co_name+"] "

	if not os.path.isfile(filename):
		raise IOError(func_name+filename+" does not exists or cannot be read.")

	with open(filename, 'rb') as binfile:
		data = np.fromfile(file=binfile, dtype=np.double, count=entries)	

		# Trim ghost cells if required
		if trim_ghosts is not None:
			trim_ghosts = int(trim_ghosts)
			 
			if shape is not None:
				data = data.reshape(shape, order='F')

				data = data[trim_ghosts:-trim_ghosts,\
							trim_ghosts:-trim_ghosts,\
							trim_ghosts:-trim_ghosts]

			else:
				data = data[trim_ghosts:-trim_ghosts]

		# Turn back to 1D
		if unshape:
			if data.ndim > 1:
				data = np.swapaxes(data, 1, 2)

			data = np.ravel(data)

	return data



def write2grid(file_, format_, dtype, var, endl=False):
	if format_ == 'binary':
		file_.write(struct.pack(dtype, var))

	elif format_ == 'ascii':
		if not endl:
			file_.write(str(var)+' ')
		else:
			file_.write(str(var)+'\n')


def zeus2polaris(nframe, output_grid, verbose=True, grid_format='binary'):

	func_name = "["+sys._getframe().f_code.co_name+"] "

	start_time = time.time()

	# Read grid cell edges
	nghost	= 3
	r_e 	= read_data("z_x1ap", trim_ghosts=nghost)
	th_e 	= read_data("z_x2ap", trim_ghosts=nghost)
	ph_e 	= read_data("z_x3ap", trim_ghosts=nghost)
	r_c 	= read_data("z_x1bp", trim_ghosts=nghost)
	th_c 	= read_data("z_x2bp", trim_ghosts=nghost)
	ph_c 	= read_data("z_x3bp", trim_ghosts=nghost)

	# Grid Dimension
	dr 		= np.diff(r_e)
	dth		= np.diff(th_e)
	dph		= np.diff(ph_e)
	n_r 	= r_e.size
	n_th 	= th_e.size
	n_ph 	= ph_e.size
	n_cells	= (n_r * n_th * n_ph)
	shape 	= (n_r+2*nghost, n_ph+2*nghost, n_th+2*nghost)
	mu 		= 2.36
	mu_CO	= 28.01
	m_H		= 1.67355751e-24 # [g]

	# Read physical quantities
	rho		= read_data("o_d__"+str(nframe), trim_ghosts=nghost, shape=shape)
	Vr 		= read_data("o_v1_"+str(nframe), trim_ghosts=nghost, shape=shape)
	Vph 	= read_data("o_v2_"+str(nframe), trim_ghosts=nghost, shape=shape)
	Vth 	= read_data("o_v3_"+str(nframe), trim_ghosts=nghost, shape=shape)
	Br 		= read_data("o_b1_"+str(nframe), trim_ghosts=nghost, shape=shape)
	Bph 	= read_data("o_b2_"+str(nframe), trim_ghosts=nghost, shape=shape)
	Bth 	= read_data("o_b3_"+str(nframe), trim_ghosts=nghost, shape=shape)
	# CO 	= read_data("o_CO_"+str(nframe), trim_ghosts=nghost, shape=shape)
	# HCO		= read_data("o_xHCO_"+str(nframe), trim_ghosts=nghost, shape=shape)

	# Replace gas density by a constant density (density in cgs)
	#rho = np.full(rho.shape, 3.94959572359999973744e-16)

	# Set null velocities to a low value to avoid problems within polaris
	minvel =  1e-5
	Vr = np.where(Vr != 0.0, Vr, minvel)
	Vph = np.where(Vph != 0.0, Vph, minvel)
	Vth = np.where(Vth != 0.0, Vth, minvel)

	# Convert from Lorentz-Heaviside to Gaussian (CGS) unit system
	Br 		= Br*np.sqrt(4*np.pi)
	Bth		= Bth*np.sqrt(4*np.pi)
	Bph		= Bph*np.sqrt(4*np.pi)

	# POLARIS grid variables
	id_grid	= 30
	IDs 	= [28, 29]
	f_r 	= 0
	f_ph 	= 0
	f_th 	= 0

	# Define output format
	if grid_format == 'binary':
		output_grid += ".dat"
		file = open(output_grid, "wb")

	elif grid_format == 'ascii':
		output_grid += ".txt"
		file = open(output_grid, "w")

	else:
		raise ValueError(func_name+'grid_format must be either "binary" or "ascii".')				

	print_(f"{func_name} Writing to file: {output_grid}", verbose=verbose)

	# Write the header
	# Note:	"H": format string for unsigned short
	# 		"d": format string for double
	write2grid(file, grid_format, "H", id_grid, endl=True)
	write2grid(file, grid_format, "H", len(IDs), endl=True)

	# Write the IDs of phys. quantities
	[write2grid(file, grid_format, "H", q , endl=True) for q in IDs]

	# Write grid parameters
	# NOTE: the original grid is uniform in phi but not in r and theta
	write2grid(file, grid_format, "d", r_e[0])
	write2grid(file, grid_format, "d", read_data('z_x1ap')[-nghost] , endl=True)
	write2grid(file, grid_format, "H", n_r)
	write2grid(file, grid_format, "H", n_ph)
	write2grid(file, grid_format, "H", n_th, endl=True)

	# Write the shape parameters f_r, f_theta and f_phi 
	write2grid(file, grid_format, "d", f_r, endl=True)   
	write2grid(file, grid_format, "d", f_ph, endl=True)  
	write2grid(file, grid_format, "d", f_th, endl=True)  

	# Write the cell borders in the r, phi and theta direction
	for r in r_e[1:]:
		write2grid(file, grid_format, "d", r)	
	if grid_format == 'ascii':
		write2grid(file, 'ascii', "d", '', endl=True)

	for phi in ph_e[1:]:
		write2grid(file, grid_format, "d", phi)
	if grid_format == 'ascii':
		write2grid(file, 'ascii', "d", '', endl=True)

	for theta in th_e[1:]:
		write2grid(file, grid_format, "d", theta)	
	if grid_format == 'ascii':
		write2grid(file, 'ascii', "d", '', endl=True)

	# Display a progress bar if verbose is enabled
	for cell in tqdm(range(n_cells)) if verbose else range(n_cells):
		# Write physical data into cells
		write2grid(file, grid_format, "d", rho	[cell])
		write2grid(file, grid_format, "d", 0.01*rho	[cell], endl=True)

	# Set the central cell to zero in all quantities.
	write2grid(file, grid_format, "d", 0)
	write2grid(file, grid_format, "d", 0, endl=True)
	
	file.close()	

	
	runtime = time.time()-start_time
	print_(f'{func_name} Elapsed time: {time.strftime("%H:%M:%S",time.gmtime(runtime))}', verbose=verbose)


if __name__ == '__main__':

	# to do: implement argparse

	if len(sys.argv) > 1:
		nframe	= sys.argv[1].zfill(5)
	else:
		nframe = input(f"Please enter the time frame number: " ).zfill(5)
	
	output_grid = f'polaris_grid/lmd2.4-1k-Slw_{nframe}'


	zeus2polaris(nframe=nframe, output_grid=output_grid, verbose=True, grid_format='ascii')
