#+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!
#                                                                       #
#                                   slf2nc.py                           #
#                                                                       #
#+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!

#
# Author : Caio Eadi Stringari
# Date   : 08 July 2015
#
# Will convert the unstructured selafin file to a regular netcdf file.
#
# Only works for 2d files at the moment. Use POSTEL to extract 2d planes
# if you want to play with 3d files.
#
# A BlueKenue regular grid (.r2s) file must be passed as the third argument.
#
# You need to install at least numpy, matplotlib, netCDF4 and Progressbar.
#
# Dont forget to define your own pytel path in line 35.
#
# Call the function as "python sel2nc.py Input.slf Output.nc Grid.r2s"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Imports
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import sys, os
from os import path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

# pytel
pytel = 'pytel/'
sys.path.append(path.join(path.dirname(sys.argv[0]),pytel))
from parsers.parserSELAFIN import SELAFIN

from progressbar import *

# netcdf4
import netCDF4

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if len(sys.argv)-1 != 3:

	print('Wrong number of Arguments, stopping now...')
	print('Usage: python slf2reg.py Input (.slf) Output (.nc) Grid (.r2s)')

	sys.exit()

else:

	sflfile = sys.argv[1]
	cdffile = sys.argv[2]
	grdfile = sys.argv[3]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Telemac:

	def __init__(self):
		pass

	@classmethod
	def iocheck(self,fname):
		io = os.path.isfile(fname)
		if io:
			pass
		else:
			raise IOError('File {0} not found.'.format(fname))


	@classmethod
	def geometry(self,fname):

		Telemac.iocheck(fname)

		slf = SELAFIN(fname)

		varnames  = slf.VARNAMES
		npoin2    = slf.NPOIN2
		nplan2    = slf.NPLAN
		nelem2    = slf.NELEM3
		xcoord    = slf.MESHX
		ycoord    = slf.MESHY
		elem2     = np.dstack((slf.MESHX[slf.IKLE3],slf.MESHY[slf.IKLE3]))
		ikle2     = np.array(slf.IKLE3)

		return varnames,npoin2,nplan2,nelem2,xcoord,ycoord,elem2,ikle2

	@classmethod
	def getvalues(self,slf,t):
		var = {}
		for v,variable in enumerate(slf.VARNAMES):
			values     = slf.getVALUES(t)     # Get all values for the timestep
			var[variable] = values[v]
		return var

	@classmethod
	def trinterp(self,x,y,tri,xi,yi,z,method="linear"):

		if method == "linear":
			I = mtri.LinearTriInterpolator(tri, z)
			zi = I(xi, yi)
		elif method == "geom":
			I = mtri.CubicTriInterpolator(tri, z)
			zi = I(xi, yi)
		elif method == "min_E":
			I = mtri.CubicTriInterpolator(tri, z)
			zi = I(xi, yi)

		return zi

class BlueKenue():

	def __init__(self):
		pass

	@classmethod
	def regulargrid(self,fname):

		Telemac.iocheck(fname)

		f = open(fname).readlines()
		for line in f:
			if ":xOrigin" in line:
				xO=float(line.split()[1].strip("\r\n"))
			if ":yOrigin" in line:
				yO=float(line.split()[1].strip("\r\n"))
			if "xCount" in line:
				xC=float(line.split()[1].strip("\r\n"))
			if "yCount" in line:
				yC=float(line.split()[1].strip("\r\n"))
			if "xDelta" in line:
				xD=float(line.split()[1].strip("\r\n"))
			if "yDelta" in line:
				yD=float(line.split()[1].strip("\r\n"))
			if "Angle" in line:
				aN=float(line.split()[1].strip("\r\n"))
				if aN > 0:
					ValueError('Rotated grids not supported yet.')
			if ":EndHeader" in line: break
		return xO,yO,xC,yC,xD,yD,aN


class NetCdf:

	def __init__(self):
		pass

	@classmethod
	def writenc(self,fname,time,xm,ym,ikle,X,Y,ncnames,slfnames,units):

		# Names
		names = zip(slfnames,ncnames)

		# Dimensions
		nt = len(time)
		nx = len(xg)
		ny = len(yg)
		nv = len(ncnames)

		# ProgressBar
		widgets = [" ", Percentage(), ' ', Bar(marker=RotatingMarker()),
		       ' ', ETA(), ' ', FileTransferSpeed()]
		pbar = ProgressBar(widgets=widgets, maxval=nt*nv).start()

		ncout = netCDF4.Dataset(fname, 'w', format='NETCDF4')

		ncout.createDimension('times',  nt)    # Temporal  dimension
		ncout.createDimension('x',      nx)    # X dimension
		ncout.createDimension('y',      ny)    # Y dimension
		#
		t = ncout.createVariable('times',   'i4', ('times',))
		x = ncout.createVariable('x-coord', 'f4', ('x',))
		y = ncout.createVariable('y-coord', 'f4', ('y',))
		#
		x.units = "m (Eastings)"
		y.units = "m (Northings)"
		t.units = "s"
		#
		t[:] = time
		x[:] = xg
		y[:] = yg
		#
		# Triangulation
		trgl = mtri.Triangulation(xm, ym, ikle)
		#
		# Interpolating and Writing
		k = 0
		for i,name in enumerate(names):
			print ("   Writing {}".format(name[1]))
			# Creating netcdf variables
			V = ncout.createVariable(name[1], 'f8', ('times', 'x', 'y',))
			for t, time in enumerate(times):
				z = slf.getVALUES(t)[i]
				zi = Telemac.trinterp(xm,ym,trgl,X,Y,z)
				V[t,:,:] = zi.T
				k+=1
				pbar.update(k)
		pbar.finish()
		ncout.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main Calls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == "__main__":

	print("Script called : {}".format(sys.argv[0]))

	print (" Opening .slf file")

	slf = SELAFIN(sflfile)
	vnames,npt2,npl2,nle2,xm,ym,elem,ikle = Telemac.geometry(sflfile)

	# Modifying variable names:
	names = []
	for name in vnames:
		nname = name.strip(" ").replace(" ","_").lower()
		names.append(nname)

	# Units
	units = []
	for i,u in enumerate(slf.VARUNITS):
		un = u.strip(" ").lower()
		units.append(un)

	# Times
	times = slf.tags["times"]

	print (" Opening .r2s file")
	xO,yO,xC,yC,xD,yD,aN = BlueKenue.regulargrid(grdfile)

	xg     = np.linspace(xO,(xC*xD)+xO,xC,endpoint=True)
	yg     = np.linspace(yO,(yC*yD)+yO,yC,endpoint=True)
	X, Y   = np.meshgrid(xg,yg)

	print " Writing .nc file"

	NetCdf.writenc(cdffile,times,xm,ym,ikle,X,Y,names,vnames,units)

	print("Converted files sucessfully ")
