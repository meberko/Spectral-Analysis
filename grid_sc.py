#!/usr/bin/env python

# global script for automation
# this is for sherpa in ciao 4.3

# examples
# grid.py \
#        --rmf=example.rmf \
#        --arf=example.arf \
#        --fits=output.qgrid.fits \
#        --model="xsphabs.A*xspowerlaw.B" \
#        --par1=A.nH:0.01,0.1,0.4,1,4,10:200 \
#        --par2=B.PhoIndex:0,1,2,3,4:200 \
#        --etcp=B.norm:1.0

import subprocess as sp
import sys
import re
import numpy
from time import time, localtime
#from scipy.stats import norm

from optparse import OptionParser

usage = "usage: %prog --model model ..."
parser = OptionParser(usage=usage)
parser.add_option("--debug", 	dest="debug", 	help="debug", default=0)
parser.add_option("--model", 	dest="model", 	help="model", default="xsphabs.A*xspowerlaw.B:xsphabs.A")
parser.add_option("--par1", 	dest="par1", 	help="par1",  default="A.nH:0.01,0.1,0.4,1,4,10:200")
parser.add_option("--par2", 	dest="par2", 	help="par2",  default="B.PhoIndex:0,1,2,3,4:200")
parser.add_option("--frac", 	dest="frac", 	help="frac",  default="0.25,0.5,0.75")
parser.add_option("--range", 	dest="range", 	help="range", default="0.3,8.0")
parser.add_option("--etcp", 	dest="etcp", 	help="etcp",  default="B.norm:1.0")
parser.add_option("--abund", 	dest="abund", 	help="abund", default="angr")
parser.add_option("--rmf", 	dest="rmf", 	help="rmf",   default="sample.rmf")
parser.add_option("--arf", 	dest="arf", 	help="arf",   default="sample.arf")
parser.add_option("--script", 	dest="script", 	help="script",default="sherpa_grid.sph")
parser.add_option("--full", 	dest="full", 	help="full",  default=False, action="store_true")
parser.add_option("--fits", 	dest="fits", 	help="fits",  default="")

(options, args) = parser.parse_args()

debug=options.debug

par1name, par1val, par1pt = format(options.par1).split(":")
par2name, par2val, par2pt = format(options.par2).split(":")

#print par1name, par1val, par1pt
#print par2name, par2val, par2pt

elo, ehi = format(options.range).split(",")

print elo, ehi

etcp = format(options.etcp)
etcp = etcp.replace(':','=')
etcp = etcp.replace(',','\n')

abund=format(options.abund)
#print etcp

def grid(minx,maxx,miny,maxy,bin=1,phase=0.5):

	binx=(maxx-minx+1)/bin
	biny=(maxy-miny+1)/bin

	xgrid = "%.1f:%.1f:#%d" % (minx+phase-1,maxx+phase,binx)
	ygrid = "%.1f:%.1f:#%d" % (miny+phase-1,maxy+phase,biny)

	return xgrid, ygrid

# ----------------------------------------------------------------

somevar = {
	'frac'   : options.frac,
	'model'  : options.model,
	'etcp'	 : etcp,
	'par1name': par1name,
	'par2name': par2name,
	'par1val' : par1val,
	'par2val' : par2val,
	'par1pt'  : par1pt,
	'par2pt'  : par2pt,
	'arf'	  : options.arf,
	'rmf'	  : options.rmf,
	'elo'	  : elo,
	'ehi'	  : ehi,
	}


sherpa_input=r"""

#import scipy.interpolate as interpol
import numpy as np

frac=[%(frac)s]
set_source("faked", "%(model)s")

par1s = [%(par1val)s]
par2s = [%(par2val)s]

%(etcp)s

set_xsabund("%(abund)s")

cpar1s = []
npar1s = len(par1s)
for ii in range(0, npar1s-1):
	for par in np.arange(par1s[ii],par1s[ii+1],1.*npar1s*(par1s[ii+1]-par1s[ii])/%(par1pt)s):
		cpar1s.append(par)

cpar2s = []
npar2s = len(par2s)
for ii in range(0, npar2s-1):
	for par in np.arange(par2s[ii],par2s[ii+1],1.*npar2s*(par2s[ii+1]-par2s[ii])/%(par2pt)s):
		cpar2s.append(par)

print "\t".join(["gridID","par1","par2","E25","E50","E75"])
print "\t".join(["S","NF","NF","N","N","N"])


arfdata = unpack_arf("%(arf)s")
rmfdata = unpack_rmf("%(rmf)s")

""" % {
	'frac'   : options.frac,
	'model'  : options.model,
	'etcp'	 : etcp,
	'abund'	 : abund,
	'par1name': par1name,
	'par2name': par2name,
	'par1val' : par1val,
	'par2val' : par2val,
	'par1pt'  : par1pt,
	'par2pt'  : par2pt,
	'arf'	  : options.arf,
	'rmf'	  : options.rmf,
	'elo'	  : elo,
	'ehi'	  : ehi,
	}




# ----------------------------------------------------------------
# for fits output
# generate dummy table first


fits="""

import pyfits as pf
import modfits as mf
import struct as st


col1=pf.Column(name='gridID', format='7A', array=np.array([' ']))
col2=pf.Column(name='par1',   format='E',  array=np.array([0.0]))
col3=pf.Column(name='par2',   format='E',  array=np.array([0.0]))
col4=pf.Column(name='E25',    format='E',  array=np.array([0.0]), unit='keV')
col5=pf.Column(name='E50',    format='E',  array=np.array([0.0]), unit='keV')
col6=pf.Column(name='E75',    format='E',  array=np.array([0.0]), unit='keV')

cols = pf.ColDefs([col1,col2,col3,col4,col5,col6])
tbhdu = pf.new_table(cols)
hdu = pf.PrimaryHDU(np.arange(1))
thdulist = pf.HDUList([hdu, tbhdu])

thdulist[1].header.update('source',"%(prog)s","creator")
thdulist[1].header.update('frac',  "%(frac)s","quantile fractions")
thdulist[1].header.update('model', "%(model)s","spectral models")
thdulist[1].header.update('par1',  "%(par1)s","par1")
thdulist[1].header.update('par2',  "%(par2)s","par2")
thdulist[1].header.update('etcp',  "%(etcp)s","etc par")
thdulist[1].header.update('arf',   "%(arf)s","arf")
thdulist[1].header.update('rmf',   "%(rmf)s","rmf")
#thdulist[1].header.update('range', "%(elo)s,%(ehi)s","energy range")

thdulist.writeto('%(fits)s',clobber=True)

fits=mf.modfits('%(fits)s')
fits.rewind(1)

""" % {'fits': options.fits,
	'prog': sys.argv[0],
 	'frac'   : options.frac,
 	'model'  : options.model,
 	'etcp'	 : options.etcp,
 	'par1'	 : options.par1,
 	'par2'	 : options.par2,
 	'arf'	 : options.arf,
 	'rmf'	 : options.rmf,
  	'elo'	 : elo,
  	'ehi'	 : ehi,
	}

# ----------------------------------------------------------------
# scipy interpolation
#lin = interpol.interp1d(acy,x)
#qtil=lin(frac) __OUTMODE1

gmode="""
for idx, par1 in enumerate(par1s):
	for par2 in cpar2s:
		%(par1name)s=par1
		%(par2name)s=par2
		fake_pha("faked", arf=arfdata, rmf=rmfdata, exposure=1.e8, grouped=False)
		data=get_data("faked")
		x=data.get_x()
		y=data.counts
		w=((x >= %(elo)s) & (x<= %(ehi)s)).nonzero()[0]
		x=x[w]
		y=y[w]
		acy = np.add.accumulate(y)
		acy = acy/max(acy)
		qtil = np.interp(frac,acy,x) __OUTMODE1


for idx, par2 in enumerate(par2s):
	for par1 in cpar1s:
		%(par1name)s=par1
		%(par2name)s=par2
		fake_pha("faked", arf=arfdata, rmf=rmfdata, exposure=1.e8, grouped=False)
		data=get_data("faked")
		x=data.get_x()
		y=data.counts
		w=((x >= %(elo)s) & (x<= %(ehi)s)).nonzero()[0]
		x=x[w]
		y=y[w]
		acy = np.add.accumulate(y)
		acy = acy/max(acy)
		qtil = np.interp(frac,acy,x) __OUTMODE2


""" % {
	'frac'   : options.frac,
	'model'  : options.model,
	'etcp'	 : etcp,
	'par1name': par1name,
	'par2name': par2name,
	'par1val' : par1val,
	'par2val' : par2val,
	'par1pt'  : par1pt,
	'par2pt'  : par2pt,
	'arf'	  : options.arf,
	'rmf'	  : options.rmf,
	'elo'	  : elo,
	'ehi'	  : ehi,
	}

fitsoutmode1= """
		fits.append("gx%-5d" % idx + st.pack('>fffff',par1,par2,qtil[0],qtil[1],qtil[2]))
"""


fitsoutmode2= """
		fits.append("gy%-5d" % idx  + st.pack('>fffff',par1,par2,qtil[0],qtil[1],qtil[2]))

fits.close()
"""

rdboutmode1= '' """
 		print "\t".join(["gx"+str(idx),str(par1),str(par2),
			str(round(qtil[0],4)),str(round(qtil[1],4)),str(round(qtil[2],4))])
"""

rdboutmode2= '' """
		print "\t".join(["gy"+str(idx),str(par1),str(par2),
			str(round(qtil[0],4)),str(round(qtil[1],4)),str(round(qtil[2],4))])
"""

if options.fits != "" :
	gmode=gmode.replace("__OUTMODE1",fitsoutmode1)
	gmode=gmode.replace("__OUTMODE2",fitsoutmode2)
else:
	gmode=gmode.replace("__OUTMODE1",rdboutmode1)
	gmode=gmode.replace("__OUTMODE2",rdboutmode2)


# ----------------------------------------------------------------
#		lin = interpol.interp1d(acy,x)
#		qtil=lin(frac) __OUTMODE1

mmode="""

for idx, par1 in enumerate(cpar1s):
	for par2 in cpar2s:
		%(par1name)s=par1
		%(par2name)s=par2
		fake_pha("faked", arf=arfdata, rmf=rmfdata, exposure=1.e8, grouped=False)
		data=get_data("faked")
		x=data.get_x()
		y=data.counts
		w=((x >= %(elo)s) & (x<= %(ehi)s)).nonzero()[0]
		x=x[w]
		y=y[w]
		acy = np.add.accumulate(y)
		acy = acy/max(acy)
		qtil = np.interp(frac,acy,x) __OUTMODE1

""" % {
	'frac'   : options.frac,
	'model'  : options.model,
	'etcp'	 : etcp,
	'par1name': par1name,
	'par2name': par2name,
	'par1val' : par1val,
	'par2val' : par2val,
	'par1pt'  : par1pt,
	'par2pt'  : par2pt,
	'arf'	  : options.arf,
	'rmf'	  : options.rmf,
	'elo'	  : elo,
	'ehi'	  : ehi,
	}

fitsoutmode1= """
		fits.append("gx%-5d" % idx  + st.pack('>fffff',par1,par2,qtil[0],qtil[1],qtil[2]))

fits.close()
"""

if options.fits != "" :
	mmode=mmode.replace("__OUTMODE1",fitsoutmode1)
else:
	mmode=mmode.replace("__OUTMODE1",rdboutmode1)

# ----------------------------------------------------------------


header=r"""#
# %(prog)s \
#	--model "%(model)s" \
#	--frac "%(frac)s" \
#	--par1 "%(par1)s" \
#	--par2 "%(par2)s" \
#	--etcp "%(etcp)s" \
#	--rmf "%(rmf)s" \
#	--arf "%(arf)s" \
#	--range "%(elo)s,%(ehi)s"
#""" % {
	'prog': sys.argv[0],
	'frac'   : options.frac,
	'model'  : options.model,
	'etcp'	 : options.etcp,
	'par1'	 : options.par1,
	'par2'	 : options.par2,
	'arf'	 : options.arf,
	'rmf'	 : options.rmf,
	'elo'	 : elo,
	'ehi'	 : ehi,
	}

scname=format(options.script)
fsc=open(scname, "w")
fsc.write(sherpa_input)
if options.fits != "" : fsc.write(fits)
if options.full == True: fsc.write(mmode)
else: fsc.write(gmode)
fsc.close()

#coutput="sherpa_grid.out"
cmd="sherpa"
cmd=["ASCDS_OVERRIDE=1;. /Users/kaya/ciao-4.8/bin/ciao.sh;", cmd, " -b "+scname,
#	"|grep -v sherpa ",
#	"> "+coutput,
	]

p=sp.Popen(" ".join(cmd), stdout=sp.PIPE, shell=True)
out, err = p.communicate()

if options.fits == "":
	print header
	print out.rstrip("\n\r")
