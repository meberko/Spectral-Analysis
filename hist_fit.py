import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from scipy.optimize import curve_fit
from matplotlib import ticker

def func(x,k,a):
	return k*(x**-a)

class HistFit:
	def __init__(self, ip):
		self.inputfile_padilla=ip

	def histPadilla(self):
		lum1=[]
		with open(self.inputfile_padilla) as f:
			content=f.readlines()
		content=[x.strip() for x in content]
		content=np.array(map(float,content))

		for src in content:
			lum1.append(src)

		n1,bins1,patches1=plt.hist(lum1,np.linspace(0.1,3,20),histtype='step',cumulative=-1,color='r')

		bin_centers1=bins1[:-1]+0.5*(bins1[1:]-bins1[:-1])

		popt1, pcov1 = curve_fit(func,bin_centers1,n1,p0=[1,1.5])

		plt.plot(bin_centers1, func(bin_centers1, *popt1), 'r-')

		k1=popt1[0]
		a1=popt1[1]

		print "Padilla pcov = %s"%pcov1
		print "Padilla K = %s"%k1
		print "Padilla alpha = %s"%a1
		perr1=np.sqrt(np.diag(pcov1))
		print "Padilla perr = %s"%perr1
		print

	def hist(self, ifname):
		with open(ifname) as f:
			content=f.readlines()
		content=[x.strip() for x in content]
		content=np.array(map(float,content))

		lum=[]

		for src in content:
		    lum.append(src)

		n, bins, patches=plt.hist(lum,np.linspace(2.0,6.0,20),histtype='step',cumulative=-1)

		bin_centers=bins[:-1]+0.5*(bins[1:]-bins[:-1])


		popt, pcov = curve_fit(func,bin_centers,n,p0=[1.0,2.0])

		plt.plot(bin_centers, func(bin_centers, *popt),'b-')
		plt.arrow(0.19,1.5,0,-0.4, head_width=0.01, head_length=0.1, fc='k', ec='k')

		k=popt[0]
		a=popt[1]

		print "pcov = %s"%pcov
		print "K = %s"%k
		print "alpha = %s"%a
		perr=np.sqrt(np.diag(pcov))
		print "perr = %s"%perr
		font = { 'size': 25 }
		#plt.text(6, 8, r'$\alpha = $%.2f$\pm$%.2f'%(float(a),float(perr[1])), **font)

def log_10_product(x, pos):
	"""The two args are the value and tick position.
	Label ticks with the product of the exponentiation"""
	return '%1i' % (x)

if __name__=='__main__':
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	h = HistFit('padilla_flux.txt')
	h.hist('flux_100_nhfree.txt')
	h.histPadilla()
	plt.axis([0.1,1.5,0.0,16.0])
	#x = np.linspace(0.1,2,100)
	#y = 43.5652331909*x**(-1.49548499323)
	#plt.plot(x,y, 'b--')
	#plt.plot([0.19,0.19],[1,1000], 'b--')
	axis_font = { 'size': 40 }
	plt.tick_params(axis='both', which='major', labelsize=30)
	plt.tick_params(axis='both', which='minor', labelsize=30)
	plt.ylim((10**0,2*10**1))
	plt.xlim((1.5*10**0,8*10**0))
	plt.gca().set_xscale("log")
	plt.gca().set_yscale("log")
	formatter = FuncFormatter(log_10_product)
	ax.yaxis.set_major_formatter(formatter)
	ax.xaxis.set_minor_formatter(formatter)
	plt.xlabel('2-8 keV Flux [$10^{-15}$ erg cm$^{-2}$ s$^{-1}$]', **axis_font)
	plt.ylabel('N(>S)', **axis_font)
	plt.show()
