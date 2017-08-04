import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

		n1,bins1,patches1=plt.hist(lum1,np.linspace(0.1,100,50),histtype='step',cumulative=-1,color='r')

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

		n, bins, patches=plt.hist(lum,np.linspace(2.0,6.0,5),histtype='step',cumulative=-1)

		bin_centers=bins[:-1]+0.5*(bins[1:]-bins[:-1])


		popt, pcov = curve_fit(func,bin_centers,n,p0=[1.0,2.0])

		plt.plot(bin_centers, func(bin_centers, *popt),'b-')

		k=popt[0]
		a=popt[1]

		print "pcov = %s"%pcov
		print "K = %s"%k
		print "alpha = %s"%a
		perr=np.sqrt(np.diag(pcov))
		print "perr = %s"%perr
		print

if __name__=='__main__':
	titlefont = { "size": 30 }
	axisfont = { "size": 20 }
	h = HistFit('padilla_flux.txt')
	h.hist('flux_100_nhfree.txt')
	h.histPadilla()
	plt.plot([0.19,0.19],[0,100], 'b--')
	#plt.axis([0.1,1.5,0.0,16.0])
	plt.ylim((10**0,10**2))
	plt.xlim((1.0*10**(-2),1.0*10**3))
	plt.gca().set_xscale("log")
	plt.gca().set_yscale("log")
	plt.xlabel('Flux [e-15 ergs / (s*cm^2)]', **axisfont)
	plt.ylabel('# of Sources >F', **axisfont)
	plt.title('LogN-LogF of Soft Sources with net cts>100 Inside 1 pc', **titlefont)
	plt.show()
