import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

inputfile='lum_100_nh17.txt'
inputfile1='padilla_lum.txt'

with open(inputfile) as f:
	input=f.readlines()
input=[x.strip() for x in input]
input=np.array(map(float,input))

lum=[]
lum1=[]

for src in input:
    lum.append(src)

#n, bins, patches=plt.hist(lum,np.linspace(0.2,1.1,16),histtype='step',cumulative=-1)
n, bins, patches=plt.hist(lum,np.linspace(0.5,2.0,5000),histtype='step',cumulative=-1)

bin_centers=bins[:-1]+0.5*(bins[1:]-bins[:-1])

def func(x,k,a):
	return k*(x**-a)

popt, pcov = curve_fit(func,bin_centers,n,p0=[1,1.5])

plt.plot(bin_centers, func(bin_centers, *popt),'b-')

k=popt[0]
a=popt[1]

print "pcov = %s"%pcov
print "K = %s"%k
print "alpha = %s"%a

perr=np.sqrt(np.diag(pcov))
print "perr = %s"%perr

with open(inputfile1) as f:
	input=f.readlines()
input=[x.strip() for x in input]
input=np.array(map(float,input))

for src in input:
	lum1.append(src)

n1,bins1,patches1=plt.hist(lum1,np.linspace(0.3,1.5,150),histtype='step',cumulative=-1,color='r')

bin_centers1=bins1[:-1]+0.5*(bins1[1:]-bins1[:-1])

popt1, pcov1 = curve_fit(func,bin_centers1,n1,p0=[1,1.5])

plt.plot(bin_centers1, func(bin_centers1, *popt1), 'r-')

k1=popt1[0]
a1=popt1[1]

print "Padilla pcov = %s"%pcov1
print "Padilla K = %s"%k1
print "Padilla alpha = %s"%a1
perr1=np.sqrt(np.diag(pcov1))
print "Padilla perr = %s"%perr

#plt.axis([0.1,1.5,0.0,16.0])
plt.ylim((10**0,10**2))
plt.xlim((2*10**-1,2*10**0))
plt.gca().set_xscale("log")
plt.gca().set_yscale("log")
plt.xlabel('Luminosity (x10^32 ergs/s)')
plt.ylabel('# of Sources >L')
plt.title('LogN-LogL of Soft Sources with net cts>100 Inside 1 pc with nH=6')
plt.show()
