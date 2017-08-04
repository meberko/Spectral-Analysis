import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys,math

class Analyzer():
    def __init__(self):
        self.data = {}
        self.data['names'] = []
        self.data['ra'] = []
        self.data['dec'] = []
        self.data['ph1'] = []
        self.data['ph1_erru'] = []
        self.data['kT1'] = []
        self.data['nh1'] = []
        self.data['hr2'] = []
        self.data['hr2_err'] = []
        self.data['hr5'] = []
        self.data['hr5_err'] = []
        self.data['kT1_erru'] = []
        self.data['kT1_errl'] = []
        self.data['src_sig'] = []
        self.data['chi_sqr'] = []
        self.data['net_cts'] = []
        self.data['net_cts_sigma_up'] = []
        self.data['net_cts_sigma_low'] = []
        self.data['r'] = []
        self.data['chi_sqr_kT1'] = []

        # Cleaned PLI and kT
        self.data['pli'] = []
        self.data['kT'] = []

        self.data_gt_50 = {}
        self.data_gt_50['names'] = []
        self.data_gt_50['ra'] = []
        self.data_gt_50['dec'] = []
        self.data_gt_50['net_cts'] = []
        self.data_gt_50['net_cts_err'] = []
        self.data_gt_50['pli'] = []
        self.data_gt_50['r'] = []
        self.data_gt_50['src_sig'] = []
        self.data_gt_50['nh1'] = []
        self.data_gt_50['hr2'] = []
        self.data_gt_50['hr5'] = []
        self.data_gt_50['hr2_r_lt_25'] = []
        self.data_gt_50['hr2_r_gt_25'] = []
        self.data_gt_50['hr2_err'] = []
        self.data_gt_50['hr5_err'] = []

        self.data_gt_100 = {}
        self.data_gt_100['names'] = []
        self.data_gt_100['ra'] = []
        self.data_gt_100['dec'] = []
        self.data_gt_100['net_cts'] = []
        self.data_gt_100['net_cts_err'] = []
        self.data_gt_100['pli'] = []
        self.data_gt_100['pli_err'] = []
        self.data_gt_100['r'] = []
        self.data_gt_100['src_sig'] = []
        self.data_gt_100['nh1'] = []
        self.data_gt_100['hr2'] = []
        self.data_gt_100['hr2_r_lt_25'] = []
        self.data_gt_100['hr2_r_gt_25'] = []
        self.data_gt_100['hr2_err'] = []
        self.data_gt_100['hr5'] = []
        self.data_gt_100['hr5_err'] = []
        self.data_gt_100['chi_diff'] = []

        self.data_gt_200 = {}
        self.data_gt_200['ra'] = []
        self.data_gt_200['dec'] = []
        self.data_gt_200['net_cts'] = []
        self.data_gt_200['net_cts_err'] = []
        self.data_gt_200['pli'] = []
        self.data_gt_200['r'] = []
        self.data_gt_200['nh1'] = []
        self.data_gt_200['hr2'] = []
        self.data_gt_200['hr2_err'] = []

        self.data_lt_100 = {}
        self.data_lt_100['ra'] = []
        self.data_lt_100['dec'] = []
        self.data_lt_100['net_cts'] = []
        self.data_lt_100['net_cts_err'] = []
        self.data_lt_100['pli'] = []
        self.data_lt_100['r'] = []
        self.data_lt_100['nh1'] = []
        self.data_lt_100['hr2'] = []
        self.data_lt_100['hr2_err'] = []

        self.data_gt_100_lt_200 = {}
        self.data_gt_100_lt_200['ra'] = []
        self.data_gt_100_lt_200['dec'] = []
        self.data_gt_100_lt_200['net_cts'] = []
        self.data_gt_100_lt_200['net_cts_err'] = []
        self.data_gt_100_lt_200['pli'] = []
        self.data_gt_100_lt_200['r'] = []
        self.data_gt_100_lt_200['nh1'] = []
        self.data_gt_100_lt_200['hr2'] = []
        self.data_gt_100_lt_200['hr2_err'] = []

        self.soft_ratios = []
        self.hard_ratios = []
        self.r_ratios = []
        self.soft_ratios_eq_area = []
        self.hard_ratios_eq_area = []
        self.r_ratios_eq_area = []
        self.soft_area_normalized = []
        self.hard_area_normalized = []
        self.pli_r = [22.5,37.5,52.5,67.5,82.5]
        self.pli_r_errs = [7.5,7.5,7.5,7.5,7.5]
        self.pli_avgs = []
        self.pli_avg_errs = []

        # Collect from power law self.data
        err_content = self.getErrData()

        # Collect from thermal self.data
        vapec_content = self.getVapecData()

        for arr in err_content:
            if len(arr)==17 and arr[1]!='CATALOG_NAME':
                self.data['names'].append(arr[1])
                self.data['ra'].append(float(arr[2]))
                self.data['dec'].append(float(arr[3]))
                self.data['ph1'].append(float(arr[4]))
                self.data['ph1_erru'].append(float(arr[5]))
                self.data['src_sig'].append(float(arr[6]))
                self.data['net_cts'].append(float(arr[7]))
                self.data['nh1'].append(float(arr[8]))
                self.data['chi_sqr'].append(float(arr[9]))
                self.data['net_cts_sigma_up'].append(float(arr[10]))
                self.data['net_cts_sigma_low'].append(float(arr[11]))
                self.data['r'].append(float(arr[12]))
                self.data['hr2'].append(float(arr[13]))
                self.data['hr2_err'].append(float(arr[14]))
                self.data['hr5'].append(float(arr[15]))
                self.data['hr5_err'].append(float(arr[16]))

        for arr in vapec_content:
            if len(arr)==9 and arr[1]!='CATALOG_NAME':
                self.data['kT1'].append(float(arr[4]))
                self.data['kT1_erru'].append(float(arr[7]))
                self.data['kT1_errl'].append(float(arr[8]))
                self.data['chi_sqr_kT1'].append(float(arr[6]))

        # Clean kTs
        for k in self.data['kT1']:
            if not np.isnan(k):
                self.data['kT'].append(k)
            else:
                self.data['kT'].append(-1)

        # Clean PLIs
        for p in self.data['ph1']:
            if not np.isnan(p):
                self.data['pli'].append(p)
            else:
                self.data['pli'].append(-10)

        # Find sources with various net counts ranges
        for i in range(0,len(self.data['chi_sqr_kT1'])):
            if self.data['net_cts'][i] > 50:
                self.data_gt_50['names'].append(self.data['names'][i])
                self.data_gt_50['ra'].append(self.data['ra'][i])
                self.data_gt_50['dec'].append(self.data['dec'][i])
                self.data_gt_50['pli'].append(self.data['pli'][i])
                self.data_gt_50['nh1'].append(self.data['nh1'][i])
                self.data_gt_50['r'].append(self.data['r'][i])
                self.data_gt_50['src_sig'].append(self.data['src_sig'][i])
                self.data_gt_50['net_cts'].append(self.data['net_cts'][i])
                self.data_gt_50['net_cts_err'].append(self.data['net_cts_sigma_up'][i])
                self.data_gt_50['hr2'].append(self.data['hr2'][i])
                self.data_gt_50['hr5'].append(self.data['hr5'][i])
                if self.data['r'][i] < 25:
                    self.data_gt_50['hr2_r_lt_25'].append(self.data['hr2'][i])
                if self.data['r'][i] > 25:
                    self.data_gt_50['hr2_r_gt_25'].append(self.data['hr2'][i])
                self.data_gt_50['hr2_err'].append(self.data['hr2_err'][i])
                self.data_gt_50['hr5_err'].append(self.data['hr5_err'][i])
            if self.data['net_cts'][i] > 100:
                self.data_gt_100['names'].append(self.data['names'][i])
                self.data_gt_100['ra'].append(self.data['ra'][i])
                self.data_gt_100['dec'].append(self.data['dec'][i])
                self.data_gt_100['pli'].append(self.data['pli'][i])
                self.data_gt_100['pli_err'].append(self.data['ph1_erru'][i])
                self.data_gt_100['nh1'].append(self.data['nh1'][i])
                self.data_gt_100['r'].append(self.data['r'][i])
                self.data_gt_100['src_sig'].append(self.data['src_sig'][i])
                self.data_gt_100['net_cts'].append(self.data['net_cts'][i])
                self.data_gt_100['net_cts_err'].append(self.data['net_cts_sigma_up'][i])
                self.data_gt_100['hr2'].append(self.data['hr2'][i])
                if self.data['r'][i] < 25:
                    self.data_gt_100['hr2_r_lt_25'].append(self.data['hr2'][i])
                if self.data['r'][i] > 25:
                    self.data_gt_100['hr2_r_gt_25'].append(self.data['hr2'][i])
                self.data_gt_100['hr2_err'].append(self.data['hr2_err'][i])
                self.data_gt_100['hr5'].append(self.data['hr5'][i])
                self.data_gt_100['hr5_err'].append(self.data['hr5_err'][i])
                self.data_gt_100['chi_diff'].append(self.data['chi_sqr'][i] - self.data['chi_sqr_kT1'][i])
            if self.data['net_cts'][i] > 200:
                self.data_gt_200['ra'].append(self.data['ra'][i])
                self.data_gt_200['dec'].append(self.data['dec'][i])
                self.data_gt_200['pli'].append(self.data['pli'][i])
                self.data_gt_200['nh1'].append(self.data['nh1'][i])
                self.data_gt_200['r'].append(self.data['r'][i])
                self.data_gt_200['net_cts'].append(self.data['net_cts'][i])
                self.data_gt_200['net_cts_err'].append(self.data['net_cts_sigma_up'][i])
                self.data_gt_200['hr2'].append(self.data['hr2'][i])
                self.data_gt_200['hr2_err'].append(self.data['hr2_err'][i])
            if self.data['net_cts'][i] < 100:
                self.data_lt_100['ra'].append(self.data['ra'][i])
                self.data_lt_100['dec'].append(self.data['dec'][i])
                self.data_lt_100['pli'].append(self.data['pli'][i])
                self.data_lt_100['nh1'].append(self.data['nh1'][i])
                self.data_lt_100['r'].append(self.data['r'][i])
                self.data_lt_100['net_cts'].append(self.data['net_cts'][i])
                self.data_lt_100['net_cts_err'].append(self.data['net_cts_sigma_up'][i])
                self.data_lt_100['hr2'].append(self.data['hr2'][i])
                self.data_lt_100['hr2_err'].append(self.data['hr2_err'][i])
            if self.data['net_cts'][i] > 100 and self.data['net_cts'][i] < 200:
                self.data_gt_100_lt_200['ra'].append(self.data['ra'][i])
                self.data_gt_100_lt_200['dec'].append(self.data['dec'][i])
                self.data_gt_100_lt_200['pli'].append(self.data['pli'][i])
                self.data_gt_100_lt_200['nh1'].append(self.data['nh1'][i])
                self.data_gt_100_lt_200['r'].append(self.data['r'][i])
                self.data_gt_100_lt_200['net_cts'].append(self.data['net_cts'][i])
                self.data_gt_100_lt_200['net_cts_err'].append(self.data['net_cts_sigma_up'][i])
                self.data_gt_100_lt_200['hr2'].append(self.data['hr2'][i])
                self.data_gt_100_lt_200['hr2_err'].append(self.data['hr2_err'][i])
        print('Num Sources Net Counts > 50: %d out of %d' % (len(self.data_gt_50['hr2']), len(self.data['hr2'])))
        print('Num Sources Net Counts > 100: %d out of %d' % (len(self.data_gt_100['hr2']), len(self.data['hr2'])))
        print('Num Sources Net Counts > 200: %d out of %d' % (len(self.data_gt_200['hr2']), len(self.data['hr2'])))
        print('Num Sources Net Counts < 100: %d out of %d' % (len(self.data_lt_100['hr2']), len(self.data['hr2'])))
        print('Num Sources 100 < Net Counts < 200: %d out of %d' % (len(self.data_gt_100_lt_200['hr2']), len(self.data['hr2'])))

    def getErrData(self):
        with open('data_err.txt') as f:
            content = f.readlines()
        content = [x.strip() for x in content]
        content = [x.split() for x in content]
        return content

    def getVapecData(self):
        with open('data_vapec.txt') as fv:
            content = fv.readlines()
        content = [x.strip() for x in content]
        content = [x.split() for x in content]
        return content

    def plot(self):
        # LOTS of plotting
        plt.figure()

        """
        # Histogram plots
        plt.title('Source Significance Histogram')
        plt.hist(self.data['src_sig'],np.linspace(-10,25,51))
        plt.xlabel('Source Significance (sigma)')
        plt.ylabel('Frequency')
        plt.axis([0,10,0,50])

        plt.title('Net Counts Histogram')
        plt.hist(self.data['net_cts'],np.linspace(0,500,11))
        plt.axis([0,500,0,150])
        plt.xlabel('Net Counts')
        plt.ylabel('Frequency')

        plt.title('kT Histogram')
        plt.hist(self.data['kT'],np.linspace(0,30,16))
        plt.xlabel('kT')
        plt.ylabel('Frequency')

        plt.title('ChiSqr Differential Histogram (Net Counts > 100)')
        plt.hist(self.data_gt_100['chi_diff'], np.linspace(-2,2,20))
        plt.xlabel('ChiSqr Differential (Power Law ChiSqr - Thermal ChiSqr)')
        plt.ylabel('Frequency')

        plt.title('HR2 Histogram for Net Counts > 100')
        plt.hist(self.data_gt_100['hr2'], np.linspace(-1,1,41))
        plt.xlabel('HR2')
        plt.ylabel('Frequency')

        plt.title('PLI Histogram (Net Counts > 100)')
        plt.hist(self.data_gt_100['pli'],np.linspace(-3,3,19))
        plt.xlabel('PLI')
        plt.ylabel('Frequency')
        #plt.axis([0,10,0,50])
        """
        axis_font = { 'size': 20}
        title_font = { 'size': 25}

        plt.subplot(1,2,1)
        plt.title('Sources < 1 pc from Sgr A*', **title_font)
        plt.hist(self.data_gt_100['hr2_r_lt_25'], np.linspace(-1,1,11), rwidth=0.9)
        plt.xlabel('HR2', **axis_font)
        plt.ylabel('Frequency', **axis_font)

        plt.subplot(1,2,2)
        plt.title('Sources > 1 pc from Sgr A*', **title_font)
        plt.hist(self.data_gt_100['hr2_r_gt_25'], np.linspace(-1,1,11), rwidth=0.9)
        plt.xlabel('HR2', **axis_font)
        plt.ylabel('Frequency', **axis_font)

        """
        # Errorbar plots
        i=0
        axis_font = {'size': 20}
        title_font = {'size': 30}
        plt.title('HR3 as a function of HR2 (Net Counts > 100) (Green: R<25", Red: R>25")', **title_font)
        for h in self.data_gt_100['hr2']:
            if self.data_gt_100['r'][i] < 25:
                plt.errorbar(h, self.data_gt_100['hr5'][i],xerr=self.data_gt_100['hr2_err'][i], yerr=self.data_gt_100['hr5_err'][i], ls='None',ecolor='g')
                plt.scatter(h,  self.data_gt_100['hr5'][i], color='g')
            else:
                plt.errorbar(h,self.data_gt_100['hr5'][i],xerr=self.data_gt_100['hr2_err'][i], yerr=self.data_gt_100['hr5_err'][i], ls='None',ecolor='r')
                plt.scatter(h, self.data_gt_100['hr5'][i], color='r')
            i+=1
        plt.scatter(0.65, 0.91, color='black', marker='D', s=50, zorder=2)
        plt.scatter(0.66, 0.92, color='black', marker='D', s=50, zorder=2)
        plt.scatter(0.67, 0.92, color='black', marker='D', s=50, zorder=2)
        plt.scatter(0.66, 0.92, color='black', marker='D', s=50, zorder=2)
        #idx1 = self.data_gt_100['names'].index('174540.04-290030.9')
        #idx2 = self.data_gt_100['names'].index('174541.02-290017.6')
        #idx3 = self.data_gt_100['names'].index('174540.07-290005.7')
        #idx4 = self.data_gt_100['names'].index('174538.07-290022.4')
        #plt.errorbar(self.data_gt_100['hr2'][idx1], self.data_gt_100['hr5'][idx1],xerr=self.data_gt_100['hr2_err'][idx1], yerr=self.data_gt_100['hr5_err'][idx1], ls='None',ecolor='b')
        #plt.scatter(self.data_gt_100['hr2'][idx1],  self.data_gt_100['hr5'][idx1], color='black', marker='D', s=50, zorder=2)
        #plt.errorbar(self.data_gt_100['hr5'][idx2], self.data_gt_100['hr5'][idx2],xerr=self.data_gt_100['hr2_err'][idx2], yerr=self.data_gt_100['hr5_err'][idx2], ls='None',ecolor='black')
        #plt.scatter(self.data_gt_100['hr5'][idx2],  self.data_gt_100['hr5'][idx2], color='black', marker='D', s=50, zorder=2)
        #plt.errorbar(self.data_gt_100['hr5'][idx3], self.data_gt_100['hr5'][idx3],xerr=self.data_gt_100['hr2_err'][idx3], yerr=self.data_gt_100['hr5_err'][idx3], ls='None',ecolor='black')
        #plt.scatter(self.data_gt_100['hr5'][idx3],  self.data_gt_100['hr5'][idx3], color='black', marker='D', s=50, zorder=2)
        #plt.errorbar(self.data_gt_100['hr5'][idx4], self.data_gt_100['hr5'][idx4],xerr=self.data_gt_100['hr2_err'][idx4], yerr=self.data_gt_100['hr5_err'][idx4], ls='None',ecolor='black')
        #plt.scatter(self.data_gt_100['hr5'][idx4],  self.data_gt_100['hr5'][idx4], color='black', marker='D', s=50, zorder=2)
        plt.axis([-0.2,1.2,-0.2,1.6])
        plt.xlabel('HR2', **axis_font)
        plt.ylabel('HR3', **axis_font)

        plt.title('HR2 as a function of Net Counts (Net Counts > 100)')
        plt.errorbar(self.data_gt_100['net_cts'],self.data_gt_100['hr2'],xerr=self.data_gt_100['net_cts_err'], yerr=self.data_gt_100['hr2_err'], ls='None')
        plt.axis([90,1500,-0.4,1.2])
        plt.xlabel('Net Counts')
        plt.ylabel('HR2')
        plt.xscale('log')

        print self.pli_avgs
        print self.pli_avg_errs
        plt.title('PLI Weighted Averages by Radius Bins (Net Counts > 100)')
        plt.errorbar(self.pli_r, self.pli_avgs, xerr=self.pli_r_errs, yerr=self.pli_avg_errs, ls='None', capsize=1)
        plt.axis([0,90,-1.0,2.0])
        plt.xlabel('Radius (")')
        plt.ylabel('PLI')
        axis_font = { 'size': 20}
        title_font = { 'size': 30}
        i=0
        p_r=[]
        p_h=[]
        p_he = []
        for x in self.data_gt_100['r']:
            if x > 5:
                p_r.append(x/25.8)
                p_h.append(self.data_gt_100['hr2'][i])
                p_he.append(self.data_gt_100['hr2_err'][i])
            i+=1
        plt.title('HR2 as a function of Radial Distance from Sgr A* (Net Counts > 100)', **title_font)
        plt.errorbar(p_r, p_h, yerr=p_he, ls='None', capsize=2)
        plt.axis([0,3.2,-0.4,1.2])
        plt.plot([0,3.2], [0.48,0.48])
        plt.xlabel('Radial Distance from Sgr A* (pc)', **axis_font)
        plt.ylabel('HR2', **axis_font)

        # Scatter plots
        plt.title('Detection Threshold As Function of Net Counts')
        plt.scatter(self.data['net_cts'], self.data['src_sig'])
        plt.axis([0, 3100,0,50])
        plt.xlabel('Net Counts')
        plt.ylabel('Source Significance (sigma)')

        plt.title('Percent PLI > 1 as a function of Radius')
        plt.plot(self.r_ratios, self.soft_ratios)
        plt.xlabel('Radius (")')
        plt.ylabel('Percent PLI > 1')

        plt.title('Percent PLI > 1 as a function of Radius (Each Point Derived From Equal Area)')
        plt.plot(self.r_ratios_eq_area, self.soft_ratios_eq_area)
        plt.xlabel('Radius (")')
        plt.ylabel('Percent PLI > 1')

        plt.title('# of Soft Sources Normalized by Annulus Area as a function of Radius (Net Counts > 50)')
        plt.scatter([0,5,10,15,25,30,40,45,50,55,60,70,80], self.soft_area_normalized)
        plt.axis([-0.75,81,-0.001,0.060])
        plt.xlabel('Outer Radius (\')')
        plt.ylabel('# of Soft Sources / Area of Annulus')

        plt.title('# of Hard Sources Normalized by Annulus Area as a function of Radius (Net Counts > 50)')
        plt.scatter([0,5,10,15,25,30,40,45,50,55,60,70,80], self.hard_area_normalized)
        plt.axis([-0.75,81,-0.001,0.060])
        plt.xlabel('Outer Radius (\')')
        plt.ylabel('# of Hard Sources / Area of Annulus')

        plt.title('HR2 as a function of Radius (Net Counts > 100)')
        plt.scatter(self.data_gt_100['r'], self.data_gt_100['hr2'])
        plt.axis([0,80,-1.0,1.2])
        plt.xlabel('Radius (")')
        plt.ylabel('HR2')

        plt.title('PLI as a function of Radius (Net Counts > 100)')
        plt.scatter(self.data_gt_100['r'], self.data_gt_100['pli'])
        plt.axis([0,80,-3,6])
        plt.xlabel('Radius (")')
        plt.ylabel('PLI')

        plt.title('PLI as a function of Net Counts')
        plt.scatter(self.data['net_cts'],self.data['pli'])
        plt.axis([0,500,-3,6])
        plt.xlabel('Net Counts')
        plt.ylabel('PLI')

        plt.title('PLI as a function of nH (Net Counts > 100)')
        plt.scatter(self.data_gt_100['nh1'],self.data_gt_100['pli'])
        plt.xlabel('nH')
        plt.ylabel('PLI')

        plt.title('kT as a function of Net Counts')
        plt.errorbar(self.data['net_cts'], self.data['kT'],xerr=self.data['net_cts_sigma_up'],yerr=self.data['kT1_erru'], ls='None')
        plt.axis([0,500,0,30])
        plt.xlabel('Net Counts')
        plt.ylabel('kT')

        plt.title('HR2 as a function of kT')
        plt.errorbar(self.data['kT'], self.data['hr2'], xerr=self.data['kT1_erru'], yerr=self.data['hr2_err'], ls='None')
        plt.axis([0,30,-1,1.2])
        plt.xlabel('kT')
        plt.ylabel('HR2')

        plt.title('HR5 as a function of HR2 (Net Counts > 100)')
        plt.errorbar(self.data_gt_100['hr2'],self.data_gt_100['hr5'],xerr=self.data_gt_100['hr2_err'],yerr=self.data_gt_100['hr5_err'], ls='None')
        plt.axis([0,1,0,1.5])
        plt.xlabel('HR2')
        plt.ylabel('HR5')

        plt.title('Source Significance as a function of Net Counts')
        plt.scatter(self.data['net_cts'],self.data['src_sig'])
        #plt.axis([0,500,-2,10])
        plt.xlabel('Net Counts')
        plt.ylabel('Source Signficance (sigma)')

        plt.title('KT as a function of PLI')
        plt.scatter(self.data['pli'], self.data['kT'])
        plt.axis([-3,6,0,10])
        plt.xlabel('PLI')
        plt.ylabel('kT')

        plt.title('Chi Sqr as a function of Thermal Chi Sqr')
        plt.scatter(self.data['chi_sqr'],self.data['chi_sqr_kT1'])
        plt.plot(range(0,3))
        plt.axis([0,1.5,0,1.5])
        plt.xlabel('Thermal Chi Sqr')
        plt.ylabel('Chi Sqr')
        """

        plt.show()

    def makeReg(self):
        i=0
        with open('pli.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for p in self.data['pli']:
                if p != -10:
                    if p > 1:
                        f.write(('circle(%f,%f,1")' % (self.data['ra'][i],self.data['dec'][i])))
                        f.write('\n')
                    else:
                        f.write(('circle(%f,%f,1") # color=red' % (self.data['ra'][i],self.data['dec'][i])))
                        f.write('\n')
                i+=1
        i=0
        with open('pli_gt_100.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for p in self.data_gt_100['pli']:
                if p!=-10:
                    if p>1:
                        f.write(('circle(%f,%f,1")' % (self.data_gt_100['ra'][i],self.data_gt_100['dec'][i])))
                        f.write('\n')
                    else:
                        f.write(('circle(%f,%f,1") # color=red' % (self.data_gt_100['ra'][i],self.data_gt_100['dec'][i])))
                        f.write('\n')
                i+=1
        i=0
        with open('pli_gt_200.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for p in self.data_gt_200['pli']:
                if p!=-10:
                    if p>1:
                        f.write(('circle(%f,%f,1")' % (self.data_gt_200['ra'][i],self.data_gt_200['dec'][i])))
                        f.write('\n')
                    else:
                        f.write(('circle(%f,%f,1") # color=red' % (self.data_gt_200['ra'][i],self.data_gt_200['dec'][i])))
                        f.write('\n')
                i+=1
        i=0
        with open('pli_lt_100.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for p in self.data_lt_100['pli']:
                if p!=-10:
                    if p>1:
                        f.write(('circle(%f,%f,1")' % (self.data_lt_100['ra'][i],self.data_lt_100['dec'][i])))
                        f.write('\n')
                    else:
                        f.write(('circle(%f,%f,1") # color=red' % (self.data_lt_100['ra'][i],self.data_lt_100['dec'][i])))
                        f.write('\n')
                i+=1
        i=0
        with open('pli_gt_100_lt_200.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for p in self.data_gt_100_lt_200['pli']:
                if p!=-10:
                    if p>1:
                        f.write(('circle(%f,%f,1")' % (self.data_gt_100_lt_200['ra'][i],self.data_gt_100_lt_200['dec'][i])))
                        f.write('\n')
                    else:
                        f.write(('circle(%f,%f,1") # color=red' % (self.data_gt_100_lt_200['ra'][i],self.data_gt_100_lt_200['dec'][i])))
                        f.write('\n')
                i+=1
        i=0
        with open('hr2.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for h in self.data['hr2']:
                if h<=0.3:
                    f.write(('circle(%f,%f,1")' % (self.data['ra'][i],self.data['dec'][i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (self.data['ra'][i],self.data['dec'][i])))
                    f.write('\n')
                    i+=1
        i=0
        with open('hr2_gt_50.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for h in self.data_gt_50['hr2']:
                if h<=0.3:
                    f.write(('circle(%f,%f,1")' % (self.data_gt_50['ra'][i],self.data_gt_50['dec'][i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (self.data_gt_50['ra'][i],self.data_gt_50['dec'][i])))
                    f.write('\n')
                i+=1
        i=0
        with open('hr2_gt_100.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for h in self.data_gt_100['hr2']:
                if h<=0.3:
                    f.write(('circle(%f,%f,1")' % (self.data_gt_100['ra'][i],self.data_gt_100['dec'][i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (self.data_gt_100['ra'][i],self.data_gt_100['dec'][i])))
                    f.write('\n')
                i+=1
        i=0
        with open('hr2_gt_200.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for h in self.data_gt_200['hr2']:
                if h<=0.3:
                    f.write(('circle(%f,%f,1")' % (self.data_gt_200['ra'][i],self.data_gt_200['dec'][i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (self.data_gt_200['ra'][i],self.data_gt_200['dec'][i])))
                    f.write('\n')
                i+=1
        i=0
        with open('hr2_lt_100.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for h in self.data_lt_100['hr2']:
                if h<=0.3:
                    f.write(('circle(%f,%f,1")' % (self.data_lt_100['ra'][i],self.data_lt_100['dec'][i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (self.data_lt_100['ra'][i],self.data_lt_100['dec'][i])))
                    f.write('\n')
                i+=1
        i=0
        with open('hr2_gt_100_lt_200.reg', 'w') as f:
            f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
            for h in self.data_gt_100_lt_200['hr2']:
                if h<=0.3:
                    f.write(('circle(%f,%f,1")' % (self.data_gt_100_lt_200['ra'][i],self.data_gt_100_lt_200['dec'][i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (self.data_gt_100_lt_200['ra'][i],self.data_gt_100_lt_200['dec'][i])))
                    f.write('\n')
                i+=1

    def runKSTest(self):
        ks1 = []
        ks2 = []
        i=0
        for h in self.data_gt_50['hr2']:
            if self.data_gt_50['r'][i] < 25:
                ks1.append(h)
            else:
                ks2.append(h)
            i+=1
        print stats.ks_2samp(ks1,ks2)

    def runFlatChiSqrTest(self):
        plt.figure()
        for n in np.linspace(0.1,2,191):
            chi = 0
            j=0
            for h in self.data['hr2']:
                chi+=(h-n)**2/(388*float(self.data['hr2_err'][j])**2)
                j+=1
            plt.scatter(n,chi)
            print(('Constant HR2: %f, Chi Sqr: %f')%(n,chi))
        plt.show()

    def printGT50(self):
        with open('sources_nc_gt_50.txt', 'w') as f:
            i=0
            for r in self.data_gt_50['r']:
                if self.data_gt_50['src_sig'][i] > 7:
                    f.write(self.data_gt_50['names'][i]+' '+str(self.data_gt_50['src_sig'][i])+'\n')
                i+=1

    def printGT100(self):
        with open('sources_nc_gt_100.txt', 'w') as f:
            i=0
            for r in self.data_gt_100['r']:
                if self.data_gt_100['src_sig'][i] > 7:
                    f.write(self.data_gt_100['names'][i]+' '+str(self.data_gt_100['src_sig'][i])+'\n')
                i+=1

    def printInside25GT100(self):
        with open('sources_lt_25_nc_gt_100.txt', 'w') as f:
            i=0
            for r in self.data_gt_100['r']:
                if r < 25:
                    f.write(self.data_gt_100['names'][i]+' '+str(self.data_gt_100['hr2'][i])+'\n')
                i+=1

    def printOutside25GT100(self):
        with open('sources_gt_25_nc_gt_100.txt', 'w') as f:
            i=0
            for r in self.data_gt_100['r']:
                if r > 25:
                    f.write(self.data_gt_100['names'][i]+'\n')
                i+=1

    def printHardGT50(self):
        with open('hard_sources_nc_gt_50','w') as f:
            i=0
            for h in self.data_gt_50['hr2']:
                if h > 0.3:
                    f.write(('%s\t%f\t%f\t%f\t%f\t%f\n')%(  self.data_gt_50['names'][i],
                                                        self.data_gt_50['r'][i],
                                                        self.data_gt_50['ra'][i],
                                                        self.data_gt_50['dec'][i],
                                                        self.data_gt_50['net_cts'][i],
                                                        self.data_gt_50['hr2'][i]))
                i+=1

    def printSoftGT50(self):
        i=0
        for h in self.data_gt_50['hr2']:
            if h < 0.3:
                print(('%s\t%f\t%f\t%f\t%f\t%f')%(  self.data_gt_50['names'][i],
                                                    self.data_gt_50['r'][i],
                                                    self.data_gt_50['ra'][i],
                                                    self.data_gt_50['dec'][i],
                                                    self.data_gt_50['net_cts'][i],
                                                    self.data_gt_50['hr2'][i]))
            i+=1

    def printSoftGT100(self):
        i=0
        soft_ct=0
        for h in self.data_gt_100['hr2']:
            if self.data_gt_100['r'][i]<25:
                if h < 0.3:
                    soft_ct+=1
                    print(('%s\t%f\t%f\t%f\t%f\t%f')%(  self.data_gt_100['names'][i],
                                                        self.data_gt_100['r'][i],
                                                        self.data_gt_100['ra'][i],
                                                        self.data_gt_100['dec'][i],
                                                        self.data_gt_100['net_cts'][i],
                                                        self.data_gt_100['hr2'][i]))
                i+=1
        print soft_ct
        print i
        print float(soft_ct)/float(i)

    def softHardCounting(self):
        soft = hard = 0
        for i in range(0,388):
            if i%25==0 or i==387:
                self.soft_ratios.append(float(soft)/float(25))
                self.hard_ratios.append(float(hard)/float(25))
                self.r_ratios.append(self.data['r'][i])
                soft = hard = 0
            if self.data['hr2'][i] < 0:
                soft+=1
            else:
                hard+=1

        soft = hard = 0
        for i in range(1,30):
            self.r_ratios_eq_area.append(22*np.sqrt(i))
            total = soft = hard = 0
            for j in range(0,388):
                if 22*np.sqrt(i-1) < self.data['r'][j] and self.data['r'][j] < 22*np.sqrt(i):
                    total+=1
                    if self.data['hr2'][j]<0:
                        soft+=1
                    else:
                        hard+=1
            if total!=0:
                self.soft_ratios_eq_area.append(float(soft)/float(total))
                self.hard_ratios_eq_area.append(float(hard)/float(total))
            else:
                self.soft_ratios_eq_area.append(0)
                self.hard_ratios_eq_area.append(0)

    def normalizedCounting(self):
        i=0
        ann_arr5 = math.pi*5**2
        ann_arr10 = math.pi*(10**2-5**2)
        ann_arr15 = math.pi*(15**2-10**2)
        ann_arr25 = math.pi*(25**2-15**2)
        ann_arr30 = math.pi*(30**2-25**2)
        ann_arr40 = math.pi*(40**2-30**2)
        ann_arr45 = math.pi*(45**2-40**2)
        ann_arr50 = math.pi*(50**2-45**2)
        ann_arr55 = math.pi*(55**2-50**2)
        ann_arr60 = math.pi*(60**2-55**2)
        ann_arr70 = math.pi*(70**2-60**2)
        ann_arr80 = math.pi*(80**2-70**2)
        n5=n10=n15=n25=n30=n40=n45=n50=n55=n60=n70=n80=0
        h5=h10=h15=h25=h30=h40=h45=h50=h55=h60=h70=h80=0

        for h in self.data['hr2']:
            if 0 <= self.data['r'][i] and self.data['r'][i] < 5:
                if h < 0.3:
                    n5 += 1
                else:
                    h5 += 1
            if 5 <= self.data['r'][i] and self.data['r'][i] < 10:
                if h < 0.3:
                    n10 += 1
                else:
                    h10 += 1
            if 10 <= self.data['r'][i] and self.data['r'][i] < 15:
                if h < 0.3:
                    n15 += 1
                else:
                    h15 += 1
            if 15 <= self.data['r'][i] and self.data['r'][i] < 25:
                if h < 0.3:
                    n25 += 1
                else:
                    h25 += 1
            if 25 <= self.data['r'][i] and self.data['r'][i] < 30:
                if h < 0.3:
                    n30 += 1
                else:
                    h30 += 1
            if 30 <= self.data['r'][i] and self.data['r'][i] < 40:
                if h < 0.3:
                    n40 += 1
                else:
                    h40 += 1
            if 40 <= self.data['r'][i] and self.data['r'][i] < 45:
                if h < 0.3:
                    n45 += 1
                else:
                    h45 += 1
            if 45 <= self.data['r'][i] and self.data['r'][i] < 50:
                if h < 0.3:
                    n50 += 1
                else:
                    h50 += 1
            if 50 <= self.data['r'][i] and self.data['r'][i] < 55:
                if h < 0.3:
                    n55 += 1
                else:
                    h55 += 1
            if 55 <= self.data['r'][i] and self.data['r'][i] < 60:
                if h < 0.3:
                    n60 += 1
                else:
                    h60 += 1
            if 60 <= self.data['r'][i] and self.data['r'][i] < 70:
                if h < 0.3:
                    n70 += 1
                else:
                    h70 += 1
            if 70 <= self.data['r'][i] and self.data['r'][i] < 80:
                if h < 0.3:
                    n80 += 1
                else:
                    h80 += 1
            i+=1
        self.soft_area_normalized.append(0)
        self.soft_area_normalized.append((float(n5)-3)/float(ann_arr5))
        self.soft_area_normalized.append(float(n10)/float(ann_arr10))
        self.soft_area_normalized.append(float(n15)/float(ann_arr15))
        self.soft_area_normalized.append(float(n25)/float(ann_arr25))
        self.soft_area_normalized.append(float(n30)/float(ann_arr30))
        self.soft_area_normalized.append(float(n40)/float(ann_arr40))
        self.soft_area_normalized.append(float(n60)/float(ann_arr60))
        self.soft_area_normalized.append(float(n45)/float(ann_arr45))
        self.soft_area_normalized.append(float(n50)/float(ann_arr50))
        self.soft_area_normalized.append(float(n55)/float(ann_arr55))
        self.soft_area_normalized.append(float(n70)/float(ann_arr70))
        self.soft_area_normalized.append(float(n80)/float(ann_arr80))
        self.hard_area_normalized.append(0)
        self.hard_area_normalized.append(float(h5)/float(ann_arr5))
        self.hard_area_normalized.append(float(h10)/float(ann_arr10))
        self.hard_area_normalized.append(float(h15)/float(ann_arr15))
        self.hard_area_normalized.append(float(h25)/float(ann_arr25))
        self.hard_area_normalized.append(float(h30)/float(ann_arr30))
        self.hard_area_normalized.append(float(h40)/float(ann_arr40))
        self.hard_area_normalized.append(float(h45)/float(ann_arr45))
        self.hard_area_normalized.append(float(h50)/float(ann_arr50))
        self.hard_area_normalized.append(float(h55)/float(ann_arr55))
        self.hard_area_normalized.append(float(h60)/float(ann_arr60))
        self.hard_area_normalized.append(float(h70)/float(ann_arr70))
        self.hard_area_normalized.append(float(h80)/float(ann_arr80))

    def pliAvgCalculations(self):
        i=0
        n15=n30=n45=n60=n75=0
        sum15=sum30=sum45=sum60=sum75=0
        squerr15=squerr30=squerr45=squerr60=squerr75=0
        for p in self.data_gt_100['pli']:
            if 15 < self.data_gt_100['r'][i] < 30:
                n15+=1
                sum15+=p
                squerr15+=self.data_gt_100['pli_err'][i]**2/p**2
            if 30 < self.data_gt_100['r'][i] < 45:
                n30+=1
                sum30+=p
                squerr30+=self.data_gt_100['pli_err'][i]**2/p**2
            if 45 < self.data_gt_100['r'][i] < 60:
                n45+=1
                sum45+=p
                squerr45+=self.data_gt_100['pli_err'][i]**2/p**2
            if 60 < self.data_gt_100['r'][i] < 75:
                n60+=1
                sum60+=p
                squerr60+=self.data_gt_100['pli_err'][i]**2/p**2
            if 75 < self.data_gt_100['r'][i] < 90:
                n75+=1
                sum75+=p
                squerr75+=self.data_gt_100['pli_err'][i]**2/p**2
            i+=1
        avg15 = sum15/n15
        err15 = avg15*np.sqrt(squerr15)/n15**2
        self.pli_avgs.append(avg15)
        self.pli_avg_errs.append(err15)

        avg30 = sum30/n30
        err30 = avg30*np.sqrt(squerr30)/n30**2
        self.pli_avgs.append(avg30)
        self.pli_avg_errs.append(err30)

        avg45 = sum45/n45
        err45 = avg45*np.sqrt(squerr45)/n45**2
        self.pli_avgs.append(avg45)
        self.pli_avg_errs.append(err45)

        avg60 = sum60/n60
        err60 = avg60*np.sqrt(squerr60)/n60**2
        self.pli_avgs.append(avg60)
        self.pli_avg_errs.append(err60)

        avg75 = sum75/n75
        err75 = avg75*np.sqrt(squerr75)/n75**2
        self.pli_avgs.append(avg75)
        self.pli_avg_errs.append(err75)

if __name__ == '__main__':
    a = Analyzer()
    a.softHardCounting()
    a.normalizedCounting()
    a.printHardGT50()
    #a.printGT50()
    #a.printGT100()
    #a.printInside25GT100()
    #a.printOutside25GT100()
    #a.pliAvgCalculations()
    #a.makeReg()
    #a.printSoftGT100()
    #a.runFlatChiSqrTest()
    #a.runKSTest()
    a.plot()
