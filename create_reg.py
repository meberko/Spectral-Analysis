import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys,math

def main():
    names = []
    ra = []
    dec = []

    ph1 = []
    kT1 = []
    nh1 = []
    hr2 = []
    hr2_err = []
    hr5 = []
    hr5_err = []
    kT1_erru = []
    kT1_errl = []

    # Cleaned PLI and kT
    pli = []
    kT = []

    names_gt_50 = []
    ra_gt_50 = []
    dec_gt_50 = []
    net_cts_gt_50 = []
    net_cts_err_gt_50 = []
    pli_gt_50 = []
    r_gt_50 = []
    nh1_gt_50 = []
    hr2_gt_50 = []
    hr2_gt_50_r_lt_25 = []
    hr2_gt_50_r_gt_25 = []
    hr2_err_gt_50 = []

    names_gt_100 = []
    ra_gt_100 = []
    dec_gt_100 = []
    net_cts_gt_100 = []
    net_cts_err_gt_100 = []
    pli_gt_100 = []
    r_gt_100 = []
    nh1_gt_100 = []
    hr2_gt_100 = []
    hr2_gt_100_r_lt_25 = []
    hr2_gt_100_r_gt_25 = []
    hr2_err_gt_100 = []
    hr5_gt_100 = []
    hr5_err_gt_100 = []

    ra_gt_200 = []
    dec_gt_200 = []
    net_cts_gt_200 = []
    net_cts_err_gt_200 = []
    pli_gt_200 = []
    r_gt_200 = []
    nh1_gt_200 = []
    hr2_gt_200 = []
    hr2_err_gt_200 = []

    ra_lt_100 = []
    dec_lt_100 = []
    net_cts_lt_100 = []
    net_cts_err_lt_100 = []
    pli_lt_100 = []
    r_lt_100 = []
    nh1_lt_100 = []
    hr2_lt_100 = []
    hr2_err_lt_100 = []

    ra_gt_100_lt_200 = []
    dec_gt_100_lt_200 = []
    net_cts_gt_100_lt_200 = []
    net_cts_err_gt_100_lt_200 = []
    pli_gt_100_lt_200 = []
    r_gt_100_lt_200 = []
    nh1_gt_100_lt_200 = []
    hr2_gt_100_lt_200 = []
    hr2_err_gt_100_lt_200 = []

    soft_ratios = []
    hard_ratios = []
    r_ratios = []
    soft_ratios_eq_area = []
    hard_ratios_eq_area = []
    r_ratios_eq_area = []
    soft_area_normalized = []
    hard_area_normalized = []

    src_sig = []
    chi_sqr = []
    net_cts = []
    net_cts_sigma_up = []
    net_cts_sigma_low = []
    r = []
    chi_sqr_kT1 = []
    chi_diff = []

    # Collect from power law data
    with open('data_err.txt') as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    content = [x.split() for x in content]

    for arr in content:
        if len(arr)==16 and arr[1]!='CATALOG_NAME':
            names.append(arr[1])
            ra.append(float(arr[2]))
            dec.append(float(arr[3]))
            ph1.append(float(arr[4]))
            src_sig.append(float(arr[5]))
            net_cts.append(float(arr[6]))
            nh1.append(float(arr[7]))
            chi_sqr.append(float(arr[8]))
            net_cts_sigma_up.append(float(arr[9]))
            net_cts_sigma_low.append(float(arr[10]))
            r.append(float(arr[11]))
            hr2.append(float(arr[12]))
            hr2_err.append(float(arr[13]))
            hr5.append(float(arr[14]))
            hr5_err.append(float(arr[15]))

    # Collect from thermal data
    with open('data_vapec.txt') as fv:
        content = fv.readlines()
    content = [x.strip() for x in content]
    content = [x.split() for x in content]
    i=0
    for arr in content:
        if len(arr)==9 and arr[1]!='CATALOG_NAME':
            kT1.append(float(arr[4]))
            kT1_erru.append(float(arr[7]))
            kT1_errl.append(float(arr[8]))
            chi_sqr_kT1.append(float(arr[6]))
            #if np.isnan(kT1_erru[i]):
            #    print hr2[i]
            i+=1

    # Clean kTs
    for k in kT1:
        if not np.isnan(k):
            kT.append(k)
        else:
            kT.append(-1)

    # Clean PLIs
    for p in ph1:
        if not np.isnan(p):
            pli.append(p)
        else:
            pli.append(-10)
    i=0
    """
    # Flat chi square tests
    for n in np.linspace(0.1,2,191):
        chi = 0
        j=0
        for h in hr2:
            chi+=(h-n)**2/(388*float(hr2_err[j])**2)
            j+=1
        print(('Constant HR2: %f, Chi Sqr: %f')%(n,chi))

    soft = hard = 0
    for i in range(0,317):
        if i%25==0 or i==316:
            soft_ratios.append(float(soft)/float(25))
            hard_ratios.append(float(hard)/float(25))
            r_ratios.append(r[i])
            soft = hard = 0
        if hr2[i] < 0:
            soft+=1
        else:
            hard+=1

    soft = hard = 0
    for i in range(1,30):
        r_ratios_eq_area.append(22*np.sqrt(i))
        total = soft = hard = 0
        for j in range(0,317):
            if 22*np.sqrt(i-1) < r[j] and r[j] < 22*np.sqrt(i):
                total+=1
                if hr2[j]<0:
                    soft+=1
                else:
                    hard+=1
        if total!=0:
            soft_ratios_eq_area.append(float(soft)/float(total))
            hard_ratios_eq_area.append(float(hard)/float(total))
        else:
            soft_ratios_eq_area.append(0)
            hard_ratios_eq_area.append(0)
    #print r_ratios_eq_area
    #print soft_ratios_eq_area
    """

    # Find sources with various net counts ranges
    i=0
    for i in range(0,len(chi_sqr_kT1)):
        if net_cts[i] > 50:
            names_gt_50.append(names[i])
            ra_gt_50.append(ra[i])
            dec_gt_50.append(dec[i])
            pli_gt_50.append(pli[i])
            nh1_gt_50.append(nh1[i])
            r_gt_50.append(r[i])
            net_cts_gt_50.append(net_cts[i])
            net_cts_err_gt_50.append(net_cts_sigma_up[i])
            hr2_gt_50.append(hr2[i])
            if r[i] < 25:
                hr2_gt_50_r_lt_25.append(hr2[i])
            if r[i] > 25:
                hr2_gt_50_r_gt_25.append(hr2[i])
            hr2_err_gt_50.append(hr2_err[i])
        if net_cts[i] > 100:
            names_gt_100.append(names[i])
            ra_gt_100.append(ra[i])
            dec_gt_100.append(dec[i])
            pli_gt_100.append(pli[i])
            nh1_gt_100.append(nh1[i])
            r_gt_100.append(r[i])
            net_cts_gt_100.append(net_cts[i])
            net_cts_err_gt_100.append(net_cts_sigma_up[i])
            hr2_gt_100.append(hr2[i])
            if r[i] < 25:
                hr2_gt_100_r_lt_25.append(hr2[i])
            if r[i] > 25:
                hr2_gt_100_r_gt_25.append(hr2[i])
            hr2_err_gt_100.append(hr2_err[i])
            hr5_gt_100.append(hr5[i])
            hr5_err_gt_100.append(hr5_err[i])
            chi_diff.append(chi_sqr[i] - chi_sqr_kT1[i])
        if net_cts[i] > 200:
            ra_gt_200.append(ra[i])
            dec_gt_200.append(dec[i])
            pli_gt_200.append(pli[i])
            nh1_gt_200.append(nh1[i])
            r_gt_200.append(r[i])
            net_cts_gt_200.append(net_cts[i])
            net_cts_err_gt_200.append(net_cts_sigma_up[i])
            hr2_gt_200.append(hr2[i])
            hr2_err_gt_200.append(hr2_err[i])
        if net_cts[i] < 100:
            ra_lt_100.append(ra[i])
            dec_lt_100.append(dec[i])
            pli_lt_100.append(pli[i])
            nh1_lt_100.append(nh1[i])
            r_lt_100.append(r[i])
            net_cts_lt_100.append(net_cts[i])
            net_cts_err_lt_100.append(net_cts_sigma_up[i])
            hr2_lt_100.append(hr2[i])
            hr2_err_lt_100.append(hr2_err[i])
        if net_cts[i] > 100 and net_cts[i] < 200:
            ra_gt_100_lt_200.append(ra[i])
            dec_gt_100_lt_200.append(dec[i])
            pli_gt_100_lt_200.append(pli[i])
            nh1_gt_100_lt_200.append(nh1[i])
            r_gt_100_lt_200.append(r[i])
            net_cts_gt_100_lt_200.append(net_cts[i])
            net_cts_err_gt_100_lt_200.append(net_cts_sigma_up[i])
            hr2_gt_100_lt_200.append(hr2[i])
            hr2_err_gt_100_lt_200.append(hr2_err[i])
    print('Num Sources Net Counts > 50: %d out of %d' % (len(hr2_gt_50), len(hr2)))
    print('Num Sources Net Counts > 100: %d out of %d' % (len(hr2_gt_100), len(hr2)))
    print('Num Sources Net Counts > 200: %d out of %d' % (len(hr2_gt_200), len(hr2)))
    print('Num Sources Net Counts < 100: %d out of %d' % (len(hr2_lt_100), len(hr2)))
    print('Num Sources 100 < Net Counts < 200: %d out of %d' % (len(hr2_gt_100_lt_200), len(hr2)))

    # Counting sources in various regions for number density plots
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

    for h in hr2_gt_50:
        if 0 <= r_gt_50[i] and r_gt_50[i] < 5:
            if h < 0.3:
                n5 += 1
            else:
                h5 += 1
        if 5 <= r_gt_50[i] and r_gt_50[i] < 10:
            if h < 0.3:
                n10 += 1
            else:
                h10 += 1
        if 10 <= r_gt_50[i] and r_gt_50[i] < 15:
            if h < 0.3:
                n15 += 1
            else:
                h15 += 1
        if 15 <= r_gt_50[i] and r_gt_50[i] < 25:
            if h < 0.3:
                n25 += 1
            else:
                h25 += 1
        if 25 <= r_gt_50[i] and r_gt_50[i] < 30:
            if h < 0.3:
                n30 += 1
            else:
                h30 += 1
        if 30 <= r_gt_50[i] and r_gt_50[i] < 40:
            if h < 0.3:
                n40 += 1
            else:
                h40 += 1
        if 40 <= r_gt_50[i] and r_gt_50[i] < 45:
            if h < 0.3:
                n45 += 1
            else:
                h45 += 1
        if 45 <= r_gt_50[i] and r_gt_50[i] < 50:
            if h < 0.3:
                n50 += 1
            else:
                h50 += 1
        if 50 <= r_gt_50[i] and r_gt_50[i] < 55:
            if h < 0.3:
                n55 += 1
            else:
                h55 += 1
        if 55 <= r_gt_50[i] and r_gt_50[i] < 60:
            if h < 0.3:
                n60 += 1
            else:
                h60 += 1
        if 60 <= r_gt_50[i] and r_gt_50[i] < 70:
            if h < 0.3:
                n70 += 1
            else:
                h70 += 1
        if 70 <= r_gt_50[i] and r_gt_50[i] < 80:
            if h < 0.3:
                n80 += 1
            else:
                h80 += 1
        i+=1
    soft_area_normalized.append(0)
    soft_area_normalized.append((float(n5)-3)/float(ann_arr5))
    soft_area_normalized.append(float(n10)/float(ann_arr10))
    soft_area_normalized.append(float(n15)/float(ann_arr15))
    soft_area_normalized.append(float(n25)/float(ann_arr25))
    soft_area_normalized.append(float(n30)/float(ann_arr30))
    soft_area_normalized.append(float(n40)/float(ann_arr40))
    soft_area_normalized.append(float(n45)/float(ann_arr45))
    soft_area_normalized.append(float(n50)/float(ann_arr50))
    soft_area_normalized.append(float(n55)/float(ann_arr55))
    soft_area_normalized.append(float(n60)/float(ann_arr60))
    soft_area_normalized.append(float(n70)/float(ann_arr70))
    soft_area_normalized.append(float(n80)/float(ann_arr80))
    hard_area_normalized.append(0)
    hard_area_normalized.append(float(h5)/float(ann_arr5))
    hard_area_normalized.append(float(h10)/float(ann_arr10))
    hard_area_normalized.append(float(h15)/float(ann_arr15))
    hard_area_normalized.append(float(h25)/float(ann_arr25))
    hard_area_normalized.append(float(h30)/float(ann_arr30))
    hard_area_normalized.append(float(h40)/float(ann_arr40))
    hard_area_normalized.append(float(h45)/float(ann_arr45))
    hard_area_normalized.append(float(h50)/float(ann_arr50))
    hard_area_normalized.append(float(h55)/float(ann_arr55))
    hard_area_normalized.append(float(h60)/float(ann_arr60))
    hard_area_normalized.append(float(h70)/float(ann_arr70))
    hard_area_normalized.append(float(h80)/float(ann_arr80))

    # Make hard v soft region file
    i=0
    with open('pli.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for p in pli:
            if p!=-10:
                if p>1:
                    f.write(('circle(%f,%f,1")' % (ra[i],dec[i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (ra[i],dec[i])))
                    f.write('\n')
            i+=1
    i=0
    with open('pli_gt_100.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for p in pli_gt_100:
            if p!=-10:
                if p>1:
                    f.write(('circle(%f,%f,1")' % (ra_gt_100[i],dec_gt_100[i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (ra_gt_100[i],dec_gt_100[i])))
                    f.write('\n')
            i+=1
    i=0
    with open('pli_gt_200.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for p in pli_gt_200:
            if p!=-10:
                if p>1:
                    f.write(('circle(%f,%f,1")' % (ra_gt_200[i],dec_gt_200[i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (ra_gt_200[i],dec_gt_200[i])))
                    f.write('\n')
            i+=1
    i=0
    with open('pli_lt_100.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for p in pli_lt_100:
            if p!=-10:
                if p>1:
                    f.write(('circle(%f,%f,1")' % (ra_lt_100[i],dec_lt_100[i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (ra_lt_100[i],dec_lt_100[i])))
                    f.write('\n')
            i+=1
    i=0
    with open('pli_gt_100_lt_200.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for p in pli_gt_100_lt_200:
            if p!=-10:
                if p>1:
                    f.write(('circle(%f,%f,1")' % (ra_gt_100_lt_200[i],dec_gt_100_lt_200[i])))
                    f.write('\n')
                else:
                    f.write(('circle(%f,%f,1") # color=red' % (ra_gt_100_lt_200[i],dec_gt_100_lt_200[i])))
                    f.write('\n')
            i+=1
    i=0
    with open('hr2.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for h in hr2:
            if h<=0.3:
                f.write(('circle(%f,%f,1")' % (ra[i],dec[i])))
                f.write('\n')
            else:
                f.write(('circle(%f,%f,1") # color=red' % (ra[i],dec[i])))
                f.write('\n')
                i+=1
    i=0
    with open('hr2_gt_100.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for h in hr2_gt_100:
            if h<=0.3:
                f.write(('circle(%f,%f,1")' % (ra_gt_100[i],dec_gt_100[i])))
                f.write('\n')
            else:
                f.write(('circle(%f,%f,1") # color=red' % (ra_gt_100[i],dec_gt_100[i])))
                f.write('\n')
            i+=1
    i=0
    with open('hr2_gt_200.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for h in hr2_gt_200:
            if h<=0.3:
                f.write(('circle(%f,%f,1")' % (ra_gt_200[i],dec_gt_200[i])))
                f.write('\n')
            else:
                f.write(('circle(%f,%f,1") # color=red' % (ra_gt_200[i],dec_gt_200[i])))
                f.write('\n')
            i+=1
    i=0
    with open('hr2_lt_100.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for h in hr2_lt_100:
            if h<=0.3:
                f.write(('circle(%f,%f,1")' % (ra_lt_100[i],dec_lt_100[i])))
                f.write('\n')
            else:
                f.write(('circle(%f,%f,1") # color=red' % (ra_lt_100[i],dec_lt_100[i])))
                f.write('\n')
            i+=1
    i=0
    with open('hr2_gt_100_lt_200.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for h in hr2_gt_100_lt_200:
            if h<=0.3:
                f.write(('circle(%f,%f,1")' % (ra_gt_100_lt_200[i],dec_gt_100_lt_200[i])))
                f.write('\n')
            else:
                f.write(('circle(%f,%f,1") # color=red' % (ra_gt_100_lt_200[i],dec_gt_100_lt_200[i])))
                f.write('\n')
            i+=1

    # Print soft sources
    print('Catalog Name\t\tRadius (")\tRA\t\tDEC\t\tNet Cts\t\tHR2')

    i=0
    s=0
    ks1 = []
    ks2 = []
    for h in hr2_gt_50:
        if r_gt_50[i] < 25:
            ks1.append(h)
        else:
            ks2.append(h)
        if h < 0.3:
            print(('%s\t%f\t%f\t%f\t%f\t%f')%(names_gt_50[i],r_gt_50[i],ra_gt_50[i],dec_gt_50[i],net_cts_gt_50[i],hr2_gt_50[i]))
            s+=1
        i+=1
    print s
    print stats.ks_2samp(ks1,ks2)

    # LOTS of plotting
    plt.figure()

    """
    plt.hist(src_sig,np.linspace(-10,25,51))
    plt.title('Source Significance Histogram')
    plt.xlabel('Source Significance (sigma)')
    plt.ylabel('Frequency')
    plt.axis([0,10,0,50])

    plt.hist(kT,np.linspace(0,30,16))
    plt.title('kT Histogram')
    plt.xlabel('kT')
    plt.ylabel('Frequency')

    plt.hist(chi_diff, np.linspace(-2,2,20))
    plt.title('ChiSqr Differential Histogram (Net Counts > 100)')
    plt.xlabel('ChiSqr Differential (Power Law ChiSqr - Thermal ChiSqr)')
    plt.ylabel('Frequency')

    plt.hist(hr2, np.linspace(-1,1,41))
    plt.title('HR2 Histogram for Net Counts > 100')
    plt.xlabel('HR2')
    plt.ylabel('Frequency')

    plt.hist(pli_gt_100,np.linspace(-3,3,19))
    plt.title('PLI Histogram (Net Counts > 100)')
    plt.xlabel('PLI')
    plt.ylabel('Frequency')
    #plt.axis([0,10,0,50])

    plt.hist(hr2_gt_50_r_gt_25, np.linspace(-1,1,11))
    plt.title('HR2 Histogram for Net Counts > 50, R > 25"')
    plt.xlabel('HR2')
    plt.ylabel('Frequency')
    """

    i=0
    for h in hr2_gt_100:
        if r[i] < 25:
            plt.errorbar(h,hr5_gt_100[i],xerr=hr2_err_gt_100[i], yerr=hr5_err_gt_100[i], ls='None',ecolor='g')
            plt.scatter(h, hr5_gt_100[i], color='g')
        else:
            plt.errorbar(h,hr5_gt_100[i],xerr=hr2_err_gt_100[i], yerr=hr5_err_gt_100[i], ls='None',ecolor='r')
            plt.scatter(h, hr5_gt_100[i], color='r')
        i+=1
    plt.axis([0,1,0,1.5])
    plt.title('HR5 as a function of HR2 (Net Counts > 100) (Green: R<25", Red: R>25")')
    plt.xlabel('HR2')
    plt.ylabel('HR5')

    """
    plt.errorbar(net_cts_gt_100,hr2_gt_100,xerr=net_cts_err_gt_100, yerr=hr2_err_gt_100, ls='None')
    plt.axis([90,1500,-0.4,1.2])
    plt.title('HR2 as a function of Net Counts (Net Counts > 100)')
    plt.xlabel('Net Counts')
    plt.ylabel('HR2')
    plt.xscale('log')

    plt.errorbar(r_gt_100,hr2_gt_100, yerr=hr2_err_gt_100, ls='None')
    plt.axis([0,80,-0.4,1.2])
    plt.title('HR2 as a function of Radius (Net Counts > 100)')
    plt.xlabel('Radius (")')
    plt.ylabel('HR2')

    plt.plot(r_ratios,soft_ratios)
    plt.title('Percent PLI > 1 as a function of Radius')
    plt.xlabel('Radius (")')
    plt.ylabel('Percent PLI > 1')

    plt.plot(r_ratios_eq_area,soft_ratios_eq_area)
    plt.title('Percent PLI > 1 as a function of Radius (Each Point Derived From Equal Area)')
    plt.xlabel('Radius (")')
    plt.ylabel('Percent PLI > 1')

    plt.scatter([0,5,10,15,25,30,40,45,50,55,60,70,80],soft_area_normalized)
    plt.title('# of Soft Sources Normalized by Annulus Area as a function of Radius (Net Counts > 50)')
    plt.axis([-0.75,81,-0.001,0.060])
    plt.xlabel('Outer Radius (\')')
    plt.ylabel('# of Soft Sources / Area of Annulus')

    plt.scatter([0,5,10,15,25,30,40,45,50,55,60,70,80],hard_area_normalized)
    plt.title('# of Hard Sources Normalized by Annulus Area as a function of Radius (Net Counts > 50)')
    plt.axis([-0.75,81,-0.001,0.060])
    plt.xlabel('Outer Radius (\')')
    plt.ylabel('# of Hard Sources / Area of Annulus')

    plt.scatter(r_gt_100,hr2_gt_100)
    plt.axis([0,80,-1.0,1.2])
    plt.title('HR2 as a function of Radius (Net Counts > 100)')
    plt.xlabel('Radius (")')
    plt.ylabel('HR2')

    plt.scatter(r_gt_100,pli_gt_100)
    plt.axis([0,80,-3,6])
    plt.title('PLI as a function of Radius (Net Counts > 100)')
    plt.xlabel('Radius (")')
    plt.ylabel('PLI')

    plt.scatter(net_cts,pli)
    plt.axis([0,500,-3,6])
    plt.title('PLI as a function of Net Counts')
    plt.xlabel('Net Counts')
    plt.ylabel('PLI')

    plt.scatter(nh1_gt_100,pli_gt_100)
    plt.title('PLI as a function of nH (Net Counts > 100)')
    plt.xlabel('nH')
    plt.ylabel('PLI')
    plt.errorbar(net_cts,kT,xerr=net_cts_sigma_up,yerr=kT1_erru, ls='None')
    plt.axis([0,500,0,30])
    plt.title('kT as a function of Net Counts')
    plt.xlabel('Net Counts')
    plt.ylabel('kT')

    plt.errorbar(kT,hr2,xerr=kT1_erru,yerr=hr2_err, ls='None')
    plt.axis([0,30,-1,1.2])
    plt.title('HR2 as a function of kT')
    plt.xlabel('kT')
    plt.ylabel('HR2')

    plt.errorbar(hr2_gt_100,hr5_gt_100,xerr=hr2_err_gt_100,yerr=hr5_err_gt_100, ls='None')
    plt.axis([0,1,0,1.5])
    plt.title('HR5 as a function of HR2 (Net Counts > 100)')
    plt.xlabel('HR2')
    plt.ylabel('HR5')

    plt.scatter(net_cts,src_sig)
    #plt.axis([0,500,-2,10])
    plt.title('Source Significance as a function of Net Counts')
    plt.xlabel('Net Counts')
    plt.ylabel('Source Signficance (sigma)')

    plt.scatter(pli,kT)
    plt.axis([-3,6,0,10])
    plt.title('KT as a function of PLI')
    plt.xlabel('PLI')
    plt.ylabel('kT')

    plt.scatter(chi_sqr,chi_sqr_kT1)
    plt.plot(range(0,3))
    plt.axis([0,1.5,0,1.5])
    plt.title('Chi Sqr as a function of Thermal Chi Sqr')
    plt.xlabel('Thermal Chi Sqr')
    plt.ylabel('Chi Sqr')
    """

    plt.show()

if __name__ == '__main__':
    main()
