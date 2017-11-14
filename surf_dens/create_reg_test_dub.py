import matplotlib.pyplot as plt
import numpy as np
import sys,math
from scipy import optimize,integrate
import scipy.odr as odr

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
    soft_area_normalized_a = []
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
    # Flat chi square tests
    for n in np.linspace(0.1,2,191):
        chi = 0
        j=0
        for h in hr2:
            chi+=(h-n)**2/(388*float(hr2_err[j])**2)
            j+=1
        print(('Constant HR2: %f, Chi Sqr: %f')%(n,chi))

    """
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

#C>50 soft sources
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
        if 25 <= r_gt_50[i] and r_gt_50[i] < 40:
            if h < 0.3:
                n40 += 1
            else:
                h40 += 1
        """
        if 30 <= r_gt_50[i] and r_gt_50[i] < 40:
            if h < 0.3:
                n40 += 1
            else:
                h40 += 1
        """
        if 40 <= r_gt_50[i] and r_gt_50[i] < 55:
            if h < 0.3:
                n55 += 1
            else:
                h55 += 1
        """
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
        """
        if 55 <= r_gt_50[i] and r_gt_50[i] < 70:
            if h < 0.3:
                n70 += 1
            else:
                h70 += 1
        """
        if 60 <= r_gt_50[i] and r_gt_50[i] < 70:
            if h < 0.3:
                n70 += 1
            else:
                h70 += 1
        """
        if 70 <= r_gt_50[i] and r_gt_50[i] < 80:
            if h < 0.3:
                n80 += 1
            else:
                h80 += 1
        i+=1
    #soft_area_normalized.append(0.1)
    #soft_area_normalized.append((float(n5)-3)/float(ann_arr5))
    soft_area_normalized.append(float(n10-1)/float(ann_arr10)) #-1 for 50
    soft_area_normalized.append(float(n15-3)/float(ann_arr15)) #-2 for 100nmsp, -1 100msp, -3 for 50msp, -4 for 50nmsp
    soft_area_normalized.append(float(n25-1)/float(ann_arr25)) #-3 for 100nmsp, -1for 50msp, -5 for 50nmsp
#soft_area_normalized.append(float(n30-1)/float(ann_arr30)) #-1 for 50msp, -2 for 50 nmsp
    soft_area_normalized.append(float(n40-5)/float(math.pi*(40**2-25**2))) #-4 for 50
#soft_area_normalized.append(float(n45-1)/float(ann_arr45)) #-1 for 50
#soft_area_normalized.append(float(n50)/float(ann_arr50))
    soft_area_normalized.append(float(n55-1)/float(math.pi*(55**2-40**2)))
#soft_area_normalized.append(float(n60)/float(ann_arr60)) #-1 for 50nmsp
    soft_area_normalized.append(float(n70)/float(math.pi*(70**2-55**2)))
    soft_area_normalized.append(float(n80-1)/float(ann_arr80)) #-2 for 100, -1 for 50
    hard_area_normalized.append(0.1)
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

#C>100 soft sources
    na5=na10=na15=na25=na30=na40=na45=na50=na55=na60=na70=na80=0
    ha5=ha10=ha15=ha25=ha30=ha40=ha45=ha50=ha55=ha60=ha70=ha80=0

    print hr2_gt_100
    print len(hr2_gt_100)
    print r_gt_100
    print len(r_gt_100)
    i=0
    for ha in hr2_gt_100:
        if 0 <= r_gt_100[i] and r_gt_100[i] < 5:
            if ha < 0.3:
                na5 += 1
            else:
                ha5 += 1
        if 5 <= r_gt_100[i] and r_gt_100[i] < 10:
            if ha < 0.3:
                na10 += 1
            else:
                ha10 += 1
        if 10 <= r_gt_100[i] and r_gt_100[i] < 15:
            if ha < 0.3:
                na15 += 1
            else:
                ha15 += 1
        if 15 <= r_gt_100[i] and r_gt_100[i] < 25:
            if ha < 0.3:
                na25 += 1
            else:
                ha25 += 1
        if 25 <= r_gt_100[i] and r_gt_100[i] < 30:
            if ha < 0.3:
                na30 += 1
            else:
                ha30 += 1
        if 30 <= r_gt_100[i] and r_gt_100[i] < 40:
            if ha < 0.3:
                na40 += 1
            else:
                ha40 += 1
        if 40 <= r_gt_100[i] and r_gt_100[i] < 45:
            if ha < 0.3:
                na45 += 1
            else:
                ha45 += 1
        if 45 <= r_gt_100[i] and r_gt_100[i] < 50:
            if ha < 0.3:
                na50 += 1
            else:
                ha50 += 1
        if 50 <= r_gt_100[i] and r_gt_100[i] < 55:
            if ha < 0.3:
                na55 += 1
            else:
                ha55 += 1
        if 55 <= r_gt_100[i] and r_gt_100[i] < 60:
            if ha < 0.3:
                na60 += 1
            else:
                ha60 += 1
        if 60 <= r_gt_100[i] and r_gt_100[i] < 70:
            if ha < 0.3:
                na70 += 1
            else:
                ha70 += 1
        if 70 <= r_gt_100[i] and r_gt_100[i] < 80:
            if ha < 0.3:
                na80 += 1
            else:
                ha80 += 1
        i+=1
    #soft_area_normalized.append(0.1)
    #soft_area_normalized.append((float(n5)-3)/float(ann_arr5))
    soft_area_normalized_a.append(float(na10)/float(ann_arr10)) #-1 for 50
    soft_area_normalized_a.append(float(na15-1)/float(ann_arr15)) #-2 for 100nmsp, -1 100msp, -3 for 50msp, -4 for 50nmsp
    soft_area_normalized_a.append(float(na25)/float(ann_arr25)) #-3 for 100nmsp, -1for 50msp, -5 for 50nmsp
    soft_area_normalized_a.append(float(na30)/float(ann_arr30)) #-1 for 50msp, -2 for 50 nmsp
    soft_area_normalized_a.append(float(na40)/float(ann_arr40)) #-4 for 50
    soft_area_normalized_a.append(float(na45)/float(ann_arr45)) #-1 for 50
    soft_area_normalized_a.append(float(na50)/float(ann_arr50))
    soft_area_normalized_a.append(float(na55)/float(ann_arr55))
    soft_area_normalized_a.append(float(na60)/float(ann_arr60)) #-1 for 50nmsp
    soft_area_normalized_a.append(float(na70)/float(ann_arr70))
    soft_area_normalized_a.append(float(na80-2)/float(ann_arr80)) #-2 for 100, -1 for 50

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
    for h in hr2_gt_100:
        if h < 0.3:
            print(('%s\t%f\t%f\t%f\t%f\t%f')%(names_gt_100[i],r_gt_100[i],ra_gt_100[i],dec_gt_100[i],net_cts_gt_100[i],hr2_gt_100[i]))
        i+=1

    # LOTS of plotting
    plt.figure(figsize=(3.5,2.5), dpi=300)

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
    """

#x=np.array([0.1,5,10,15,25,30,40,45,50,55,60,70,80])
#    x=np.array([5,10,15,25,30,40,45,50,55,60,70,80])
#x=np.array([2.5,7.5,12.5,20,27.5,35,42.5,47.5,52.5,57.5,65,75,])
    x=np.array([7.5,12.5,20,27.5,35,42.5,47.5,52.5,57.5,65,75,])
    x_50=np.array([7.5,12.5,20,32.5,47.5,62.5,75])
    y=np.array(soft_area_normalized)
    y_a=np.array(soft_area_normalized_a)
    #    x=np.array([0.2,0.4,0.6,1.0,1.2,1.6,1.8,2.0,2.2,2.4,2.8,3.2])
    #    y=np.array([8.48827,11.31769,11.88357,4.24413,3.08664,1.81891,0.99862,0.89350,0.40420,0.73811,0.65294,1.27324])

    """Running Scenario 1"""
    """
    def fitfunc(x, r_b, k):
        return np.piecewise(x, [x < r_b], [lambda x: k*(x**0), lambda x: k*((r_b**(0+2))*(x**(-2)))])

    popt, pcov = optimize.curve_fit(fitfunc, x, y, p0=[19, 0.1])
    xd = np.linspace(0.1, 80, 13)
    plt.loglog(xd, fitfunc(xd, *popt))
    plt.xlim(5, 90)

    r_b=popt[0]
    k=popt[1]
#g=popt[2]

    print "r_b ="
    print r_b
    print "K ="
    print k
    print "gamma ="
#print g
    residuals=y-fitfunc(xd, *popt)
    chisq=sum((residuals**2)/fitfunc(xd, *popt))
    print "chisq ="
    print chisq

    n=len(x)
    rchisq=chisq/n
    print "rX^2= "
    print rchisq

    print pcov
#    print residuals
#    print residuals**2
#    print fitfunc(xd, *popt)
    """

#########

    """Running Scenario 2"""

    """
    def fitfunc(x, r_b, k, g1, g2):
        return np.piecewise(x, [x < r_b], [lambda x: k*(x**g1), lambda x: k*((r_b**(g1-g2))*(x**(g2)))])

    popt, pcov = optimize.curve_fit(fitfunc, x, y, p0=[16, 0.016, 0, -1])
    xd = np.linspace(0.1, 100, 1000)
    plt.loglog(xd, fitfunc(xd, *popt))
    plt.minorticks_off()
    #plt.xlim(5, 90)
    label=[3,4,7,3,2,1,2,1,2,4,8]
    #label=np.array(len())
    #for label, x, y in zip(label, x, y):
#  plt.annotate(label,xy=(x,y),xytext=(5,10),textcoords='offset points')
#   x=np.array([7.5,12.5,20,27.5,35,42.5,47.5,52.5,57.5,65,75,])
#   y=np.array(soft_area_normalized)

    r_b=popt[0]
    k=popt[1]
    g1=popt[2]
    g2=popt[3]

    print "pcov ="
    print pcov

    print "r_b ="
    print r_b
    print "K ="
    print k
    print "gamma_1 ="
    print g1
    print "gamma_2 ="
    print g2

    perr=np.sqrt(abs(np.diag(pcov)))
    print "perr = %s"%perr


    #alt fitting for errors

    def fitfunc_a(x, k_a, g1_a):
        return np.piecewise(x, [x < r_b], [lambda x: k_a*(x**g1_a), lambda x: k_a*((r_b**(g1_a-g2))*(x**(g2)))])

    popt_a, pcov_a = optimize.curve_fit(fitfunc_a, x, y_alt, p0=[0.016, 0])
    xp=np.linspace(0.1,r_b,1000)
    plt.loglog(xp, fitfunc_a(xp, *popt_a),'k--')
    plt.minorticks_off()
    plt.xlim(5, 90)

#    r_b_a=popt_a[0]
    k_a=popt_a[0]
    g1_a=popt_a[1]
#    g2_a=popt_a[3]

    print "pcov ="
    print pcov_a

#    print "r_b ="
#    print r_b_a
    print "K ="
    print k_a
    print "gamma_1 ="
    print g1_a
#    print "gamma_2 ="
#    print g2_a

    perr_a=np.sqrt(abs(np.diag(pcov_a)))
    print "perr = %s"%perr_a

    #alt fitting with +2
    def fitfunc_a2(x, k_a2, g1_a2):
        return np.piecewise(x, [x < r_b], [lambda x: k_a2*(x**g1_a2), lambda x: k_a2*((r_b**(g1_a2-g2))*(x**(g2)))])

    popt_a2, pcov_a2 = optimize.curve_fit(fitfunc_a2, x, y_alt2, p0=[0.016, 0])
    xp=np.linspace(0.1,r_b,1000)
    plt.loglog(xp, fitfunc_a2(xp, *popt_a2),'k--')
    plt.minorticks_off()
    plt.xlim(5, 90)

#    r_b_a=popt_a2[0]
    k_a2=popt_a2[0]
    g1_a2=popt_a2[1]
#    g2_a=popt_a2[3]

    print "gamma_1 (+2 src) = %s"%g1_a2
    perr_a2=np.sqrt(abs(np.diag(pcov_a2)))
    print "perr = %s"%perr_a2

    #alt fitting with +5
    def fitfunc_a5(x, k_a5, g1_a5):
        return np.piecewise(x, [x < r_b], [lambda x: k_a5*(x**g1_a5), lambda x: k_a5*((r_b**(g1_a5-g2))*(x**(g2)))])

    popt_a5, pcov_a5 = optimize.curve_fit(fitfunc_a5, x, y_alt5, p0=[0.016, 0])
    xp=np.linspace(0.1,r_b,1000)
    plt.loglog(xp, fitfunc_a5(xp, *popt_a5),'k--')
    plt.minorticks_off()
    plt.xlim(5, 90)

#    r_b_a=popt_a5[0]
    k_a5=popt_a5[0]
    g1_a5=popt_a5[1]
#    g2_a=popt_a5[3]

    perr_a5=np.sqrt(abs(np.diag(pcov_a5)))
    print "perr = %s"%perr_a5
    print "gamma_1 (+5 src) = %s +/- %s"%(g1_a5,perr_a5[1])
    print r_b
    print g2

#plt.fill_between(xp,fitfunc(xp,*popt),fitfunc_a(xp,*popt_a),color="none",hatch="///",edgecolor="k",linewidth=0.0)
    plt.arrow(x[0],y[0],0,fitfunc_a(x[0],*popt_a)-y[0],length_includes_head="True",head_width=0.3,head_length=0.001,ec='k')
    plt.arrow(x[0],fitfunc_a(x[0],*popt_a),0,fitfunc_a2(x[0],*popt_a2)-fitfunc_a(x[0],*popt_a),length_includes_head="True",head_width=0.3,head_length=0.001,ec='k')
    plt.arrow(x[0],fitfunc_a2(x[0],*popt_a2),0,fitfunc_a5(x[0],*popt_a5)-fitfunc_a2(x[0],*popt_a2),length_includes_head="True",head_width=0.3,head_length=0.001,ec='k')
    #print x
    #print y[0]
    #print fitfunc_a(x[0],*popt_a)

#residuals=y-fitfunc(xd, *popt)
#   chisq=sum((residuals**2)/fitfunc(xd, *popt))
#    print "chisq ="
#    print chisq
#    r_sq=sum((residuals)**2)
#    print "r_sq ="
#    print r_sq

#    n=len(x)
#    rchisq=chisq/n
#    print "rX^2= "
#    print rchisq

#    print residuals
#    print residuals**2
#   print fitfunc(xd, *popt)
    """

#######

    """Running Scenario 3"""
    """
    def fitfunc(x, r_b, gam, B, nr_b):
        return nr_b*2**((B-gam)/10)*((x/r_b)**(-gam))*((1+(x/r_b)**10)**((gam-B)/10))

    popt, pcov = optimize.curve_fit(fitfunc, x, y, p0=[19, 2, 1, 1])
    xd = np.linspace(0.1, 100, 13)
    plt.loglog(xd, fitfunc(xd, *popt))
    plt.xlim(5, 90)

    r_b=popt[0]
    gam=popt[1]
    B=popt[2]
    nr_b=popt[3]

    print "pcov ="
    print pcov

    print "r_b ="
    print r_b
    print "gamma ="
    print gam
    print "Beta ="
    print B
    print "n(r_b) ="
    print nr_b
    residuals=y-fitfunc(xd, *popt)
    chisq=sum((residuals**2)/fitfunc(xd, *popt))
#    print "chisq ="
#    print chisq
    r_sq=sum((residuals)**2)
    print "r_sq ="
    print r_sq

    n=len(x)
    rchisq=chisq/n
    print "rX^2= "
    print rchisq

    print residuals
    print residuals**2
    print fitfunc(xd, *popt)
    """


    """Running Scenario 3a"""
    """
    def fitfunc(x, k, z, gam):
        return k*((np.sqrt(x**2+z**2))**-gam)

    popt, pcov = optimize.curve_fit(fitfunc, x, y, p0=[19, 2, 1, 1])
    xd = np.linspace(0.1, 100, 13)
    plt.loglog(xd, fitfunc(xd, *popt))
    plt.xlim(5, 90)

    r_b=popt[0]
    gam=popt[1]
    B=popt[2]
    nr_b=popt[3]

    print "pcov ="
    print pcov

    print "r_b ="
    print r_b
    print "gamma ="
    print gam
    print "Beta ="
    print B
    print "n(r_b) ="
    print nr_b
    residuals=y-fitfunc(xd, *popt)
    chisq=sum((residuals**2)/fitfunc(xd, *popt))
    #    print "chisq ="
    #    print chisq
    r_sq=sum((residuals)**2)
    print "r_sq ="
    print r_sq

    n=len(x)
    rchisq=chisq/n
    print "rX^2= "
    print rchisq

    print residuals
    print residuals**2
    print fitfunc(xd, *popt)
    """


    """Running Scenario 4"""
    """
    def integrand(t, args):
        w, p, j = args
        return math.sin(t * w)/t + p*j

    def curve(w, p, j):
        res = integrate.quad(integrand, 0.0, math.pi, args=(w, p, j,)) #why adding a variable causes an error?
        return res[0]

    vcurve = np.vectorize(curve, excluded=set([1]))

    truexdata = np.asarray([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    trueydata = vcurve(truexdata, 1.0)

    xdata = truexdata + 0.1 * np.random.randn(8)
    ydata = trueydata + 0.1 * np.random.randn(8)

    popt, pcov = optimize.curve_fit(vcurve, xdata, ydata, p0=[2.0, 2.0])

    print popt

    plt.plot(xdata, vcurve(xdata, *popt))
    """
#####
    """Complex 3D Model"""

    """
    def integrand(r,x,r_b,gamma,B,n_rb):
        return ((r*(n_rb*2**((B-gamma)/10)*((r/r_b)**-gamma)*((1+(r/r_b)**10)**(gamma-B)/10)))/(np.sqrt(r**2-x**2)))

    def curve(x, n_rb, B, gamma, r_b):
        y=(x,n_rb,B,gamma,r_b)
        res = 2*integrate.quad(integrand, (x+2.5), np.inf, args=y, epsabs=0, limit=100)
        return res[0]

    vcurve = np.vectorize(curve, excluded=set([1]))

#truexdata = np.asarray([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
#trueydata = vcurve(truexdata, 1.0)

#xdata = truexdata + 0.1 * np.random.randn(8)
#ydata = trueydata + 0.1 * np.random.randn(8)

    popt, pcov = optimize.curve_fit(vcurve, x, y, p0=[1.0, 1.0, 1.15, 1.0])

    print pcov
    print popt

    perr=np.sqrt(abs(np.diag(pcov)))
    print perr

    print "nr_b = %s +/- %s"%(popt[0],perr[0])
    print "beta = %s +/- %s"%(popt[1],perr[1])
    print "gamma = %s +/- %s"%(popt[2],perr[2])
    print "r_b = %s +/- %s"%(popt[3],perr[3])

    print x
    print y

    xd = np.linspace(0.1, 100, 130)
    plt.loglog(xd, vcurve(xd, *popt))
    plt.minorticks_off()
    """

########

    """Simple 3D Model"""
    """
    def integrand(z, x_50, k, gam):
        return k*((np.sqrt(x**2+z**2))**(-gam))

    def curve(x_50, k, gam):
        y=(x_50,k,gam)
        res = integrate.quad(integrand, -np.sqrt(125+x_50**2), np.sqrt(125+x_50**2), args=y)
        return res[0]

    vcurve = np.vectorize(curve, excluded=set([1]))

    popt, pcov = optimize.curve_fit(vcurve, x_50, y, p0=[1.0,1.5])

    perr=np.sqrt(np.diag(pcov))
    print "gamma = %s +/- %s" % (popt[1],perr[1])
    print"k = %s +/- %s" % (popt[0],perr[0])

    print popt
    print pcov

    print curve
    print vcurve

#   xd = np.linspace(0.1, 100, 5000)
    xd = np.linspace(0.1, 1000, 5000)
    plt.loglog(xd, vcurve(xd, *popt))
    #plt.axis(0.1,10,0.001)
    print "here"

    #C>100 source fitting
    def integrand_a(z_a, x, k_a, gam_a):
        return k_a*((np.sqrt(x**2+z_a**2))**(-gam_a))

    def curve_a(x, k_a, gam_a):
        y=(x,k_a,gam_a)
        res = integrate.quad(integrand, -np.sqrt(125+x**2), np.sqrt(125+x**2), args=y)
        return res[0]

    vcurve_a = np.vectorize(curve_a, excluded=set([1]))
    popt_a, pcov_a = optimize.curve_fit(vcurve_a, x, y_a, p0=[1.0,1.5])
    perr_a=np.sqrt(np.diag(pcov_a))
    print "gamma_a = %s +/- %s" % (popt_a[1],perr_a[1])
    print"k_a = %s +/- %s" % (popt_a[0],perr_a[0])

    plt.loglog(xd, vcurve_a(xd, *popt_a))
    """

    """Broken PL Simple 3D Model - not done yet"""
    """
    def integrand(z,x,g1,g2,k):
        return np.piecewise(x, [x < r_b], [lambda x: k*(np.sqrt(x**2+z**2))**(-g1), lambda x: k*((r_b**(g1-g2))*(x**(g2)))])

    def curve(x, g, k):
        y=(x,g,k)
        res = integrate.quad(integrand, -(125+x**2)**(1/2), (125+x**2)**(1/2), args=y)
        return res[0]

    vcurve = np.vectorize(curve, excluded=set([1]))

    truexdata = np.asarray([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    def trueydata(x,g,k):
        return vcurve(x, g, k)

    popt, pcov = optimize.curve_fit(vcurve, x, y, p0=[2.0, 1.0])

    print('gamma = %s' % popt[0])
    print('k = %s' % popt[1])

    print popt
    print pcov

    print curve
    print vcurve

    plt.loglog(x, vcurve(x, *popt))
    """


    """
    def integrand(r, x, r_b):
        return x*2*r+r_b

    def curve(r, x, r_b):
        res=integrate.quad(integrand, x+0.258, np.inf, [r, x, r_b])
        return res[0]

    fitfunc=np.vectorize(curve, excluded=set([1]))

    popt, pcov = optimize.curve_fit(fitfunc, x, y, p0=[19, 1])
    xd = np.linspace(0.1, 100, 13)
    plt.loglog(xd, fitfunc(xd, *popt))
    plt.xlim(5, 90)

    r_b=popt[0]
    # gamma=popt[1]
    # B=popt[2]
    #n_rb=popt[3]

    print "pcov ="
    print pcov

    print "r_b ="
    print r_b
    print "gamma ="
#print gamma
    print "Beta ="
    # print B
    print "n(r_b) ="
    #print n_rb
    residuals=y-fitfunc(xd, *popt)
    chisq=sum((residuals**2)/fitfunc(xd, *popt))
    r_sq=sum((residuals)**2)
    print "r_sq ="
    print r_sq

    n=len(x)
    rchisq=chisq/n
    print "rX^2= "
    print rchisq

    print residuals
    print residuals**2
    print fitfunc(xd, *popt)
    """

    """Individual Powerlaw"""
    """
    def fitfunc(x_50, k, gam):
        return k*(x_50**(-gam))

    popt, pcov = optimize.curve_fit(fitfunc, x_50, y, p0=[0.016, 1.5])
    xd = np.linspace(0.1, 100, 1000)
    plt.loglog(xd, fitfunc(xd, *popt),c='r')
    plt.minorticks_off()
    #plt.xlim(5, 90)
    #label=[3,4,7,3,2,1,2,1,2,4,8]
    #label=np.array(len())
    #for label, x, y in zip(label, x, y):
    #    plt.annotate(label,xy=(x,y),xytext=(5,10),textcoords='offset points')
    #x=np.array([7.5,12.5,20,27.5,35,42.5,47.5,52.5,57.5,65,75,])
    #y=np.array(soft_area_normalized)

    #r_b=popt[0]
    k=popt[0]
    gam=popt[1]
    #g2=popt[3]

    print "pcov ="
    print pcov

    #print "r_b ="
    #print r_b
    print "K ="
    print k
    print "gamma ="
    print gam
    #print "gamma_2 ="
    #print g2

    perr=np.sqrt(abs(np.diag(pcov)))
    print "perr = %s"%perr
    """
    #start from here
    #using ODR
    def fitfunc(B, x_50):
        return B[0]*(x_50**(-B[1]))

    Model=odr.Model(fitfunc)
    Data=odr.RealData(x_50,y)
    Odr=odr.ODR(Data, Model, beta0=[0.016,1.5], maxit=10000)
    Odr.set_job(fit_type=2)
    output=Odr.run()
    beta=output.beta
    betastd=output.sd_beta
    output.pprint()

    xd = np.linspace(0.1, 100, 1000)
    plt.loglog(xd,fitfunc(beta,xd),'r')

    #label=[3,4,7,5,4,6,8] #with MSPs
    #label=[3,3,3,4,4,5,8] #without MSPs
        #for label, x_50, y in zip(label, x_50, y):
        #plt.annotate(label,xy=(x_50,y),xytext=(-5,7),textcoords='offset points',color='r')
        #x_50=np.array([7.5,12.5,20,32.5,47.5,62.5,75])
        #y=np.array(soft_area_normalized)

    freqs=np.array(soft_area_normalized)
    bins=np.array([5,10,15,25,40,55,70,80])
    widths=bins[1:]-bins[:-1]
    heights=freqs.astype(np.float)
    plt.fill_between(bins.repeat(2)[1:-1],np.linspace(0.0004,0.0004,14),heights.repeat(2),edgecolor='r',facecolor='None')

#fit C>100 sources
    def fitfunc_a(x, k_a, gam_a):
        return k_a*(x**(-gam_a))

    popt_a, pcov_a = optimize.curve_fit(fitfunc_a, x, y_a, p0=[0.016, 1.5])
    xp=np.linspace(0.1,25,1000)
    plt.loglog(xp, fitfunc_a(xp, *popt_a),c='k')
    plt.minorticks_off()

    #label=[3,4,5] #with MSPs
    #label=[3,3,2] #without MSPs
    #x=np.array([7.5,12.5,20])
    #for label, x, y_a in zip(label, x, y_a):
    #plt.annotate(label,xy=(x,y_a),xytext=(-5,-12),textcoords='offset points',color='k')
    #x=np.array([7.5,12.5,20,27.5,35,42.5,47.5,52.5,57.5,65,75,])
    #y_a=np.array(soft_area_normalized_a)

    k_a=popt_a[0]
    gam_a=popt_a[1]

    print "pcov_a ="
    print pcov_a
    perr_a=np.sqrt(np.diag(pcov_a))
    print "K = %s +/- %s" % (k_a,perr_a[0])
    print "gamma_a = %s +/- %s" % (gam_a,perr_a[1])

    freqs_a=np.array(soft_area_normalized_a[0:3])
    bins_a=np.array([5,10,15,25])
    widths_a=bins[1:]-bins[:-1]
    heights_a=freqs_a.astype(np.float)
    plt.fill_between(bins_a.repeat(2)[1:-1],np.linspace(0.0004,0.0004,float(len(bins_a.repeat(2)[1:-1]))),heights_a.repeat(2),edgecolor='k',facecolor='None')




    #xbins=np.array([2.5,2.5,2.5,5,2.5,5,2.5,2.5,2.5,2.5,5,5])
    #xbins=np.array([2.5,2.5,5,2.5,5,2.5,2.5,2.5,2.5,5,5])
    #xbins2=np.array(xbins)
    #xbins2[xbins>=x] = x[xbins>=x]*.999999
    #plt.errorbar(x,y,xerr=[xbins2,xbins],linestyle="None",ecolor="r")
    #plt.errorbar(x,y_a,xerr=[xbins2,xbins],linestyle="None",ecolor="k")
    #plt.errorbar(x,y,xerr=xbins,linestyle="None",ecolor="b")
    #plt.scatter([0,5,10,15,25,30,40,45,50,55,60,70,80],soft_area_normalized)
#plt.scatter(x,y_a,c='k',marker='o')
#plt.scatter(x_50,y,c='r',marker='o',edgecolor='face')
    #plt.minorticks_off()
    #plt.title('# of Soft Sources Normalized by Annulus Projected Area as a function of Radius (Net Counts > 50)')
#plt.title('# of Soft Sources Normalized by Surface Area as a function of Radius (Net Counts > 100)')
    #plt.title('# of Soft sources (with MSPs) Normalized by Surface Area as a function of Radius')
#plt.axis([1,81,-0.01,3.5])
#plt.axis([-0.75,81,-0.001,0.060])
#plt.axis([-0.1,81,-0.001,0.060])
    plt.gca().set_xticks([0.100,0.125,0.150,0.175,0.200,0.225,0.258,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.58,5.0,7.5,10.0,12.5,15.0,17.5,20.0,22.5,25.8,50.0,75.0,100.0,125.0])
    plt.gca().set_xticklabels(["","","","","","","0.01","","","","","","","","","0.1","","","","","","","","","1","","","",""], size=7)
    plt.gca().set_yticks([0.075,0.060,0.045,0.030,0.015,0.0135,0.0120,0.0105,0.0090,0.0075,0.0060,0.0045,0.0030,0.0015,0.00135,0.00120,0.00105,0.00090,0.00075,0.00060,0.00045,0.00030,0.00015])
    plt.gca().set_yticklabels(["","","","","10.0","","","","","","","","","1.0","","","","","","","","","0.1"],size=7)
    plt.axis([5.0,81,-0.001,0.060]) #net cts>50
#plt.axis([5.0,81,0.00025,0.05]) #net cts>100
    #plt.gca().set_yticks([0.01, 0.001,])
    #plt.gca().set_yticklabels(["6.6","0.66"])
#plt.xlabel('Radius from Sgr A* (pc)') #plot 3,4
    plt.xlabel('Projected Radius from Sgr A* (pc)',size=7) #plot 2
    #plt.ylabel('# of Soft Sources / Surface Area (source/sq. arcsec)') #plot 2 arcsec
    #plt.ylabel(r'# of Soft Sources / Surface Area (pc$^{-2}$)') #plot 2
    plt.ylabel(r'Surface Density (pc$^{-2}$)',size=7) #plot 2
    #plt.ylabel(r'Soft Source Density (sources/pc$^{2}$)') #plot 3,4

    """
    plt.scatter([0,5,10,15,25,30,40,45,50,55,60,70,80],hard_area_normalized)
    plt.title('# of Hard Sources Normalized by Annulus Area as a function of Radius (Net Counts > 50)')
    plt.axis([-0.75,81,-0.001,0.060])
    plt.xlabel('Outer Radius (\")')
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
    plt.savefig('Figure_4.eps', dpi=300, format='eps', bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main()
