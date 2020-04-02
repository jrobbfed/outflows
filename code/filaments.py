from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, Column
from glob import glob
import os
import matplotlib
from scipy import stats
# from astropy.visualization.mpl_style import astropy_mpl_style_1
def main():
    matplotlib.style.use('paper')
    # gamma_random = mc_gamma(1e7, project_2d=True, select_gamma3d=[0,90])
    # np.save("mc_gamma_random_1e7", gamma_random)
    # gamma_parallel = mc_gamma(1e7, project_2d=True, select_gamma3d=[0,20])
    # np.save("mc_gamma_parallel_1e7", gamma_parallel)
    # gamma_perpendicular = mc_gamma(1e7, project_2d=True, select_gamma3d=[70,90])
    # np.save("mc_gamma_perpendicular_1e7", gamma_perpendicular)

    gamma_random = np.load("mc_gamma_random_1e7.npy")
    gamma_parallel = np.load("mc_gamma_parallel_1e7.npy")
    gamma_perpendicular = np.load("mc_gamma_perpendicular_1e7.npy")
    print(len(gamma_random))

    # gamma_3d = mc_gamma(1e7, project_2d=False, select_gamma3d=[0,90]) 
    # plt.hist(gamma_3d, bins=90, range=[0,90], histtype='step', color='k')
    # plt.hist(gamma_3d[gamma_3d <= 20], bins=90, range=[0,90], histtype='stepfilled', color='grey')
    # plt.hist(gamma_3d[gamma_3d >= 70], bins=90, range=[0,90], histtype='stepfilled', color='black')
    # plt.xlim(0,90)
    # plt.ylim(0,100000)
    # plt.gca().xaxis.set_ticks(np.arange(0, 100, 10))
    # plt.xlabel(r"3D Angle Between Vectors, $\gamma_{3D}$ (deg)")
    # plt.ylabel("Number of Angles")
    # plt.tight_layout()
    # plt.savefig("hist_gamma3d.pdf", bbox_inches='tight')

    # plt.hist(mc_gamma(1e6, project_2d=True, select_gamma3d=[0,45]), cumulative=True, bins=90, density=True, alpha=0.5)
    # plt.hist(mc_gamma(1e6, project_2d=True, select_gamma3d=[45,90]), cumulative=True, bins=90, density=True, alpha=0.5)
    # plt.hist(mc_gamma(1e6, project_2d=True, select_gamma3d=[0,90]), cumulative=True, bins=90, density=True, alpha=0.5)
    # plt.show()
    # raise
    ##########
    ###Plotting the CDF of gamma.
    max_dist = 0.05

    bins=90

    fig, ax = plt.subplots()
    hist_kwargs=dict(bins=bins, range=[0,90], linewidth=3)
    gamma_all = gamma_cdf("physics_outflows.txt", ax=ax,
     label='Full', hist_kwargs=hist_kwargs)
    # gamma_def = gamma_cdf("physics_outflows.txt", ax=ax,
    #  label=r'Definite', use_confidence='D', hist_kwargs=hist_kwargs)
    # gamma_near = gamma_cdf("physics_outflows.txt", ax=ax,
    #  max_dist=max_dist, label=r'd$_{{\rm fil}}$ $\leq$ {} pc'.format(max_dist),
    #  hist_kwargs=hist_kwargs)
    gamma_both = gamma_cdf("physics_outflows.txt", ax=ax,
     label=r'Def & d$_{{\rm fil}}$ $\leq$ {} pc'.format(max_dist), use_confidence='D', max_dist=max_dist, hist_kwargs=hist_kwargs) 

    # print("AD Test for the all sample vs. random", stats.anderson_ksamp([gamma_all, gamma_random]))
    # print("AD Test for the all sample vs. perpendicular", stats.anderson_ksamp([gamma_all, gamma_perpendicular]))
    # print("AD Test for the def sample vs. random", stats.anderson_ksamp([gamma_def, gamma_random]))
    # print("AD Test for the def sample vs. perpendicular", stats.anderson_ksamp([gamma_def, gamma_perpendicular]))
    # print("AD Test for the near sample vs. random", stats.anderson_ksamp([gamma_near, gamma_random]))
    # print("AD Test for the near sample vs. perpendicular", stats.anderson_ksamp([gamma_near, gamma_perpendicular]))
    # print("AD Test for the both sample vs. random", stats.anderson_ksamp([gamma_both, gamma_random]))
    # print("AD Test for the both sample vs. perpendicular", stats.anderson_ksamp([gamma_both, gamma_perpendicular]))

    # ax.plot([0,90], [0, 1], '-', color='gray')
    count, bins = np.histogram(gamma_random, bins=bins, density=True, range=[0,90])

    bin_centers = (bins[:-1] + bins[1:]) / 2.
    ax.plot(bin_centers, np.cumsum(count), color='k', ls='-', zorder=3, lw=1.5)
    ax.text(40, 0.7, "Random")
    # ax.hist(gamma_random, normed=True, cumulative=True, histtype='step', **hist_kwargs)

    # perp_xy = ascii.read("stephens17_perpendicular.csv")
    count, bins = np.histogram(gamma_parallel, bins=bins, density=True)
    # ax.plot(perp_xy['col1'], perp_xy['col2'], ':', color='gray')
    bin_centers = (bins[:-1] + bins[1:]) / 2.
    ax.plot(bin_centers, np.cumsum(count), color='k', ls=':', zorder=3, lw=1.5)
    # ax.text(7, 0.9, r"0-20$^\circ$")
    ax.text(7, 0.9, "Parallel")
    # ax.hist(gamma_parallel, normed=True, cumulative=True, histtype='step', **hist_kwargs)

    # para_xy = ascii.read("stephens17_parallel.csv")
    count, bins = np.histogram(gamma_perpendicular, bins=bins, density=True)
    # ax.plot(para_xy['col1'], para_xy['col2'], ':')
    bin_centers = (bins[:-1] + bins[1:]) / 2.
    ax.plot(bin_centers, np.cumsum(count), color='k', ls='--', zorder=3, lw=1.5)
    # ax.hist(gamma_perpendicular, normed=True, cumulative=True, histtype='step', **hist_kwargs)
    ax.text(66, 0.45, "Perpendicular")
    # ax.text(68, 0.45, r"70-90$^\circ$"))
    plt.xlim(-2,92)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(0, 105, 15))

    ax.set_xlabel(r"$\gamma$ (deg)")
    ax.set_ylabel("Cumulative Distribution Function")
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig("gamma_cdf_all_both.pdf", bbox_inches='tight')
    plt.clf()


    # fig, ax = plt.subplots()
    # Tbol_cut = 70
    # hist_kwargs=dict(bins=bins, range=[0,90], linewidth=3)
    # gamma_Tbol_hi = gamma_cdf("physics_outflows.txt", ax=ax,
    #  label=r'$T_{{\rm bol}}$ > {} K'.format(Tbol_cut), Tbol_cut=Tbol_cut,
    #  Tbol_comp=np.greater,
    #  hist_kwargs=hist_kwargs)
    # gamma_Tbol_low = gamma_cdf("physics_outflows.txt", ax=ax,
    #  label=r'$T_{{\rm bol}} \leq$ {} K'.format(Tbol_cut), Tbol_cut=Tbol_cut,
    #  Tbol_comp=np.less_equal,
    #  hist_kwargs=hist_kwargs) 

    # print("AD Test for the low- and hi-Tbol samples (with 70 K cutoff)", stats.anderson_ksamp([gamma_Tbol_hi, gamma_Tbol_low]))

    # print("AD Test for the Tbol_low sample vs. random", stats.anderson_ksamp([gamma_Tbol_low, gamma_random]))
    # print("AD Test for the Tbol_low sample vs. perpendicular", stats.anderson_ksamp([gamma_Tbol_low, gamma_perpendicular]))
    # print("AD Test for the Tbol_hi sample vs. random", stats.anderson_ksamp([gamma_Tbol_hi, gamma_random]))
    # print("AD Test for the Tbol_hi sample vs. perpendicular", stats.anderson_ksamp([gamma_Tbol_hi, gamma_perpendicular]))

    # count, bins = np.histogram(gamma_random, bins=bins, density=True, range=[0,90])

    # bin_centers = (bins[:-1] + bins[1:]) / 2.
    # ax.plot(bin_centers, np.cumsum(count), color='k', ls='-', zorder=3, lw=1.5)
    # ax.text(40, 0.7, "Random")
    # # ax.hist(gamma_random, normed=True, cumulative=True, histtype='step', **hist_kwargs)

    # # perp_xy = ascii.read("stephens17_perpendicular.csv")
    # count, bins = np.histogram(gamma_parallel, bins=bins, density=True)
    # # ax.plot(perp_xy['col1'], perp_xy['col2'], ':', color='gray')
    # bin_centers = (bins[:-1] + bins[1:]) / 2.
    # ax.plot(bin_centers, np.cumsum(count), color='k', ls=':', zorder=3, lw=1.5)
    # # ax.text(7, 0.9, r"0-20$^\circ$")
    # ax.text(7, 0.9, "Parallel")
    # # ax.hist(gamma_parallel, normed=True, cumulative=True, histtype='step', **hist_kwargs)

    # # para_xy = ascii.read("stephens17_parallel.csv")
    # count, bins = np.histogram(gamma_perpendicular, bins=bins, density=True)
    # # ax.plot(para_xy['col1'], para_xy['col2'], ':')
    # bin_centers = (bins[:-1] + bins[1:]) / 2.
    # ax.plot(bin_centers, np.cumsum(count), color='k', ls='--', zorder=3, lw=1.5)
    # # ax.hist(gamma_perpendicular, normed=True, cumulative=True, histtype='step', **hist_kwargs)
    # ax.text(66, 0.45, "Perpendicular")
    # # ax.text(68, 0.45, r"70-90$^\circ$"))
    # plt.xlim(-2,92)
    # start, end = ax.get_xlim()
    # ax.xaxis.set_ticks(np.arange(0, 105, 15))
    # plt.tight_layout()
    # plt.savefig("gamma_cdf_Tbol.pdf", bbox_inches='tight')
    # plt.clf()


    # fig, ax = plt.subplots()
    # Lbol_cut = 7
    # hist_kwargs=dict(bins=bins, range=[0,90], linewidth=3)
    # gamma_Lbol_hi = gamma_cdf("physics_outflows.txt", ax=ax,
    #  label=r'$L_{{\rm bol}}$ > {} $L_{{\odot}}$'.format(Lbol_cut), Lbol_cut=Lbol_cut,
    #  Lbol_comp=np.greater,
    #  hist_kwargs=hist_kwargs)
    # gamma_Lbol_low = gamma_cdf("physics_outflows.txt", ax=ax,
    #  label=r'$L_{{\rm bol}} \leq$ {} $L_{{\odot}}$'.format(Lbol_cut), Lbol_cut=Lbol_cut,
    #  Lbol_comp=np.less_equal,
    #  hist_kwargs=hist_kwargs) 

    # print("AD Test for the low- and hi-Lbol samples ()", stats.anderson_ksamp([gamma_Lbol_hi, gamma_Lbol_low]))

    # print("AD Test for the Lbol_low sample vs. random", stats.anderson_ksamp([gamma_Lbol_low, gamma_random]))
    # print("AD Test for the Lbol_low sample vs. perpendicular", stats.anderson_ksamp([gamma_Lbol_low, gamma_perpendicular]))
    # print("AD Test for the Lbol_hi sample vs. random", stats.anderson_ksamp([gamma_Lbol_hi, gamma_random]))
    # print("AD Test for the Lbol_hi sample vs. perpendicular", stats.anderson_ksamp([gamma_Lbol_hi, gamma_perpendicular]))

    # count, bins = np.histogram(gamma_random, bins=bins, density=True, range=[0,90])

    # bin_centers = (bins[:-1] + bins[1:]) / 2.
    # ax.plot(bin_centers, np.cumsum(count), color='k', ls='-', zorder=3, lw=1.5)
    # ax.text(40, 0.7, "Random")
    # # ax.hist(gamma_random, normed=True, cumulative=True, histtype='step', **hist_kwargs)

    # # perp_xy = ascii.read("stephens17_perpendicular.csv")
    # count, bins = np.histogram(gamma_parallel, bins=bins, density=True)
    # # ax.plot(perp_xy['col1'], perp_xy['col2'], ':', color='gray')
    # bin_centers = (bins[:-1] + bins[1:]) / 2.
    # ax.plot(bin_centers, np.cumsum(count), color='k', ls=':', zorder=3, lw=1.5)
    # # ax.text(7, 0.9, r"0-20$^\circ$")
    # ax.text(7, 0.9, "Parallel")
    # # ax.hist(gamma_parallel, normed=True, cumulative=True, histtype='step', **hist_kwargs)

    # # para_xy = ascii.read("stephens17_parallel.csv")
    # count, bins = np.histogram(gamma_perpendicular, bins=bins, density=True)
    # # ax.plot(para_xy['col1'], para_xy['col2'], ':')
    # bin_centers = (bins[:-1] + bins[1:]) / 2.
    # ax.plot(bin_centers, np.cumsum(count), color='k', ls='--', zorder=3, lw=1.5)
    # # ax.hist(gamma_perpendicular, normed=True, cumulative=True, histtype='step', **hist_kwargs)
    # ax.text(66, 0.45, "Perpendicular")
    # # ax.text(68, 0.45, r"70-90$^\circ$"))
    # plt.xlim(-2,92)
    # start, end = ax.get_xlim()
    # ax.xaxis.set_ticks(np.arange(0, 105, 15))
    # plt.tight_layout()
    # plt.savefig("gamma_cdf_Lbol.pdf", bbox_inches='tight')
    # plt.clf()



    # ###########
    # #plotting the histogram of gamma.
    # fig, ax = plt.subplots()
    # hist_kwargs = dict(bins=9, range=[0,90])
    # ylabel='Number'
    # gamma_cdf("physics_outflows.txt", ax=ax, stackhist=True, max_dist=max_dist, use_confidence='D',
    #  label=True, cdf=False, hist=False, hist_kwargs=hist_kwargs,
    #  ylabel=ylabel, return_gamma=False)
    # # ax.set_ylim(0,40)
    # ax.set_xlim(-5,95)
    # ax.legend(loc='upper left', fontsize=12)
    # # gamma_cdf("physics_outflows.txt", ax=ax,
    # #  label=r'Definite Outflows', use_confidence='D',
    # #   cdf=False, hist=True, hist_kwargs=hist_kwargs)
    # # gamma_cdf("physics_outflows.txt", ax=ax,
    # #  max_dist=max_dist, label=r'd $\leq$ {} pc'.format(max_dist),
    # #   cdf=False, hist=True, hist_kwargs=hist_kwargs)
    # # gamma_cdf("physics_outflows.txt", ax=ax,
    # #  label=r'both', use_confidence='D', max_dist=max_dist,
    # #   cdf=False, hist=True, hist_kwargs=hist_kwargs) 
    # plt.tight_layout()
    # plt.savefig("gamma_hist.pdf")
    # plt.clf()



    # #Writing out a table with each outflow source, the closest filament, the minimum distance to that filament in pixels, the
    # #position angle of the filament as tabulated in fil_pa.txt, and all the other filaments within 10 arcmin of the outflow source.

    # # write_table()




    # #Calculate the position angle of each filament by finding the closest point to the outflow source and doing spline interpolation
    # #on the filament and calculating the derivative of the spline interpolation.

    # l_fil_pa = []
    # t_out = ascii.read("physics_outflows.txt")
    # t_fil = ascii.read("closest_filament.txt")
    # wcs = WCS("../cubes/mask_imfit_c18o_pix_2_Tmb.fits")

    # smooth_factor = 0.1 #sprep smoothing factor as a fraction of the distance to source..
    # npix = 20 #Number of pixels on either side of the closest point to use for spline interpolation.

    # for row_fil in t_fil:
    #     source = row_fil["source"]
    #     if source != 'hops 68':
    #         continue 
    #     out_row = t_out[t_out["source"] == source][0]
    #     source_coord = SkyCoord(out_row["RAJ2000"], out_row["DEJ2000"], unit=u.deg)
    #     source_xy = source_coord.to_pixel(wcs)
    #     print(source_coord, source_xy)
        

    #     closest_fil = row_fil["filament"]
    #     fil_file =  "../filaments/filaments_fullmap/{}.txt".format(closest_fil)
    #     fil_coords = np.loadtxt(fil_file)
    #     x, y = fil_coords[:,0], fil_coords[:,1]

        
    #     if source == "hops 75":
    #         shift_closest_point = -1
    #         smooth_factor = 0.1
    #     elif source == "hops 177":
    #         pa = 0
    #         l_fil_pa.append(pa)
    #         continue
    #     elif source == "hops 355":
    #         smooth_factor = 0.01
    #     else:
    #         shift_closest_point = 0
    #         smooth_factor = 0.1

    #     pa = calc_pa(x, y, point_x=source_xy[0], point_y=source_xy[1], smooth_factor=smooth_factor,
    #         npix=npix, shift_closest_point=shift_closest_point, plot=True, plot_arcsec=True, zero_at_source=True)[0]
    #     plt.gca().set_aspect('equal')
    #     plt.xlim(30,-95)
    #     plt.savefig("/Users/jesse/repos/outflows/filaments/plots/{}_{}.pdf".format(source.strip(r' '), closest_fil),
    #         bbox_inches='tight')
    #     plt.clf()
    #     print(pa)
    #     l_fil_pa.append(pa)
    # print(l_fil_pa)
    # print(Column(l_fil_pa))
    # ascii.write([t_fil["source"], t_fil['filament'], Column(l_fil_pa)], "fil_pa.txt", names=['source', 'filament', 'pa'], overwrite=True)
    pass

    # t = write_gamma(Table())
    # print(t)


# def ad_test():
#     pass

def mc_gamma(n=1e7, smallest_angle=True, project_2d=True,
    select_gamma3d=None):
    """
    ### Sample a uniform distribution of 
    ### two random vectors and return the distribution
    ### of the angle between them (gamma). 
    Optionally project to 2D and get the projected angle distribution.
    smallest_angle constrains gamma to be between 0 and 90 degrees.
    if select_gamma3d: only retain those angles before projecting.
    From Stephens et al. 2017 and Kong et al. 2019
    """
    n = int(n)
    z1 = np.random.uniform(-1,1,n) #  
    theta1 = np.random.uniform(0,2*np.pi,n) # 
    x1 = (1. - z1**2)**0.5 * np.cos(theta1) 
    y1 = (1. - z1**2)**0.5 * np.sin(theta1) 

    # second vector 
    z2 = np.random.uniform(-1,1,n) #  
    theta2 = np.random.uniform(0,2*np.pi,n) # 
    x2 = (1. - z2**2)**0.5 * np.cos(theta2) 
    y2 = (1. - z2**2)**0.5 * np.sin(theta2)  

    gamma3d = np.arccos(x1 * x2 + y1 * y2 + z1 * z2) * 180. / np.pi 

    if project_2d:
        #Project onto y-z plane, assuming x-axis is along line of sight,
        #thus the x-coordinate of the vectors is set to 0. Since they
        #are no longer unit vectors, must divide by the their lengths inside
        #the arccos to get angle.
        r1 = (y1**2. + z1**2.)**0.5
        r2 = (y2**2. + z2**2.)**0.5
        gamma = np.arccos((y1 * y2 + z1 * z2) / ((r1) * (r2))) * 180. / np.pi
    else:
        # gamma_3D in Stephens et al. (2017) 
        gamma = gamma3d

    if smallest_angle:
        gamma = np.array([g if g <= 90. else 180. - g for g in gamma])

    if select_gamma3d:
        print(gamma3d.shape, gamma.shape)
        gamma = gamma[(gamma3d >= select_gamma3d[0]) & (gamma3d <= select_gamma3d[1])]

    return gamma



def calc_gamma(pa1, pa2, force_positive=True):
    """
    Calculate the smallest angular difference between two position angles.
    """
    if force_positive:
        if pa1 < -180:
            pa1 += 360
        elif pa1 < 0:
            pa1 += 180
        if pa2 < -180:
            pa2 += 360
        elif pa2 < 0:
            pa2 += 180

    gamma = min([abs(pa1 - pa2), 180 - abs(pa1 - pa2)])

    return gamma


def write_gamma(table, fmt='ascii', fil_table="fil_pa.txt",
 outflow_table="physics_outflows.txt", averagelobes=False, force_positive=True):

    t_fil = ascii.read(fil_table)
    t_outflow = ascii.read(outflow_table)
    g = np.array([])
    s = np.array([])

    for source in t_fil['source']: 
        pa_fil = t_fil[t_fil['source'] == source]['pa'][0] 
        #print(pa_lobes[t_outflow["source"] == source]) 
        pa_lobes = t_outflow[t_outflow["source"] == source]["pa"]
        if averagelobes:
            gamma = calc_gamma(pa_fil, np.mean(pa_lobes), force_positive=force_positive)
            g = np.append(g, gamma)
            s = np.append(s, source)
        else:
            for pa_lobe in pa_lobes:
                gamma = calc_gamma(pa_fil, pa_lobe, force_positive=force_positive)
                g = np.append(g, gamma)
                s = np.append(s, source)
        print(source, gamma) 
    try:
        table = Table.read(table, format=fmt)
    except:
        pass
    table['source'] = s
    table['gamma'] = g

    table.write("gamma.txt", format='ascii', overwrite=True)
    return table


def gamma_cdf(table, plotfile=None,
 gamma_colname='gamma', dist_colname='fil_dist', dist_unit='pc', confidence_colname='confidence',
 use_confidence=False, Tbol_cut=False, Tbol_comp=np.greater, source_table="../catalogs/hops.fits", Tbol_colname="Tbol",
 Lbol_cut=False, Lbol_comp=np.greater, Lbol_colname="Lbol",
 ax=None, cdf=True, hist=False, stackhist=False,
 max_dist=False, xlabel=r"$\gamma$ (deg)", ylabel="Cumulative Distribution Function",
 hist_kwargs=dict(bins=18, range=[0,90]),
 step_kwargs=dict(),return_gamma=True,
 label=None):
    from astropy.stats import histogram

    try:
        table = ascii.read(table)
    except:
        pass

    if stackhist:
        l_gamma = []
        l_label = []
        gamma_full = table[gamma_colname]
        label_full = "Full"
        l_gamma.append(gamma_full)
        l_label.append(label_full)
        if use_confidence:
            gamma_def = table[gamma_colname][table[confidence_colname] == use_confidence]
            l_gamma.append(gamma_def)
            label_def = 'Definite'
            l_label.append(label_def)
        if max_dist:
            gamma_near = table[gamma_colname][table[dist_colname] <= max_dist]
            l_gamma.append(gamma_near)
            label_near = r'd$_{{\rm fil}}$ $\leq$ {} pc'.format(max_dist)
            l_label.append(label_near)
        if (use_confidence is not False) & (max_dist is not False):
            gamma_both = table[gamma_colname][(table[dist_colname] <= max_dist) & (table[confidence_colname] == use_confidence)]
            l_gamma.append(gamma_both)
            label_both = 'Both'
            l_label.append(label_both)

    else:


        if max_dist:
            dist = table[dist_colname]
            table = table[dist <= max_dist]
            # gamma = gamma[dist <= max_dist]

        if use_confidence:
            confidence = table[confidence_colname]
            table = table[confidence == use_confidence]

        if Tbol_cut:
            t_hops = Table.read(source_table)
            only_hops = ['hops' in source for source in table['source']]
            Tbol = np.array([t_hops[Tbol_colname][t_hops['HOPS'] == int(source.split()[1])][0] for source in table['source'] if 'hops' in source])
            table = table[only_hops]
            print(Tbol_comp(Tbol,Tbol_cut))
            table = table[Tbol_comp(Tbol,Tbol_cut)]

        if Lbol_cut:
            t_hops = Table.read(source_table)
            only_hops = ['hops' in source for source in table['source']]
            Lbol = np.array([t_hops[Lbol_colname][t_hops['HOPS'] == int(source.split()[1])][0] for source in table['source'] if 'hops' in source])
            table = table[only_hops]
            print(Lbol_comp(Lbol,Lbol_cut))
            table = table[Lbol_comp(Lbol,Lbol_cut)]

        gamma = table[gamma_colname]
        print(len(gamma))
    if not ax:
        fig, ax = plt.subplots()

    # patches = ax.step(bins[:-1], cdf.cumsum()/max(cdf.cumsum()), label=label, **step_kwargs)
    # if step:
        # patches = ax.step(np.sort(gamma), np.sort(gamma).cumsum()/max(np.sort(gamma).cumsum()),
        #  label=label, **step_kwargs)
    if cdf:
        # bins, count = np.histogram(gamma, )
        count, bins = np.histogram(gamma, bins=hist_kwargs['bins'], range=hist_kwargs['range'], density=True)
        # count, bins, patches = ax.hist(gamma, density=True, cumulative=True, histtype='step',
        #     label=label, alpha=0, **hist_kwargs)
        # count, bins, patches = ax.hist(gamma, density=True, cumulative=True, histtype='step',
        #     label=label, alpha=1, **hist_kwargs)
        # print(patches)
        bin_centers = (bins[:-1] + bins[1:]) / 2.
        ax.plot(bin_centers, np.cumsum(count), label=label)
        # patches[0].set_xy(patches[0].get_xy()[:-1]) #Remove vertical right edge.
    elif hist:
        counts, bins, patches = ax.hist(gamma, label=label,
            **hist_kwargs)

    elif stackhist:
        counts, bins, patches = ax.hist(l_gamma, stacked=True, label=l_label, **hist_kwargs)
        print(bins)

    # patches[0].set_xy(patches[0].get_xy()[:-1]) #Remove vertical right edge.

    if label:
        ax.legend(loc='best')


    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if plotfile:
        plt.savefig(plotfile)
    if return_gamma:
        return gamma
    # return count, bins
    # pass

def calc_pa(x, y, point_x=None, point_y=None, from_north=True, smooth_factor=0.1,
    shift_closest_point=0, try_flip=True, try_npix=True, npix=20, try_sortx=True, plot=True,
    return_mindist=False, return_slope=False, force_flip = False, plot_arcsec=False, zero_at_source=False):
    """
    Calculate the position angle of a filament with pixel coordinates
    x, y. Use spline interpolation to find the tangent
    along the filament at the point closest to (point_x, point_y).
    Assumes x -> ra and y -> dec (i.e. increasing X goes W and increasing Y goes N)

    from_north toggles the astronomical convention of measuring angle E of North. 

    smooth_factor is the fraction of the distance to use for the spline
    interpolation smoothing parameter s. 0.1 works pretty well by experiment.

    shift_closest_point optionally nudges the closest filament point over by some number
    of places. For the HOPS 75 filament (91_north), this should be set to -1 or the tangent
    will be wonky.
    """


    dist = ((point_x - x)**2. + (point_y - y)**2.)**0.5
    min_dist = min(dist)
    i_closest = np.argmin(dist)

    i_closest += int(shift_closest_point)

    try:
        tck = interpolate.splrep(x,y,s=smooth_factor*min_dist)
        flipped=False

    except ValueError as ve:
        print(ve)
        print("Filament has non-unique x values, flipping the x and y coordinates to use spline interpolation.")
        try:
            x,y = y,x
    #         plt.plot(y,x)
            tck = interpolate.splrep(x,y,s=smooth_factor*min_dist)
            flipped = True
        
        except ValueError as ve:
            print("Flipping x and y still didn't help! Unflipping and only interpolating {} points on either\
            side of the minimum distance point.".format(npix))
            try:
                x,y = y,x
                flipped = False
                x = x[i_closest - npix : i_closest + npix+1]
                y = y[i_closest - npix : i_closest + npix+1]
                i_closest = len(x) // 2
                if force_flip:
                    x,y = y,x
                    flipped = True
                tck = interpolate.splrep(x,y,s=smooth_factor*min_dist)
            except ValueError as ve:
                print("This also didn't work! Last resort: sorting by x-value.")
    #             x, y, dist = x[::3], y[::3], dist[::3]
    #             i_closest = np.argmin(dist)
                i_sort = np.argsort(x)
                x = x[i_sort]
                y = y[i_sort]
                dist = ((point_x - x)**2. + (point_y - y)**2.)**0.5
                i_closest = np.argmin(dist)
                if force_flip:
                    x,y = y,x
                    flipped = True
                tck = interpolate.splrep(x,y,s=smooth_factor*min_dist)
          
    print("Smoothing parameter s:", smooth_factor*min_dist)
    # tck_s15 = interpolate.splrep(x,y,s=15)
    # print(tck)

    x0=x[i_closest]
    print("X-coordinate of closest position:", x0)
    y0 = interpolate.splev(x0,tck)
    dydx = interpolate.splev(x0,tck,der=1)


    if flipped:
        slope = 1/dydx
    else:
        slope = dydx

    if from_north:
        pa = np.arctan(slope)*180/(np.pi) + 90.
    else:
        pa = np.arctan(slope)*180/(np.pi)

    tngnt = lambda x: dydx*x + (y0-dydx*x0)

    if plot:
        if plot_arcsec:
            scalex = -2. #arcsec per pix this should change for diff datasets.
            scaley = 2. #arcsec per pix this should change for diff datasets.
        else:
            scalex = scaley = 1.

        if zero_at_source:
            shiftx, shifty = point_x, point_y 
        else:
            shiftx = shifty = 0


        if flipped:
            # shiftx, shifty = shifty, shift
            plt.plot((y[::2]-shiftx)*scalex,(x[::2]-shifty)*scaley, '.', color='k', ms=8, zorder=0)
            plt.plot((y0-shiftx)*scalex, (x0-shifty)*scaley, "o", color='tab:blue', zorder=3)
            plt.plot((tngnt(x) - shiftx)*scalex,(x-shifty)*scaley, label="tangent", color='k', ls='--')
            plt.plot((interpolate.splev(x,tck) - shiftx)*scalex, (x-shifty)*scaley, color='red', lw=1.5, zorder=1)
        else:
            plt.plot((x-shiftx)*scalex,(y-shifty)*scaley, 'k.')
            plt.plot((x0-shiftx)*scalex, (y0-shifty)*scaley, "or")
            plt.plot((x-shiftx)*scalex, (tngnt(x)-shifty)*scaley, label="tangent")
            plt.plot((x-shiftx)*scalex, (interpolate.splev(x,tck) - shifty)*scaley)

        if plot_arcsec:
            plt.xlabel("RA (arcsec)")
            plt.ylabel("DEC (arcsec)")
        else:
            plt.xlabel("X [pixel]")
            plt.ylabel("Y [pixel]")

        plt.plot((point_x-shiftx)*scalex, (point_y-shifty)*scaley, "k*", label='source')
        plt.gca().invert_xaxis()
        # plt.legend()

        # plt.plot()

    return_list = [pa]

    if return_mindist:
        return_list.append(min_dist*u.pix)

    if return_slope:
        return_list.append(slope)

    return return_list


# This code is used to write out the closest_filament table.

def write_table(out="closest_filament_new.txt", out_names=["source", "filament", "dist_min", "pa", "within_10arcmin"],
 min_length=15, wcs=WCS("../cubes/mask_imfit_c18o_pix_2_Tmb.fits"),
 fil_glob="../filaments/filaments_fullmap/*.txt", outflow_table="physics_outflows.txt",
 max_dist=300.,
 pa_table="fil_pa.txt"):

    outflow_table = ascii.read(outflow_table)
    pa_table = ascii.read(pa_table)

    files = glob(fil_glob)
    print(files)

    l_min_dist = []
    l_fil_min_dist = []
    l_source = []
    l_fil_near = []
    l_pa = []
    for outflow in outflow_table:
        if outflow["Source"] in l_source:
            continue
    #     if outflow["Source"] != 'davis 11':
    #         continue
        outflow_xy = SkyCoord(outflow["RAJ2000"], outflow["DEJ2000"], unit=u.deg).to_pixel(wcs, 0)

        min_dist = np.inf
        fil_min_dist = None
        fils_near = []
        for file in files:
            filament = np.loadtxt(file)
            fil_name = os.path.splitext(os.path.basename(file))[0]
            # print(filament.shape)
            if len(filament) < min_length:
                continue
            dist = ((outflow_xy[0] - filament[:,0])**2. + (outflow_xy[1] - filament[:,1])**2.)**0.5
            # print(dist)
            if min(dist) < min_dist:
                min_dist = min(dist)
                fil_min_dist = fil_name
            if min(dist) <= max_dist:
                fils_near.append(fil_name)

        #Save the filament which is closest to the outflow and the minimum distance from the source to that filament.
        print("The closest filament to {} is {} with a minimum distance of {:.2g} pixel".format(
            outflow["Source"], fil_min_dist, min_dist))
        l_source.append(outflow["Source"])
        l_min_dist.append(min_dist)
        l_fil_min_dist.append(fil_min_dist)
        l_fil_near.append(" ".join(fils_near))

        pa = pa_table['pa'][pa_table['source'] == outflow["Source"]][0]
        print(pa)
        l_pa.append(pa)

    ascii.write([l_source, l_fil_min_dist, l_min_dist, l_pa, l_fil_near],
                out, names=out_names, overwrite=True)
    #     print(outflow_table["DEJ2000"])

if __name__ == '__main__':
    main()