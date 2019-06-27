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
def main():

    #Writing out a table with each outflow source, the closest filament, the minimum distance to that filament in pixels, the
    #position angle of the filament as tabulated in fil_pa.txt, and all the other filaments within 10 arcmin of the outflow source.

    # write_table()




    #Calculate the position angle of each filament by finding the closest point to the outflow source and doing spline interpolation
    #on the filament and calculating the derivative of the spline interpolation.

    # l_fil_pa = []
    # t_out = ascii.read("physics_outflows.txt")
    # t_fil = ascii.read("closest_filament.txt")
    # wcs = WCS("../cubes/mask_imfit_c18o_pix_2_Tmb.fits")

    # smooth_factor = 0.1 #sprep smoothing factor as a fraction of the distance to source..
    # npix = 20 #Number of pixels on either side of the closest point to use for spline interpolation.

    # for row_fil in t_fil:
    #     source = row_fil["source"]
        
    #     out_row = t_out[t_out["Source"] == source][0]
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
    #         npix=npix, shift_closest_point=shift_closest_point, plot=True)[0]
    #     plt.savefig("/Users/jesse/repos/outflows/filaments/plots/{}_{}.pdf".format(source.strip(r' '), closest_fil))
    #     plt.clf()
    #     print(pa)
    #     l_fil_pa.append(pa)
    # print(l_fil_pa)
    # print(Column(l_fil_pa))
    # ascii.write([t_fil["source"], t_fil['filament'], Column(l_fil_pa)], "fil_pa.txt", names=['source', 'filament', 'pa'], overwrite=True)
    pass

    t = write_gamma(Table())
    print(t)

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


def gamma_cdf(table, plotfile="gamma_cdf.pdf",
 gamma_colname='gamma', dist_colname='',
 hist_kwargs=dict()):

    pass




def calc_pa(x, y, point_x=None, point_y=None, from_north=True, smooth_factor=0.1,
    shift_closest_point=0, try_flip=True, try_npix=True, npix=20, try_sortx=True, plot=True,
    return_mindist=False, return_slope=False, force_flip = False):
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
        if flipped:
            plt.plot(y,x, '.')
            plt.plot(y0, x0, "or")
            plt.plot(tngnt(x),x, label="tangent")
            plt.plot(interpolate.splev(x,tck), x)
        else:
            plt.plot(x,y, '.')
            plt.plot(x0, y0, "or")
            plt.plot(x,tngnt(x), label="tangent")
            plt.plot(x, interpolate.splev(x,tck))

        plt.xlabel("X [pixel]")
        plt.ylabel("Y [pixel]")

        plt.plot(point_x, point_y, "k*", label='source')

        plt.legend()

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