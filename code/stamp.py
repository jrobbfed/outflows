#Here will be code to make postage stamps of 
#outflows. I should include functions to:
# - Make red/blue lobe contours
# - Find and plot young stars in vicinity
# DONE - Fit gaussian to spectrum?
# - Extract (mean?/median?/sum?) CO spectrum around
# (radius?) a star/core.
# 

import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from regions import CircleSkyRegion, RectangleSkyRegion
from spectral_cube import SpectralCube
from astropy.wcs import WCS
from astropy.table import Table
import matplotlib
from matplotlib.patches import Circle,FancyArrow
import math
import matplotlib.animation as animation

from matplotlib.widgets import SpanSelector

channel_rms = 1.*u.K

def main():
    matplotlib.style.use("presentation")

    tanabe_t = Table.read("../catalogs/outflows_nro45m_new.csv")

    t = Table.read("../catalogs/hops.fits")
    coords = SkyCoord(t['RAJ2000'], t['DEJ2000'])
    cube13 = SpectralCube.read("../cubes/mask_imfit_13co_pix_2_Tmb.fits")
    cube12 = SpectralCube.read("../cubes/mask_imfit_12co_pix_2_Tmb.fits")

    # for outflow in tanabe_t:
    #     if not outflow['new']:
    #         do_fit = 1
    #         nsigma_vel = 2.5
    #         fit_radius = 15*u.arcsec
    #         plot_finder(cube12, coord=SkyCoord(outflow['RA_J2000'], outflow['DEC_J2000'], unit=u.deg),
    #                     fit_cube=cube12, fit_radius=fit_radius,
    #                     nsigma_vel=nsigma_vel, blue_vel=9.*u.km/u.s, red_vel=13.7*u.km/u.s,
    #                     fit_spectrum=do_fit,
    #                     figname="finders/tanabe_old/finder_tanabe_{:02d}_fit12co_{}sigma_fit{:g}arcsec.pdf".format(
    #                         outflow["Number"], "{:g}".format(nsigma_vel).replace(".","p"), fit_radius.value),
    #                     region_width=10*u.arcmin, region_height=10*u.arcmin,
    #                     blue_levels=np.arange(5., 50, 5), red_levels=np.arange(5., 50, 5),
    #                     show_contours=True, show_catalogs=True, show_spectrum=True, show_fit=do_fit,
    #                     show_fitcircle=False, show_vrange=True, interactive=True
    #                     )
    #         break


    t_davis_nohops = Table.read("../catalogs/davis09_jets_without_hops.txt", format='ascii.tab')
    # print(t_davis_nohops["RAJ2000"] )
    c_davis_nohops = SkyCoord(t_davis_nohops["_RAJ2000"], t_davis_nohops["_DEJ2000"], unit='deg')

    for c,davis_id in zip(c_davis_nohops, t_davis_nohops['SMZ']):
        if davis_id == 124:
        # if True:
            try:
                blue_vel = 5*u.km/u.s
                red_vel = 9.5*u.km/u.s
                dofit=1
                width=height=20*u.arcmin
                show_name=0
                # print(hops_id)
                plot_finder(cube12, coord=c,
                        fit_cube=cube12, fit_radius=15*u.arcsec, fit_spectrum=dofit,
                        nsigma_vel=2, blue_vel=blue_vel, red_vel=red_vel,
                        figname="finders/davis_nohops/finder_{}{}_{}.pdf".format("davis",davis_id,"fit12co_2sigma_fit15arcsec"),
                        region_width=width, region_height=height,
                        blue_levels=np.arange(5, 50, 5), red_levels=np.arange(5., 50, 5),
                        show_contours=True, show_catalogs=True, show_spectrum=True, show_fit=dofit,
                        show_fitcircle=False, show_vrange=True, show_outflows=True, show_name=show_name,
                        interactive=False
                        ) 
            except ValueError as ve:
                print("Region is outside of cube, moving to next source.")
            except IndexError as ie:
                print("Region is outside of cube, moving to next source.")
            else:
                pass
                # break

    # for c,hops_id in zip(coords, t['HOPS']):
    #     if hops_id == 409:
    #         try:
    #             blue_vel = 8*u.km/u.s
    #             red_vel = 14.7*u.km/u.s
    #             dofit=0
    #             width=height=5*u.arcmin
    #             show_name=0
    #             print(hops_id)
    #             plot_finder(cube12, coord=c,
    #                     fit_cube=cube12, fit_radius=15*u.arcsec, fit_spectrum=dofit,
    #                     nsigma_vel=2, blue_vel=blue_vel, red_vel=red_vel,
    #                     figname="finders/all_hops/finder_{}{}_{}.pdf".format("hops",hops_id,"fit12co_2sigma_fit15arcsec"),
    #                     region_width=width, region_height=height,
    #                     blue_levels=np.arange(5, 50, 5), red_levels=np.arange(5., 50, 5),
    #                     show_contours=True, show_catalogs=True, show_spectrum=True, show_fit=dofit,
    #                     show_fitcircle=False, show_vrange=True, show_outflows=True, show_name=show_name,
    #                     interactive=False
    #                     ) 
    #         except ValueError as ve:
    #             print("Region is outside of cube, moving to next source.")
    #         except IndexError as ie:
    #             print("Region is outside of cube, moving to next source.")
    #         else:
    #             break
            

# def plot_stamps():

# def integrate_rms():

def calc_linewings(cube, spec_method=SpectralCube.mean,
        autoguess=True, gaussian_kwargs=dict(), fit_kwargs=dict(),
        nsigma_vel=2.):

    spec = extract_spectrum(cube, method=spec_method)
    gauss = fit_gaussian(spec.spectral_axis.value, spec.value, autoguess=autoguess,
            gaussian_kwargs=gaussian_kwargs, fit_kwargs=fit_kwargs)

    spec_unit = spec.spectral_axis.unit

    blue_vel = (gauss.mean - gauss.stddev*nsigma_vel)*spec_unit
    red_vel = (gauss.mean + gauss.stddev*nsigma_vel)*spec_unit
    return blue_vel, red_vel

def plot_finder(cube,
        #extract_subcube arguments
        coord=SkyCoord("5h35m22.43s", "-5d01m14.1s"),
        region_class=RectangleSkyRegion,
        region_width=10*u.arcmin, region_height=10*u.arcmin,
        
        #Fitting arguments
        fit_spectrum=True, fit_cube="../cubes/12co_jrf_north.fits",
        fit_radius=30*u.arcsec, 

        #More Arguments passed to fitting and spectral extraction function:
        spec_method=SpectralCube.mean, autoguess=True, gaussian_kwargs=dict(), fit_kwargs=dict(),
        nsigma_vel=2., spectral_unit=u.km/u.s,

        #Contour arguments
        blue_vel=8.2*u.km/u.s, red_vel=13.8*u.km/u.s,
        blue_levels=np.arange(5., 50, 5), red_levels=np.arange(5., 50, 5),
        channel_sigma=1.*u.K, auto_sigma=True, sigma_contours=True,
        red_sigma=None, blue_sigma=None,
        blue_contour_kwargs={"colors":'blue', "linewidths":1},
        red_contour_kwargs={'colors':'red', "linewidths":1},
        annotate=False,
        # xlabel="RA [J2000]", ylabel="DEC [J2000]",

        #Arguments to choose what to show.
        show_contours=True, show_catalogs=True, show_fitcircle=True, show_outflows=True,
        show_spectrum=True, show_fit=True, show_vrange=True,

        #Catalog arguments.
        catalogs=["../catalogs/hops.fits", "../catalogs/davis09_h2jets.fits"],
        catalog_kwargs=[dict(marker="o", s=10, color='black', lw=1, zorder=3),
                        dict(marker="+", s=60, color='black', lw=1, zorder=3)],

        #Outflow arguments.
        outflow_file="../catalogs/outflows_nro45m_new.csv", 
        outflow_kwargs=dict(linewidth=2, linestyle=":",width=0,zorder=2,color='gray'),
        outflow_length=0.05, show_name=False, 

        #Spectrum plotting arguments.
        plot_spec_kwargs=dict(ls='-', color='black'), plot_fit_kwargs=dict(ls=':', color='tab:orange'),
        plot_spec_xlabel="Velocity [km/s]", plot_spec_ylabel=r"Mean T$_{MB}$",
        spec_label=r"$^{12}$CO",

        #Arguments for showing the red/blue velocity range on the spectrum.
        blue_axvspan_kwargs=dict(facecolor='blue', alpha=0.2, ymin=0, ymax=1),
        red_axvspan_kwargs=dict(facecolor='red', alpha=0.2, ymin=0, ymax=1),
        
        #Figure arguments.
        # ncols=2, nrows=1,
        fig=None, figname="finder.pdf",
        
        interactive=False):
    """
    fit_spectrum: Fit a gaussian to the spectrum collapsed within fit_radius by
        spec_method.
    nsigma_vel: The number of standard deviations away from the fit mean velocity to
        included in the red/blue contours.
    """
    if fig is None:
        fig = plt.figure(1, clear=True, figsize=(10,5))

    #make sure subplots go where they should
    nrows = 1
    
    if (show_contours or show_catalogs) and show_spectrum:
        ncols = 2
        i_cont = 1
        i_cat = 1
        i_spec = 2
    else:
        ncols = 1
        i_cont = 1
        i_cat = 1
        i_spec = 1

    
    cube = extract_subcube(cube, region_class,
            region_kwargs=dict(center=coord, width=region_width, height=region_height))

    if fit_spectrum:
        subcube_fit = extract_subcube(fit_cube, CircleSkyRegion, 
            region_kwargs={'center':coord,
                           'radius':fit_radius})
        print(coord, fit_radius)
        print(subcube_fit)
        spec = extract_spectrum(subcube_fit, spectral_unit=spectral_unit, method=spec_method)
        gauss = fit_gaussian(spec.spectral_axis.value, spec.value, autoguess=autoguess,
            gaussian_kwargs=gaussian_kwargs, fit_kwargs=fit_kwargs)

        spec_unit = spec.spectral_axis.unit

        blue_vel = (gauss.mean - gauss.stddev*nsigma_vel)*spec_unit
        red_vel = (gauss.mean + gauss.stddev*nsigma_vel)*spec_unit

    if show_contours: 
        ax_cont = fig.add_subplot(nrows, ncols, i_cont, projection=WCS(cube.header).celestial)
        ax_cont.set_aspect('equal')
        ax_cont, bluecont, redcont = plot_redblue(cube, ax=ax_cont, plot_file="plot_redblue.pdf",
            blue_vel=blue_vel, red_vel=red_vel,
            blue_levels=blue_levels, red_levels=red_levels,
            channel_sigma=channel_sigma, auto_sigma=auto_sigma, sigma_contours=sigma_contours,
            red_sigma=red_sigma, blue_sigma=blue_sigma,
            blue_contour_kwargs=blue_contour_kwargs,
            red_contour_kwargs=red_contour_kwargs,
            # xlabel="RA [J2000]", ylabel="DEC [J2000]",
            annotate=annotate)
        ax_cont.set_xlabel("RA [J2000]")
        ax_cont.set_ylabel("DEC [J2000]")

        if show_fitcircle:
            c = Circle((coord.ra.to(u.deg).value, coord.dec.to(u.deg).value),
                    fit_radius.to(u.deg).value, edgecolor='black', lw=2,
                    linestyle='solid', facecolor='none', zorder=3,
                    transform=ax_cont.get_transform('fk5'))
            ax_cont.add_patch(c)
            

    if show_catalogs:
        try:
            ax_cat = ax_cont
        except NameError as ne:
            ax_cat = fig.add_subplot(nrows, ncols, i_cont)
        try:
            z = zip(catalogs, catalog_kwargs)
        except TypeError:
            z = zip([catalogs], [catalog_kwargs])

        for catalog, kwargs in z:
            ax_cat = plot_catalog(catalog, ax=ax_cat,
                    scatter_kwargs=kwargs,
                    autoscale=False)

    
    if show_outflows:
        try:
            ax_out = ax_cont
        except NameError as ne:
            ax_out = fig.add_subplot(nrows, ncols, i_cont)

        ax_out = plot_arrows(outflow_file, ax=ax_out,
                arrow_kwargs=outflow_kwargs, length=outflow_length, show_name=show_name)


    if show_spectrum:
        ax_spec = fig.add_subplot(nrows, ncols, i_spec)
        try:
            spec.spectral_axis
        except NameError as ne:
            print("Spectrum has not yet been extracted, doing so now...") 
            print("""
            Using fit_cube and fitting parameters to extract the spectrum to show,
            if this is not what you want, you need to adjust the fit_cube, fit_radius,
            and spec_method (and maybe coord) arguments!
            """)
            subcube = extract_subcube(fit_cube, CircleSkyRegion, 
                region_kwargs={'center':coord,
                           'radius':fit_radius})
            spec = extract_spectrum(subcube, spectral_unit=spectral_unit, method=spec_method)

        ax_spec.plot(spec.spectral_axis.value, spec.value,
                label=spec_label, **plot_spec_kwargs)
        if fit_spectrum:
            ax_spec.set_xlim(gauss.mean - gauss.stddev*6, gauss.mean + gauss.stddev*6)
        else:
            ax_spec.set_xlim(spec.spectral_axis[0].value, spec.spectral_axis[-1].value)
        asp = np.diff(ax_spec.get_xlim())[0] / np.diff(ax_spec.get_ylim())[0]
        ax_spec.set_aspect(asp)
        ax_spec.set_xlabel(plot_spec_xlabel)
        ax_spec.set_ylabel(plot_spec_ylabel)
        print(blue_vel, red_vel)
        
    
    if show_fit:
        try:
            ax_fit = ax_spec
        except NameError as ne:
            ax_spec = fig.add_subplot(nrows, ncols, i_spec)
            ax_fit = ax_spec

        try:
            gauss.mean
        except NameError as ne:
            print("Spectrum has not been fitted, probably because custom velocity ranges were used")
            print("Fitting spectrum now...")
            subcube_fit = extract_subcube(fit_cube, CircleSkyRegion, 
            region_kwargs={'center':coord,
                           'radius':fit_radius})
            spec = extract_spectrum(fit_cube, method=spec_method)
            gauss = fit_gaussian(spec.spectral_axis.value, spec.value, autoguess=autoguess,
                gaussian_kwargs=gaussian_kwargs, fit_kwargs=fit_kwargs)

            spec_unit = spec.spectral_axis.unit

            blue_vel = (gauss.mean - gauss.stddev*nsigma_vel)*spec_unit
            red_vel = (gauss.mean + gauss.stddev*nsigma_vel)*spec_unit


        ax_fit.plot(spec.spectral_axis.value, gauss(spec.spectral_axis.value),
                label="Fit", **plot_fit_kwargs)
    ax_spec.legend()

    if show_vrange:
        if type(blue_vel.value) is np.ndarray:
            bluevspan = ax_spec.axvspan(blue_vel[0].value, blue_vel[1].value, **blue_axvspan_kwargs)
        else:
            bluevspan = ax_spec.axvspan(ax_spec.get_xlim()[0], blue_vel.value, **blue_axvspan_kwargs)

        if type(red_vel.value) is np.ndarray:
            redvspan = ax_spec.axvspan(red_vel[0].value, red_vel[1].value, **red_axvspan_kwargs)
        else:
            redvspan = ax_spec.axvspan(red_vel.value, ax_spec.get_xlim()[1], **red_axvspan_kwargs)
    
            
    fig.subplots_adjust(wspace=0.3)

    if interactive:
        
        # slider_ax = fig.add_axes([0.25, 0.15, 0.65, 0.03])


        # def onclick(event):
        #     print(event.xdata)
        #     if event.xdata <= gauss.mean:
        #         change_vspan = bluevspan
        #         change_cont = bluecont
        #     else:
        #         change_vspan = redvspan
        #         change_cont = redcont

        print('blue_levels', blue_levels)
        span = SpanSelector(ax_spec,
                lambda xmin, xmax: onselect_vrange(xmin, xmax, vsplit=gauss.mean,
                    bluevspan=bluevspan, redvspan=redvspan,
                    ax_cont=ax_cont, bluecont=bluecont, redcont=redcont, cube=cube,
                    spec_unit=spec_unit, blue_levels=blue_levels, red_levels=red_levels,
                    channel_sigma=channel_sigma, auto_sigma=auto_sigma, sigma_contours=sigma_contours,
                    red_sigma=red_sigma, blue_sigma=blue_sigma, blue_contour_kwargs=blue_contour_kwargs,
                    red_contour_kwargs=red_contour_kwargs),
                'horizontal', useblit=True) 

        # cid = fig.canvas.mpl_connect('button_press_event', onclick)
        

        # span = SpanSelector(ax_spec,
        #         lambda xmin, xmax: onselect_vrange(xmin, xmax, vspan=change_vspan,
        #             ax_cont=ax_cont, cont=change_cont, cube=cube),
        #         'horizontal', useblit=True) 
        # fig.canvas.draw()

        fig.show()

    # raise
    else: 
        fig.savefig(figname)
    # plt.clf()



# def onclick(event):
#     print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
#           ('double' if event.dblclick else 'single', event.button,
#            event.x, event.y, event.xdata, event.ydata))


def onselect_vrange(vmin, vmax, bluevspan=None, redvspan=None,
        vsplit=None, ax_cont=None, bluecont=None, redcont=None, cube=None,
        spec_unit=None, channel_sigma=1.*u.K, blue_levels=None, red_levels=None,
        auto_sigma=None, sigma_contours=None,
        red_sigma=None, blue_sigma=None, blue_contour_kwargs=None,
        red_contour_kwargs=None, dummy=0):
    
    print(blue_levels)
    if vmax < vsplit:
        cont = bluecont
        vspan = bluevspan
        levels = np.copy(blue_levels)
        contour_kwargs = blue_contour_kwargs
    else:
        cont = redcont
        vspan = redvspan
        levels = np.copy(red_levels)
        contour_kwargs = red_contour_kwargs
    print(levels)

    xy = vspan.get_xy()
    xy[:,0][0:2] = vmin
    xy[:,0][2:4] = vmax
    xy[:,0][4] = vmin
    vspan.set_xy(xy)


    color = contour_kwargs['colors']
    hexcolor = matplotlib.colors.to_hex(color)
    
    ncoll = len(ax_cont.collections)
    ncoll_new = 0
    
    # while len(ax_cont.collections) > 0:
    # print("All collections in axes: ", ax_cont.collections)
    # bluecont.remove()
    for coll in ax_cont.collections:
        print(type(coll))
        if type(coll) == matplotlib.collections.LineCollection:
            print("This is one contour.")
            coll.remove()
            # del(coll)

            # print(coll.get_color()[0])
            # # print("Removing all contours with color ", color)
            # if matplotlib.colors.rgb2hex(coll.get_color()[0]) == hexcolor:
            #     print("removing contour")
            #     coll.remove()
        
    for coll in ax_cont.collections:
        print(type(coll))
        # if type(coll) == matplotlib.collections.LineCollection:
        #     print("This is one contour.")
        #     print(coll.get_color())
        
            
    print("After deleting contours: ",
            [coll.get_color()[0] for coll in ax_cont.collections if type(coll) == matplotlib.collections.LineCollection])
    print("Number of collections: ", len(ax_cont.collections))
    # ax_cont.redraw_in_frame()

            # try:
            #     coll.remove()

                
                # while ncoll_new < ncoll:
                #     for coll in ax_cont.collections:
                #         if matplotlib.colors.rgb2hex(coll.get_color()[0]) == hexcolor:
                #             coll.remove()
                #             ncoll_new = len(ax_cont.collections)

                #         ncoll = len(ax_cont.collections)
            # except ValueError:
            #     pass

    print("vmin, vmax = ", vmin, vmax)
    mom0, n_channels = integrate_cube(cube, vmin*spec_unit, vmax*spec_unit, return_nchan=True)
    print(mom0.shape)


    channel_width = cube.spectral_axis[1] - cube.spectral_axis[0]
    print(channel_width)
    sigma = channel_width * np.sqrt(n_channels * channel_sigma ** 2.)

    if sigma_contours:
        # newlevels = levels*sigma.value
        newlevels = levels*sigma.value
        # levels *= sigma.value


    print(sigma, channel_width, n_channels, channel_sigma, levels)

    # print("mom0 mean = {}, stddev = {}, min = {}, max = {}".format(
    #     np.nanmean(mom0).value, np.nanstd(mom0).value, np.nanmin(mom0).value, np.nanmax(mom0).value))
    
    cont = ax_cont.contour(mom0.data, levels=[8000,100000], **contour_kwargs)
    bluecont = cont

    # plt.savefig(plot_file) 








def plot_spectrum(spec, ax=None, plot_file="plot_spectrum.pdf",
        plot_kwargs=dict()):
    if ax is None:
        ax = plt.subplot()
    ax.plot(spec.spectral_axis.value, spec.value, **plot_kwargs)
    return ax


def integrate_cube(cube, v_low, v_hi, return_nchan=False):
    """
    Returns a Projection which is the integrated intensity of cube between
    v_low and v_hi
    """
    print("integrating {} from {} to {}".format(cube, v_low, v_hi))
    cube = read_cube(cube)
    slab = cube.spectral_slab(v_low, v_hi)
    if return_nchan:
        return slab.moment0(), slab.shape[0]

    return slab.moment0()



# def plot_patches():

def rotate_vector(v, angle, anchor=[0,0]):
    """Rotate a vector `v` by the given angle in radians, relative to the anchor point."""
    x, y = v

    x = x - anchor[0]
    y = y - anchor[1]
    # Here is a compiler optimization; inplace operators are slower than
    # non-inplace operators like above. This function gets used a lot, so
    # performance is critical.

    cos_theta = math.cos(angle)
    sin_theta = math.sin(angle)

    nx = x*cos_theta - y*sin_theta
    ny = x*sin_theta + y*cos_theta

    nx = nx + anchor[0]
    ny = ny + anchor[1]
    return [nx, ny]


# def rotate_arrow(x, y, dx, dy, angle=0):
#     """
#     rotate an arrow going from x,y to x+dx,y+dy by angle in degrees.
#     """
    

def plot_arrows(catalog="../catalogs/outflows_nro45m_new.csv", wcs=None, ax=None, plot_file="plot_catalog_arrows.pdf",
        ra_colname="RA_J2000", dec_colname="DEC_J2000", name_colname="Number", length=0.1, pa_colname="PA",
        center=True, visible_only=True, show_name=False,
        arrow_kwargs=dict(), autoscale=False):
    if ax is None:
        ax = plt.subplot(projection=wcs)
    ax.autoscale(enable=autoscale)
    cat = Table.read(catalog)
    ra,dec,pa,name = cat[ra_colname], cat[dec_colname], cat[pa_colname], cat[name_colname]
    print(name)

    dx = 0
    dy = length

    for x,y,angle,n in zip(ra,dec,pa,name):
        
        dx,dy = rotate_vector([0,length], -angle*np.pi/180)

        if center:
            x,y = x-dx/2., y-dy/2.
        # print("Adding arrow")
        # coord = SkyCoord(x,y,unit='deg')
        # xpix, ypix = coord.to_pixel(wcs)
        # xlim = ax.get_xlim()
        # ylim = ax.get_ylim()

        # if xpix < xlim[1] and xpix > xlim[0] and 

        # if ax.contains_point((x,y)) or ax.contains_point((x+dx, y+dy)):
        arrow = FancyArrow(x,y,dx,dy, transform=ax.get_transform('world'), **arrow_kwargs)
        ax.add_patch(arrow)

        if show_name:
            print("showing name.")
            ax.text(x+dx/2+0.01, y+dy/2, n, transform=ax.get_transform('world'))

    plt.savefig(plot_file)
    # plt.clf()

    return ax


        



    

    
    

def plot_catalog(catalog, wcs=None, ax=None, plot_file="plot_catalog.pdf",
        ra_colname="RAJ2000", dec_colname="DEJ2000",
        scatter_kwargs=dict(), autoscale=False):
    if ax is None:
        ax = plt.subplot(projection=wcs)

    ax.autoscale(enable=autoscale)

    cat = Table.read(catalog)

    ra, dec = cat[ra_colname], cat[dec_colname]

    try:
        ax.scatter(ra, dec, transform=ax.get_transform('world'), **scatter_kwargs)
    except TypeError as te:
        print("Axes needs to be WCSAxes with world transform.")
    
    # plt.savefig(plot_file)

    return ax
        

def plot_redblue(cube, ax=None, plot_file="plot_redblue.pdf",
        blue_vel=0*u.km/u.s, red_vel=10*u.km/u.s,
        blue_levels=None, red_levels=None,
        channel_sigma=1.*u.K, auto_sigma=False, sigma_contours=False,
        red_sigma=None, blue_sigma=None,
        blue_contour_kwargs={"colors":'blue', "linewidths":1},
        red_contour_kwargs={'colors':'red', "linewidths":1},
        xlabel="RA [J2000]", ylabel="DEC [J2000]",
        annotate=False):

    cube = read_cube(cube)
    spectral_axis = cube.spectral_axis

    try:
        assert len(blue_vel) == 2
    except TypeError as te:
        print("Only one blue velocity inputted, integrating cube between the first channel at {} and {}.".format(
           cube.spectral_extrema[0], blue_vel))
        blue_vel = [cube.spectral_extrema[0], blue_vel]

    blue_mom0, n_bluechannels = integrate_cube(cube, blue_vel[0], blue_vel[1], return_nchan=True)

    try:
        assert len(red_vel) == 2
    except TypeError as te:
        print("Only one red velocity inputted, integrating cube between {} and the last channel at {}.".format(
           red_vel, cube.spectral_extrema[1]))
        red_vel = [red_vel, cube.spectral_extrema[1]]

    red_mom0, n_redchannels = integrate_cube(cube, red_vel[0], red_vel[1], return_nchan=True)

    if ax is None:
        ax = plt.subplot(projection=WCS(cube.header).celestial)


    if auto_sigma:
        channel_width = cube.spectral_axis[1] - cube.spectral_axis[0]
        blue_sigma =  channel_width * np.sqrt(n_bluechannels * channel_sigma ** 2.)
        print("blue sigma = ", blue_sigma)
        red_sigma =  channel_width * np.sqrt(n_redchannels * channel_sigma ** 2.)
        print("red sigma = ", red_sigma)

    if sigma_contours:
        red_levels = red_levels * red_sigma.value
        blue_levels = blue_levels * blue_sigma.value

    if annotate:
        ax.annotate("{} to {}".format(red_vel[0], red_vel[1]),
            (0.1, 0.05), color="red",
            xycoords="axes fraction")
        ax.annotate("{} to {}".format(blue_vel[0], blue_vel[1]),
            (0.1, 0.1), color="blue",
            xycoords="axes fraction")
        
    print("Plotting blue contours at ", blue_levels)
    bluecont = ax.contour(blue_mom0.data, levels=blue_levels,
            **blue_contour_kwargs)
    print("Plotting red contours at ", red_levels)
    redcont = ax.contour(red_mom0.data, levels=red_levels,
            **red_contour_kwargs)

    # plt.savefig(plot_file) 

    return ax, bluecont, redcont

    
     
def fit_gaussian(x, y, autoguess=False, gaussian_kwargs={"amplitude":1., "mean":0, "stddev":1.},
        fit_kwargs={}):
    if autoguess:
        gaussian_kwargs={"amplitude":y.max() - y.min(),
                "mean":x[y.argmax()],
                "stddev":abs(x[y.argmax()] - x[y.argmin()])}
        print("Guessing these parameters: ", gaussian_kwargs)
    g_init = models.Gaussian1D(**gaussian_kwargs)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y, **fit_kwargs) 
    return g

def extract_subcube(cube, region_class=CircleSkyRegion, region_kwargs={}):
    """
    Using regions package, extract a region of the shape denoted by 
    region_class. The region_class constructor will take the parameters
    in region_kwargs:
    The most commonly needed region classes are probably: 

    CircleSkyRegion, with parameters:
        center: SkyCoord
        radius: Quantity in angular units

    RectangleSkyRegion, with parameters:
        center: SkyCoord 
        width: Quantity in angular units
        height Quantity in angular units
        angle: Quantity in angular units, optional
            Rotation angle of the rectangle, measured anti-clockwise.

    For the rest of the region options, see the regions documentation at
    https://astropy-regions.readthedocs.io/en/stable/shapes.html
    """
    cube = read_cube(cube)
    region = region_class(**region_kwargs)
    cube = cube.subcube_from_regions([region])

    return cube


def extract_spectrum(cube, spectral_unit=None, method=SpectralCube.mean, axis=(1,2)):
    """
    Extract a spectrum from a spectralcube object. Choose from a number of methods
    to do this, from SpectralCube. 
    method may be SpectralCube.[mean, max, median, std, sum, mad_std]

    returns a OneDSpectrum object with units.
    """
    cube = read_cube(cube)
    if spectral_unit:
        cube = cube.with_spectral_unit(spectral_unit)
    
    return method(cube, axis=axis)

def read_cube(cube):
    if type(cube) == str:
        cube = SpectralCube.read(cube)
    elif type(cube) == SpectralCube:
        pass
    else:
        raise
    return cube


if __name__ == "__main__":
    main()
    