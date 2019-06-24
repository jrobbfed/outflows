### Adapted from code for expanding shells paper.
###
from spectral_cube import SpectralCube
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion, RectangleSkyRegion
from astropy.modeling import fitting, models
from scipy.optimize import curve_fit
#Distance to Orion (Menten et al.)
dist = 414*u.pc
#CARMA-NRO Orion Map Noise Values
sig12, sig13, sig18 = 0.86*u.K, 0.64*u.K, 0.47*u.K #From Kong et al. 2018a

path_12 = "../cubes/mask_imfit_13co_pix_2_Tmb.fits"
path_13 = "../cubes/mask_imfit_12co_pix_2_Tmb.fits"
path_13_regrid_to_12 = "../cubes/mask_imfit_13co_pix_2_Tmb_regrid12co.fits"

def main():

    from stamp import extract_subcube
    c12 = SpectralCube.read(path_12)
    c13 = SpectralCube.read(path_13)
    c13_regrid = SpectralCube.read(path_13_regrid_to_12)
    t_hops = Table.read("../../catalogs/hops.fits")
    sig12, sig13 = 0.86*u.K, 0.64*u.K

    hops_169 = t_hops[t_hops["HOPS"] == 169][0]
    coord = SkyCoord(hops_169["RAJ2000"], hops_169["DEJ2000"], unit=u.deg)
    width=height=4*u.arcmin
    vmin, vmax = -2*u.km/u.s, 4.7*u.km/u.s
    sub12 = extract_subcube(c12, region_class=RectangleSkyRegion,
                          region_kwargs=dict(center=coord, width=width, height=height))
    sub13 = extract_subcube(c13, region_class=RectangleSkyRegion,
                            region_kwargs=dict(center=coord, width=width, height=height))
    sub13_regrid = extract_subcube(c13_regrid, region_class=RectangleSkyRegion,
                            region_kwargs=dict(center=coord, width=width, height=height))
    rms12 = rms(sub12, velocity_range=[[-2,0]*u.km/u.s, [18,20]*u.km/u.s])
    print(rms12, sig12)
    rms13 = rms(sub13, velocity_range=[[0, 2]*u.km/u.s, [15, 17]*u.km/u.s])
    print(rms13, sig13)
    rms13_regrid = rms(sub13_regrid, velocity_range=[[0, 2]*u.km/u.s, [15, 17]*u.km/u.s])
    print(rms13_regrid)


def calc_physics(ra=None, dec=None, r=0.17*u.pc, dr=0.05*u.pc,
 vexp=4*u.km/u.s, v0=14*u.km/u.s, cube_12co=None, cube_13co=None,
 dist=414*u.pc, snr_cutoff=5., shell_snr_cutoff=3., shell=True,
 plot=False, linewidth_mode='fwhm',
 average_Tex=True):
    """
    shell_snr_cutoff is the sigma cutoff for extracting the shell
    voxels.
    """
    import shell_model
    #print(cube_12co.header['CDELT3'])
    from spectral_cube.lower_dimensional_structures import Projection


    pix_size = (cube_12co.header['CDELT2']*u.deg).to(u.arcsec)
    vstep = (cube_12co.header['CDELT3']*(u.m/u.s)).to(u.km/u.s)

    if shell:
        model_pars = {
            'dist':dist, # pc
            'pix_size':pix_size, # arcsec
            'vstep':vstep, # km/s
            'acen':ra.to(u.deg), # deg
            'dcen':dec.to(u.deg), # deg
            'thickness':0.0, # pc
            'fwhm':0.0, # km/s
            'beta':0.0, # spectral index
            'R':r, # pc
            'dr':dr, # pc
            'vexp':vexp, # km/s
            'depth_offset':0.0, # pc
            'vel_offset':0.0, # km/s
            'v0':v0, # km/s
            'ignore_cloud':1, #Ignore cloud.
            'method':'sample',
            'write_fits':False,
            'samples_per_voxel':27}

        #Extract 12co shell voxels using model.
        model_cube = SpectralCube.read(shell_model.ppv_model(**model_pars))
        #shell_masked, shell_mask = extract_shell(
        #    cube_file=cube_12co, model_pars=model_pars, return_mask=True)
        #Extract subcubes with same ra/dec range as shell voxel cube, but 
        #full velocity range.
        subcube_shell_12co = cube_12co.subcube(
            model_cube.longitude_extrema[1],
            model_cube.longitude_extrema[0],
            model_cube.latitude_extrema[0],
            model_cube.latitude_extrema[1])
        subcube_shell_13co = cube_13co.subcube(
            model_cube.longitude_extrema[1],
            model_cube.longitude_extrema[0],
            model_cube.latitude_extrema[0],
            model_cube.latitude_extrema[1])
        #print(subcube_shell_12co, subcube_shell_13co, model_cube)
        if plot:
            plt.figure()
            plt.subplot(131)
            plt.imshow(subcube_shell_12co.moment0().data)
            plt.subplot(132)
            plt.imshow(subcube_shell_13co.moment0().data)
            plt.subplot(133)
            plt.imshow(model_cube.moment0().data)
            plt.show()
    else:
        subcube_shell_12co = cube_12co
        subcube_shell_13co = cube_13co

    rms_12co = rms_map(cube=subcube_shell_12co)
    rms_13co = rms_map(cube=subcube_shell_13co)

    ###
    ### Use 13co if 3sigma detected, otherwise use corrected 12co if 3sigma detected.
    mask_use_13co = (subcube_shell_13co >= shell_snr_cutoff * rms_13co)
    mask_use_12co = (~mask_use_13co) & (subcube_shell_12co >= shell_snr_cutoff * rms_12co)   

    ### Excitation Temperature
    Tex = cube_Tex(subcube_shell_12co, average_first=True, plot=False,
        average=average_Tex)
    print("Tex is {}".format(Tex))
    ### Correct 12co for opacity.
    subcube_shell_12co_correct = opacity_correct(
        subcube_shell_12co, cube_thin=subcube_shell_13co,
        snr_cutoff=snr_cutoff, plot_ratio='ratio.png')
    #print(subcube_shell_12co_correct)

    subcube_shell_12co_correct = subcube_shell_12co_correct.with_mask(mask_use_12co)
    subcube_shell_13co = subcube_shell_13co.with_mask(mask_use_13co)
    
    # if plot:
    #     plt.figure()
    #     plt.subplot(121)
    #     plt.imshow(subcube_shell_12co_correct.moment0().data)
    #     plt.subplot(122)
    #     plt.imshow(subcube_shell_13co.moment0().data)
    #     plt.show()


    if shell:
        ### Extract shell voxels from opacity-corrected 12co
        shell_12co_correct = extract_shell(
            subcube_shell_12co_correct, keep_latlon=True, model_cube=model_cube)
        shell_13co = extract_shell(subcube_shell_13co,
            keep_latlon=True, model_cube=model_cube)
    else:
        shell_12co_correct = subcube_shell_12co_correct
        shell_13co = subcube_shell_13co

    # if plot:
    #     plt.figure()
    #     plt.subplot(121)
    #     plt.imshow(shell_12co_correct.moment0().data)
    #     plt.subplot(122)
    #     plt.imshow(shell_13co.moment0().data)
    #     plt.show()
    ###
    ### Use 13co if 3sigma detected, otherwise use corrected 12co if 3sigma detected.
    # mask_use_13co = (shell_13co >= shell_snr_cutoff * rms_13co)
    # mask_use_12co = (~mask_use_13co) & (subcube_shell_12co >= shell_snr_cutoff * rms_12co)


    ### Calculate column density of H2 from 13co where 13co is 3sig,
    ### otherwise use opacity-corrected 12co.

    
    shell_nH2_12co = column_density_H2(shell_12co_correct, Tex=Tex, molecule="12co")
    shell_nH2_13co = column_density_H2(shell_13co, Tex=Tex, molecule="13co")

    #print(shell_nH2_12co.shape, shell_nH2_13co.shape)
    print(np.nansum([shell_nH2_12co, shell_nH2_13co], axis=0), shell_nH2_12co.unit)
    shell_nH2 = Projection(np.nansum([shell_nH2_12co, shell_nH2_13co], axis=0),
     header=shell_nH2_12co.header, wcs=shell_nH2_12co.wcs, unit=shell_nH2_12co.unit,
     dtype=shell_nH2_12co.dtype, meta=shell_nH2_12co.meta, mask=shell_nH2_12co.mask)

    #print(shell_nH2.shape)

    # if plot:
    #     plt.figure()
    #     plt.subplot(131)
    #     plt.title("H2 from 12co")
    #     plt.imshow(shell_nH2_12co.data)
    #     plt.subplot(132)
    #     plt.title("H2 from 13co")
    #     plt.imshow(shell_nH2_13co.data)
    #     plt.subplot(133)
    #     plt.title("Total H2")
    #     plt.imshow(shell_nH2)


    ### Calculate Mass, Momentum, and Energy of Shell!
    if shell:
        shell_mass = mass(shell_nH2, distance=414*u.pc, molecule='H2',
         mass_unit=u.Msun)
        shell_momentum = momentum(shell_mass, vexp)
        shell_energy = energy(shell_mass, vexp, fwhm_to_3Dsigma=False)
        shell_luminosity = (shell_energy / (r / vexp)).to(u.erg/u.s)
    else:

        mass_map = mass(shell_nH2, distance=414*u.pc, molecule='H2',
            mass_unit=u.Msun, return_map=True)
        if linewidth_mode == "fwhm":
            vel_map = subcube_shell_13co.linewidth_fwhm()
        elif linewidth_mode == "sigma":
            vel_map = subcube_shell_13co.linewidth_sigma()
        elif linewidth_mode == "sigma3D":
            vel_map = 3.**0.5 * subcube_shell_13co.linewidth_sigma()

        # plt.figure()
        # plt.imshow(vel_map)
        # plt.show()
        #vel_average = np.nanmean(vel_map.value)
        shell_mass = u.Quantity(np.nansum(mass_map))
        shell_momentum = u.Quantity(np.nansum(mass_map * vel_map))
        #shell_momentum = vel
        shell_energy = 0.5 * u.Quantity(np.nansum(mass_map * vel_map * vel_map))
        #shell_energy = 0.5 * shell_mass
        #print(u.Quantity(0.5*shell_mass*u.Quantity(np.nanmean(vel_map))**2.).to(u.erg))
        if plot:
            from aplpy import FITSFigure
            from astropy.io import fits
            hdu = shell_nH2.hdu
            hdu.data = np.log10(shell_nH2.data)
            fig = FITSFigure(hdu.data)
            fig.show_colorscale(vmin=21.5, vmax=23.65, interpolation='none')
            fig.add_colorbar()
            # plt.imshow(np.log10(shell_nH2.value), interpolation='none',
            #     vmin=21.5, vmax=23.648)
            # plt.colorbar()
            # plt.title("log(n(H2) [cm^-2])")

            #plt.show()
            plt.savefig("../subregions/berneregion_lognH2.png")

            hdu = fits.open("../berne_Nh_2.fits")[0]
            hdu.data = np.log10(hdu.data)
            fig = FITSFigure(hdu.data)
            fig.show_colorscale(vmin=21.5, vmax=23.65, interpolation='none')
            fig.add_colorbar()

            #plt.show()
            plt.savefig("../subregions/berne_lognH2.png")

            # plt.figure()
            # plt.imshow(vel_map.to(u.km/u.s).value, interpolation='none')
            # #energy_map = (0.5*mass_map*vel_map*vel_map).to(u.erg).value
            # vels = vel_map.to(u.km/u.s).value.flatten()
            # plt.figure()
            # plt.hist(vels[vels>0], normed=True, log=True, bins=40)
            # plt.xlabel("Velocity FWHM (km/s)")
            # plt.ylabel("PDF")
            # plt.title('Velocity FWHM PDF')
            # plt.show()
        # plt.figure()
        # plt.imshow(vel_map.data)
        # plt.colorbar()
        # plt.show()

    return shell_mass, shell_momentum, shell_energy
    
#Equation 1 in Arce+ 2011, from Rohlfs and Wilson 1996
def Tex(Tpeak, thick=True):
    """
    Find the excitation temperature given peak temperature
    of an optically thick line. 
    """
    from astropy.units.core import UnitsError
    if thick:
        try:
            Tex = 5.53 / np.log(1 + (5.53)/(Tpeak+0.82))
        except UnitsError:
            Tex = 5.53*u.K / np.log(1 + (5.53*u.K)/(Tpeak+0.82*u.K))
    else:
        raise("Line must be optically thick.")
    return Tex

def cube_Tex(cube, thick=True, snr_cutoff=0,
 empty_velocity_range=[[-3.,-0.1], [19.,20.]]*u.km/u.s,
 average=True, average_first=False, plot=None):
    
    if snr_cutoff:
        rms = rms_map(cube, empty_velocity_range)
        cube = snr_mask(cube, rms, snr_cutoff)
    if average:
        if average_first:
            average_spec = average_spectrum(cube, return_std=False)
            Tpeak = average_spec.max()
        else:
            Tpeak = np.nanmean(cube.max(axis=0).filled_data)*cube.unit
    else:
        Tpeak = cube.max(axis=0)

    if plot:
        plt.figure()
        if average_first:
            vel = cube.spectral_axis.to(u.km/u.s).value
            plt.plot(vel, average_spec, label='Average Spectrum')
            plt.plot([vel[0], vel[-1]], [Tpeak.value, Tpeak.value], '--', label='Tpeak')
            plt.xlabel("velocity [km/s]")
            plt.ylabel("T [K]")
        else:
            plt.imshow(cube.max(axis=0).data, interpolation='none')
            plt.title("Tpeak Map: Average is {}".format(Tpeak))
            plt.colorbar(label="K")
        plt.savefig(plot)
    #print(Tpeak)
    return Tex(Tpeak, thick=thick)

def rms_mom0(cube, channel_rms=0.86*u.K):
        channel_width = cube.spectral_axis[1] - cube.spectral_axis[0]
        rms_mom0 =  channel_width * np.sqrt(cube.spectral_axis.size * channel_rms ** 2.)
        return rms_mom0 

def rms(cube=None, velocity_range=[[-2.,-0.], [18.,20.]]*u.km/u.s, sigma_clipped=False,
    sigma_clip_kwargs={},
    return_map=False):
    """
    Returns 2D array (or one value) of the standard deviation of a spectral cube,
    calculated only in the emission-free channels.
    """
    from astropy.stats import sigma_clipped_stats
    if velocity_range.ndim == 1:
        channel_range = [cube.closest_spectral_channel(v) for v in velocity_range]
        emissionless_channels = np.arange(channel_range[0], channel_range[1]+1)
    else:
        channel_range = [[cube.closest_spectral_channel(vpair[0]),
                      cube.closest_spectral_channel(vpair[1])]
                     for vpair in velocity_range]
        emissionless_channels = np.concatenate(
            [np.arange(c[0], c[1]+1) for c in channel_range])

    emissionless_cube = cube.filled_data[emissionless_channels,:,:]


    if return_map:
        rms_map = np.nanstd(emissionless_cube, axis=0)
        return rms_map

    elif sigma_clipped:
        mean, median, rms = sigma_clipped_stats(emissionless_cube, **sigma_clip_kwargs)
        return rms

    else:
        rms = np.nanstd(emissionless_cube)
        return rms

def regrid(cube=None, new_axis=None, smooth=True):
    """
    Regrids a SpectralCube `cube` to a new spectral axis.
    To preserve Nyquist sampling, the new spectral axis must be coarser
    than the old one, and smooth must be set to True. new_axis should have 
    units of velocity.
    See http://spectral-cube.readthedocs.io/en/latest/smoothing.html
    """
    from astropy.convolution import Gaussian1DKernel
    if smooth:
        fwhm_factor = np.sqrt(8*np.log(2)) #(FWHM/std. dev.) of gaussian
        current_resolution = cube.spectral_axis[1] - cube.spectral_axis[0]
        target_resolution = new_axis[1] - new_axis[0]
        smooth_factor = target_resolution / current_resolution

        cube = cube.spectral_smooth(Gaussian1DKernel(smooth_factor/fwhm_factor))

    cube = cube.spectral_interpolate(new_axis, 
        suppress_smooth_warning=smooth)
    return cube

def sigma_mom0(cube, channel_sigma=sig12):
    channel_width = cube.spectral_axis[1] - cube.spectral_axis[0]
    sigma_mom0 =  channel_width * np.sqrt(cube.spectral_axis.size * channel_sigma ** 2.)
    return sigma_mom0 

def average_spectrum(cube, axis=(1,2), weights=None, ignore_nan=True,
                    return_std=True):
    """
    Calculate the (weighted) average spectrum in a spectral cube
    Optionally calculate and return the (weighted) standard deviation
    spectrum. `weights` should be a cube/numpy array with the same shape as `cube`
    `axis` denotes the axis/axes to average over. For a standard SpectralCube,
    the two spatial axes are (1,2).
    Returns
    -------
    average_spectrum : 1D array_like
    std_spectrum: 1D array_like, optional
    """
    
    if weights:
        spec_mean = (cube * weights).sum(axis) / weights.sum(axis)
    else:
        spec_mean = cube.mean(axis)
#     except np.AxisError:
#         #Weight is a 2D array
#         spec_mean = weighted_cube.sum(axis) / weights.nansum()
#     except AttributeError:
#         #Weight is a constant value.
#         spec_mean = weighted_cube.sum(axis) / weights
    
    if return_std:
        if weights:
            resids = (cube - spec_mean[:, np.newaxis, np.newaxis])**2.
            spec_std = ((resids * weights).sum(axis) / weights.sum(axis)) ** 0.5
        else:
            spec_std = cube.std(axis)
        
        return spec_mean, spec_std
    else:
        return spec_mean




# def average_spectrum(cube=None, weights=1., axis=(1,2), return_std=True,
#     ignore_nan=True):
#     """
#     Calculate the (weighted) average spectrum in a spectral cube
#     Optionally calculate and return the (weighted) standard deviation
#     spectrum. `weights` should be a numpy array with the same shape as `cube`
#     `axis` denotes the axis/axes to average over. For a standard SpectralCube,
#     the two spatial axes are (1,2).
#     Returns
#     -------
#     average_spectrum : 1D array_like
#     std_spectrum: 1D array_like, optional
#     """
#     weighted_data = cube.filled_data[:,:,:] * weights
#     average_spectrum = np.nanmean(weighted_data, axis=axis)
#     #average_spectrum = np.average(cube.filled_data[:,:,:],
#     # weights=weights, axis=axis)

#     if return_std:
#         resids_squared = (cube.filled_data[:,:,:] - average_spectrum[:,np.newaxis,np.newaxis])**2. * weights
#         std_spectrum = np.sqrt(
#             np.nanmean(resids_squared, axis=axis))
#         return average_spectrum, std_spectrum
#     else:
#         return average_spectrum

# def mask_snr(return_mask=, channel_sigma=sig12):

def opacity_correct_mom0(cube_thick, cube_thin, abundance_ratio=62., return_opacity=False,
    np_func=np.nanmean):
    """
    Correct an optically thick emission-line cube for optical depth using
    an optically thin emission-line cube. This method used the ratios of the mean
    integrated intensity of the two cubes to calculate the opacity correction. To 
    find a velocity-dependent opacity correction, use opacity_correct_vel. 

    Returns the opacity correction factor to multiply cube_thick.

    cube_thick: SpectralCube or Projection
        if Projection, then treat as moment0 map
    cube_thin: SpectralCube or Projection
        if Projection, then treat as moment0 map
    """
    from stamp import read_cube
    try:
        mom0_thick = read_cube(cube_thick).moment0()
    except:
        mom0_thick = cube_thick
    try:
        mom0_thin = read_cube(cube_thin).moment0()
    except:
        mom0_thin = cube_thin
    correct_factor = 1./(
        (np_func(mom0_thick.data) / np_func(mom0_thin.data)) / abundance_ratio)
    if return_opacity:
        from scipy.optimize import minimize
        costfunc = lambda tau: np.abs(correct_factor - (tau / (1-np.exp(-tau))))
        opacity = minimize(costfunc, correct_factor).x[0]
        return (correct_factor, opacity)
    else:
        return correct_factor


def fit_parabola(x, y, x_shift=0., fit_xrange=None, autoguess=True,
                poly_kwargs={}, fix_minima=True,
                fitfunc=fitting.LinearLSQFitter,
                weights=None, ignore_nan=True,
                fit_kwargs={},
                shift_back=True):
    if autoguess:
        poly_kwargs.update(
            dict(c0 = np.interp(x_shift, x, y), c1 = 0, c2 = 1, fixed={'c1':fix_minima}))

    p_init = models.Polynomial1D(2, **poly_kwargs)
    fitter = fitfunc()
    if ignore_nan:
        i_notnan = ~np.isnan(y)
        x = x[i_notnan]
        y = y[i_notnan]
        try:
            weights = weights[i_notnan]
        except:
            pass
    p_fit = fitter(p_init, x - x_shift, y, weights=weights, **fit_kwargs)
    # print(p_init, p_fit)
    if shift_back:
        p_shiftback = p_fit.copy()
        p_shiftback.c0 = p_fit.c0 - p_fit.c1*x_shift + p_fit.c2*x_shift**2.
        p_shiftback.c1 = p_fit.c1 - 2*p_fit.c2*x_shift
        p_shiftback.c2 = p_fit.c2
        return p_shiftback
    else:   
        return p_fit

def fit_poly4(x, y, sigma=None, p0=None,
    fit_xrange=None, fix_minima=False,
    fit_kwargs=dict(), return_params=False
    ):

    if np.array(fix_minima).all():
        x_min1, x_min2 = fix_minima[0], fix_minima[1]
        def p4_2minima_fixed(x, c0, c1, c4):
            return p4_2minima(x, c0, c1, c4, x_min1, x_min2)
        popt, pcov = curve_fit(p4_2minima_fixed, x, y,
         sigma=sigma, p0=p0, **fit_kwargs)

        p4_fit = lambda l: p4_2minima(l, popt[0], popt[1], popt[2], x_min1, x_min2)

    else:
        popt, pcov = curve_fit(p4, x, y, p0=p0, sigma=sigma, **fit_kwargs)
        p4_fit = lambda l: p4(l, *popt)
    if return_params:
        return p4_fit, popt
    else:
        return p4_fit



def p4(x, c0, c1, c2, c3, c4):
    return c0 + c1*x + c2*x**2 + c3*x**3 + c4*x**4

def p4_2minima(x, c0, c1, c4, x_min1, x_min2):
    """
    This is a 4th order polynomial with two minima, at x_min1 and x_min2.
    I derived this by setting the first derivative of a general 4th-order polynomial
    equal to zero and factoring out the roots, then expressing the general 
    polynomical coefficients in terms of the x-values of the minima. The x-value of the
    local maximum is also shown.
    """
    x_max = -c1/(4*c4*x_min1*x_min2)
#     print(x_max)
    c2 = 2*c4*(x_min1*x_max + x_min2*x_max + x_min1*x_min2)
    c3 = -(4./3)*c4*(x_max + x_min1 + x_min2)
    return c0 + c1*x + c2*x**2 + c3*x**3 + c4*x**4

def cube_ratio(cubes=[None, None], rms=[None, None], return_ratiorms=True):
    """
    Returns spectral cube with ratio (uncertainties on the ratio) between
    two cubes. Uncertainity calculated with error propagation assuming
    the error on each cube is given by the rms in emission-free channels.
    """
    cube_ratio = cubes[0] / cubes[1]
    if return_ratiorms:
        cube_ratio_rms = cube_ratio * (
            (cubes[0] / rms[0]) ** -2. +
            (cubes[1] / rms[1]) ** -2.) ** 0.5
        return cube_ratio, cube_ratio_rms
    else:
        return cube_ratio

def opacity_correct(cube_thick, cube_thin=None, abundance_ratio=62., min_pix=0,
    snr_cutoff=0., rms_thick=0.86*u.K, rms_thin=0.39*u.K, calc_rms=False,
    empty_velocity_range=[[-2.,0], [18.,20.]]*u.km/u.s, vsys=7*u.km/u.s,
    weighted_average=True,
    plot_ratio=None, plot_xlim=[0,18], plot_ylim=[0,20], errorbar_kwargs=dict(marker='s', ls=''),
    plot_kwargs=dict(),
    weighted_fit=True, fixed_fit=True, return_ratio=False,
    fit_range=[6, 8]*u.km/u.s, auto_range=False, dv_range = 2*u.km/u.s,
    fit_kwargs=dict(autoguess=True, shift_back=True),
    fit=True, fit_func="parabola", return_factor=False):
    """
    Correct an optically thick emission line cube using an (assumed) optically
    thin emission line cube. The abundance ratio betweeen the cube and cube_thin
    isotopologues is given as `abundance_ratio`.
    
    `regrid_cube`: Optionally regrid `cube_thin` to the velocity grid of `cube`, preserving
    Nyquist sampling.
    
    Uses method detailed in Zhang et al 2016 (c.f. Dunham+14, Arce+2001)
    """

    if calc_rms:
        rms_thick = rms(cube_thick, empty_velocity_range)
        rms_thin = rms(cube_thin, empty_velocity_range)

    if snr_cutoff:
        cube_thick_ratiomask = cube_thick.with_mask(cube_thick > snr_cutoff*rms_thick).with_mask(cube_thin > snr_cutoff*rms_thin)
        cube_thin_ratiomask = cube_thin.with_mask(cube_thin > snr_cutoff*rms_thin).with_mask(cube_thick > snr_cutoff*rms_thick)
    else:
        cube_thick_ratiomask = cube_thick
        cube_thin_ratiomask = cube_thin
    print("thick min: ", cube_thick_ratiomask.min(), "thin min: ", cube_thin_ratiomask.min())
    ratio, ratio_rms = cube_ratio(
        cubes=[cube_thick_ratiomask, cube_thin_ratiomask],
        rms=[rms_thick, rms_thin],
        return_ratiorms=True)
    print("ratio min: ", ratio.min())

    weights = ratio_rms ** -2

    if weighted_average:
        ratiospec, ratiospec_rms = average_spectrum(ratio, weights=weights)
    else:
        ratiospec, ratiospec_rms = average_spectrum(ratio)

    if min_pix:
        count_pix = np.count_nonzero(~ratio.apply_numpy_function(np.isnan), axis=(1,2))
        print(count_pix)
        lt_min_pix = count_pix < min_pix
        print(ratiospec)
        ratiospec[lt_min_pix] = np.nan
        ratiospec_rms[lt_min_pix] = np.nan

    if fit:
        #Fit with quadratic
        if auto_range:
            vmin = ratiospec.spectral_axis[np.nanargmin(ratiospec)]
            fit_range = u.Quantity([vmin - dv_range, vmin + dv_range])
        ii_specfit = (ratiospec.spectral_axis > fit_range[0]) & (ratiospec.spectral_axis < fit_range[1])

        if fit_func == "parabola":
            if weighted_fit:
                weights = (ratiospec_rms[ii_specfit])**-1
            else:
                weights = None
            p_fit = fit_parabola(ratiospec.spectral_axis[ii_specfit], ratiospec[ii_specfit],
                x_shift=vsys, weights=weights, fix_minima=fixed_fit)
            ratiospec_fit_clip = p_fit(cube_thick.spectral_axis).clip(0, abundance_ratio)

        elif fit_func == "poly4":
            print(vsys.to(ratiospec.spectral_axis.unit).value)
            if weighted_fit:
                sigma = ratiospec_rms[ii_specfit].value
            else:
                sigma = None
            if fixed_fit:
                fix_minima = vsys.to(ratiospec.spectral_axis.unit).value
            else:
                fix_minima = None
            p_fit = fit_poly4(ratiospec.spectral_axis[ii_specfit].value, ratiospec[ii_specfit].value,
                sigma=sigma,
                fix_minima=fix_minima)
            ratiospec_fit_clip = p_fit(cube_thick.spectral_axis.value).clip(0, abundance_ratio)
        # vel = ratio.spectral_axis.to(u.km/u.s).value
        #print(vel[notnan], average_ratio[notnan])

        factor = abundance_ratio * (ratiospec_fit_clip[:,np.newaxis,np.newaxis])**-1.
        cube_correct = cube_thick * factor #X_12,13 / (T_12/T_13)

    if plot_ratio:
        ax = plt.subplot()
        ax.errorbar(ratiospec.spectral_axis.value/1000, ratiospec.value, yerr=ratiospec_rms.value,
         **errorbar_kwargs)
        # plt.plot(np.linspace(-1000,1000)+7000, p_fit2(np.linspace(-1000,1000)*u.m/u.s))
        # ax.axhline(abundance_ratio)
        ax.set_xlim(plot_xlim)
        ax.set_ylim(plot_ylim)

        if fit:
            if fit_func == "parabola":
                ax.plot(ratiospec.spectral_axis.value/1000, p_fit(ratiospec.spectral_axis).clip(0,abundance_ratio), ':',
                        ratiospec.spectral_axis.value[ii_specfit]/1000, p_fit(ratiospec.spectral_axis[ii_specfit]).clip(0,abundance_ratio), '-',
                    **plot_kwargs)
            elif fit_func == "poly4":
                ax.plot(ratiospec.spectral_axis.value/1000, p_fit(ratiospec.spectral_axis.value).clip(0,abundance_ratio), ':',
                        ratiospec.spectral_axis.value[ii_specfit]/1000, p_fit(ratiospec.spectral_axis[ii_specfit].value).clip(0,abundance_ratio), '-',
                    **plot_kwargs)

        ax.set_ylabel(r"T$_{12}$ / T$_{13}$", size=16)
        ax.set_xlabel("Velocity [km / s]", size=16)
        plt.tight_layout()
        try:
            plt.savefig(plot_ratio)
        except:
            pass
    if return_factor:
        return cube_correct, factor

    elif return_ratio:
        return ratiospec, ratiospec_rms
    else:
        return cube_correct

def column_density_H2(cube, Tex,
    molecule="12co", transition="1-0",
    opacity_correction=True, Qrot_order=100,
    beam_filling_factor=1., moment0=True):
    """
    Calculate the column density of a molecule from an 
    optically thin emission-line cube.

    For optically thin lines only! Must correct 
    optically thick line for opacity FIRST.
    
    Can return the column density per channel (dn/dv) or 
    the total column density with a moment0 of dn/dv.
    """
    import astropy.constants as const
    import astropy.units as u
    if molecule == "12co":
        if transition == "1-0":
            #From Zhang et al. 2016
            nu_ul = 115.271 * u.GHz
            A_ul = 7.203E-8 * (1/u.s) 
            g_u = 3 # 2*J_u + 1
            E_u_k = 5.53 * u.K
            B0_k = 2.765 * u.K
            X_factor = 1e-4 # 12CO/H2
    if molecule == "13co":
        if transition == "1-0":
            #From Zhang et al. 2016
            nu_ul = 110.201 * u.GHz
            A_ul = 6.294E-8 * (1/u.s) 
            g_u = 3 # 2*J_u + 1
            E_u_k = 5.29 * u.K
            B0_k = 2.645 * u.K
            X_factor = 1e-4 / 62. # 13CO/H2 = 1.61e-6
            #X_factor = 1.43e-6 #From Berné 2014

    factor = (
        (8*np.pi*const.k_B*nu_ul**2.)\
        / (const.h*const.c**3.*A_ul*g_u)
             )
    ### APPLY THE ABUNDANCE RATIO BETWEEN CO/H2. CHECK THIS
    factor = factor * Qrot_partial(Tex, B0_k, N=Qrot_order)\
     * np.exp(E_u_k/Tex) / beam_filling_factor / X_factor

    # print(factor)
    if moment0:
        return (cube.moment0() * factor).to(u.cm**-2)
    else:
        return (cube * factor).to(u.cm**-2 / (u.m/u.s))

def Qrot_partial(Tex, B0_k=2.765*u.K, N=100):
    """
    Calculate Partial sum of partition function
    at a given excitation temperature. 
    B_0/k depends on the transition:
        12CO (1-0): B_0/k = 2.765 K
    """
    #Tex = u.Quantity(Tex, u.K)
    Qrot = np.zeros(Tex.shape)
    for J in range(N+1):
        Qrot = Qrot + (2*J + 1) * np.exp(
            -B0_k*J*(J+1)/Tex)
    
    return Qrot


def dmdv(cube, molecule='12co', Tex=30*u.K, distance=414*u.pc,
    spectral_unit=u.km/u.s, mass_unit=u.Msun,
    return_cube=False):
    """
    Compute mass spectrum (sum of mass per channel)
    of a given, optionally masked, SpectralCube. For now, 
    assume optically thin, but should implement an optional
    opacity correction.
    """
    cube = cube.with_spectral_unit(spectral_unit)
    nH2 = column_density_H2(cube, Tex=Tex, moment0=False,
        molecule=molecule)
    mass_cube = mass(nH2, distance=distance, molecule='H2',
        return_map=True, mass_unit=mass_unit).to(mass_unit/spectral_unit)
    dmdv = mass_cube.sum((1,2))
    if return_cube:
        return dmdv, mass_cube
    else:
        return dmdv


def ambient_correct(dmdv, n_models=1, autoguess=True, 
    find_peaks_kwargs=dict(height=0.1, width=2, rel_height=0.5),
    gaussian_kwargs=dict(), fit_func=fitting.LevMarLSQFitter(),
    plot=False, return_fit=False):
    from stamp import fit_gaussian
    from spectral_cube import OneDSpectrum

    x, y = dmdv.spectral_axis, dmdv

    dmdv_fit = fit_gaussian(x, y, n_models=n_models, autoguess=autoguess,
                            gaussian_kwargs=gaussian_kwargs, find_peaks_kwargs=find_peaks_kwargs,
                            fit_func=fit_func)

    try:
        dmdv_out = (y.value - dmdv_fit(x.value))*dmdv.unit
        # print(dmdv_out.spectral_axis)
    except u.UnitsError:
        dmdv_out = y - dmdv_fit(x)
    dmdv_out = dmdv_out.clip(0., np.inf) 

    dmdv_out = OneDSpectrum(dmdv_out, wcs=dmdv.wcs, spectral_unit=dmdv.spectral_axis.unit)


    if plot:
        plt.plot(x, y, label="Total")
        try:
            plt.plot(np.linspace(x[0], x[-1], 1000), dmdv_fit(np.linspace(x[0], x[-1], 1000)), label="Cloud Fit")
        except u.UnitsError:
            plt.plot(np.linspace(x[0], x[-1], 1000).value, dmdv_fit(np.linspace(x[0], x[-1], 1000).value), label="Cloud Fit")
        plt.plot(x, dmdv_out, label="Outflow")
        plt.semilogy()
        plt.ylim(0.001, 10)
        plt.xlim(-2, 20)
        plt.xlabel("v [km/s]")
        plt.ylabel(r"dM/dv [M$_\odot$ (km/s)$^{-1}$]")

    if return_fit:
        return dmdv_out, dmdv_fit
    else:
        return dmdv_out


def sum_dmdv(dmdv, vsys=None, vrange='auto', dv=1*u.km/u.s,
    return_lobes=False):
    """
    Integrate an input mass spectrum. Integrate to within +/- dv of vsys,
    or if there are multiple components, within +dv of the max vsys and -dv of the min vsys.

    Alternative, can set a specific velocity range to integrate over.

    dmdv must be astro.units.quantity.Quantity object with units of mass/velocity.
    """

    try:
        dmdv.unit.to(u.Msun/(u.km/u.s))
    except u.UnitConversionError:
        print("dmdv must be in units of mass/velocity.")
        raise

    channel_width = dmdv.spectral_axis[1] - dmdv.spectral_axis[0]

    if vrange == 'auto':
        if vsys is None:
            print("Must set vsys in order to use auto vrange.")
            raise
        elif np.size(vsys) > 1:
            i_blue = dmdv.spectral_axis < (min(vsys) - dv)
            i_red = dmdv.spectral_axis > (max(vsys) + dv)
        else:
            i_blue = dmdv.spectral_axis < vsys - dv
            i_red = dmdv.spectral_axis > vsys + dv

        mass_blue = dmdv[i_blue].nansum()*channel_width
        mass_red = dmdv[i_red].nansum()*channel_width
        mass = mass_blue + mass_red

    elif np.ndim(vrange) > 2:
        mass = 0
        for vra in vrange:
            i_vra = (dmdv.spectral_axis >= min(vra)) & (dmdv.spectral_axis <= max(vra))
            mass += dmdv[i_vra].nansum()*channel_width
    else:
        i_vra = (dmdv.spectral_axis >= vrange[0]) & (dmdv.spectral_axis <= vrange[1])
        mass = dmdv[i_vra].nansum()*channel_width

    if return_lobes:
        return mass_blue, mass_red
    else:
        return mass

def momentum_dmdv(dmdv, vsys=None, vrange='auto', dv=1*u.km/u.s, 
    return_lobes=False):

    try:
        dmdv.unit.to(u.Msun/(u.km/u.s))
    except u.UnitConversionError:
        print("dmdv must be in units of mass/velocity.")
        raise

    channel_width = dmdv.spectral_axis[1] - dmdv.spectral_axis[0]

    if vrange == 'auto':
        if vsys is None:
            print("Must set vsys in order to use auto vrange.")
            raise
        elif np.size(vsys) > 1:
            i_blue = dmdv.spectral_axis <= (min(vsys) - dv)
            i_red = dmdv.spectral_axis >= (max(vsys) + dv)
            i_between = (dmdv.spectral_axis > (min(vsys) - dv)) & (dmdv.spectral_axis < (max(vsys) + dv))
            v_out_blue = abs(min(vsys) - dmdv.spectral_axis)
            # print(v_out_blue)
            v_out_red = abs(max(vsys) - dmdv.spectral_axis)
        else:
            i_blue = dmdv.spectral_axis < vsys - dv
            i_red = dmdv.spectral_axis > vsys + dv
            i_between = (dmdv.spectral_axis > vsys - dv) & (dmdv.spectral_axis < vsys + dv)
            v_out_blue = abs(vsys - dmdv.spectral_axis)
            v_out_red = abs(vsys - dmdv.spectral_axis)

        v_out = v_out_blue
        # print(v_out)
        v_out[i_blue] = v_out_blue[i_blue]
        v_out[i_red] = v_out_red[i_red]
        v_out[i_between] = 0
        # print(v_out)
        mom_blue = (dmdv[i_blue]*v_out[i_blue]).nansum()*channel_width
        mom_red = (dmdv[i_red]*v_out[i_red]).nansum()*channel_width
        mom = mom_blue + mom_red

    elif np.ndim(vrange) > 2:
        mom = 0
        for vra in vrange:
            i_vra = (dmdv.spectral_axis >= min(vra)) & (dmdv.spectral_axis <= max(vra))
            v_out = abs(vsys - dmdv.spectral_axis)
            mom += (dmdv[i_vra]*v_out[i_vra]).nansum()*channel_width

    else:
        i_vra = (dmdv.spectral_axis >= vrange[0]) & (dmdv.spectral_axis <= vrange[1])
        v_out = abs(vsys - dmdv.spectral_axis)
        mom = (dmdv[i_vra]*v_out[i_vra]).nansum()*channel_width

    if return_lobes:
        return mom_blue, mom_red
    else:
        return mom



def energy_dmdv(dmdv, vsys=None, vrange='auto', dv=1*u.km/u.s, 
    return_lobes=False):

    try:
        dmdv.unit.to(u.Msun/(u.km/u.s))
    except u.UnitConversionError:
        print("dmdv must be in units of mass/velocity.")
        raise

    channel_width = dmdv.spectral_axis[1] - dmdv.spectral_axis[0]

    if vrange == 'auto':
        if vsys is None:
            print("Must set vsys in order to use auto vrange.")
            raise
        elif np.size(vsys) > 1:
            i_blue = dmdv.spectral_axis <= (min(vsys) - dv)
            i_red = dmdv.spectral_axis >= (max(vsys) + dv)
            i_between = (dmdv.spectral_axis > (min(vsys) - dv)) & (dmdv.spectral_axis < (max(vsys) + dv))
            v_out_blue = abs(min(vsys) - dmdv.spectral_axis)
            v_out_red = abs(max(vsys) - dmdv.spectral_axis)
        else:
            i_blue = dmdv.spectral_axis < vsys - dv
            i_red = dmdv.spectral_axis > vsys + dv
            i_between = (dmdv.spectral_axis > vsys - dv) & (dmdv.spectral_axis < vsys + dv)
            v_out_blue = abs(vsys - dmdv.spectral_axis)
            v_out_red = abs(vsys - dmdv.spectral_axis)

        v_out = v_out_blue

        v_out[i_blue] = v_out_blue[i_blue]
        v_out[i_red] = v_out_red[i_red]
        v_out[i_between] = 0

        energy_blue = 0.5*(dmdv[i_blue]*v_out[i_blue]**2.).nansum()*channel_width
        energy_red = 0.5*(dmdv[i_red]*v_out[i_red]**2.).nansum()*channel_width
        energy = energy_blue + energy_red

    elif np.ndim(vrange) > 2:
        energy = 0
        for vra in vrange:
            i_vra = (dmdv.spectral_axis >= min(vra)) & (dmdv.spectral_axis <= max(vra))
            v_out = abs(vsys - dmdv.spectral_axis)
            energy += 0.5*(dmdv[i_vra]*v_out[i_vra]**2.).nansum()*channel_width
    else:
        i_vra = (dmdv.spectral_axis >= vrange[0]) & (dmdv.spectral_axis <= vrange[1])
        v_out = abs(vsys - dmdv.spectral_axis)
        energy = 0.5*(dmdv[i_vra]*v_out[i_vra]**2.).nansum()*channel_width

    if return_lobes:
        return energy_blue, energy_red
    else:
        return energy

def vmax_cube(cube, vsys=0*u.km/u.s, snr=0, min_pixels=1,
         empty_velocity_range=[-2,0]*u.km/u.s,
         return_lobes=True, plot=False):
    """
    Returns the maximum outflow velocity on either side of vsys in the given cube.
    """
    if snr:
        rms = rms(cube, velocity_range=empty_velocity_range)
        cube = cube.with_mask(snr*rms)

    if np.size(vsys) > 1:
            v_out_blue = min(vsys) - cube.spectral_axis
            v_out_red = -max(vsys) + cube.spectral_axis
    else:
            v_out_blue = vsys - cube.spectral_axis
            v_out_red = -vsys + cube.spectral_axis
    #Selects channels where at least min_pixels are detected to the given snr threshold. 
    chans_detected = np.count_nonzero(~cube.apply_numpy_function(np.isnan), axis=(1,2)) >= min_pixels
    print(chans_detected)
    print(v_out_red)

    try:
        vmax_blue = min(v_out_blue[(~chans_detected) & (v_out_blue > 0)]) - cube.spectral_axis.diff()[0]
    except ValueError:
        print("outflow is detected at all blue channels, returning lower limit of vmax_blue.")
        vmax_blue = max(v_out_blue)

    try:
        vmax_red = min(v_out_red[(~chans_detected) & (v_out_red > 0)]) - cube.spectral_axis.diff()[0]
    except ValueError:
        print("outflow is detected at all red channels, returning lower limit of vmax_red.")
        vmax_red = max(v_out_red)
    # vmax_blue = cube.spectral_axis[]

    if plot:
        plt.plot(cube.spectral_axis.to(u.km/u.s), cube.mean((1,2)))
        plt.axvline((vsys - vmax_blue).value, label=r"Blue v$_{max}$", ls=':', color='tab:blue')
        plt.axvline((vsys + vmax_red).value,  label=r"Red v$_{max}$", ls=':', color='tab:red')

    return vmax_blue.to(u.km/u.s), vmax_red.to(u.km/u.s)

def rmax_cube(cube, coord_zero=None, dist=414*u.pc, plot=False):
    """
    Calculates the maximum distance from a specific coordinate in a cube. If coord_zero
    is None, defaults to calculating the distanc from the central pixel. The cube should
    be masked to only calculate the maximum distance inside the outflow mask. 
    """
    from astropy.wcs.utils import proj_plane_pixel_scales

    if coord_zero is None:
        coord_zero = SkyCoord(cube.longitude_extrema.mean(), cube.latitude_extrema.mean()) 

    dec, ra = cube.spatial_coordinate_map
    dec = dec - coord_zero.dec
    ra = ra - coord_zero.ra
    dist_arr = (ra**2 + dec**2)**0.5 * (2*np.pi / (360*u.deg)) * dist

    isnan = np.isnan(cube[0])
    dist_arr[isnan] = np.nan
    rmax = np.nanmax(dist_arr)

    if plot:
        dist_arr[isnan] = np.nan
        plt.imshow(dist_arr.value, origin='lower')
        plt.colorbar()

    return rmax

def outflow_angles(cube, pa_guess=0*u.deg, autoguess=False,
    coord_zero=None, plot=False, bins='knuth',
    plothist_kwargs=dict(),
    plotfit_kwargs=dict(),
    fit_kwargs=dict(autoguess=True),
    return_error=False,
    count_voxels=False):
    """
    Input a cube with outflow lobe mask.
    """
    from stamp import fit_gaussian
    from astropy.visualization import hist as plot_hist
    from astropy.stats import histogram as hist

    dec, ra = cube.spatial_coordinate_map

    if coord_zero is None:
        coord_zero = SkyCoord(np.mean(cube.longitude_extrema), np.mean(cube.latitude_extrema))

    ra = ra - coord_zero.ra
    dec = dec - coord_zero.dec

    #Ranges from [-pi, pi] radians, East of North
    pa_arr = np.arctan2(ra, dec)

    if autoguess:
        mom0 = cube.moment0()
        # print("In outflow_angle autoguess: ", mom0)
        ra_peak = ra[mom0 == np.nanmax(mom0)][0]
        dec_peak = dec[mom0 == np.nanmax(mom0)][0]
        pa_guess = np.arctan2(ra_peak.value, dec_peak.value)*180*u.deg/(np.pi)

    print("Guessing position angle of {}".format(pa_guess))

    #Rotate the positional angle array by the inital guess value for the outflow pa.
    c, s = np.cos(pa_guess), np.sin(pa_guess)
    R = np.array([[c,-s], [s, c]])
    rot = R.dot(np.array([ra.flatten(),dec.flatten()]))
    ra_rot, dec_rot = rot.reshape(np.array([ra,dec]).shape)
    pa_arr_rot = np.arctan2(ra_rot, dec_rot)

    pa_arr_rot = pa_arr_rot[~np.isnan(cube[0])]

    if plot:
        hist, bin_edges, patches = plot_hist(pa_arr_rot.flatten()*180/np.pi, bins=bins, **plothist_kwargs)
    else:
        hist, bin_edges = hist(pa_arr_rot.flatten()*180/np.pi, bins=bins)

    bin_centers = np.array([0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(len(bin_edges)-1)])

    if return_error:
        g, cov = fit_gaussian(bin_centers, hist, return_cov=True, **fit_kwargs)
    else:  
        g = fit_gaussian(bin_centers, hist, **fit_kwargs)


    position_angle = pa_guess + g.mean*u.deg
    if return_error:
        e_position_angle = np.sqrt(cov[1,1])*u.deg
    fwqm = 3.3302*g.stddev
    opening_angle = fwqm*u.deg
    if return_error:
        e_opening_angle = 3.3302*np.sqrt(cov[2,2])*u.deg

    if plot:
        x = np.linspace(bin_centers[0], bin_centers[-1], 100)
        y = g(x)
        plt.plot(x, y, **plotfit_kwargs)
        plt.xlabel("PA - {: .2g} [deg]".format(pa_guess))
        plt.ylabel("Pixels")


    if return_error:
        return position_angle, opening_angle, e_position_angle, e_opening_angle
    else:
        return position_angle, opening_angle

# def parse_vels():
#     """
#     Read table and parse the columns, calculating the red/blue velocity
#     """

def mass(column_density, distance=414*u.pc, molecule='H2',
    return_map=False, mass_unit=u.Msun):
    """
    If column_density is a cube, return_map will return mass cube.
    """
    if molecule == 'H2':
        mass_per_molecule = 2.34e-24*u.gram 

    pixel_angle = abs(column_density.header['CDELT2']) * u.deg
    pixel_area = (pixel_angle.to(u.radian).value * distance)**2.
    #n_pixels = nH2[~np.isnan(nH2)].size
    try:
        mass_map = (column_density * mass_per_molecule * pixel_area).to(mass_unit)
    except:
        mass_map = (column_density * mass_per_molecule * pixel_area)


    if return_map:
        return mass_map
    else:
        return u.Quantity(mass_map.nansum())

def momentum(mass, velocity, unit=u.Msun*(u.km/u.s)):
    return (mass * velocity).to(unit)
def energy(mass, velocity, unit=u.erg, fwhm_to_3Dsigma = False):
    if fwhm_to_3Dsigma:
        return ((np.sqrt(3) / np.sqrt(8 * np.log(2))) * 0.5 * mass * velocity ** 2.).to(unit)
    else:
        return (0.5 * mass * velocity ** 2.).to(unit)

def mass_loss_rate(momentum, wind_velocity=200.*u.km/u.s, wind_timescale=1.*u.Myr):
    """
    From Arce et al 2011:
    The wind velocity is typically assumed to be close to the escape velocity
    of the star, which is about (1–4) × 102 km s−1 for
    low- and intermediate-mass stars (Lamers & Cassinelli 1999).
    Here we assume vw = 200 km s−1.

    For the most part, the candidate sources appear to be in the Class II
    or Class III stage, which implies that these have ages of about
    1–3 Myr (Evans et al. 2009). Assuming that the wind has been active
    for most of the lifetime of the star, then we can assume that τw ∼ 1 Myr.
    """
    return (momentum / (wind_velocity * wind_timescale)).to(u.Msun/u.yr)

def wind_energy_rate(momentum, wind_velocity=200.*u.km/u.s, wind_timescale=1.*u.Myr, sigma_3D=2.9*u.km/u.s):
    return ((1./2.) * mass_loss_rate(momentum) * wind_velocity * sigma_3D).to(u.erg/u.s)

def dissipation_time(cube):
    return 0.5 * d / cube.linewidth_sigma().mean()


def dissipation_rate(E_turb):
    """Estimates of the value of η from numerical
     simulations of clouds range between ∼1 and 10 (McKee 1989; Mac Low 1999).
      Assuming a gas density averaged over the entire cloud complex of∼103 cm−3 andη=5
      ,weobtain tdiss ∼5×10^6 yr,which results in a turbulent energy dissipation rate 
      of approximately 10^33 erg/s
    """
    return E_turb / 5e6*u.yr

def dpdt(M, R, sigma_los):
    return (6.4E-4*u.Msun*u.km/(u.s*u.yr)) * (M / (500*u.Msun)) * (R / (0.5*u.pc))**(-1.) * (sigma_los / (u.km/u.s))**2.


if __name__ == '__main__':
    main()