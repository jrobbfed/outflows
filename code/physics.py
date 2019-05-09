from spectral_cube import SpectralCube
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.io import ascii
import shells
#Calculate various physical quantities from
#spectral cubes, spectra, and shell parameters.
nro_12co = "../nro_maps/12CO_20170514_FOREST-BEARS_spheroidal_grid7.5_dV0.099kms_xyb_YS_regrid0.11kms_reproj.fits"
nro_13co = "../nro_maps/13CO_BEARS-FOREST_20170913_7.5grid_Spheroidal_Tmb_0.11kms_xy_YS.fits" 
nro_13co_divided = "../nro_maps/13CO_20170518_FOREST-BEARS_spheroidal_grid7.5_dV0.11kms_xyb_YS_regridto12CO_divide1.4.fits" 



#Old numbering
# best_shells = [3,6,9,11,17,18,21,24,25,30,36,37]
# north_shells = [18,19,20,21,22,23,24,29,40]
# central_shells = [16,17,26,30,36,38,39]
# south_shells = [3,4,5,6,7,15,28,33,34,35]
# l1641_shells = [1,2,8,9,10,11,12,13,14,25,27,31,32,37,41,42]
#New N-S ordering
best_shells = [1,5,6,7,11,13,25,28,32,37,40,42]
north_shells = [1,2,3,4,5,6,7,8,9]
central_shells = [10,11,12,13,14,15,16]
south_shells = [17,18,19,20,21,22,23,24,25,26]
l1641_shells = [27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42]

berne_dec = [-5.8*u.deg, -4.95*u.deg]
berne_ra = [5.62*15*u.deg, 5.554*15*u.deg]
def main():
    import shell_model
    import shells
    dist = 414*u.pc
    # north_energy, central_energy, south_energy, l1641n_energy = [
    # 7.8e46, 2e47, 1.4e47, 1.6e47]
    # north_m, central_m, south_m, l1641n_m = np.array([4048, 3736, 5001, 5196])*u.Msun
    # north_r, central_r, south_r, l1641n_r = np.array([2, 1, 2.5, 2])*u.pc
    # north_sigma, central_sigma, south_sigma, l1641n_sigma = np.array([1.6, 1.7, 1.6, 1.6])*u.km/u.s
    # north_dpdt, central_dpdt, south_dpdt, l1641n_dpdt = dpdt(north_m, north_r, north_sigma), dpdt(central_m, central_r, central_sigma), dpdt(south_m, south_r, south_sigma), dpdt(l1641n_m, l1641n_r, l1641n_sigma)

    # turbs = [north_energy, central_energy, south_energy, l1641n_energy]
    # turb_dpdt = [north_dpdt, central_dpdt, south_dpdt, l1641n_dpdt]
    # print(turbs, turb_dpdt)
    # table_shell_parameters(param_file="shell_parameters_full_NtoS.txt", best_n=best_shells,
    #     table_name="NEW_shell_parameters_NtoS_tex_texp.txt", usecols=[0,1,2,3,4,5,6,7,8], show_texp=True)
    # table_shell_physics(param_file="shell_parameters_full_NtoS.txt", scale_energy=1e44, scale_mdot=1e-7, scale_L=1e31, scale_Edot=1e31,
    #   scale_dpdt=1e-4, low_name="_properties_low_1213co5sig_NtoSorder",
    #     mid_name="_properties_mid_1213co5sig_NtoSorder", hi_name="_properties_hi_1213co5sig_NtoSorder",
    #     table_name="shell_physics_all_5sig_NtoS.txt")
    table_subregions(scale_energy=1e46, table_name='subregions.tex')


    plot_physicsrange(column=1, plotname="../paper/massrange_all.pdf", low_name="_properties_low_1213co5sig_NtoSorder",
 mid_name="_properties_mid_1213co5sig_NtoSorder", hi_name="_properties_hi_1213co5sig_NtoSorder")
    plot_physicsrange(column=2, plotname="../paper/momentumrange_all.pdf", low_name="_properties_low_1213co5sig_NtoSorder",
 mid_name="_properties_mid_1213co5sig_NtoSorder", hi_name="_properties_hi_1213co5sig_NtoSorder")
    plot_physicsrange(column=3, plotname="../paper/energyrange_all.pdf", low_name="_properties_low_1213co5sig_NtoSorder",
 mid_name="_properties_mid_1213co5sig_NtoSorder", hi_name="_properties_hi_1213co5sig_NtoSorder",
 scale=1e46)

    # plot_physicsrange(column=1, plotname=None, best_n=north_shells)
    # plot_physicsrange(column=2, plotname=None, best_n=north_shells)
    # plot_physicsrange(column=3, plotname=None, best_n=north_shells)
    # plot_physicsrange(column=3, plotname=None, best_n=set(north_shells).intersection(best_shells))
    # # plot_physicsrange(column=1, plotname=None, best_n=central_shells)
    # # plot_physicsrange(column=2, plotname=None, best_n=central_shells)
    # plot_physicsrange(column=3, plotname=None, best_n=central_shells)
    # plot_physicsrange(column=3, plotname=None, best_n=set(central_shells).intersection(best_shells))    
    # # plot_physicsrange(column=1, plotname=None, best_n=south_shells)
    # # plot_physicsrange(column=2, plotname=None, best_n=south_shells)
    # plot_physicsrange(column=3, plotname=None, best_n=south_shells)
    # plot_physicsrange(column=3, plotname=None, best_n=set(south_shells).intersection(best_shells))    
    # # plot_physicsrange(column=1, plotname=None, best_n=south_shells)
    # print(plot_physicsrange(column=1, plotname=None, best_n=l1641_shells))
    # print(plot_physicsrange(column=2, plotname=None, best_n=l1641_shells))
    # print("Energy: ", plot_physicsrange(column=3, plotname=None, best_n=l1641_shells))
    # plot_physicsrange(column=3, plotname=None, best_n=set(l1641_shells).intersection(best_shells))    

    # all_low_tables = glob.glob("shell*properties*low_13co_1.4.txt")
    # all_hi_tables = glob.glob("shell*properties*hi_13co_1.4.txt")
    # robust_low_tables = glob.glob("shell[369]_properties_low_13co_1.4.txt")\
    #  + glob.glob("shell1[178]_properties_low_13co_1.4.txt")\
    #  + glob.glob("shell2[145]_properties_low_13co_1.4.txt")\
    #  + glob.glob("shell3[078]_properties_low_13co_1.4.txt")
    # robust_hi_tables = glob.glob("shell[369]_properties_hi_13co_1.4.txt")\
    #  + glob.glob("shell1[178]_properties_hi_13co_1.4.txt")\
    #  + glob.glob("shell2[145]_properties_hi_13co_1.4.txt")\
    #  + glob.glob("shell3[078]_properties_hi_13co_1.4.txt")

    # all_low_tables = glob.glob("shell*properties*low.txt")
    # all_hi_tables = glob.glob("shell*properties*hi.txt")
    # robust_low_tables = glob.glob("shell[369]_properties_low.txt")\
    #  + glob.glob("shell1[178]_properties_low.txt")\
    #  + glob.glob("shell2[145]_properties_low.txt")\
    #  + glob.glob("shell3[078]_properties_low.txt")
    # robust_hi_tables = glob.glob("shell[369]_properties_hi.txt")\
    #  + glob.glob("shell1[178]_properties_hi.txt")\
    #  + glob.glob("shell2[145]_properties_hi.txt")\
    #  + glob.glob("shell3[078]_properties_hi.txt")
    # mass_all = hist_physics(table_list=all_low_tables+all_hi_tables,
    #     column=1, table_list_shaded=robust_low_tables+robust_hi_tables,
    #     plotname='hist_mass_all.png')
    # momentum_all = hist_physics(table_list=all_low_tables+all_hi_tables,
    #     column=2, table_list_shaded=robust_low_tables+robust_hi_tables,
    #     plotname='hist_momentum_all.png')
    # energy_all = hist_physics(table_list=all_low_tables+all_hi_tables,
    #     column=3, table_list_shaded=robust_low_tables+robust_hi_tables,
    #     plotname='hist_energy_all.png')

    # mass_all_low = hist_physics(table_list=all_low_tables, column=1,
    #     plotname='hist_mass_low_all_new13co.png', table_list_shaded=robust_low_tables
    #     #xlim=[-50,300]
    #     )
    # mass_all_hi = hist_physics(table_list=all_hi_tables, column=1,
    #     plotname='hist_mass_hi_all_new13co.png', table_list_shaded=robust_hi_tables
    #     #xlim=[-50,300]
    #     )
    # momentum_all_low = hist_physics(table_list=all_low_tables, column=2,
    #     plotname='hist_momentum_low_all_new13co.png', table_list_shaded=robust_low_tables)
    # momentum_all_hi = hist_physics(table_list=all_hi_tables, column=2,
    #     plotname='hist_momentum_hi_all_new13co.png', table_list_shaded=robust_hi_tables)
    # energy_all_low = hist_physics(table_list=all_low_tables, column=3,
    #     plotname='hist_energy_low_all_new13co.png', table_list_shaded=robust_low_tables)
    # energy_all_hi = hist_physics(table_list=all_hi_tables, column=3,
    #     plotname='hist_energy_hi_all_new13co.png', table_list_shaded=robust_hi_tables)

    # print("Total for all shells:\n Mass - {} to {} Msun\n\
    #  Momentum - {} to {} Msun km/s\n\
    #  Energy - {} to {} erg\n".format(
    #     mass_all_low,mass_all_hi,
    #     momentum_all_low,momentum_all_hi,
    #     energy_all_low,energy_all_hi))


    # mass_robust_low = hist_physics(table_list=robust_low_tables, column=1,
    #     plotname='hist_mass_low_robust_new13co.png')
    # mass_robust_hi = hist_physics(table_list=robust_hi_tables, column=1,
    #     plotname='hist_mass_hi_robust_new13co.png')
    # momentum_robust_low = hist_physics(table_list=robust_low_tables, column=2,
    #     plotname='hist_momentum_low_robust_new13co.png')
    # momentum_robust_hi = hist_physics(table_list=robust_hi_tables, column=2,
    #     plotname='hist_momentum_hi_robust_new13co.png')
    # energy_robust_low = hist_physics(table_list=robust_low_tables, column=3,
    #     plotname='hist_energy_low_robust_new13co.png')
    # energy_robust_hi = hist_physics(table_list=robust_hi_tables, column=3,
    #     plotname='hist_energy_hi_robust_new13co.png')

    # print("Total for most confident 12 shells:\n Mass - {} to {} Msun\n\
    #  Momentum - {} to {} Msun km/s\n\
    #  Energy - {} to {} erg\n".format(
    #     mass_robust_low,mass_robust_hi,
    #     momentum_robust_low,momentum_robust_hi,


    #     energy_robust_low,energy_robust_hi))
    # cube_12co = SpectralCube.read(nro_12co)
    # cube_13co = SpectralCube.read(nro_13co)
    #cube_13co = SpectralCube.read(nro_13co_divided)
    #return

    # cube_12co = SpectralCube.read(nro_12co)
    # cube_13co = SpectralCube.read(nro_13co)

    # shell_list = shells.get_shells(region_file="../shell_candidates/AllShells_NtoS_old.reg",
    # velocity_file="../shell_candidates/AllShells_vrange_NtoS.txt")

    # for n in range(12,13): #shell = shell_list[n] print("Doing Shell {}".format(n)) shell = shell_list[n-1] print(n, shell.ra,shell.dec) # break #     params = np.loadtxt("shell_parameters_full_NtoS.txt") params = params[params[:,0] == 1.*n, 1:][0] r_best, r_sig = params[0], params[1] dr_best, dr_sig = params[2], params[3] vexp_best, vexp_sig = params[4], params[5] v0_best, v0_sig = params[6], params[7] #     dv = 0.1 v0_sample = np.arange(v0_best-v0_sig, v0_best+v0_sig, dv) print(v0_sample) N = v0_sample.size print(N) #     ## Loop over several v0 values for the minimum r, dr, vexp.  properties_low = np.empty((N, 4)) r = (r_best - r_sig) * u.pc dr = (dr_best - dr_sig) * u.pc vexp = (vexp_best - vexp_sig) * u.km/u.s #     for i,v0 in enumerate(v0_sample): #         v0 = v0 * u.km/u.s   # v0 = 14.25*u.km/u.s print(shell.ra, shell.dec, r, dr, vexp, v0) try: mass, momentum, energy = calc_physics( ra=shell.ra, dec=shell.dec, r=r, dr=dr, vexp=vexp, v0=v0, dist=dist, cube_12co=cube_12co, cube_13co=cube_13co, shell=True, plot=False, shell_snr_cutoff=5.) properties_low[i] = [v0.value, mass.value, momentum.value, energy.value] #             print("Shell Physical Properties:") print("------------------------------") print("Mass = {}".format(mass)) print("Expansion Velocity = {}".format(vexp)) print("Momentum = {}".format(momentum)) print("Energy = {}".format(energy)) except ValueError: print("Shell {} failed due to mismatched data shape.".format(n)) #     np.savetxt("shell{}_properties_low_1213co5sig_NtoSorder.txt".format(n), properties_low) #     ### Loop over several v0 values for the mid r, dr, vexp.  properties_mid = np.empty((N, 4)) r = (r_best) * u.pc dr = (dr_best) * u.pc vexp = (vexp_best) * u.km/u.s #     for i,v0 in enumerate(v0_sample): #         v0 = v0 * u.km/u.s   # v0 = 14.25*u.km/u.s print(shell.ra, shell.dec, r, dr, vexp, v0) try: mass, momentum, energy = calc_physics( ra=shell.ra, dec=shell.dec, r=r, dr=dr, vexp=vexp, v0=v0, dist=dist, cube_12co=cube_12co, cube_13co=cube_13co, shell_snr_cutoff=5.) properties_mid[i] = [v0.value, mass.value, momentum.value, energy.value] #             print("Shell Physical Properties:") print("------------------------------") print("Mass = {}".format(mass)) print("Expansion Velocity = {}".format(vexp)) print("Momentum = {}".format(momentum)) print("Energy = {}".format(energy)) except ValueError: print("Shell {} failed due to mismatched data shape.".format(n)) #     np.savetxt("shell{}_properties_mid_1213co5sig_NtoSorder.txt".format(n), properties_mid) #     ## Loop over v0 values with the maximum r, dr, vexp.  properties_hi = np.empty((N, 4)) r = (r_best + r_sig) * u.pc dr = (dr_best + dr_sig) * u.pc vexp = (vexp_best + vexp_sig) * u.km/u.s for i,v0 in enumerate(v0_sample): #         v0 = v0 * u.km/u.s   print(shell.ra, shell.dec, r, dr, vexp, v0) try: mass, momentum, energy = calc_physics( ra=shell.ra, dec=shell.dec, r=r, dr=dr, vexp=vexp, v0=v0, dist=dist, cube_12co=cube_12co, cube_13co=cube_13co, shell_snr_cutoff=5.) properties_hi[i] = [v0.value, mass.value, momentum.value, energy.value] #             print("Shell Physical Properties:") print("------------------------------") print("Mass = {}".format(mass)) print("Expansion Velocity = {}".format(vexp)) print("Momentum = {}".format(momentum)) print("Energy = {}".format(energy)) except ValueError: print("Shell {} failed due to mismatched data shape.".format(n)) #     np.savetxt("shell{}_properties_hi_1213co5sig_NtoSorder.txt".format(n), properties_hi)

###########################################-----------------------------------

    # # plt.figure()
    # plt.imshow(subcube_shell_12co_correct.moment0().data, interpolation='none')
    # plt.colorbar()
    # plt.title("Opacity-Corrected Integrated 12CO in K*m/s")
    # plt.savefig("corrected_12co.png")


    ## Calculate physics of cube subregions.

    # cube_12co = SpectralCube.read(nro_12co)
    # cube_13co = SpectralCube.read(nro_13co)
    # north_12co = subcube_region(cube_12co, region_index=2)
    # north_13co = subcube_region(cube_13co, region_index=2)
    # north_mass, north_momentum, north_energy = calc_physics(
    #     cube_12co=north_12co, cube_13co=north_13co, shell=False
    #     ,shell_snr_cutoff=5., average_Tex=False, linewidth_mode='sigma3D'
    #     )

    # central_12co = subcube_region(cube_12co, region_index=9)
    # central_13co = subcube_region(cube_13co, region_index=9)
    # central_mass, central_momentum, central_energy = calc_physics(
    #     cube_12co=central_12co, cube_13co=central_13co, shell=False
    #     ,shell_snr_cutoff=5., average_Tex=False, linewidth_mode='sigma3D'
    #     )

    # south_12co = subcube_region(cube_12co, region_index=7)
    # south_13co = subcube_region(cube_13co, region_index=7)
    # south_mass, south_momentum, south_energy = calc_physics(
    #     cube_12co=south_12co, cube_13co=south_13co, shell=False
    #     ,shell_snr_cutoff=5., average_Tex=False, linewidth_mode='sigma3D'
    #     )

    # l1641n_12co = subcube_region(cube_12co, region_index=8)
    # l1641n_13co = subcube_region(cube_13co, region_index=8)
    # l1641n_mass, l1641n_momentum, l1641n_energy = calc_physics(
    #     cube_12co=l1641n_12co, cube_13co=l1641n_13co, shell=False
    #     ,shell_snr_cutoff=5., average_Tex=False, linewidth_mode='sigma3D'
    #     )

    # berne_12co = cube_12co.subcube(berne_ra[0], berne_ra[1], berne_dec[0], berne_dec[1])
    # berne_13co = cube_13co.subcube(berne_ra[0], berne_ra[1], berne_dec[0], berne_dec[1])
    # # berne_mass, berne_momentum, berne_energy = calc_physics(
    # #     cube_12co=berne_12co, cube_13co=berne_13co, shell=False
    # #     ,shell_snr_cutoff=5., linewidth_mode='sigma'
    # #     )
    # # print("Using sigma: {} {} {}".format(berne_mass.to(u.Msun), berne_momentum.to(u.Msun*u.km/u.s), berne_energy.to(u.erg)))

    # berne_mass, berne_momentum, berne_energy = calc_physics(
    #     cube_12co=berne_12co, cube_13co=berne_13co, shell=False
    #     ,shell_snr_cutoff=5., linewidth_mode='sigma3D', plot=False,
    #     average_Tex=False)
    # print("Using FWHM: {} {} {}".format(berne_mass.to(u.Msun), berne_momentum.to(u.Msun*u.km/u.s), berne_energy.to(u.erg)))


    # for mass in [north_mass, central_mass, south_mass, l1641n_mass]:
    #     print(u.Quantity(mass).to(u.Msun))
    # for momentum in [north_mass, central_mass, south_mass, l1641n_mass]:
    #     print(u.Quantity(mass).to(u.Msun))
    # for energy in [north_energy, central_energy, south_energy, l1641n_energy]:
    #      print(u.Quantity(energy).to(u.erg))

# def table_shell_physics


    cube_12co = SpectralCube.read(nro_12co)
    cube_13co = SpectralCube.read(nro_13co)

    # shell_list = shells.get_shells()
    shell_list = shells.get_shells(velocity_file="../shell_candidates/AllShells_vrange_NtoS.txt",
            region_file="../shell_candidates/AllShells_NtoS.reg")

    # for n in range(39,40):
    #     #shell = shell_list[n]
    #     print("Doing Shell {}".format(n))
    #     shell = shell_list[n-1]


    #     params = np.loadtxt("shell_parameters_full_NtoS_Shell39Updated.txt")
    #     params = params[params[:,0] == 1.*n, 1:][0]
    #     r_best, r_sig = params[0], params[1]
    #     dr_best, dr_sig = params[2], params[3]
    #     vexp_best, vexp_sig = params[4], params[5]
    #     v0_best, v0_sig = params[6], params[7]

    #     dv = 0.1
    #     v0_sample = np.arange(v0_best-v0_sig, v0_best+v0_sig+dv, dv)
    #     N = v0_sample.size
    #     print(N)

    #     ### Loop over several v0 values for the minimum r, dr, vexp.
    #     properties_low = np.empty((N, 4))
    #     r = (r_best - r_sig) * u.pc
    #     dr = (dr_best - dr_sig) * u.pc
    #     vexp = (vexp_best - vexp_sig) * u.km/u.s
        
    #     for i,v0 in enumerate(v0_sample):

    #         v0 = v0 * u.km/u.s   
    #         # v0 = 14.25*u.km/u.s 
    #         print(shell.ra, shell.dec, r, dr, vexp, v0) 
    #         try:
    #             mass, momentum, energy = calc_physics(
    #             ra=shell.ra, dec=shell.dec, r=r, dr=dr, vexp=vexp, v0=v0, dist=dist,
    #             cube_12co=cube_12co, cube_13co=cube_13co, shell=True, plot=False, shell_snr_cutoff=5.)
    #             properties_low[i] = [v0.value, mass.value, momentum.value, energy.value]

    #             print("Shell Physical Properties:")
    #             print("------------------------------")
    #             print("Mass = {}".format(mass))
    #             print("Expansion Velocity = {}".format(vexp))
    #             print("Momentum = {}".format(momentum))
    #             print("Energy = {}".format(energy))
    #         except ValueError:
    #             print("Shell {} failed due to mismatched data shape.".format(n))
        
    #     np.savetxt("shell{}_properties_low_1213co5sig_NtoS_NewShell39.txt".format(n), properties_low)

    #     ### Loop over several v0 values for the mid r, dr, vexp.
    #     properties_mid = np.empty((N, 4))
    #     r = (r_best) * u.pc
    #     dr = (dr_best) * u.pc
    #     vexp = (vexp_best) * u.km/u.s
   
    #     for i,v0 in enumerate(v0_sample):

    #         v0 = v0 * u.km/u.s   
    #         # v0 = 14.25*u.km/u.s 
    #         print(shell.ra, shell.dec, r, dr, vexp, v0) 
    #         try:
    #             mass, momentum, energy = calc_physics(
    #             ra=shell.ra, dec=shell.dec, r=r, dr=dr, vexp=vexp, v0=v0, dist=dist,
    #             cube_12co=cube_12co, cube_13co=cube_13co, shell_snr_cutoff=5.)
    #             properties_mid[i] = [v0.value, mass.value, momentum.value, energy.value]

    #             print("Shell Physical Properties:")
    #             print("------------------------------")
    #             print("Mass = {}".format(mass))
    #             print("Expansion Velocity = {}".format(vexp))
    #             print("Momentum = {}".format(momentum))
    #             print("Energy = {}".format(energy))
    #         except ValueError:
    #             print("Shell {} failed due to mismatched data shape.".format(n))
        
    #     np.savetxt("shell{}_properties_mid_1213co5sig_NtoS_NewShell39.txt".format(n), properties_mid)



    #     ### Loop over v0 values with the maximum r, dr, vexp.
    #     properties_hi = np.empty((N, 4))
    #     r = (r_best + r_sig) * u.pc
    #     dr = (dr_best + dr_sig) * u.pc
    #     vexp = (vexp_best + vexp_sig) * u.km/u.s
    #     for i,v0 in enumerate(v0_sample):

    #         v0 = v0 * u.km/u.s   
    #         print(shell.ra, shell.dec, r, dr, vexp, v0) 
    #         try:
    #             mass, momentum, energy = calc_physics(
    #             ra=shell.ra, dec=shell.dec, r=r, dr=dr, vexp=vexp, v0=v0, dist=dist,
    #             cube_12co=cube_12co, cube_13co=cube_13co, shell_snr_cutoff=5.)
    #             properties_hi[i] = [v0.value, mass.value, momentum.value, energy.value]

    #             print("Shell Physical Properties:")
    #             print("------------------------------")
    #             print("Mass = {}".format(mass))
    #             print("Expansion Velocity = {}".format(vexp))
    #             print("Momentum = {}".format(momentum))
    #             print("Energy = {}".format(energy))
    #         except ValueError:
    #             print("Shell {} failed due to mismatched data shape.".format(n))
        


    #     np.savetxt("shell{}_properties_hi_1213co5sig_NtoS_NewShell39.txt".format(n), properties_hi)




def subcube_region(cube=None, region_file='../subregions/subregions.reg', region_index=2,
    region_unit = u.degree):
    import pyregion
    shape = pyregion.open(region_file)[region_index]
    coords = shape.coord_list * region_unit
    ra, dec, width, height = coords[0], coords[1], coords[2], coords[3]
    return cube.subcube(ra + width/2., ra - width/2., dec - height/2., dec + height/2.)


def hist_physics(table_list=None, table_list_shaded=None, mode='median', column=1, plotname="hist_mass_low_all.png",
    return_total=True, bins='auto', xlim=None):
    """
    table_list_shaded gives shell property tables that I want to shade in the histogram.
    """
    #print(table_list)
    x = []
    x_shaded = []
    for t in table_list:
        if table_list_shaded and t in table_list_shaded:
            x_shaded_sample = np.loadtxt(t)[:,column]
            if mode == 'median':
                x_shaded = np.append(x_shaded, np.median(x_shaded_sample))
        x_sample = np.loadtxt(t)[:,column]
        if mode == 'median':
            x = np.append(x, np.median(x_sample))

    plt.figure()
    n, b, patches = plt.hist(
        [x, x_shaded], histtype='step',
        bins=bins,# stacked=True,
        label=["All shells", "Best 12 shells"],
        facecolor='black', edgecolor='black')
    hatches = ['', '']
    fills = [False,True]
    for patch_set, hatch, fill in zip(patches, hatches, fills):
        plt.setp(patch_set, hatch=hatch)
        plt.setp(patch_set, fill=fill)

    if column == 1:
        plt.xlabel(r"Mass [$M_\odot$]")
    if column == 2:
        plt.xlabel(r"Momentum [$M_\odot$ km/s]")
    if column == 3:
        plt.xlabel(r"Kinetic Energy [erg]")
    if xlim:
        plt.xlim(xlim)

    plt.ylabel("count")
    plt.legend()

    plt.savefig(plotname)

    if return_total:
        return np.sum(x)

def table_shell_physics(param_file="shell_parameters_full_NtoS.txt", low_name="_properties_low_1213co5sig", mid_name="_properties_mid_1213co5sig",
    hi_name="_properties_hi_1213co5sig", name_tail=".txt", all_n=np.arange(1,43), best_n=best_shells, np_func=np.median,
    table_name="shell_physics_all_5sig.txt", usecols=[0,1,2,3], colnames=["v_exp", "mass", "momentum", "energy"],
    scale_energy=1., scale_mdot=1., scale_Edot=1., scale_L=1., scale_momentum=1., scale_dpdt=1.):
    
    #from astropy.table import Table
    
    shell_list = shells.get_shells()
    with open(table_name, 'w') as f:

        print("Shell&"
              "Mass&"
              "Momentum&"
              "Energy&"
              "Mechanical Luminosity&"
              "Momentum Injection Rate&"
              "Wind Mass Loss Rate&"
              "Wind Energy Injection Rate\\\\", file=f)

        print("Name&"
              "[M$_\\odot$]&"
              "[M$_\\odot$ km s$^{{-1}}$]&"
              "[$10^{{{}}}$ erg]&"
              "[$10^{{{}}}$ erg s$^{{-1}}$]&"
              "[$10^{{{}}}$ M$_\\odot$ km s$^{{-1}}$ yr$^{{-1}}$]&"
              "[$10^{{{}}}$ M$_\\odot$ yr$^{{-1}}$]&"
              "[$10^{{{}}}$ erg s$^{{-1}}$]\\\\".format(
                int(np.log10(scale_energy)),
                int(np.log10(scale_L)),
                int(np.log10(scale_dpdt)),
                int(np.log10(scale_mdot)),
                int(np.log10(scale_Edot))
                ), file=f)
        n_dpdt = np.zeros(3)
        c_dpdt = np.zeros(3)
        s_dpdt = np.zeros(3)
        l_dpdt = np.zeros(3)
        for n in all_n:
            params = np.loadtxt(param_file)
            params = params[params[:,0] == 1.*n, 1:][0]
            r_best, r_sig = params[0], params[1]
            dr_best, dr_sig = params[2], params[3]
            vexp_best, vexp_sig = params[4], params[5]
            v0_best, v0_sig = params[6], params[7]

            t_dyn = (r_best*u.pc / (vexp_best*u.km/u.s)).to(u.s).value
            t_dyn_yr = t_dyn / 3.154e7

            low = np_func(np.loadtxt("shell{}{}{}".format(n, low_name, name_tail), usecols=usecols), axis=0)
            mid = np_func(np.loadtxt("shell{}{}{}".format(n, mid_name, name_tail), usecols=usecols), axis=0)
            hi = np_func(np.loadtxt("shell{}{}{}".format(n, hi_name, name_tail), usecols=usecols), axis=0)
            #wind velocity of 200 km/s and wind timescale of 1e6 yrs, and sigma_3D of 1.7 km/s

            mass_three = np.array([low[1], mid[1], hi[1]])
            momentum_three = np.array([low[2], mid[2], hi[2]])
            energy_three = np.array([low[3], mid[3], hi[3]]) / scale_energy
            L_three = np.array([low[3], mid[3], hi[3]]) / t_dyn / scale_L
            dpdt_three = momentum_three / t_dyn_yr / scale_dpdt
            mdot_three = mass_loss_rate(momentum_three*u.Msun*u.km/u.s).value / scale_mdot
            Edot_three = wind_energy_rate(momentum_three*u.Msun*u.km/u.s).value / scale_Edot

            mid_mass, hi_mass, low_mass = np.median(mass_three), np.max(mass_three), np.min(mass_three)
            mid_momentum, hi_momentum, low_momentum = np.median(momentum_three)/scale_momentum, np.max(momentum_three)/scale_momentum, np.min(momentum_three)/scale_momentum
            mid_energy, hi_energy, low_energy = np.median(energy_three), np.max(energy_three), np.min(energy_three)
            mid_L, hi_L, low_L = np.median(L_three), np.max(L_three), np.min(L_three)
            mid_dpdt, hi_dpdt, low_dpdt = np.median(dpdt_three), np.max(dpdt_three), np.min(dpdt_three)
            mid_mdot, hi_mdot, low_mdot = np.median(mdot_three), np.max(mdot_three), np.min(mdot_three)
            mid_Edot, hi_Edot, low_Edot = np.median(Edot_three), np.max(Edot_three), np.min(Edot_three)

            if n in north_shells: 
                n_dpdt += np.array([low_dpdt, mid_dpdt, hi_dpdt])
            elif n in central_shells:
                c_dpdt += np.array([low_dpdt, mid_dpdt, hi_dpdt])
            elif n in south_shells:
                s_dpdt += np.array([low_dpdt, mid_dpdt, hi_dpdt])
            elif n in l1641_shells:
                l_dpdt += np.array([low_dpdt, mid_dpdt, hi_dpdt])

            def nice_round(x, max_decimal=1):
                if x < 1:
                    x = round(x, max_decimal)
                else:
                    x = round(x)
                return x


            print("${:.4g}$&${:.4g}~[{:.4g}, {:.4g}]$"
                  "&${:.4g}~[{:.4g}, {:.4g}]$"
                  "&${:.4g}~[{:.4g}, {:.4g}]$"
                  "&${:.4g}~[{:.4g}, {:.4g}]$"
                  "&${:.4g}~[{:.4g}, {:.4g}]$"
                  "&${:.4g}~[{:.4g}, {:.4g}]$"
                  "&${:.4g}~[{:.4g}, {:.4g}]$\\\\".format(
                n, nice_round(mid_mass), nice_round(low_mass), nice_round(hi_mass),
                   nice_round(mid_momentum), nice_round(low_momentum), nice_round(hi_momentum),
                   nice_round(mid_energy), nice_round(low_energy), nice_round(hi_energy),
                   nice_round(mid_L), nice_round(low_L), nice_round(hi_L),
                   nice_round(mid_dpdt), nice_round(low_dpdt), nice_round(hi_dpdt),
                   nice_round(mid_mdot), nice_round(low_mdot), nice_round(hi_mdot),
                   nice_round(mid_Edot), nice_round(low_Edot), nice_round(hi_Edot)),
                  file=f)

        print("dpdt North: ", n_dpdt)
        print("dpdt Central: ", c_dpdt)
        print("dpdt South: ", s_dpdt)
        print("dpdt L1641: ", l_dpdt)
        print("dpdt TOT: ", n_dpdt + c_dpdt + s_dpdt + l_dpdt)

def table_subregions(cube_12co="../nro_maps/12CO_20170514_FOREST-BEARS_spheroidal_grid7.5_dV0.099kms_xyb_YS_regrid0.11kms_reproj.fits",
    cube_13co="../nro_maps/13CO_BEARS-FOREST_20170913_7.5grid_Spheroidal_Tmb_0.11kms_xy_YS.fits",
    region_file="subregions.reg", outflow_file="outflows.txt", all_n=np.arange(1,43), best_n=best_shells,
    table_name="../paper/tables/subregions.tex", usecols=[0,1,2,3,4,5,6,7,8],
    shell_snr_cutoff=5., scale_energy=1e46):
    
    #from astropy.table import Table
    
    # cube_12co = SpectralCube.read(nro_12co)
    # cube_13co = SpectralCube.read(nro_13co)
    # north_12co = subcube_region(cube_12co, region_index=2)
    # north_13co = subcube_region(cube_13co, region_index=2)
    # north_mass, north_momentum, north_energy = calc_physics(
    #     cube_12co=north_12co, cube_13co=north_13co, shell=False
    #     ,shell_snr_cutoff=5., average_Tex=False
    #     )
    cube_12co = SpectralCube.read(nro_12co)
    cube_13co = SpectralCube.read(nro_13co)
    ### ENERGIES
    north_low, north_mid, north_hi = plot_physicsrange(column=3, plotname=None, best_n=north_shells)
    best_north_low, best_north_mid, best_north_hi = plot_physicsrange(column=3, plotname=None, best_n=set(north_shells).intersection(best_shells))
    # # plot_physicsrange(column=1, plotname=None, best_n=central_shells)
    # # plot_physicsrange(column=2, plotname=None, best_n=central_shells)
    central_low, central_mid, central_hi = plot_physicsrange(column=3, plotname=None, best_n=central_shells)
    best_central_low, best_central_mid, best_central_hi = plot_physicsrange(column=3, plotname=None, best_n=set(central_shells).intersection(best_shells))
    # # plot_physicsrange(column=1, plotname=None, best_n=south_shells)
    # # plot_physicsrange(column=2, plotname=None, best_n=south_shells)
    south_low, south_mid, south_hi = plot_physicsrange(column=3, plotname=None, best_n=south_shells)
    best_south_low, best_south_mid, best_south_hi = plot_physicsrange(column=3, plotname=None, best_n=set(south_shells).intersection(best_shells))
    # # plot_physicsrange(column=1, plotname=None, best_n=south_shells)
    # # plot_physicsrange(column=1, plotname=None, best_n=l1641_shells)
    l1641n_low, l1641n_mid, l1641n_hi = plot_physicsrange(column=3, plotname=None, best_n=l1641_shells)
    best_l1641n_low, best_l1641n_mid, best_l1641n_hi = plot_physicsrange(column=3, plotname=None, best_n=set(l1641_shells).intersection(best_shells))

    lows = [north_low, central_low, south_low, l1641n_low]
    mids = [north_mid, central_mid, south_mid, l1641n_mid]
    his = [north_hi, central_hi, south_hi, l1641n_hi]

    print("Lows: ", lows)
    print("Mids: ", mids)
    print("His: ", his)

    # north_12co = subcube_region(cube_12co, region_index=2)
    # north_13co = subcube_region(cube_13co, region_index=2)
    # north_mass, north_momentum, north_energy = calc_physics(
    #     cube_12co=north_12co, cube_13co=north_13co, shell=False
    #     ,shell_snr_cutoff=shell_snr_cutoff, average_Tex=False
    #     )


    # central_12co = subcube_region(cube_12co, region_index=9)
    # central_13co = subcube_region(cube_13co, region_index=9)
    # central_mass, central_momentum, central_energy = calc_physics(
    #     cube_12co=central_12co, cube_13co=central_13co, shell=False
    #     ,shell_snr_cutoff=shell_snr_cutoff, average_Tex=False
    #     )

    # south_12co = subcube_region(cube_12co, region_index=7)
    # south_13co = subcube_region(cube_13co, region_index=7)
    # south_mass, south_momentum, south_energy = calc_physics(
    #     cube_12co=south_12co, cube_13co=south_13co, shell=False
    #     ,shell_snr_cutoff=shell_snr_cutoff, average_Tex=False
    #     )

    # l1641n_12co = subcube_region(cube_12co, region_index=8)
    # l1641n_13co = subcube_region(cube_13co, region_index=8)
    # l1641n_mass, l1641n_momentum, l1641n_energy = calc_physics(
    #     cube_12co=l1641n_12co, cube_13co=l1641n_13co, shell=False
    #     ,shell_snr_cutoff=shell_snr_cutoff, average_Tex=False
    #     )
    north_energy, central_energy, south_energy, l1641n_energy = [
    7.8e46, 2e47, 1.4e47, 1.6e47]
    north_m, central_m, south_m, l1641n_m = np.array([4048, 3736, 5001, 5196])*u.Msun
    north_r, central_r, south_r, l1641n_r = np.array([2, 1, 2.5, 2])*u.pc
    north_sigma, central_sigma, south_sigma, l1641n_sigma = np.array([1.6, 1.7, 1.6, 1.6])*u.km/u.s
    north_dpdt, central_dpdt, south_dpdt, l1641n_dpdt = dpdt(north_m, north_r, north_sigma), dpdt(central_m, central_r, central_sigma), dpdt(south_m, south_r, south_sigma), dpdt(l1641n_m, l1641n_r, l1641n_sigma)

    turbs = [north_energy, central_energy, south_energy, l1641n_energy]
    turb_dpdt = [north_dpdt, central_dpdt, south_dpdt, l1641n_dpdt]
    print(turbs, turb_dpdt)
    outflow_table = ascii.read(outflow_file)
    print(outflow_table['energy'])

    with open(table_name, 'w') as f:

        print("\\begin{table*}\n"
        "\\centering"
        "\\caption{\\label{tab:impact} Impact by Cloud Region}"
        "\\begin{tabular}{cccc}"
        "\\hline"
        "\\hline", file=f)

        print("Subregion&"
              "$E_{{\\rm shells}}$ [$10^{{{}}}$ erg]&"
              "$E_{{\\rm outflows}}$ [$10^{{{}}}$ erg]&"
              "$E_{{\\rm turb}}$ [$10^{{{}}}$ erg]&".format(
                int(np.log10(scale_energy)), int(np.log10(scale_energy)), int(np.log10(scale_energy))),
                file=f)

        

        for i in range(len(lows)):

            subregion = outflow_table['subregion'][i]
            outflow = outflow_table['energy'][i]
            shell_low = lows[i]
            shell_mid = mids[i]
            shell_hi = his[i]
            turb = turbs[i]
            print("{}"
                  "&${:.2g}_{{-{:.2g}}}^{{+{:.2g}}}$"
                  "&${:.2g}$"
                  "&${:.2g}$\\\\".format(
                subregion, shell_mid/scale_energy, (shell_mid-shell_low)/scale_energy, (shell_hi-shell_mid)/scale_energy,
                outflow/scale_energy, turb/scale_energy),
                file=f)

        print("{}"
                  "&${:.2g}_{{-{:.2g}}}^{{+{:.2g}}}$"
                  "&${:.2g}$"
                  "&${:.2g}$\\\\".format(
                "Total", np.sum(mids)/scale_energy, (np.sum(mids) - np.sum(lows))/scale_energy,
                (np.sum(his)-np.sum(mids))/scale_energy,
                np.sum(outflow_table['energy'])/scale_energy, np.sum(turbs)/scale_energy),
                file=f)        

        print("\\\n"
        "\\hline"
        "\\end{tabular}"
        "\\end{table*}", file=f)

def table_shell_parameters(param_file="shell_parameters_full_NtoS.txt", all_n=np.arange(1,43), best_n=best_shells,
    table_name="shell_parameters_NtoS_tex.txt", usecols=[0,1,2,3,4,5,6,7,8], show_texp=False):
    from astropy.coordinates import SkyCoord
    #from astropy.table import Table
    kmspc_toMyr = 0.9778
    
    with open(table_name, 'w') as f:
        if show_texp:
            print("Shell&"
                  "$\\alpha (J2000)/\\delta (J2000)$&"
                  "$R$ [pc]&"
                  "$dr$ [pc]&"
                  "$v_{\\rm exp}$ km s$^{-1}$&"
                  "$v_{\\rm 0}$ km s$^{-1}$&"
                  "$t_{\\tm exp}$ Myr&", file=f)
            
        else:
            print("Shell&"
                  "$R$ [pc]&"
                  "$dr$ [pc]&"
                  "$v_{\\rm exp}$ km s$^{-1}$&"
                  "$v_{\\rm 0}$ km s$^{-1}$&", file=f)

        data = np.loadtxt(param_file)
        shell_list = shells.get_shells(velocity_file="../shell_candidates/AllShells_vrange_NtoS.txt",
            region_file="../shell_candidates/AllShells_NtoS.reg")
        for n in all_n:
            i = n-1
            line = data[i]

            shell = shell_list[i]
            hms = SkyCoord(shell.ra, shell.dec).ra.hms
            dms = SkyCoord(shell.ra, shell.dec).dec.dms
            ra_str = "${}^{{\\rm h}}{}^{{\\rm m}}{}^{{\\rm s}}.{}$".format(
                int(hms.h), int(hms.m), int(round(hms.s,1)), str(round(hms.s,1))[-1])
            dec_str = "${}\\arcdeg{}\\arcmin{}\\arcsec$".format(
                int(dms.d), int(dms.m), int(round(dms.s)))

            if show_texp:
                texp = (line[1]/line[5]) * kmspc_toMyr #Myr
                e_texp = texp * np.sqrt((line[2]/line[1])**2. + (line[6]/line[5])**2.)
                print("${}$&{}/{}&${:.3f}\pm{:.3f}$&${:.3f}\pm{:.3f}$&${:.2f}\pm{:.2f}$&${:.2f}\pm{:.2f}$&${:.2f}\pm{:.2f}$\\\\".format(
                    int(line[0]), ra_str, dec_str, line[1], line[2], line[3], line[4], line[5],
                    line[6], line[7], line[8], texp, e_texp),
                    file=f)
            else:
                print("${}$&${:.3f}\pm{:.3f}$&${:.3f}\pm{:.3f}$&${:.2f}\pm{:.2f}$&${:.2f}\pm{:.2f}$\\\\".format(
                    int(line[0]), line[1], line[2], line[3], line[4], line[5],
                    line[6], line[7], line[8]),
                    file=f)




def plot_physicsrange(low_name="_properties_low_1213co5sig_NtoSorder",
 mid_name="_properties_mid_1213co5sig_NtoSorder", hi_name="_properties_hi_1213co5sig_NtoSorder", name_tail=".txt",
    all_n=np.arange(1,43), best_n=best_shells, mode='median',
    column=1, plotname='massrange_all.png', lw=0.7, ms=4., mew=0.7, scale=1.):
    import matplotlib.lines as mlines
    plt.figure()

    total_low = 0
    total_mid = 0
    total_hi = 0
    L_mid = 0
    L_hi = 0
    L_low = 0
    for n in all_n:
        params = np.loadtxt("shell_parameters_full.txt")
        params = params[params[:,0] == 1.*n, 1:][0]
        r_best, r_sig = params[0], params[1]
        dr_best, dr_sig = params[2], params[3]
        vexp_best, vexp_sig = params[4], params[5]
        v0_best, v0_sig = params[6], params[7]
        tdyn = (r_best*u.pc / (vexp_best*(u.km/u.s))).to(u.s).value

        print(n, "shell{}{}{}".format(n, low_name, name_tail))
        low = np.loadtxt("shell{}{}{}".format(n, low_name, name_tail))[:,column]
        mid = np.loadtxt("shell{}{}{}".format(n, mid_name, name_tail))[:,column]
        hi = np.loadtxt("shell{}{}{}".format(n, hi_name, name_tail))[:,column]
        if mode == 'median':
            low = np.median(low)
            mid = np.median(mid)
            hi = np.median(hi)
        
        if n in best_n:
            #print(low,mid,hi)
            L_low += np.min([low,mid,hi]) / tdyn
            L_mid += np.median([low,mid,hi]) / tdyn
            L_hi += np.max([low,mid,hi]) / tdyn
            total_low += np.min([low,mid,hi])
            total_mid += np.median([low,mid,hi])
            total_hi += np.max([low,mid,hi])


            plt.plot([np.min([low,mid,hi])/scale, np.max([low,mid,hi])/scale], [n,n],
                color='k', ls='-', lw=lw)
            plt.plot(np.median([low,mid,hi])/scale, n, marker='o', color='k', ms=ms)
            

        else:
            plt.plot([np.min([low,mid,hi])/scale, np.max([low,mid,hi])/scale], [n,n],
                color='k', ls=':', lw=lw)
            plt.plot(np.median([low,mid,hi])/scale, n, marker='o', markerfacecolor='white', color='k',
                ms=ms, mew=mew)
            

    if column == 1:
        plt.xlim([-10.,505.])
        plt.xlabel(r"Mass [$M_\odot$]")
    if column == 2:
        plt.xlim([-50., 2700])
        plt.xlabel(r"Momentum [$M_\odot$ km/s]")
    if column == 3:
        plt.xlim([-0.2, 16.])
        plt.xlabel(r"Kinetic Energy [$10^{46}$ erg]")

    plt.ylabel("Shell Number")
    best_line = mlines.Line2D([], [],
     color='k', marker='o', linestyle='solid', label='Best 12 shells',
      lw=lw, ms=ms, mew=mew)
    other_line = mlines.Line2D([], [],
     color='k',marker='o', markerfacecolor='white', linestyle='dotted', label='Other shells',
      lw=lw, ms=ms, mew=mew)
    plt.legend(handles=[best_line, other_line], loc='best')
    try:
        plt.savefig(plotname, dpi=200)
    except TypeError:
        pass
    #   plt.show()
    print("Mechanical Luminosity of Shells [erg/s]: ", L_low, L_mid, L_hi)
    return total_low, total_mid, total_hi

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
            Tpeak = np.max(average_spec)
        else:
            Tpeak = np.average(cube.max(axis=0))
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


def extract_shell(cube_file=nro_12co, model_pars=None,
 mask_minimum=0.00001, return_mask=False, model_cube=None,
 keep_latlon=False, snr_cutoff=None, rms=0.):
    """
    Return a masked cube from an observed spectral cube file,
    where the mask is True wherever a model shell cube
    with parameters given by `model_pars` dictionary is > mask_minimum.

    """
    if not model_cube:
        model_cube = SpectralCube.read(shell_model.ppv_model(**model_pars))
    if type(cube_file) == str:
        if keep_latlon:
            obs_cube = SpectralCube.read(cube_file).spectral_slab(
                model_cube.spectral_extrema[0],
                model_cube.spectral_extrema[1])
        else:
            obs_cube = SpectralCube.read(cube_file).subcube(
                                model_cube.longitude_extrema[1],
                                model_cube.longitude_extrema[0],
                                model_cube.latitude_extrema[0],
                                model_cube.latitude_extrema[1],
                                model_cube.spectral_extrema[0],
                                model_cube.spectral_extrema[1])
    else:
        if keep_latlon:
            obs_cube = cube_file.spectral_slab(
                model_cube.spectral_extrema[0],
                model_cube.spectral_extrema[1])
        else:
            obs_cube = cube_file.subcube(
                                model_cube.longitude_extrema[1],
                                model_cube.longitude_extrema[0],
                                model_cube.latitude_extrema[0],
                                model_cube.latitude_extrema[1],
                                model_cube.spectral_extrema[0],
                                model_cube.spectral_extrema[1])
    print("Before changing wcs, obs_cube shape: ", obs_cube.shape)
    #Reset the cube wcs to the values corresponding to the subcube.
    obs_cube = SpectralCube(obs_cube.hdu.data, wcs=model_cube.wcs) * u.K
    shell_mask = model_cube > mask_minimum*u.dimensionless_unscaled
    print(model_cube.shape, obs_cube.shape)
    obs_cube_masked = obs_cube.with_mask(shell_mask)
    if snr_cutoff:
        obs_cube_masked = obs_cube_masked.with_mask(obs_cube > snr_cutoff * rms)
    #obs_array_masked = obs_cube_masked.filled_data[:,:,:]
    if return_mask:
        return obs_cube_masked, shell_mask
    else:
        return obs_cube_masked

def rms_map(cube=None, velocity_range=[[-3.,-0.1], [19.,20.]]*u.km/u.s):
    """
    Returns 2D array of the standard deviation of a spectral cube,
    calculated only in the emission-free channels.
    """

    channel_range = [[cube.closest_spectral_channel(vpair[0]),
                  cube.closest_spectral_channel(vpair[1])]
                 for vpair in velocity_range]

    emissionless_channels = np.concatenate(
        [np.arange(c[0], c[1]+1) for c in channel_range])
    emissionless_cube = cube.unmasked_data[emissionless_channels,:,:]
    rms_map = np.nanstd(emissionless_cube, axis=0)
    return rms_map

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

def mask_snr(cube=None, rms_map=None, snr=5., return_mask=False):
    """
    Returns the spectral cube with low significance voxels
    masked out, using a map of rms calculated in emission-free
    channels.
    """
    if return_mask:
        return (cube > snr * rms_map)
    else:
        return cube.with_mask(cube > snr * rms_map)

def cube_ratio(cubes=[None, None], rms_maps=[None, None], return_uncertainty=True):
    """
    Returns spectral cube with ratio (uncertainties on the ratio) between
    two cubes. Uncertainity calculated with error propagation assuming
    the error on each cube is given by the rms in emission-free channels.
    """
    cube_ratio = cubes[0] / cubes[1]
    if return_uncertainty:
        cube_ratio_uncertainty = cube_ratio * np.sqrt(
            (rms_maps[0] / cubes[0]) ** 2. +
            (rms_maps[1] / cubes[1]) ** 2.)
        return cube_ratio, cube_ratio_uncertainty
    else:
        return cube_ratio

def average_spectrum(cube=None, weights=1., axis=(1,2), return_std=True,
    ignore_nan=True):
    """
    Calculate the (weighted) average spectrum in a spectral cube
    Optionally calculate and return the (weighted) standard deviation
    spectrum. `weights` should be a numpy array with the same shape as `cube`
    `axis` denotes the axis/axes to average over. For a standard SpectralCube,
    the two spatial axes are (1,2).
    Returns
    -------
    average_spectrum : 1D array_like
    std_spectrum: 1D array_like, optional
    """
    weighted_data = cube.filled_data[:,:,:] * weights
    average_spectrum = np.nanmean(weighted_data, axis=axis)
    #average_spectrum = np.average(cube.filled_data[:,:,:],
    # weights=weights, axis=axis)

    if return_std:
        resids_squared = (cube.filled_data[:,:,:] - average_spectrum[:,np.newaxis,np.newaxis])**2. * weights
        std_spectrum = np.sqrt(
            np.nanmean(resids_squared, axis=axis))
        return average_spectrum, std_spectrum
    else:
        return average_spectrum

def opacity_correct(cube_thick, cube_thin=None, abundance_ratio=62.,
    snr_cutoff=5., empty_velocity_range=[[-3.,-0.1], [19.,20.]]*u.km/u.s,
    regrid_cube=False, plot_ratio=None,
    fit=True, fit_order=2, **kwargs):
    """
    Correct an optically thick emission line cube using an (assumed) optically
    thin emission line cube. The abundance ratio betweeen the cube and cube_thin
    isotopologues is given as `abundance_ratio`.
    
    `regrid_cube`: Optionally regrid `cube_thin` to the velocity grid of `cube`, preserving
    Nyquist sampling.
    
    Uses method detailed in Zhang et al 2016 (c.f. Dunham+14, Arce+2001)
    """
    if regrid_cube:
        cube_thin = regrid(cube_thin, cube_thick.spectral_axis)

    rms_thick = rms_map(cube_thick, empty_velocity_range)
    rms_thin = rms_map(cube_thin, empty_velocity_range)

    cube_thick_masked = mask_snr(cube_thick, rms_thick,
     snr=snr_cutoff)
    cube_thin_masked = mask_snr(cube_thin, rms_thin,
     snr=snr_cutoff)

    ratio, sigma_ratio = cube_ratio(
        cubes=[cube_thick_masked, cube_thin_masked],
        rms_maps=[rms_thick, rms_thin],
        return_uncertainty=True)

    #print(ratio.size, ratio.flattened().size)
    weights = 1. / (sigma_ratio.filled_data[:,:,:]**2.)
    #print(weights)
    average_ratio, std_ratio = average_spectrum(
        cube=ratio, weights=weights)
    #print(average_ratio, std_ratio, weights)

    if fit:
        #Fit with quadratic
        notnan = ~np.isnan(average_ratio)
        vel = ratio.spectral_axis.to(u.km/u.s).value
        #print(vel[notnan], average_ratio[notnan])
        fit_coeff, fit_cov = np.polyfit(vel[notnan], average_ratio[notnan],
            fit_order, w=(1/std_ratio[notnan]), cov=True)

        fit_func = np.poly1d(fit_coeff)
        ratio_spectrum_fit = fit_func(vel)
        cube_correct = (
        cube_thick * abundance_ratio / ratio_spectrum_fit[:,np.newaxis,np.newaxis]) #X_12,13 / (T_12/T_13)

    else:
        cube_correct = cube_thick * abundance_ratio / average_ratio[:,np.newaxis,np.newaxis]

    if plot_ratio: 
        plt.plot(vel, average_ratio, 'o')

        if fit:
            xfit = np.linspace(vel[notnan][0], vel[notnan][-1], 1000)
            yfit = fit_func(xfit)
            plt.plot(xfit, yfit)
        plt.fill_between(vel,
            average_ratio+std_ratio, average_ratio-std_ratio,
                alpha=0.3)
        plt.ylabel(r"T$_{12}$ / T$_{13}$", size=16)
        plt.xlabel("Velocity [km / s]", size=16)
        #plt.ylim(1.3,3.9)
        plt.title("Weighted Average and Std. Dev. of 5sigma 12CO/13CO Ratio")
        plt.savefig(plot_ratio)

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


    if moment0:
        return (factor * cube.moment0()).to(1/(u.cm**2.))
    else:
        return (factor * cube).decompose()

def Qrot_partial(Tex, B0_k=2.765*u.K, N=20):
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



def mass(column_density, distance=414*u.pc, molecule='H2',
    return_map=False, mass_unit=u.Msun):
    """
    
    """
    if molecule == 'H2':
        mass_per_molecule = 2.34e-24*u.gram

    pixel_angle = abs(column_density.header['CDELT2']) * u.deg
    pixel_area = (pixel_angle.to(u.radian).value * distance)**2.
    #n_pixels = nH2[~np.isnan(nH2)].size
    mass_map = (column_density * mass_per_molecule * pixel_area).to(mass_unit)
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