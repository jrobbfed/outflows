from astropy.io import ascii
import numpy as np
import astropy.units as u
from astropy.table import Table, Column
def main():

    filament_table(infile="physics_outflows.txt",
     outfile="filament_table_dist400pc.tex", dist_conversion=400./414.)
    # outflow_table()

    #Optionally convert distance using this ratio.
    #M,P,E will be multiplied by this value^2. Rmax and tdyn will be multiplied by this value.
    #Mdot, Pdot, Edot will be multiplied by this value. 
    # dist_conversion = 1
    # dist_conversion = 400./414.

    physics_table(infile="physics_outflows.txt", outfile="physics_table_dist400pc.tex",
        dist_conversion=400./414.)

def physics_table(infile="physics_outflows.txt",
 source_colname='source', source_outname="Source",
 lobe_colname='lobe', lobe_outname="Lobe",
 M_lowvel_colname='M_lowvel', M_outname=r"$M$",
 M_hivel_colname='M_hivel',
 M_scale=1,
 P_lowvel_colname='P_lowvel', P_outname=r"$P$",
 P_hivel_colname='P_hivel',
 P_scale=1,
 E_lowvel_colname='E_lowvel', E_outname=r"$E$",
 E_hivel_colname='E_hivel',
 E_scale=1e43,
 rmax_colname='rmax', rmax_outname=r"$R_{\rm max}$",
 rmax_scale=1,
 vmax_colname='vmax', vmax_outname=r"$v_{\rm max}$",
 vmax_scale=1,
 tdyn_colname='tdyn', tdyn_outname=r"$t_{\rm dyn}$",
 tdyn_scale=1e4,
 Mdot_lowvel_colname='Mdot_lowvel', Mdot_outname=r"$\dot M$",
 Mdot_hivel_colname='Mdot_hivel',
 Mdot_scale=1e-6,
 Pdot_lowvel_colname='Pdot_lowvel', Pdot_outname=r"$\dot P$",
 Pdot_hivel_colname='Pdot_hivel',
 Pdot_scale=1e-6,
 Edot_lowvel_colname='Edot_lowvel', Edot_outname=r"$\dot E$",
 Edot_hivel_colname='Edot_hivel',
 Edot_scale=1e30,

 pa_colname='pa', e_pa_colname='e_pa', pa_outname=r"PA",
 oa_colname='oa', e_oa_colname='e_oa', oa_outname=r"OA",

 gamma_colname='gamma', gamma_outname=r'$\gamma$',
 dfil_colname='fil_dist', dfil_outname=r'$d_{\rm fil}$',

 outfile="physics_table.tex", write_fmt='aastex',
 angle_unit=u.deg, mass_unit=u.Msun, vel_unit=u.km/u.s, dist_unit=u.pc,
 momentum_unit=u.Msun*u.km/u.s, energy_unit=u.erg, time_unit=u.yr,

 dist_conversion=1):
    
    t = ascii.read(infile)

    source_col = t[source_colname].copy() 
    for source in t[source_colname]:

        if np.count_nonzero(t[source_colname] == source) == 2: 
            source_col[np.where(t[source_colname] == source)[0][-1]] = '-'

    source_col = [source.replace("davis", "SMZ").replace("hops", "HOPS") for source in source_col]

    units_dict = dict(M_outname=mass_unit, P_outname=momentum_unit, E_outname=energy_unit)

    lobe_col = t[lobe_colname].copy()
    M_col = ["{:.2f} - {:.2f}".format(
        round(a*(dist_conversion**2),2), round(b*(dist_conversion**2),2)) for (a,b) in zip(
        t[M_hivel_colname], t[M_lowvel_colname] + t[M_hivel_colname])]

    P_col = ["{:.2f} - {:.2f}".format(
        round(a*(dist_conversion**2),2), round(b*(dist_conversion**2),2)) for (a,b) in zip(
        t[P_hivel_colname], t[P_lowvel_colname] + t[P_hivel_colname])]

    E_col = ["{:.1f} - {:.1f}".format(
        round(a*(dist_conversion**2)/E_scale,2), round(b*(dist_conversion**2)/E_scale,2)) for (a,b) in zip(
        t[E_hivel_colname], t[E_lowvel_colname] + t[E_hivel_colname])]

    Mdot_col = ["{:.1f} - {:.1f}".format(
        round(a*dist_conversion/Mdot_scale,2), round(b*dist_conversion/Mdot_scale,2)) for (a,b) in zip(
        t[Mdot_hivel_colname], t[Mdot_lowvel_colname] + t[Mdot_hivel_colname])]

    Pdot_col = ["{:.1f} - {:.1f}".format(
        round(a*dist_conversion/Pdot_scale,2), round(b*dist_conversion/Pdot_scale,2)) for (a,b) in zip(
        t[Pdot_hivel_colname], t[Pdot_lowvel_colname] + t[Pdot_hivel_colname])]

    Edot_col = ["{:.1f} - {:.1f}".format(
        round(a*dist_conversion/Edot_scale,2), round(b*dist_conversion/Edot_scale,2)) for (a,b) in zip(
        t[Edot_hivel_colname], t[Edot_lowvel_colname] + t[Edot_hivel_colname])]

    vmax_col = ["{:.1f}".format(round(a, 1)) for a in t[vmax_colname]]
    rmax_col = ["{:.2f}".format(round(a*dist_conversion, 2)) for a in t[rmax_colname]]
    tdyn_col = ["{:.1f}".format(round(a*dist_conversion, 1)) for a in (t[tdyn_colname]/tdyn_scale)]

    pa_col = ["{:0g} $\\pm$ {:.2g}".format(pa, e_pa) for (pa, e_pa) in zip(
        t[pa_colname].round(0), t[e_pa_colname].round(1))]
    oa_col = ["{:0g} $\\pm$ {:.2g}".format(oa, e_oa) for (oa, e_oa) in zip(
        t[oa_colname].round(0), t[e_oa_colname].round(1))]

    gamma_col = ["{:0g}".format(gamma) for gamma in
        t[gamma_colname].round(0)]

    dfil_col = ["{:.2g}".format(round(d*dist_conversion,3)) for d in t[dfil_colname]]


    # t_out.show_in_browser()
    # # raise

    # print(t_out)

    cols = [source_col, lobe_col,
        M_col, P_col, E_col,
        rmax_col, vmax_col, tdyn_col,
        Mdot_col, Pdot_col, Edot_col,
        #pa_col, oa_col
        #, gamma_col, dfil_col
        ]



    colnames = [source_outname, lobe_outname,
        M_outname, P_outname, E_outname,
        rmax_outname, vmax_outname, tdyn_outname,
        Mdot_outname, Pdot_outname, Edot_outname,
        #pa_outname, oa_outname
        # , gamma_outname, dfil_outname
        ]


    t_out = Table(cols, names=colnames)
    t_out.show_in_browser()


    ascii.write(t_out, outfile, format=write_fmt,
        latexdict = {'tabletype':'deluxetable*'})



    pass

def outflow_table(infile="physics_outflows.txt", outfile="outflow_table.tex",
 t_outflow_hops="outflow_hops.csv",
    t_outflow_davis='outflow_davisnohops.csv',
    source_colname='source', source_outname="Source",
    conf_colname='confidence'):
    from astropy.coordinates import SkyCoord
    from spectral_cube import SpectralCube
    from stamp import extract_subcube, calc_linewings
    from regions import CircleSkyRegion
    c12 = SpectralCube.read("../cubes/mask_imfit_12co_pix_2_Tmb.fits")


    t = ascii.read(infile)
    t_hops = ascii.read(t_outflow_hops)
    t_davis = ascii.read(t_outflow_davis)
    
    lookup = set()
    unique_sources = [x for x in t['source'] if x not in lookup and lookup.add(x) is None]
    source_col = [source.replace("davis", "SMZ").replace("hops", "HOPS") for source in unique_sources]

    vbluered_col = []
    tanabe_col = []
    conf_col = []
    ra_col = []
    dec_col = []
    for source in unique_sources:

        source_n = int(source.split()[1])

        c = SkyCoord(t['RAJ2000'][t['source'] == source][0], t['DEJ2000'][t['source'] == source][0], unit=u.deg)
        # c = SkyCoord(t['RAJ2000'], t['DEJ2000'], unit=u.deg)
        ra = c.ra.to(u.hourangle).to_string(format='latex', precision=2)  
        dec = c.dec.to(u.deg).to_string(format='latex', precision=2)  

        if "hops" in source:
            vblue, vred = t_hops['blue_vel'][t_hops['hops'] == source_n][0], t_hops['red_vel'][t_hops['hops'] == source_n][0]
            tanabe = str(int(t_hops['tanabe'][t_hops['hops'] == source_n]))

        elif "davis" in source:
            vblue, vred = t_davis['blue_vel'][t_davis['davis'] == source_n][0], t_davis['red_vel'][t_davis['davis'] == source_n][0]
            tanabe = str(int(t_davis['tanabe'][t_davis['davis'] == source_n]))

        if vblue == '-':
            fit_12 = extract_subcube(c12, CircleSkyRegion,
                                    region_kwargs={'center':c, 'radius':15*u.arcsec})

            vblue, vred = calc_linewings(fit_12, nsigma_vel=2.)
            vblue, vred = np.round(vblue.to(u.km/u.s).value, 1), np.round(vred.to(u.km/u.s).value, 1)
            print(vblue, vred)

        if tanabe == '-1':
            tanabe = ''
        conf = t['confidence'][t['source'] == source]

        print(conf)
        if len(conf) == 2:
            conf_str = "{}/{}".format(conf[0],conf[1])
        elif t['lobe'][t['source'] == source] == 'R':
            vblue = "-"
            conf_str = "-/{}".format(conf[0])
        elif t['lobe'][t['source'] == source] == 'B':
            vred = "-"
            conf_str = "{}/-".format(conf[0])



        vbluered = "{}/{}".format(vblue, vred)
        vbluered_col = np.append(vbluered_col, vbluered)
        tanabe_col = np.append(tanabe_col, tanabe)
        conf_col = np.append(conf_col, conf_str)
        ra_col = np.append(ra_col, ra)
        dec_col = np.append(dec_col, dec)



    cols = [source_col, ra_col, dec_col, vbluered_col, conf_col, tanabe_col]
    colnames = ["Source", "R.A.", "Decl.", r"$v_{\rm blue}$/$v_{\rm red}$", "Confidence", "Tanabe"]

    t_out = Table(cols, names=colnames)
    t_out.show_in_browser()

    ascii.write(t_out, outfile, format='aastex',
        latexdict = {'tabletype':'deluxetable*'})

    pass



def filament_table(infile="physics_outflows.txt",
 source_colname='source', source_outname="Source",
 lobe_colname='lobe', lobe_outname="Lobe",
 pa_colname='pa', e_pa_colname='e_pa', pa_outname=r"PA",
 oa_colname='oa', e_oa_colname='e_oa', oa_outname=r"OA",

 gamma_colname='gamma', gamma_outname=r'$\gamma$',
 dfil_colname='fil_dist', dfil_outname=r'$d_{\rm fil}$',

 outfile="filament_table.tex", write_fmt='aastex',
 angle_unit=u.deg, dist_unit=u.pc,
 dist_conversion=1):
    
    t = ascii.read(infile)


    lookup = set()
    unique_sources = [x for x in t['source'] if x not in lookup and lookup.add(x) is None]
    source_col = [source.replace("davis", "SMZ").replace("hops", "HOPS") for source in unique_sources]


    pa_col = []
    oa_col = []
    gamma_col = []
    dfil_col = []

    for source in unique_sources:

        source_n = int(source.split()[1])
        ii_source = t['source'] == source
        lobes = t[lobe_colname][ii_source]

        print(lobes)
        if len(lobes) == 2:
            pa_pair = t[pa_colname][ii_source]
            e_pa_pair = t[e_pa_colname][ii_source]
            oa_pair = t[oa_colname][ii_source]
            e_oa_pair = t[e_oa_colname][ii_source]
            gamma_pair = t[gamma_colname][ii_source]

            pa = "${:0g} \pm {:.2g}$/${:0g} \pm {:.2g}$".format(pa_pair[0].round(0), e_pa_pair[0].round(1), pa_pair[1].round(0), e_pa_pair[1].round(1))
            oa = "${:0g} \pm {:.2g}$/${:0g} \pm {:.2g}$".format(oa_pair[0].round(0), e_oa_pair[0].round(1), oa_pair[1].round(0), e_oa_pair[1].round(1))
            gamma = "{:0g}/{:0g}".format(gamma_pair[0].round(0), gamma_pair[1].round(0))
            pa_col = np.append(pa_col, pa)
            oa_col = np.append(oa_col, oa)
            gamma_col = np.append(gamma_col, gamma)

            # dfil_colname = t[dfil_colname][ii_source]
            # conf_str = "{}/{}".format(conf[0],conf[1])
        elif t['lobe'][t['source'] == source] == 'B':
            pa = t[pa_colname][ii_source][0]
            print(ii_source)
            e_pa = t[e_pa_colname][ii_source][0]
            oa = t[oa_colname][ii_source][0]
            e_oa = t[e_oa_colname][ii_source][0]
            gamma = t[gamma_colname][ii_source][0]

            pa = "${:0g} \pm {:.2g}$/-".format(pa.round(0), e_pa.round(1))
            oa = "${:0g} \pm {:.2g}$/-".format(oa.round(0), e_oa.round(1))
            gamma = "${:0g}$/-".format(gamma.round(0))
            pa_col = np.append(pa_col, pa)
            oa_col = np.append(oa_col, oa)
            gamma_col = np.append(gamma_col, gamma)
            
        elif t['lobe'][t['source'] == source] == 'R':
            pa = t[pa_colname][ii_source][0]
            e_pa = t[e_pa_colname][ii_source][0]
            oa = t[oa_colname][ii_source][0]
            e_oa = t[e_oa_colname][ii_source][0]
            gamma = t[gamma_colname][ii_source][0]

            pa = "-/${:0g} \pm {:.2g}$".format(pa.round(0), e_pa.round(1))
            oa = "-/${:0g} \pm {:.2g}$".format(oa.round(0), e_oa.round(1))
            gamma = "-/${:0g}$".format(gamma.round(0))
            pa_col = np.append(pa_col, pa)
            oa_col = np.append(oa_col, oa)
            gamma_col = np.append(gamma_col, gamma)


        dfil = "{:.2g}".format(round(t[dfil_colname][ii_source][0]*dist_conversion,3))
        dfil_col = np.append(dfil_col, dfil)


    # t_out.show_in_browser()
    # # raise

    # print(t_out)

    cols = [source_col,
        pa_col, oa_col
        , gamma_col, dfil_col
        ]



    colnames = [source_outname,
        pa_outname, oa_outname
        , gamma_outname, dfil_outname
        ]


    t_out = Table(cols, names=colnames)
    t_out.show_in_browser()


    ascii.write(t_out, outfile, format=write_fmt,
        latexdict = {'tabletype':'deluxetable*'})



    pass

# def filament_table():
if __name__ == '__main__':
    main()