import sys
import copy
import matplotlib.pyplot as plt
import os
import numpy as np
from astropy import units as u
from astropy.analytic_functions import blackbody_lambda, blackbody_nu
from matplotlib import rc
rc('text', usetex=True)
font = {'weight' : 'normal','size':16,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)

def outfilpa(rawpa,filpa):
    return min(abs(rawpa-filpa),180.-abs(rawpa-filpa))

import outflow_catalog
rawoutflows = outflow_catalog.outflows
outflows = copy.deepcopy(rawoutflows)
for oo in outflows.keys():
    if 'blue_1' in outflows[oo].keys():
        multioutcore = outflows.pop(oo)
        outflows[oo+'_p1'] = {} 
        outflows[oo+'_p1']['filament'] = multioutcore['filament'] 
        outflows[oo+'_p1']['blue'] = multioutcore['blue_1'] 
        outflows[oo+'_p1']['red'] = multioutcore['red_1'] 
        outflows[oo+'_p2'] = {} 
        outflows[oo+'_p2']['filament'] = multioutcore['filament'] 
        outflows[oo+'_p2']['blue'] = multioutcore['blue_2'] 
        outflows[oo+'_p2']['red'] = multioutcore['red_2'] 

palist_blue = []
palist_red = []
palist = []
fil_palist_blue = []
fil_palist_red = []
fil_palist = []
filapalist = []
blueco = 0
redco = 0
bluesio = 0
redsio = 0
pairco = 0
pairsio = 0
totalcocores = len(rawoutflows.keys())
totalsiocores = 0

for oo in rawoutflows.keys():
    if 'blue' in rawoutflows[oo].keys():
        if 'SiO(5-4)' in rawoutflows[oo]['blue']: 
            totalsiocores += 1 
            continue 
    if 'red' in rawoutflows[oo].keys():
        if 'SiO(5-4)' in rawoutflows[oo]['red']: 
            totalsiocores += 1 
            continue 
    if 'blue_1' in rawoutflows[oo].keys():
        if 'SiO(5-4)' in rawoutflows[oo]['blue_1']: 
            totalsiocores += 1 
            continue 
    if 'red_1' in rawoutflows[oo].keys():
        if 'SiO(5-4)' in rawoutflows[oo]['red_1']: 
            totalsiocores += 1 
            continue 
    if 'blue_2' in rawoutflows[oo].keys():
        if 'SiO(5-4)' in rawoutflows[oo]['blue_2']: 
            totalsiocores += 1 
            continue 
    if 'red_2' in rawoutflows[oo].keys():
        if 'SiO(5-4)' in rawoutflows[oo]['red_2']: 
            totalsiocores += 1 
            continue 

tablecontents = []
for oo in outflows.keys():
    ### single bipolar
    pa_blue = np.nan
    pa_red = np.nan
    averpa = np.nan
    fil_pa = np.nan
    filpa = np.nan
    filament = outflows[oo]['filament']
    pairsioyesorno = 0
    bluesiotag = 'N'
    redsiotag = 'N'
    table_pa_blue = '-' 
    table_pa_red = '-' 
    table_fil_pa = '-' 
    table_filpa = '-' 
    table_averpa = '-' 
    rra = '' 
    ddec = '' 
    if filament['pa'] != 'deg': 
        filpa = float(filament['pa'][:-3]) - 90.
        if filpa < 0: 
            filpa += 180. 
        if filpa < 0 or filpa > 180: 
            print 'tell me why, filpa not in (0,180)' 
        table_filpa = str(int(round(filpa)))
        filapalist.append(table_filpa)
    if 'blue' in outflows[oo].keys():
        rra = '%.5f' % outflows[oo]['blue']['12CO(2-1)']['ra']
        ddec = '%.5f' % outflows[oo]['blue']['12CO(2-1)']['dec']
        rawpa_blue = outflows[oo]['blue']['12CO(2-1)']['pa']
        pa_blue = float(rawpa_blue[:-3]) - 90. # ds9 pa - 90deg 
        if pa_blue > 180:
            pa_blue = pa_blue - 360. 
        palist_blue.append(pa_blue)
        table_pa_blue = str(int(round(pa_blue))) 
        if pa_blue < -180 or pa_blue > 180:
            print 'tell me why, pa_blue not in (-180,180)',oo
            sys.exit() 
        if filament['pa'] != 'deg': 
            if pa_blue < 0: 
                fil_pa_blue = outfilpa(pa_blue+180.,filpa)
            else:
                fil_pa_blue = outfilpa(pa_blue,filpa)
            if fil_pa_blue < 0 or fil_pa_blue > 90: 
                print 'tell me why, fil_pa_blue not in (0,90)',oo
                sys.exit() 
            fil_palist_blue.append(fil_pa_blue)
        blueco += 1 # because if there is an outflow, it must be seen in CO 
        if 'SiO(5-4)' in outflows[oo]['blue'].keys():
            bluesio += 1 
            bluesiotag = 'Y' 
            pairsioyesorno += 1 
    if 'red' in outflows[oo].keys():
        rra = '%.5f' % outflows[oo]['red']['12CO(2-1)']['ra']
        ddec = '%.5f' % outflows[oo]['red']['12CO(2-1)']['dec']
        rawpa_red = outflows[oo]['red']['12CO(2-1)']['pa']
        pa_red = float(rawpa_red[:-3]) - 90. # ds9 pa - 90deg 
        if pa_red > 180:
            pa_red = pa_red - 360. 
        palist_red.append(pa_red)
        table_pa_red = str(int(round(pa_red))) 
        if pa_red < -180 or pa_red > 180:
            print 'tell me why, pa_red not in (-180,180)',oo
            sys.exit() 
        if filament['pa'] != 'deg': 
            if pa_red < 0: 
                fil_pa_red = outfilpa(pa_red+180.,filpa)
            else: 
                fil_pa_red = outfilpa(pa_red,filpa)
            if fil_pa_red < 0 or fil_pa_red > 90: 
                print 'tell me why, fil_pa_red not in (0,90)',oo
                sys.exit() 
            fil_palist_red.append(fil_pa_red)
        redco += 1 # because if there is an outflow, it must be seen in CO 
        if 'SiO(5-4)' in outflows[oo]['red'].keys():
            redsio += 1 
            redsiotag = 'Y' 
            pairsioyesorno += 1 
    if pairsioyesorno == 2:
        pairsio += 1
    if 'blue' in outflows[oo].keys() and 'red' in outflows[oo].keys() and ~np.isnan(pa_blue) and ~np.isnan(pa_red):
        averpa = (pa_blue + pa_red + 180.)/2.
        if oo == 'core_7': averpa += 180. # core_7 is special
        if averpa > 180 or averpa < 0:  
            print 'tell me why, averpa not in (0,180)',oo 
            print 'averpa',averpa,'pa_blue',pa_blue,'pa_red',pa_red 
            sys.exit() 
        palist.append(averpa) 
        table_averpa = str(int(round(averpa))) 
        if filament['pa'] != 'deg': 
            fil_pa = outfilpa(averpa,filpa)
            if fil_pa < 0 or fil_pa > 90: 
                print 'tell me why, fil_pa not in (0,90)',oo
                print 'fil_pa',fil_pa,'averpa',averpa,'pa_blue',pa_blue,'pa_red',pa_red 
                sys.exit() 
            fil_palist.append(fil_pa)
            table_fil_pa = str(int(round(fil_pa)))
        pairco += 1 # because if there is an outflow, it must be seen in CO 
    elif 'blue' in outflows[oo].keys() and ~np.isnan(pa_blue):
        if pa_blue < 0: 
            averpa = pa_blue + 180. 
        else: 
            averpa = pa_blue 
        palist.append(averpa) 
        table_averpa = str(int(round(averpa))) 
        if filament['pa'] != 'deg': 
            fil_pa = outfilpa(averpa,filpa)
            if fil_pa < 0 or fil_pa > 90: 
                print 'tell me why, fil_pa not in (0,90)',oo
                print 'fil_pa',fil_pa,'averpa',averpa,'pa_blue',pa_blue,'pa_red',pa_red 
                sys.exit() 
            fil_palist.append(fil_pa)
            table_fil_pa = str(int(round(fil_pa)))
    elif 'red' in outflows[oo].keys() and ~np.isnan(pa_red):
        if pa_red < 0: 
            averpa = pa_red + 180. 
        else: 
            averpa = pa_red 
        palist.append(averpa) 
        table_averpa = str(int(round(averpa))) 
        if filament['pa'] != 'deg': 
            fil_pa = outfilpa(averpa,filpa)
            if fil_pa < 0 or fil_pa > 90: 
                print 'tell me why, fil_pa not in (0,90)',oo
                print 'fil_pa',fil_pa,'averpa',averpa,'pa_blue',pa_blue,'pa_red',pa_red 
                sys.exit() 
            fil_palist.append(fil_pa)
            table_fil_pa = str(int(round(fil_pa)))
    tablecontents.append([oo.split('_')[1],oo,rra,ddec,table_pa_blue,bluesiotag,table_pa_red,redsiotag,table_averpa,table_filpa,table_fil_pa])

### AD test of PA_Out
from skgof import ad_test
from scipy.stats import uniform
from scipy import stats
st, pvalue = ad_test(palist,uniform(0,180))
print 'PA_Out AD test, statistic, p-value', st, pvalue
st, cv, sl = stats.anderson_ksamp([palist,np.random.uniform(0,180.,10000000)]) 
print 'PA_Out AD test 2sample, random, st, cv, sl',st,cv,sl 
##

### PA_Out distribution
median_pa = np.nanmedian(palist)
pa_low = np.nanmin(palist)
pa_high = np.nanmax(palist)
print 'pa_low,pa_high',pa_low,pa_high
print 'median_pa',median_pa
hist, bin_edges = np.histogram(palist,bins=[0,20,40,60,80,100,120,140,160,180],range=(0,180))
bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
print 'hist',hist
print 'bincenter',bincenter
hist_blue, bin_edges_blue = np.histogram(palist_blue,bins=range(-180,180,20),range=(-180,180))
bincenter_blue = (bin_edges_blue[:-1] + bin_edges_blue[1:]) / 2.
hist_red, bin_edges_red = np.histogram(palist_red,bins=range(-180,180,20),range=(-180,180))
bincenter_red = (bin_edges_red[:-1] + bin_edges_red[1:]) / 2.
# filament 2D histogram
hist_fil, bin_edges_fil = np.histogram(filapalist,bins=range(0,181,20),range=(0,180))
bincenter_fil = (bin_edges_fil[:-1] + bin_edges_fil[1:]) / 2.

### gamma distribution
fil_median_pa = np.nanmedian(fil_palist)
fil_pa_low = np.nanmin(fil_palist)
fil_pa_high = np.nanmax(fil_palist)
print 'fil_pa_low,fil_pa_high',fil_pa_low,fil_pa_high
print 'fil_median_pa',fil_median_pa
fil_hist, fil_bin_edges = np.histogram(fil_palist,bins=[0,10,20,30,40,50,60,70,80,90],range=(0,90))
fil_bincenter = (fil_bin_edges[:-1] + fil_bin_edges[1:]) / 2.
print 'fil_hist',fil_hist
print 'fil_bincenter',fil_bincenter
fil_hist_blue, fil_bin_edges_blue = np.histogram(fil_palist_blue,bins=[0,10,20,30,40,50,60,70,80,90],range=(0,90))
fil_bincenter_blue = (fil_bin_edges_blue[:-1] + fil_bin_edges_blue[1:]) / 2.
fil_hist_red, fil_bin_edges_red = np.histogram(fil_palist_red,bins=[0,10,20,30,40,50,60,70,80,90],range=(0,90))
fil_bincenter_red = (fil_bin_edges_red[:-1] + fil_bin_edges_red[1:]) / 2.

print 'blueco',blueco
print 'redco',redco
print 'bluesio',bluesio
print 'redsio',redsio
print 'total co outflow cores',totalcocores
print 'total sio outflow cores',totalsiocores
print 'total co outflow pairs',pairco
print 'total sio outflow pairs',pairsio

outputtable = 0
plothist = 0
plotfilhist = 0
montecarlo = 1

if outputtable == 1:
    numtablecontents = []
    utablecontents = []
    for ii in tablecontents:
        if ii[0][0] == 'u': 
            ii[0] = int(ii[0][1:]) 
            utablecontents.append(ii) 
        else: 
            ii[0] = int(ii[0]) 
            numtablecontents.append(ii) 
    for ii in sorted(numtablecontents, key=lambda x: x[0]):
        print ii[1][5:].rjust(10) + ' & ' + ii[2].rjust(10) + ' & ' + ii[3].rjust(10) + ' & ' + ii[4].rjust(5) + ' & ' + ii[5].rjust(3) + ' & ' + ii[6].rjust(5) + ' & ' + ii[7].rjust(3) + ' & ' + ii[8].rjust(5) + ' & ' + ii[9].rjust(5) + ' & ' + ii[10].rjust(5) + r' \\'
    for ii in sorted(utablecontents, key=lambda x: x[0]):
        print ii[1][5:].rjust(10) + ' & ' + ii[2].rjust(10) + ' & ' + ii[3].rjust(10) + ' & ' + ii[4].rjust(5) + ' & ' + ii[5].rjust(3) + ' & ' + ii[6].rjust(5) + ' & ' + ii[7].rjust(3) + ' & ' + ii[8].rjust(5) + ' & ' + ii[9].rjust(5) + ' & ' + ii[10].rjust(5) + r' \\'

if plothist == 1:
    plotdata = {
               'panel2':{
                         'line1':{'xdata':bincenter,'ydata':hist,'color':'k','width':18,'label':r'','title':''},
                         'xlim':(0,180),'ylim':[0,20],'xscale':'linear','yscale':'linear','xlabel':r'$\rm PA_{\rm Out}~(deg)$','ylabel':r'number','loc':1, 
                        },
               'panel1':{
                         'line1':{'xdata':bincenter_blue,'ydata':hist_blue,'color':'b','width':12,'label':r'blue','title':''},
                         'line2':{'xdata':bincenter_red,'ydata':hist_red,'color':'r','width':18,'label':r'red','title':''},
                         'xlim':(-180,180),'ylim':[0,15],'xscale':'linear','yscale':'linear','xlabel':r'$\rm PA~(deg)$','ylabel':r'number','loc':1,
                         #'scatter1':{'xdata':raw_coremasslist,'ydata':tkinlist,'color':'k','label':'','title':r'${\rm T}_{\rm gas}-M$'},
                        },
               'panel4':{
                         'line1':{'xdata':fil_bincenter,'ydata':fil_hist,'color':'k','width':8,'label':r'','title':''},
                         'xlim':(0,90),'ylim':[0,18],'xscale':'linear','yscale':'linear','xlabel':r'$\rm \gamma~(deg)$','ylabel':r'number','loc':1, 
                        },
               'panel3':{
                         'line1':{'xdata':fil_bincenter_blue,'ydata':fil_hist_blue,'color':'b','width':5,'label':r'blue','title':''},
                         'line2':{'xdata':fil_bincenter_red,'ydata':fil_hist_red,'color':'r','width':8,'label':r'red','title':''},
                         'xlim':(0,90),'ylim':[0,19],'xscale':'linear','yscale':'linear','xlabel':r'$\rm PA~(deg)$','ylabel':r'number','loc':1, 
                         #'scatter1':{'xdata':raw_coremasslist,'ydata':tkinlist,'color':'k','label':'','title':r'${\rm T}_{\rm gas}-M$'},
                        },
               }
    xpanels = 2
    ypanels = 2
    
    lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    
    pdfname = 'pahist.pdf'
    p=plt.figure(figsize=(1+6*xpanels,6*ypanels))
    plt.subplots_adjust(top=0.96,bottom=0.1,left=0.15,right=0.96)
    for i in range(xpanels):
        for j in range(ypanels):
            panel = i*xpanels+j+1
            ax=p.add_subplot(ypanels,xpanels,panel)
            ax.set_xscale(plotdata['panel'+str(panel)]['xscale'])
            ax.set_yscale(plotdata['panel'+str(panel)]['yscale'])
            xmin,xmax = plotdata['panel'+str(panel)]['xlim']
            ax.set_xlim(xmin,xmax)
            if plotdata['panel'+str(panel)]['ylim'] != []:
                ymin,ymax = plotdata['panel'+str(panel)]['ylim']
                ax.set_ylim(ymin,ymax)
            ax.text(0.05, 0.95,'('+lletter[panel-1]+')',horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
            ax.set_ylabel(plotdata['panel'+str(panel)]['ylabel'])
            ax.set_xlabel(plotdata['panel'+str(panel)]['xlabel'])
            #ax.set_title(r'$\rm F_{1.3mm}-\Sigma~(MIR~source~removed)$')
            for lline in filter(lambda x: x[:4] == 'line', plotdata['panel'+str(panel)].keys()): 
                ax.bar(plotdata['panel'+str(panel)][lline]['xdata'], plotdata['panel'+str(panel)][lline]['ydata'], plotdata['panel'+str(panel)][lline]['width'], align='center', color=plotdata['panel'+str(panel)][lline]['color'], label=plotdata['panel'+str(panel)][lline]['label'])
            for sscatter in filter(lambda x: x[:7] == 'scatter', plotdata['panel'+str(panel)].keys()): 
                ax.scatter(plotdata['panel'+str(panel)][sscatter]['xdata'], plotdata['panel'+str(panel)][sscatter]['ydata'], s=50, c=plotdata['panel'+str(panel)][sscatter]['color'], edgecolors=plotdata['panel'+str(panel)][sscatter]['color'], label=plotdata['panel'+str(panel)][sscatter]['label'])
            ax.legend(frameon=False,loc=plotdata['panel'+str(panel)]['loc'])
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    plt.close(p)
    os.system('cp '+pdfname+os.path.expandvars(' ${DROPATH}/cloudc1'))

if plotfilhist == 1:
    plotdata = {
               'panel1':{
                         'line1':{'xdata':bincenter_fil,'ydata':hist_fil,'color':'k','width':18,'label':r'filament','title':''},
                         'xlim':(0,180),'ylim':[0,30],'xscale':'linear','yscale':'linear','xlabel':r'$\rm PA_{\rm Fil}~(deg)$','ylabel':r'number','loc':1,
                         #'scatter1':{'xdata':raw_coremasslist,'ydata':tkinlist,'color':'k','label':'','title':r'${\rm T}_{\rm gas}-M$'},
                        },
               }
    xpanels = 1
    ypanels = 1
    
    lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    
    pdfname = 'pafilhist.pdf'
    p=plt.figure(figsize=(1+6*xpanels,6*ypanels))
    plt.subplots_adjust(top=0.96,bottom=0.1,left=0.15,right=0.96)
    for i in range(xpanels):
        for j in range(ypanels):
            panel = i*xpanels+j+1
            ax=p.add_subplot(ypanels,xpanels,panel)
            ax.set_xscale(plotdata['panel'+str(panel)]['xscale'])
            ax.set_yscale(plotdata['panel'+str(panel)]['yscale'])
            xmin,xmax = plotdata['panel'+str(panel)]['xlim']
            ax.set_xlim(xmin,xmax)
            if plotdata['panel'+str(panel)]['ylim'] != []:
                ymin,ymax = plotdata['panel'+str(panel)]['ylim']
                ax.set_ylim(ymin,ymax)
            ax.text(0.05, 0.95,'('+lletter[panel-1]+')',horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
            ax.set_ylabel(plotdata['panel'+str(panel)]['ylabel'])
            ax.set_xlabel(plotdata['panel'+str(panel)]['xlabel'])
            #ax.set_title(r'$\rm F_{1.3mm}-\Sigma~(MIR~source~removed)$')
            for lline in filter(lambda x: x[:4] == 'line', plotdata['panel'+str(panel)].keys()): 
                ax.bar(plotdata['panel'+str(panel)][lline]['xdata'], plotdata['panel'+str(panel)][lline]['ydata'], plotdata['panel'+str(panel)][lline]['width'], align='center', color=plotdata['panel'+str(panel)][lline]['color'], label=plotdata['panel'+str(panel)][lline]['label'])
            for sscatter in filter(lambda x: x[:7] == 'scatter', plotdata['panel'+str(panel)].keys()): 
                ax.scatter(plotdata['panel'+str(panel)][sscatter]['xdata'], plotdata['panel'+str(panel)][sscatter]['ydata'], s=50, c=plotdata['panel'+str(panel)][sscatter]['color'], edgecolors=plotdata['panel'+str(panel)][sscatter]['color'], label=plotdata['panel'+str(panel)][sscatter]['label'])
            ax.legend(frameon=False,loc=plotdata['panel'+str(panel)]['loc'])
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    plt.close(p)
    os.system('cp '+pdfname+os.path.expandvars(' ${DROPATH}/cloudc1'))

if montecarlo == 1: # Follow Stephens et al. (2017)
    nn = 1000000
    # first vector
    uu1 = np.random.uniform(-1,1,nn) # 
    theta1 = np.random.uniform(0,2*np.pi,nn) #
    x1 = (1. - uu1**2)**0.5 * np.cos(theta1)
    y1 = (1. - uu1**2)**0.5 * np.sin(theta1)
    z1 = uu1
    # second vector
    uu2 = np.random.uniform(-1,1,nn) # 
    theta2 = np.random.uniform(0,2*np.pi,nn) #
    x2 = (1. - uu2**2)**0.5 * np.cos(theta2)
    y2 = (1. - uu2**2)**0.5 * np.sin(theta2)
    z2 = uu2
    # gamma_3D in Stephens et al. (2017)
    gamma3D = np.arccos(x1 * x2 + y1 * y2 + z1 * z2) / np.pi * 180. 
    new_gamma3D = [ii if ii<=90. else 180.-ii for ii in gamma3D]
    print 'min(new_gamma3D)',min(new_gamma3D),'max(new_gamma3D)',max(new_gamma3D) 
    hist, bin_edges = np.histogram(new_gamma3D,bins=np.arange(0,91,1),range=(0,90))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
    print 'hist, bin_edges, bin_centers', hist, bin_edges, bin_centers
    # gamma_2D in Stephens et al. (2017)
    gamma2D = np.arccos((y1 * y2 + z1 * z2)/(y1**2+z1**2)**0.5/(y2**2+z2**2)**0.5) / np.pi * 180. 
    new_gamma2D = [ii if ii<=90. else 180.-ii for ii in gamma2D]
    print 'min(new_gamma2D)',min(new_gamma2D),'max(new_gamma2D)',max(new_gamma2D) 
#    stephens17 = np.loadtxt('stephens_table2.txt') ### test Stephens17 table 3, his AD test is 2-sample test using the randomly generated data
#    st, cv, sl = stats.anderson_ksamp([stephens17[:,3],new_gamma2D]) ### test Stephens17 table 3, his AD test is 2-sample test using the randomly generated data
    st, cv, sl = stats.anderson_ksamp([fil_palist,new_gamma2D]) 
    print 'random, st, cv, sl',st,cv,sl 
    # gamma_2D random histogram
    hist2D, bin_edges2D = np.histogram(new_gamma2D,bins=np.arange(0,91,1),range=(0,90),normed=True)
    bin_centers2D = (bin_edges2D[:-1] + bin_edges2D[1:]) / 2.
    # gamma_2D (0,20) degree histogram
    new_gamma2D_parallel = [ii for cc,ii in enumerate(new_gamma2D) if new_gamma3D[cc]<=20.]
    print 'no need to test parallel case, it is obviously false in the plot'
#    st, cv, sl = stats.anderson_ksamp([fil_palist,new_gamma2D_parallel]) 
#    print 'parallel, st, cv, sl',st,cv,sl 
    hist2D_parallel, bin_edges2D_parallel = np.histogram(new_gamma2D_parallel,bins=np.arange(0,91,1),range=(0,90),normed=True)
    bin_centers2D_parallel = (bin_edges2D_parallel[:-1] + bin_edges2D_parallel[1:]) / 2.
    # gamma_2D (70,90) degree histogram
    new_gamma2D_perpendicular = [ii for cc,ii in enumerate(new_gamma2D) if new_gamma3D[cc]>70.]
    st, cv, sl = stats.anderson_ksamp([fil_palist,new_gamma2D_perpendicular]) 
    print 'perpendicular, st, cv, sl',st,cv,sl 
    hist2D_perpendicular, bin_edges2D_perpendicular = np.histogram(new_gamma2D_perpendicular,bins=np.arange(0,91,1),range=(0,90),normed=True)
    bin_centers2D_perpendicular = (bin_edges2D_perpendicular[:-1] + bin_edges2D_perpendicular[1:]) / 2.
#    st, cv, sl = stats.anderson_ksamp([stephens17[:,3],new_gamma2D_perpendicular]) 
#    print 'perpendicular, st, cv, sl',st,cv,sl 
#    sys.exit() ### end of AD test
    # gamma_F histogram
    outfil_hist, outfil_bin_edges = np.histogram(fil_palist,bins=np.arange(0,91,1),range=(0,90),normed=True)
    outfil_bincenters = (outfil_bin_edges[:-1] + outfil_bin_edges[1:]) / 2.

    # plot data
    plotdata = {
               'panel1':{
                         'line1':{'xdata':bin_centers,'ydata':hist,'linestyle':'k-','linewidth':1,'drawsty':'steps-mid','label':'','title':''},
                         'xlim':(0,90),'ylim':[0,18000],'xscale':'linear','yscale':'linear','xlabel':r'$\rm \gamma_{\rm 3D}~(deg)$','ylabel':r'Number of Angles','loc':2, 
                        },
               'panel2':{
                         'line1':{'xdata':bin_centers2D,'ydata':np.cumsum(hist2D),'linestyle':'k-','linewidth':1,'drawsty':'default','label':r'$\rm random$','title':''},
                         'line2':{'xdata':bin_centers2D_parallel,'ydata':np.cumsum(hist2D_parallel),'linestyle':'k--','linewidth':1,'drawsty':'default','label':r'$\rm 0-20^\circ$','title':''},
                         'line3':{'xdata':bin_centers2D_perpendicular,'ydata':np.cumsum(hist2D_perpendicular),'linestyle':'k-.','linewidth':1,'drawsty':'default','label':r'$\rm 70-90^\circ$','title':''},
                         'line4':{'xdata':outfil_bincenters,'ydata':np.cumsum(outfil_hist),'linestyle':'r-','linewidth':3,'drawsty':'default','label':r'$\rm Observation$','title':''},
                         'xlim':(0-1,90+1),'ylim':[0-0.01,1+0.01],'xscale':'linear','yscale':'linear','xlabel':r'$\rm \gamma~(deg)$','ylabel':r'Cumulative Distribution','loc':4,
                         #'scatter1':{'xdata':raw_coremasslist,'ydata':tkinlist,'color':'k','label':'','title':r'${\rm T}_{\rm gas}-M$'},
                        },
               }
    xpanels = 1
    ypanels = 2
    
    lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    
    pdfname = 'montecarlo_outfil_3D_2D_projection.pdf'
    p=plt.figure(figsize=(1+6*xpanels,6*ypanels))
    plt.subplots_adjust(top=0.96,bottom=0.1,left=0.15,right=0.96)
    for i in range(xpanels):
        for j in range(ypanels):
            panel = i*xpanels+j+1
            ax=p.add_subplot(ypanels,xpanels,panel)
            ax.set_xscale(plotdata['panel'+str(panel)]['xscale'])
            ax.set_yscale(plotdata['panel'+str(panel)]['yscale'])
            xmin,xmax = plotdata['panel'+str(panel)]['xlim']
            ax.set_xlim(xmin,xmax)
            if plotdata['panel'+str(panel)]['ylim'] != []:
                ymin,ymax = plotdata['panel'+str(panel)]['ylim']
                ax.set_ylim(ymin,ymax)
            ax.text(0.05, 0.95,'('+lletter[panel-1]+')',horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
            ax.set_ylabel(plotdata['panel'+str(panel)]['ylabel'])
            ax.set_xlabel(plotdata['panel'+str(panel)]['xlabel'])
            #ax.set_title(r'$\rm F_{1.3mm}-\Sigma~(MIR~source~removed)$')
            for lline in filter(lambda x: x[:4] == 'line', plotdata['panel'+str(panel)].keys()): 
                ax.plot(plotdata['panel'+str(panel)][lline]['xdata'], plotdata['panel'+str(panel)][lline]['ydata'], plotdata['panel'+str(panel)][lline]['linestyle'], linewidth=plotdata['panel'+str(panel)][lline]['linewidth'], label=plotdata['panel'+str(panel)][lline]['label'],drawstyle=plotdata['panel'+str(panel)][lline]['drawsty'])
            for sscatter in filter(lambda x: x[:7] == 'scatter', plotdata['panel'+str(panel)].keys()): 
                ax.scatter(plotdata['panel'+str(panel)][sscatter]['xdata'], plotdata['panel'+str(panel)][sscatter]['ydata'], s=50, c=plotdata['panel'+str(panel)][sscatter]['color'], edgecolors=plotdata['panel'+str(panel)][sscatter]['color'], label=plotdata['panel'+str(panel)][sscatter]['label'])
            ax.legend(frameon=False,loc=plotdata['panel'+str(panel)]['loc'])
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    plt.close(p)
    os.system('cp '+pdfname+os.path.expandvars(' ${DROPATH}/cloudc1'))


