#MCLINE_FUNCTIONS.py
#Keegan Marr
#Last Modified: June 24, 2020

import warnings
warnings.filterwarnings('ignore')
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
matplotlib.rcParams['font.family'] = "sans-serif"
font_color = "black"
tick_color = "black"
import scipy.stats as sp
import scipy.constants as const   
import scipy.interpolate as interp
import scipy.integrate as integrate
import pyhdust.spectools as spc
import glob
import math
import emcee
import time
import datetime
import operator
from corner import corner
from astropy.io import fits
from astropy.convolution import Gaussian1DKernel, convolve, Box1DKernel
from scipy.interpolate import griddata, interp2d
from random import randint
from gaussfold import gaussfold as gf
from multiprocessing import Pool

#Reads in the reddened fullsed file from HDUST simulation
def readInReddeningFile(filename, luminosity, distance):
    file=open(filename, "r")

    #List stores: MU, LAMBDA, FLUX
    redParam=[[],[],[]]

    #Reads in each line of redParam from the file and appends each point to the
    #proper row in the table
    i=-1
    prevMu=-20000
    for line in file:
        line=line.split()
        if (float(line[0].strip())==prevMu):
            redParam[0][i].append(float(line[0].strip()))
            redParam[1][i].append(float(line[1].strip()))
            redParam[2][i].append(float(line[2].strip())*float(line[1].strip())*\
                10**23/10**4/2.99792458/10**10*2.99792458*10**8/\
                (float(line[1].strip()))*10**6/(3.9*10**33*luminosity)*\
                (4*np.pi*(distance*3.08568*10**18)**2)*float(line[1].strip()))

        else:
            i+=1
            for k in range(3):
                redParam[k].append([])
            redParam[0][i].append(float(line[0].strip()))
            redParam[1][i].append(float(line[1].strip()))
            redParam[2][i].append(float(line[2].strip())*float(line[1].strip())*\
                (float(line[1].strip()))*10**6/(3.9*10**33*luminosity)*\
                (4*np.pi*distance)*float(line[1].strip()))

            prevMu=redParam[0][i][0]

    file.close()

    return redParam


#Reads in observed Halpha profile obtained from Chris Tycner
def readInHalphaFile(fileName):
    file=open(fileName, "r")

    observ=[[],[]]

    for line in file:
        line=line.rsplit(" ", 1)
        observ[0].append(float(line[0].strip()))
        observ[1].append(float(line[1].strip()))

    file.close()

    return observ
   

#Reads in and plots the observed Halpha profile obtained from Rivinius
def plotRiviData():
    #REDACTED
    return wav[0], flx[0]
    

#Performs calculations needed to plot the H_alpha line
def findHAlphaLine(wavelength, flux, plotWave, normFlux, subtractCont):
    #Put arrays in order from least to greatest (if they are in inverse order)
    n=len(wavelength)
    if (wavelength[0]>wavelength[-1]):
        temp_wave=wavelength[:]
        temp_flux=flux[:]
        for i in range(n):
            wavelength[n-1-i]=temp_wave[i]
            flux[n-1-i]=temp_flux[i]

    speedOfLight=299792 #speed of light in km/s
    h_alpha_wavelength=6564.61 #wavelength of H_alpha emission line (angstroms)

    #Convert from wavelength to equivalent velocity
    velocity=list(map(lambda obs:(obs-h_alpha_wavelength)*speedOfLight/\
        h_alpha_wavelength, wavelength))

    #Find the velocities within a reasonable range of the H_alpha line
    leftmost=findExtremeVelocity(velocity, -2500)
    rightmost=findExtremeVelocity(velocity, 2500)

    #Adjust the index value of the rightmost velocity (it is not included in the
    #[x:x] operation
    if (rightmost<n):
        rightmost+=1

    #Take only the velocities within the desired range
    velocity=velocity[leftmost:rightmost+1]
    flux=flux[leftmost:rightmost+1]
    wavelength=wavelength[leftmost:rightmost+1]

    n=len(velocity)

    #For low resolution observ around the H_alpha line, make the velocity plot
    #flat (all 1's)
    if(n<3):
        if plotWave:
            return [wavelength, [1]*n]
        else:
            return [velocity, [1]*n]
    else:
        if normFlux:
            #number of points of the continuum that will be used in the fitting
            ncont=6 

            #Takes the six left/rightmost points to use as the continuum
            wave_fit=wavelength[:ncont-1]+wavelength[len(wavelength)-1-ncont:]
            flux_fit=flux[:ncont-1]+flux[len(flux)-1-ncont:rightmost-1]

            #Converts from a list to a numpy array
            wave_fit=np.array(wave_fit)
            flux_fit=np.array(flux_fit)

            #Calculates the normalized linear equation from the continuum points
            normal_eqn=np.polyfit(wave_fit, flux_fit, 1)

            #New set of normalized points
            normalized_continuum=list(map(lambda x: normal_eqn[0]*\
                x+normal_eqn[1], wavelength))

            #Normalizes the flux
            flux=list(map(lambda x,y: x/y, flux, normalized_continuum))
        
        if subtractCont:
            #This line is for comparing the fluxes of halpha lines with different
            #continuum emissions
            flux = np.array(flux) - np.array(flux[leftmost])

        if plotWave:
            return [wavelength, flux]
        else:
            return [velocity, flux]

#Used to find the left/rightmost velocity in a given range
def findExtremeVelocity(velocity, extremeValue):
    n=len(velocity)

    #current upper and lower bounds of the range
    lower=-1
    upper=n

    #Loop continues until the range consists of only one value
    while (upper-lower>1):
        #Finds the middle value in the range
        middle=(upper+lower//2)
        #If the endpoint of the range is greater than the middle value, make the
        #middle value the lower bound of the new range
        if (extremeValue>velocity[middle]):
            lower=middle
        #Otherwise make the middle value the upper bound of the new range
        else:
            upper=middle

    #Note that lower contains the INDEX value of the observ point that fits the
    #above description
    #Adjust the value of lower accordingly in these special cases so that it is
    #a valid index value
    if (lower==-1):
        lower=0
    elif (lower==n):
        lower=n-1

    #Returns the index value of the observ point
    return lower

#Organizes simulation observ acquisition, interpolation and convolution
def getSim(filename, inclAngle, ranges, luminosity, distance, fac_e, v_h,
            v_e, inclFlag, plotWave, normFlux, subtractCont, *inclSpacing):
    values=[]
    flux=[]
    flow=[]
    fhigh=[]
    flux=[]

    #returns list of mu, lambda, flux in the same units as reddened fullsed
    observ=readInReddeningFile(filename, luminosity, distance)

    inclIndex=0
    
    if inclFlag: #If there are a range of inclinations
        inclSpacing = inclSpacing[0][0][0] #remove from tuple
        lowiinrange = ranges[0][0]
        #determine what inclinations were simulated around the incoming i
        ilow = lowiinrange + math.floor((inclAngle-lowiinrange)/\
            inclSpacing)*inclSpacing
        #ilow = math.floor(inclAngle*(1./inclSpacing))/(1./inclSpacing)
        ihigh = ilow+inclSpacing
        #get the idx of the mu/incl in the fullsed
        idxlow = getMuIdx(observ[0], ilow)
        idxhigh = getMuIdx(observ[0], ihigh)
        #get the flux at these inclinations
        flow.append([observ[1][idxlow], observ[2][idxlow]])
        fhigh.append([observ[1][idxhigh], observ[2][idxhigh]])
        #interpolate the Halpha line at the desired inclination
        flux.append([observ[1][idxlow], np.ndarray.tolist(inclInterp(inclAngle,\
            [ilow,ihigh], flow, fhigh))])
        #convolve the flux of the inclination interpolated line
        for values in flux:
            (flux_conv, origData) = gaussconv(values, fac_e, v_h, v_e, plotWave,\
             normFlux, subtractCont)
        
    else: #If just one inclination is given
        inclIndex = getMuIdx(observ[0], inclAngle)
        flux.append([observ[1][inclIndex],observ[2][inclIndex]])
        #convolve the Halpha line
        for values in flux:
            (flux_conv, origData) = gaussconv(values, fac_e, v_h, v_e, plotWave,\
             normFlux, subtractCont)
    
    return values[0], flux_conv, origData

#Convolving the halpha line with two Gaussians
def gaussconv(values, fac_e, v_h, v_e, plotWave, normFlux, subtractCont):
    
    xyParams=findHAlphaLine(list(map(lambda x: x*10**4, values[0])), values[1],
                                            plotWave, normFlux, subtractCont)
    values[0]=xyParams[0] #velocity
    values[1]=xyParams[1] #flux before convolution
    
    #Adding junk ones to either side of the halpha continuum
    #other gaussfold complains about x and y not being the same
    #length, as lammin and lammax
    #Start by finding the second largest and smallest velocity to extend from
    velo2ndlarg = max(n for n in values[0] if n!=max(values[0]))
    velo2ndsmal = min(n for n in values[0] if n!=min(values[0]))
    #velocity
    values[0].extend(np.linspace(velo2ndlarg+1,24000,10000))
    values[0].extend(np.linspace(-24000, velo2ndsmal-1, 10000))
    #flux
    values[1].extend(np.ones(10000))
    values[1].extend(np.ones(10000))
    #Sort the junk values so that x is in numerically increasing order
    L = sorted(zip(values[0],values[1]), key=operator.itemgetter(0))
    values[0], values[1] = zip(*L)

    notscat_flux = [i*(1-fac_e) for i in values[1]]
    scat_flux = [i*fac_e for i in values[1]]
    notscat_conv = gf(values[0], notscat_flux, v_h)
    scat_conv = gf(values[0], scat_flux, v_e)
    flux_conv = scat_conv + notscat_conv #flux after convolution
    
    return flux_conv, xyParams

#Checks if the calculated mu angle is within the accepted variance
#and returns the index of that mu as it is in the observers list
def getMuIdx(muAngles, inclAngle):
    for muAng in muAngles:
        if (inclAngle==muAng[0]):
            inclIndex=muAngles.index(muAng)
            break
        elif (muAng==muAngles[-1]):
            print("Observer angle", inclAngle, "does not exist.")
            exit()
    return inclIndex

#Interpolates the flux over at a desired inclination
def inclInterp(inclAngle, inclRange, flow, fhigh):   
    flux = interp2d(np.arange(1., len(flow[0][1])+1), inclRange,
                    np.stack((flow[0][1], fhigh[0][1])), kind='linear')
                     
    return flux(np.arange(1., len(flow[0][1])+1), inclAngle)

#Cuts the velocity range to only be within a desired window
def veloWindow(velo, flux, interpFlux, veloUBound, veloLBound):
    idxs = []
    idxs.extend(np.argwhere(velo > veloUBound))
    idxs.extend(np.argwhere(velo < veloLBound))

    velo = np.delete(velo, idxs)
    flux = np.delete(flux, idxs)
    interpFlux = np.delete(interpFlux, idxs)

    return velo, flux, interpFlux

#Cuts the spectrum to within a range of -1500 to +1500 km/s
def trimSpectrum(wavelength, flux):
    minWave=-1500/299792*6564.61+6564.61
    maxWave=1500/299792*6564.61+6564.61

    newWave=[]
    newFlux=[]

    for i in range(len(wavelength)):
        if minWave<=wavelength[i]<=maxWave:
            newWave.append(wavelength[i])
            newFlux.append(flux[i])

    return (newWave, newFlux)

#Computes the equivalent width of a spectral line
def ewcalc(velocity,flux,plotWave,normFlux,subtractCont):
    #Convert from velocity [km/s] to wavelength [angstroms]
    speedOfLight=299792#speed of light in km/s
    h_alpha_wavelength=6564.61#wavelength of H_alpha emission line (angstroms)
    wavelength=list(map(lambda obs: (obs/speedOfLight)*h_alpha_wavelength+\
        h_alpha_wavelength, velocity))

    #define list names
    lambdatable = []
    Hfluxtable = []
    polylambda = []
    polyflux = []

    #manipulates observ to find h-alpha line
    normFlux=False
    plotWave=True
    (wavelength, flux)=trimSpectrum(wavelength, flux)
    xyParams=findHAlphaLine(wavelength, flux, plotWave, normFlux, subtractCont)
    Hfluxtable=xyParams[1]
    lambdatable=xyParams[0]

    for idx, num in enumerate(lambdatable):
        if num<6555 or num>6575:
            polylambda.append(num)
            polyflux.append(Hfluxtable[idx])

    #fit linear polynomial to 'poly' lists
    a = np.polyfit(polylambda, polyflux, 1)
    #normalize every flux value to polynomial 'a'
    normfluxtable = []
    for i in range(0,len(Hfluxtable)):
        normfluxtable.append(Hfluxtable[i]/(a[1]+a[0]*lambdatable[i]))
    #uses eqn(9.59) from carroll and ostille for EW calc (F_c = 1 here)
    ewfluxtable = []
    for i in range(len(normfluxtable)):
        ewfluxtable.append(1-normfluxtable[i])
    #interpolates function of 1-flux over all lambdas
    f = interp.interp1d(lambdatable,ewfluxtable,axis=0,fill_value='extrapolate')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        #integrates f*dlambda with lambda in chris's fixed integration limits
        #(i.e integral((1-flux)dlambda), lambda = 0.653997 to 0. 658586)
        ew = integrate.quad(f,6539.77,6585.86)
    ew = ew[0]*0.1 #convert ew from angstroms to nanometers
    return ew

#Writes the flux, velocity, and other parameters to a text file for use elsewhere
def writeBestFit(filename, bestIncl, fac_e, v_h, v_e, interpFlux, obsValX,
                ew_fit, bestChi2, bestRedChi2):
    #Make the out folder
    dir0 = os.getcwd()
    savepath = dir0+'/out/'
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    filename = savepath+filename.split("/")[-1]+".txt"

    #Write to the file
    with open(filename, 'w') as file:
        file.write('%'+filename+'\n')
        file.write('%Date Created:'+datetime.datetime.today()\
            .strftime('%y-%m-%d %H:%M:%S')+'\n')
        file.write('{0:20}{1:20}{2:20}{3:20}{4:20}{5:20}{6:20} \n'\
            .format('Incl [deg]', 'f_e', 'v_h [km/s]', 'v_e [km/s]', 'EW [nm]',\
            'chi2', 'redChi2'))
        file.write('{0:20}{1:20}{2:20}{3:20}{4:20}{5:20}{6:20} \n'.format(\
            str('{0:.5f}'.format(bestIncl)),str('{0:.5f}'.format(fac_e)),\
            str('{0:.5f}'.format(v_h)),str('{0:.5f}'.format(v_e)),\
            str('{0:.5f}'.format(ew_fit)),str('{0:.5f}'.format(bestChi2)),\
            str('{0:.5f}'.format(bestRedChi2))))
        file.write('\n {0:30}{1} \n'\
            .format('%Flux [Normalized]', 'Velocity [km/s]'))
        for idx, val in enumerate(interpFlux):
            towrite='{0:30}{1}'.format(str(interpFlux[idx]), str(obsValX[idx]))
            file.write('%s \n' % towrite)


#Plots the best fitting observ set simulation observ
def pltBestFit(params, inclAngles, ranges, obsValX, obsValY, inclFlag,
                reddeningFiles, luminosity, distance, plotWave, normFlux,
                subtractCont, filename, veloUBound, veloLBound, bestChi2,
                plotOrigLine, replotting, *inclSpacing):
    if inclFlag:
        bestIncl, fac_e, v_h, v_e = params
        simX, simY, origData = getSim(reddeningFiles[0], bestIncl, ranges,
                        luminosity[0], distance, fac_e, v_h, v_e, inclFlag,
                        plotWave, normFlux, subtractCont, inclSpacing)
    else:
        fac_e, v_h, v_e = params
        bestIncl = inclAngles
        simX, simY, origData = getSim(reddeningFiles[0], bestIncl, ranges,
                        luminosity[0], distance, fac_e, v_h, v_e, inclFlag,
                        plotWave, normFlux, subtractCont)

    #Interpolate the flux to match x array of sim to obs
    interpFlux = griddata(np.array(simX),simY,np.array(obsValX),method='linear',
        fill_value=1)

    #Compute the equivalent widths of the observed and simulated lines
    ew_ori=ewcalc(origData[0],origData[1],plotWave,normFlux,subtractCont)
    ew_fit=ewcalc(obsValX,interpFlux,plotWave,normFlux,subtractCont)
    ew_obs=ewcalc(obsValX,obsValY,plotWave,normFlux,subtractCont)
    print('EW Simulation Original:', str('{0:.2f}'.format(ew_ori)))
    print('EW Simulation Fitted:', str('{0:.2f}'.format(ew_fit)))
    print('EW Observed:', str('{0:.2f}'.format(ew_obs)))
    print('Best Fitting Model - Chi-squared:', str('{0:.2f}'\
        .format(bestChi2 * -1)))

    #Getting the reduced chi-squared
    #Interpolate the flux to match x array of sim to obs
    interpFlux = griddata(np.array(simX),simY,np.array(obsValX),
        method='linear',fill_value=1)
    #All chi2 are divided equally by this value
    st_dev = 1
    #Cut the velocity range over which the chi-squared is being computed
    obsValXCut, obsValYCut, interpFluxCut = veloWindow(obsValX,obsValY,\
        interpFlux,veloUBound,veloLBound)
    #Compute the reduced chi-squared (degrees of freedom set to length of params)
    bestRedChi2 = (np.sum(((obsValYCut - interpFluxCut)**2 / (st_dev)**2.)))\
        /(len(interpFluxCut)-len(params))
    print('Best Fitting Model - Reduced Chi-squared:', str('{0:.2f}'\
        .format(bestRedChi2)))

    #Write the best fitting line's flux and wavelength to a file if not replotting
    if not replotting:
        writeBestFit(filename, bestIncl, fac_e, v_h, v_e, interpFluxCut,\
         obsValXCut, ew_fit, bestChi2, bestRedChi2)

    plt.figure(num=None, figsize=(11, 8), dpi=100, facecolor='w', edgecolor='k')
    #Plot the original halpha line before convolution if desired
    if plotOrigLine:
        L = sorted(zip(origData[0],origData[1]), key=operator.itemgetter(0))
        origData[0], origData[1] = zip(*L)
        plt.plot(origData[0],origData[1], color='lightsteelblue', linestyle='--',\
        label='Simulated\nOriginal', linewidth=1.0)
    #Plot the observed observ set to compare
    if plotOrigLine:
        plt.plot(obsValX, obsValY, color='black', label='\nObserved\n1998 (HEROS)'\
            +'\nEW='+str('{0:.2f}'.format(ew_obs))+' nm', linewidth=2.0)
    else:
        plt.plot(obsValX, obsValY, color='black', label='Observed\n1998 (HEROS)'\
            +'\nEW='+str('{0:.2f}'.format(ew_obs))+' nm', linewidth=2.0)
    #Plot the best fitting observ
    plt.plot(obsValX, interpFlux, color='royalblue', label='\nSimulated'\
        +'\n'+r'$i=$ '+str('{0:.2f}'.format(bestIncl))+r'$^{\circ}$'\
        +'\n'+r'$f_e=$'+str('{0:.2f}'.format(fac_e))\
        +'\n'+r'$v_h=$'+str('{0:.2f}'.format(v_h))+' km/s'\
        +'\n'+r'$v_e=$'+str('{0:.2f}'.format(v_e))+' km/s'\
        +'\nEW='+str('{0:.2f}'.format(ew_fit))+' nm'\
        +'\n'+r'$\chi^{2}=$'+str('{0:.2f}'.format(bestChi2*-1)))
    #Plot the boundaries of the line fitting
    plt.plot([veloUBound, veloUBound], [0.5, 1.5], color='black', linestyle='--')
    plt.plot([veloLBound, veloLBound], [0.5, 1.5], color='black', linestyle='--')
    plt.xlim(-1250, 1250)
    #plt.ylim(0.75, 7.5)
    plt.ylim(ymin=0.75)
    plt.xlabel("Velocity [km/s]",fontsize='x-large')
    plt.ylabel("Relative Flux",fontsize='x-large')
    plt.xticks(size='x-large')
    plt.yticks(size='x-large')
    plt.legend(fontsize='x-large')

    plt.savefig(filename + '.png')

#Plots the convergence of the walkers as they step for each parameter
def plotConvergence(npy, filename, ranges, inclFlag):
    converged_idx = 0
    
    if inclFlag: #If multiple inclinations are given
        linspace = [np.linspace(ranges[0][0], ranges[0][1], 5),
                    np.linspace(ranges[1][0], ranges[1][1], 5),
                    np.linspace(ranges[2][0], ranges[2][1], 5),
                    np.linspace(ranges[3][0], ranges[3][1], 5)]
        
        param_to_latex = dict(inc=r'$i [^{\circ}$]', f_e=r'$f_{e}$',
                              v_h=r'$v_{h} [km/s]$', v_e=r'$v_{e}$ [km/s]')
        params = ["inc", "f_e", "v_h", "v_e"]
        fig = plt.figure(num=None, figsize=(6, 8), dpi=100, facecolor='w',\
         edgecolor='k')
        
    else: #If only one inclination is given
        linspace = [np.linspace(ranges[0][0], ranges[0][1], 5),
                    np.linspace(ranges[1][0], ranges[1][1], 5),
                    np.linspace(ranges[2][0], ranges[2][1], 5)]
        
        param_to_latex = dict(f_e=r'$f_{e}$', v_h=r'$v_{h} [km/s]$',\
         v_e=r'$v_{e} [km/s]$')
        params = ["f_e", "v_h", "v_e"]
        fig = plt.figure(num=None, figsize=(6, 8), dpi=100, facecolor='w',\
         edgecolor='k')
    
    chain = np.load(npy)
    
    gs = gridspec.GridSpec(len(params), 3)
    gs.update(hspace=0.25)
    
    for ii, param in enumerate(params):
        these_chains = chain[:, :, ii]

        max_var = max(np.var(these_chains[:, converged_idx:], axis=1))

        ax1 = plt.subplot(gs[ii, :2])

        ax1.axvline(0, color="#67A9CF", alpha=0.7, linewidth=2)

        for walker in these_chains:
            ax1.plot(np.arange(len(walker)) - converged_idx, walker,
                     drawstyle="steps",
                     color=cm.Blues_r(np.var(walker[converged_idx:]) /max_var),
                     alpha=0.5)

        ax1.set_ylabel(param_to_latex[param], fontsize='xx-large', labelpad=18,
                       rotation="vertical", color=font_color)

        # Don't show ticks on the y-axis
        ax1.yaxis.set_ticks([])

        # For the plot on the bottom, add an x-axis label. Hide all others
        if ii == len(params) - 1:
            ax1.set_xlabel("step number", fontsize='x-large', labelpad=18,
                           color=font_color)
        else:
            ax1.xaxis.set_visible(False)

        ax2 = plt.subplot(gs[ii, 2])

        ax2.hist(np.ravel(these_chains[:, converged_idx:]),
                 bins=np.linspace(ax1.get_ylim()[0], ax1.get_ylim()[1], 20),
                 orientation='horizontal',
                 facecolor="#67A9CF",
                 edgecolor="none",
                 normed=True, histtype='barstacked')

        ax2.xaxis.set_visible(False)
        ax2.yaxis.tick_right()

        ax2.set_yticks(linspace[ii])
        ax1.set_ylim(ax2.get_ylim())

        if ii == 0:
            t = ax1.set_title("Walkers", fontsize='x-large', color=font_color)
            t.set_y(1.01)
            t = ax2.set_title("Posterior", fontsize='x-large', color=font_color)
            t.set_y(1.01)

        ax1.tick_params(axis='x', pad=2, direction='out',
                        colors=tick_color, labelsize='large')
        ax2.tick_params(axis='y', pad=2, direction='out',
                        colors=tick_color, labelsize='large')

        ax1.get_xaxis().tick_bottom()

    fig.subplots_adjust(hspace=0.0, wspace=0.0, bottom=0.10, top=0.93,
                        left=0.12, right=0.87)
    #plt.show()
    plt.savefig(filename + '.png')

#This function is what emcee runs
def lnprob(params, obsValX, obsValY, inclAngles, ranges, inclFlag, reddeningFiles,
            luminosity, distance, plotWave, normFlux, subtractCont, veloUBound,
            veloLBound, *inclSpacing):
    count = 0
    inside_ranges = True
    #This ensures that only values within the given parameter ranges are used
    while inside_ranges * (count < len(params)):
        inside_ranges = (params[count] >= ranges[count, 0]) *\
            (params[count] <= ranges[count, 1])
        count += 1
    if inside_ranges:
        #Simulation observ
        if inclFlag:
            inclAngles, fac_e, v_h, v_e = params
        else:
            fac_e, v_h, v_e = params

        simX, simY, origData = getSim(reddeningFiles[0], inclAngles, ranges,
                            luminosity[0], distance, fac_e, v_h, v_e, inclFlag,
                            plotWave, normFlux, subtractCont, inclSpacing)

        #Interpolate the flux to match x array of sim to obs
        interpFlux = griddata(np.array(simX),simY,np.array(obsValX),\
            method='linear',fill_value=1)
        
        #All chi2 are divided equally by this value
        st_dev = 1

        #Cut the velocity range over which the chi-squared is being computed
        obsValX, obsValY, interpFlux = veloWindow(obsValX,obsValY,interpFlux,\
            veloUBound,veloLBound)
        chi2 = np.sum(((obsValY - interpFlux)**2 / (st_dev)**2.))

        return -1 * chi2

    else:
        return -np.inf

    


# MAIN PART OF THE SCRIPT-------------------------------------------------------
def run(filename, a_parameter, nwalk, nstep, nburn, fac_e, v_e, v_h, inclFlag,
    inclAngles, reddeningFiles, plot_Rivi, luminosity, distance,
    plotWave, normFlux, subtractCont, veloUBound, veloLBound, plotOrigLine,
    *inclSpacing):
        
    # Getting the observed observ-----------------------------------------------
    obsValX, obsValY = plotRiviData() #x = velo y = flux


    # emcee stuff beyond here---------------------------------------------------
    if inclFlag:
        ranges = np.array([inclAngles, fac_e, v_h, v_e])
        ndim = 4
    else:
        ranges = np.array([fac_e, v_h, v_e])
        ndim=3

    #First guess of parameters
    p0 = [np.random.rand(ndim) * (ranges[:, 1] - ranges[:, 0]) +
            ranges[:, 0] for i in range(nwalk)]
    start_time = time.time()
    
    with Pool() as pool:
        if inclFlag:
            sampler = emcee.EnsembleSampler(nwalk, ndim, lnprob, args=[obsValX,
            obsValY, inclAngles, ranges, inclFlag, reddeningFiles, luminosity,
            distance, plotWave, normFlux, subtractCont, veloUBound, veloLBound,
            inclSpacing], a=a_parameter, pool=pool)
        else:
            sampler = emcee.EnsembleSampler(nwalk, ndim, lnprob, args=[obsValX,
            obsValY, inclAngles, ranges, inclFlag, reddeningFiles, luminosity,
            distance, plotWave, normFlux, subtractCont, veloUBound, veloLBound],
            a=a_parameter, pool=pool)

        # Running the Burn-In---------------------------------------------------
        print("Running burn-in with " + str(nburn) + " walkers ################")

        pos, prob, state = sampler.run_mcmc(p0, nburn, progress=True)

        print("--- Completed in %s minutes ---" % (str('{0:.2f}'.format((time\
            .time() - start_time) / 60))))
        af_bi = np.mean(sampler.acceptance_fraction)
        print("Mean acceptance fraction (Burn-In): ", str('{0:.2f}'\
            .format(af_bi)))
        sampler.reset()


        # Running MCMC----------------------------------------------------------
        print("Running MCMC ############################################")
        sampler.run_mcmc(pos, nstep, rstate0=state, progress=True)

        print("--- Completed in %s minutes ---" % (str('{0:.2f}'.format((time\
            .time() - start_time) / 60))))
        af = np.mean(sampler.acceptance_fraction)
        print("Mean acceptance fraction (MCMC): ", str('{0:.2f}'.format(af)))
        #print("Mean autocorrelation time: {0:.3f} steps"\
        #.format(np.mean(sampler.get_autocorr_time())))

    #Getting the chain and flattening
    chain = sampler.chain
    flatchain = chain.reshape((-1, ndim))

    #best fit parameters
    maxprob_index = np.argmax(prob)

    #Get the best parameters and their respective errors
    params_fit = pos[maxprob_index]
    errors_fit = [sampler.flatchain[:, i].std() for i in range(ndim)]

    # Saving the chain----------------------------------------------------------
    dir0 = os.getcwd()
    savepath = dir0+'/vars/chains/'
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    chain = sampler.chain
    chain_npy = savepath+filename+'_af'+str('{0:.2f}'.format(af))+"_chain.npy"
    np.save(chain_npy, chain)

    # Saving the prob and pos af------------------------------------------------
    savepath = dir0+'/vars/probs/'
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    probs_npy = savepath+filename+'_af'+str('{0:.2f}'.format(af))+"_prob.npy"
    np.save(probs_npy, prob)

    savepath = dir0+'/vars/pos/'
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    pos_npy = savepath+filename+'_af'+str('{0:.2f}'.format(af))+"_pos.npy"
    np.save(pos_npy, pos)

    savepath = dir0+'/vars/af/'
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    af_npy = savepath+filename+'_af'+str('{0:.2f}'.format(af))+"_af.npy"
    np.save(af_npy, af)

    # Corner Plotting-----------------------------------------------------------
    savepath = dir0+'/figures/'
    if not os.path.exists(savepath):
        os.makedirs(savepath)
        
    if inclFlag:
        labels = [r'$i$', r'$f_e$', r'$v_h$', r'$v_e$']
    else:
        labels = [r'$f_e$', r'$v_h$', r'$v_e$']
    quantiles = [0.16, 0.5, 0.84]
    figure = corner(flatchain, labels=labels, range=ranges, quantiles=quantiles,
            plot_contours=True, smooth=2., smooth1d=False, plot_datapoints=True,
            label_kwargs={'fontsize': 17},title_kwargs={'fontsize': 17},
            truths=None, show_titles=True, color_hist='black', plot_density=True,
            fill_contours=False, no_fill_contours=False, normed=True)
    figure.savefig(savepath+filename+'_af'+str('{0:.2f}'\
        .format(af))+"_corner.png")

    # Convergence Plotting-----------------------------------------------------------
    plotConvergence(chain_npy, savepath+filename+'_af'+str('{0:.2f}'.format(af))\
        +'_convergence', ranges, inclFlag)

    # Halpha Plotting----------------------------------------------------------------
    replotting=False
    if inclFlag:
        pltBestFit(params_fit, inclAngles, ranges, obsValX, obsValY, inclFlag,
                reddeningFiles, luminosity, distance, plotWave, normFlux,
                subtractCont, savepath+filename+'_af'+str('{0:.2f}'.format(af))\
                +'_halpha', veloUBound, veloLBound, np.max(prob), plotOrigLine,
                replotting, inclSpacing)
    else:
        pltBestFit(params_fit, inclAngles, ranges, obsValX, obsValY, inclFlag,
                reddeningFiles, luminosity, distance, plotWave, normFlux,
                subtractCont, savepath+filename+'_af'+str('{0:.2f}'.format(af))\
                +'_halpha', veloUBound, veloLBound, np.max(prob), plotOrigLine,
                replotting)

    print("Fitting complete ########################################")
    print("--- Completed in %s minutes ---" % (str('{0:.2f}'.format((time.time()\
     - start_time) / 60))))

