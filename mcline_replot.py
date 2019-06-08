#MCLINE_REPLOT.py
#Keegan Marr
#Last Modified: May 20, 2019

import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import mcline_functions as mcf
import numpy as np
import datetime
import time
import os
from corner import corner


# MAIN PART OF THE SCRIPT--------------------------------------------------------
def main():
    #The file name of the chain to be replotted
    filename = "19-05-20_test_af0.10"
    
    chain_npy = "vars/chains/"+filename+"_chain.npy"
    prob_npy = "vars/probs/"+filename+"_prob.npy"
    pos_npy = "vars/pos/"+filename+"_pos.npy"
    af_npy = "vars/af/"+filename+"_af.npy"
        
    
    
    #COPY INPUTS FROM mcline.py
    # INPUTS========================================================================
    date = datetime.datetime.today().strftime('%y-%m-%d')
    filename = date+'_19-05-14_mod03_a=4.0_i0-90_nwalk100_nstep500'
    filename = date+'_test'

    a_parameter = 4.0 #Walker step size (internal to emcee)
    nwalk = 100 #number of walkers
    nstep = 500 #number of steps for each walker
    nwalk = 8 #number of walkers
    nstep = 10 #number of steps for each walker
    nburn = 0.1*nstep #number of steps to burn (typically 10% of nstep)

    #PRIORS ########################################################
    #Ranges of values to run MCMC over
    fac_e = [0.2, 0.4] #Fraction of light scattered by electrons
    v_e = [500.0, 650.0] #Speed of electron motion #Expect 582 km/s
    v_h = [10.0, 20.0] #Sound speed of the disk #Expect 13.58 km/s
    ################################################################


    #Upper and lower bounding velocity at which to fit
    veloUBound = 650
    veloLBound = -650

    inclFlag = True
    if inclFlag: #Range of inclinations, linewidth=2.0
        #The lower bound of this range must correspond to some observer angle
        #in your HDUST observers file
        inclAngles = [0, 90] #Inclination angles
        
        #Spacing between grid of inclination in HDUST observers file.
        #For example I have 0-90 degrees, every 2.5 degrees for a total of 37 observers
        inclSpacing = 2.5
        
    else: #One inclination (Must be precisely what is in observers file)
        inclAngles=45


    # Sim and Obs Data==============================================================
    #Reddened fullsed file
    simuDir="../hdust/66oph/19-02-15_compare_19-02-12_best_with_19-01-28_best/"
    #simuFile="66oph/fullsed/fullsed_mod09a_red.txt" #rho=1e-11 n=2.5
    simuFile="66oph/fullsed/fullsed_mod21a_red.txt" #rho=1e-11 n=2.4
    #simuFile="66oph/fullsed/fullsed_mod16a_red.txt" #rho=1e-11 n=2.3

    simuDir="../hdust/66oph/19-05-14_increasing_Rout/"
    #simuFile="66oph/fullsed/fullsed_mod01a_red.txt"
    simuFile="66oph/fullsed/fullsed_mod03a_red.txt"
    simuFile=[simuDir+simuFile]
    reddeningFiles=simuFile

    #observedDataFiles=["66oph_20180929_withwave.txt"] #66oph
    #oph66_obs = False #FOR BESS DATA. False for 66 Oph data. Otherwise it will be off center
    plot_Rivi_66Oph = True

    #Luminosity in L_sol
    luminosity=[8042.5] #66oph

    #Distance in pc
    distance=203.24 #66oph

    #Plot the original hdust simulation line before convolution
    plotOrigLine=True

    #Don't change the options below, they are vestigial options that I need to
    #remove and keep fixed.
    #Halpha plot wavelength (true) or velocity (false)
    plotWave=False

    #Choose whether the flux should be normalized (true) or not (false) for when plotting Halpha
    normFlux=True

    #Subtract the continuum to compare lines that are not produced on the same continua. Subtract (true) don't (false)
    subtractCont=False
    # ==============================================================================





    # Getting the observed data-----------------------------------------------------
    obsValX, obsValY = mcf.plotRiviData() #x = velo y = flux
    
    # emcee stuff beyond here-------------------------------------------------------
    if inclFlag:
        ranges = np.array([inclAngles, fac_e, v_h, v_e])
        ndim = 4
    else:
        ranges = np.array([fac_e, v_h, v_e])
        ndim=3
        
    #Getting the chain and flattening and prob array for halpha plot
    chain = np.load(chain_npy)
    flatchain = chain.reshape((-1, ndim))
    prob = np.load(prob_npy)
    pos = np.load(pos_npy)
    af = np.load(af_npy)


    #best fit parameters
    maxprob_index = np.argmax(prob)

    #Get the best parameters and their respective errors
    params_fit = pos[maxprob_index]
    

    dir0 = os.getcwd()
    savepath = dir0+'/figures/'
    
    # Corner Plotting----------------------------------------------------------------
    savepath = dir0+'/figures/'
    if not os.path.exists(savepath):
        os.makedirs(savepath)
        
    if inclFlag:
        labels = [r'$i$', r'$f_e$', r'$v_h$', r'$v_e$']
    else:
        labels = [r'$f_e$', r'$v_h$', r'$v_e$']
    quantiles = [0.16, 0.5, 0.84]
    figure=corner(flatchain, labels=labels, range=ranges, quantiles=quantiles,
                plot_contours=True, smooth=2., smooth1d=False, plot_datapoints=True,
                label_kwargs={'fontsize': 17},title_kwargs={'fontsize': 17},
                truths=None, show_titles=True, color_hist='black', plot_density=True,
                fill_contours=False, no_fill_contours=False, normed=True)

    # Convergence Plotting-----------------------------------------------------------
    mcf.plotConvergence(chain_npy, savepath+filename+'_af'+str('{0:.2f}'.format(af))\
        +'_convergence', ranges, inclFlag)

    # Halpha plotting----------------------------------------------------------------
    replotting=True
    if inclFlag:
        inclSpacing = [inclSpacing]
        mcf.pltBestFit(params_fit, inclAngles, ranges, obsValX, obsValY, inclFlag, reddeningFiles,
                luminosity, distance, plotWave, normFlux, subtractCont,
                savepath+filename+'_af'+str('{0:.2f}'.format(af))+'_halpha',
                veloUBound, veloLBound, np.max(prob), plotOrigLine, replotting, inclSpacing)
    else:
        mcf.pltBestFit(params_fit, inclAngles, ranges, obsValX, obsValY, inclFlag, reddeningFiles,
                luminosity, distance, plotWave, normFlux, subtractCont,
                savepath+filename+'_af'+str('{0:.2f}'.format(af))+'_halpha',
                veloUBound, veloLBound, np.max(prob), plotOrigLine, replotting)
                
    plt.show()
main()
