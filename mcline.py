#MCLINE.py
#Keegan Marr
#Last Modified: May 20, 2019

#Convolves an input function with a Gaussian profile
#This is the MASTER file for mcline

#REQUIRES THE DEV RELEASE OF EMCEE, AT LEAST VERSION 3.00
#SEE https://emcee.readthedocs.io/en/latest/user/install/

import datetime
from mcline_functions import run

#from multiprocessing import cpu_count
#print('This computer has {0} CPUs'.format(cpu_count()))


# ==============================================================================
# TO DO (Last updated May 20, 2019)
#DONE1)  Fix inclination mapping to mu for decimal values
#DONE2)  For multiple inclAngles, make sure I can run it for 0-90 (upper bound issue with rounding down)
#DONE3)  Allow the user to define a velocity range over which to fit
#DONE4)  Make a different python script that reads the chain and replots it
#DONE5)  Add legend to plotting
#DONE6)  Have the best fit halpha line write to a file (same format as reddening fullsed)
#DONE7)  Separate out functions into separate script from main inputs
#DONE8)  Add a percentage complete tracker (Based on the walker number)
#DONE9) Make it so that chains and figures don't overwrite eachother
#DONE10)  Have figures save instead of show
# 11) Remove unneeded variables like subtractCont, normFlux, plotWave
#DONE12) Make a function that reads the chain and replots.
#DONE13) Update mcline_replot.py to replot the best fit halpha line (will need the result of 7) )
#DONE14) Have the chi-squared be reported correctly as a reduced chi-squared
#DONE15) Compute the equiv width of each line and include them on the figure
#DONE16) Allow for plotting of the original line before convolution
#DONE17) Implement multithreading
# ==============================================================================

# INPUTS========================================================================
date = datetime.datetime.today().strftime('%y-%m-%d-%H%M%S')
filename = date+''

a_parameter = 4.0 #Walker step size (internal to emcee)
nwalk = 100 #number of walkers
nstep = 500 #number of steps for each walker
nburn = 0.1*nstep #number of steps to burn (typically 10%)

#PRIORS ########################################################
#Ranges of values to run MCMC over
fac_e = [0.2, 0.4] #Fraction of light scattered by electrons
v_e = [500.0, 650.0] #Speed of electron motion
v_h = [10.0, 20.0] #Sound speed of the disk
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
simuDir=""
simuFile=""
simuFile=[simuDir+simuFile]

#Flag for observed data
plot_Rivi_66Oph = True

#Luminosity in L_sol
luminosity=[]

#Distance in pc
distance= 

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


#Run the code
if inclFlag:
    run(filename, a_parameter, nwalk, nstep, nburn, fac_e, v_e, v_h, inclFlag,
        inclAngles, simuFile, plot_Rivi_66Oph, luminosity, distance,
        plotWave, normFlux, subtractCont, veloUBound, veloLBound, plotOrigLine,
        inclSpacing)
else:
    run(filename, a_parameter, nwalk, nstep, nburn, fac_e, v_e, v_h, inclFlag,
        inclAngles, simuFile, plot_Rivi_66Oph, luminosity, distance,
        plotWave, normFlux, subtractCont, veloUBound, veloLBound, plotOrigLine)
