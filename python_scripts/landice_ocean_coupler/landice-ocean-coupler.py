#!/usr/bin/python
# Coupler for the MPAS Ocean and Land Ice models
# Matt Hoffman, March 8, 2013

import numpy
import os, sys, subprocess
import time
#import matplotlib.pyplot as plt
try:
  from Scientific.IO.NetCDF import NetCDFFile
  netCDF_module = 'Scientific.IO.NetCDF'
except ImportError:
  try:
    from netCDF4 import Dataset as NetCDFFile
    netCDF_module = 'netCDF4'
  except ImportError:
      print 'Unable to import any of the following python modules:'
      print '  Scientific.IO.NetCDF \n  netcdf4 '
      print 'One of them must be installed.'
      raise ImportError('No netCDF module found')


# In the future we may want to add command line arguments using the optParse module.

# The script assumes you have setup a directory structure:
#  ocean model is setup in ./ocean
#  land ice model is setup in  ./land_ice
#  and each has a subdirectory called 'output-files' where all the output is accumulated.


# =======================================

# Set if you are using data or running a model
runLICE = False
runOCN = True

# Setup needed executables
runLICEcmd = 'mpirun -np 4 land_ice_model.exe'
runOCNcmd = 'mpirun -np 4 ocean_model.exe'

# Setup needed files
liceInput = 'land_ice_grid.nc'
liceOutput = 'output.nc'
ocnInput = 'grid.nc'
ocnOutput = 'output.nc'

# Define number of coupling intervals (assumes you have set dt for each model accordingly, and ocean dt is evenly divisble in land ice dt)
# (We may want to do more sophisticated time-keeping of the coupling if there is a need.)
numCoupleIntervals = 8

# Some parameters
rhoi = 900.0  # land ice model's ice density (kg m^-3)
grav = 9.8101 # gravitational acceleration

exchangeVel = 0.0001 # m/s
waterSpecificHeat = 3974 # J / (kg * K)
latentHeatFusion = 334000 # J / kg


# =======================================

# clean up previous runs
os.system('rm -rf land_ice/output.*')
os.system('rm -rf land_ice/restart.*')
os.system('rm -rf land_ice/log.*')
os.system('rm -rf land_ice/output-files/*')
os.system('rm -rf ocean/output.*')
os.system('rm -rf ocean/restart.*')
os.system('rm -rf ocean/output-files/*')

for t in range(numCoupleIntervals):
    print '============================================================'
    print 'Starting coupling number ', t

    # =======================================
    # 1. Get the appropriate files and variables for the coupler
    # =======================================

    # Get the appropriate source and destination files for the input/output of the coupler
    if t==0:
        if runLICE:
            # Use I.C. from input files as Source at initial time
            liceSourceFile = NetCDFFile('./land_ice/'+liceInput,'r')
            # At initial time, we want to write to the input files
            liceDestFile = NetCDFFile('./land_ice/'+liceInput,'r+')
        else:
            # TODO Read from a file
            pass

        if runOCN:
            # TODO - set this up properly for the ocean model
            # Use I.C. from input files as Source at initial time
            ocnSourceFile =  NetCDFFile('./ocean/'+ocnInput,'r')
            # At initial time, we want to write to the input files
            ocnDestFile =  NetCDFFile('./ocean/'+ocnInput,'r+')
        else:
            # TODO Read from a file
            pass

    else:
        if runLICE:
            # Use values from output files as Source at other times.
            liceSourceFile = NetCDFFile('./land_ice/' + liceOutput, 'r')

            # At other times we want to write to the restart files as the Destination
            # Need to choose the appropriate restart file - there may be smarter ways to do this, but an easy way is to take the most recent
            # os.listdir does not have a sort ability, so use a shell command (unpythonic but works)
            p = os.popen('ls -1t ./land_ice/restart.*.nc | head -n 1')
            liceRestart = p.read().rstrip()  # get stdout
            print 'For land ice model restart data, using file: ' + liceRestart
            liceDestFile = NetCDFFile(liceRestart,'r+')
        else:
            # TODO Read from a file
            pass

        if runOCN:
            # Use values from output files as Source at other times.
            # TODO - set this up properly for the ocean model
            pass

        else:
            # TODO Read from a file
            pass

    # Get the Source fields needed by the coupler
    # (Get the -1 time level for all variables - this will be the last whether this is from input or output
    try:
        liceThickness = liceSourceFile.variables['thickness'][-1,:]
        liceBasalTemp = liceSourceFile.variables['temperature'][-1,:,-1]
    except:
        print "Problem getting needed Source fields from the land ice files"
    try:
        ocnSurfTemp = ocnSourceFile.variables['temperature'][-1,:,0]
        ocnSurfDensity = ocnSourceFile.variables['density'][-1,:,0]
    except:
        print "Problem getting needed Source fields from the ocean files"

    # Get the destination fields needed by the coupler
    try:
        liceBMB = liceDestFile.variables['marineBasalMassBalTimeSeries'][:,0]  # There should only be on time level
        #liceBHF = liceDestFile.variables['basalHeatFluxTimeSeries'][:,0]  # NEED TO RESTRICT THIS TO CELLS WITH ICE
    except:
        print "Problem getting needed Destination fields from the land ice files"
    try:
        ocnSurfPressure = ocnDestFile.variables['seaSurfacePressure'][0,:]
        ocnSurfTempFlux = ocnDestFile.variables['surfaceTemperatureFlux'][0,:]
        ocnSurfMassFlux = ocnDestFile.variables['surfaceMassFlux'][0,:]
    except:
        print "Problem getting needed Destination fields from the ocean files"


    # =======================================
    # 2. Do the coupler calculations
    # =======================================

    # Calculate the ocean surface pressure - this could be a function call
    # TODO WHAT TO DO ABOUT NON-ICE SHELF OCEAN CELLS?
#    ocnSurfPressure = rhoi * grav * liceThickness  # TODO ARE THESE THE PROPER UNITS?

    # Calculate boundary layer fluxes - this could be an internal or external function call
    #calculateBdyLyrFluxes(liceBasalTemp, ocnSurfTemp, ocnSurfSalin,   liceBHF, ocnMassFlux, ocnHeatFlux, ocnSalinFlux) #TODO

### Ice Presure Melting Point ###

### Heat Flux ####
# Q_t = rho * waterSpecificHeat * exchangeVel * ( T_water - T_ice)

### Mass Flux ####
# Q_m = Q_t / (rho * latentHeatFusion)


    # =======================================
    # 3. Write the new forcing data to the files
    # =======================================

    # Write new quantities to the destination files
    # TODO Does this need to be done explicitly?  I always forget how netCDF i/o works in python precisely.


    # =======================================
    # 4. Run the models
    # =======================================
    # TODO Do we want to do this in parallel?  Probably no reason to do so.


    # === LAND ICE MODEL ===
    try:
        if runLICE:
            print 'hello'
            os.chdir('land_ice')
            # Set the namelist file properly for a restart or no and set the start and stop times
            print 'Setting up land ice namelist file'
            if t==0:
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_do_restart.*/   config_do_restart = .false./" namelist.input' , shell=True, executable='/bin/bash')  # no restart
                print 'set no restart'
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_write_output_on_startup.*/   config_write_output_on_startup = .true./" namelist.input' , shell=True, executable='/bin/bash')  # write the initial time
                print 'set write initial time'
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_start_time.*/   config_start_time = \'0000-01-01_00:00:00\'/" namelist.input' , shell=True, executable='/bin/bash')  # set the start time
                print 'set start time'
            else:
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_do_restart.*/   config_do_restart = .true./" namelist.input' , shell=True, executable='/bin/bash') # do a restart
                print 'set restart'
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_write_output_on_startup.*/   config_write_output_on_startup = .true./" namelist.input' , shell=True, executable='/bin/bash')  # don't write the initial time
                print 'disabled initial write'
                # Get the old stop time so we can set it to be the new start time. - it seems more reliable to search for the most recent restart file rather than the most recent output file.
                p = os.popen('ls -1t restart*nc | head -n 1')
                lastrestartfilename=p.read().rstrip()
                p = os.popen('ncks -v xtime -H -s "%c" ' + lastrestartfilename + ' | tail -r -c 65')
                oldstoptime=p.read().rstrip()
                print 'For land ice, getting the restart time from file: ' + lastrestartfilename + ' which gives a time of ' + oldstoptime
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_start_time.*/   config_start_time = ' + oldstoptime + '/" namelist.input ' , shell=True, executable='/bin/bash')  # set the start time to the old stop time
                print 'set start time'

            print 'Starting land ice model.'
            subprocess.check_call(runLICEcmd, shell=True, executable='/bin/bash')
            print 'Land ice model completed.'
            # Don't need to copy output file if the namelist is setup properly
            ### Copy output file so we don't lose it when the model restarts - could use Python's shutil module.  TODO Also may want to change format of the string conversion of t
            subprocess.check_call('cp ' + liceOutput + ' ./output-files/output.' + '{0:04d}'.format(t+1) + '.nc', shell=True, executable='/bin/bash')  # the number is the ending time level (if there are more than one)
            subprocess.check_call('cp log.0000.out ./output-files/log.0000.' + '{0:04d}'.format(t+1) + '.out', shell=True, executable='/bin/bash')  # the number is the ending time level (if there are more than one)
            subprocess.check_call('cp log.0000.err ./output-files/log.0000.' + '{0:04d}'.format(t+1) + '.err', shell=True, executable='/bin/bash')  # the number is the ending time level (if there are more than one)

            os.chdir('..')
    except:
        sys.exit('Land Ice model failed!')


    # === OCEAN MODEL ===
    try:
        if runOCN:
            os.chdir('ocean')
            # Set the namelist file properly for a restart or no
            print 'Setting up ocean namelist file'
            if t==0:
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_do_restart.*/   config_do_restart = .false./" namelist.input' , shell=True, executable='/bin/bash')  # no restart
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_write_output_on_startup.*/   config_write_output_on_startup = .true./" namelist.input' , shell=True, executable='/bin/bash')  # write the initial time
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_start_time.*/   config_start_time = \'0000-01-01_00:00:00\'/" namelist.input' , shell=True, executable='/bin/bash')  # set the start time
            else:
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_do_restart.*/   config_do_restart = .true./" namelist.input' , shell=True, executable='/bin/bash') # do a restart
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_write_output_on_startup.*/   config_write_output_on_startup = .false./" namelist.input' , shell=True, executable='/bin/bash')  # don't write the initial time
                # Get the old stop time so we can set it to be the new start time.
                p = os.popen('ls -1t output*nc | head -n 1')
                oldoutputfile=p.read().rstrip()
                p = os.popen('ncks -v xtime -H -s "%c" ' + oldoutputfile + ' | tail -r -c 65')
                oldstoptime=p.read().rstrip()
                subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_start_time.*/   config_start_time = ' + oldstoptime + '/" namelist.input ' , shell=True, executable='/bin/bash')  # set the start time to the old stop time


            print 'Starting ocean model.'
            subprocess.check_call(runOCNcmd, shell=True, executable='/bin/bash')
            print 'Ocean model completed.'
            # Don't need to copy output file if the namelist is setup properly
            ### Copy output file so we don't lose it when the model restarts - could use Python's shutil module.  TODO Also may want to change format of the string conversion of t
            ##subprocess.check_call('cp ' + ocnOutput + ' output.' + str(t) + '.nc', shell=True, executable='/bin/bash')
            os.chdir('..')
    except:
        sys.exit('Ocean model failed!')


# =======================================

print 'Coupled model run complete!'


# =======================================
# 5. Post-processing
# =======================================

# After the coupled run is completed, we may want to go in and use NCO to patch together all of the output files for each model.
if runLICE:
    try:
        os.chdir('land_ice/output-files')
        subprocess.check_call('ncrcat output.*.nc output-all.nc', shell=True, executable='/bin/bash')
        os.chdir('../..')
    except:
        sys.exit('Failed combining land ice output files!')

# TODO Does the ocean want to do this too?  If not, delete these lines.
if runOCN:
    try:
        os.chdir('ocean/output-files')
        subprocess.check_call('ncrcat output.*.nc output-all.nc', shell=True, executable='/bin/bash')
        os.chdir('../..')
    except:
        sys.exit('Failed combining ocean output files!')

print 'Combined output files.'
print 'Successful competion.'



