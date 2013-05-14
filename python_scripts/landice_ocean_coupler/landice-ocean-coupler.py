#!/usr/bin/python
# Coupler for the MPAS Ocean and Land Ice models
# Matt Hoffman, March 8, 2013

import numpy as np
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
runOCNcmd = 'mpirun -np 4 ocean_model'

# Setup needed files
liceInput = 'land_ice_grid.nc'
liceOutput = 'output.nc'
ocnInput = 'input.nc'
ocnOutput = 'output.nc'

# Define number of coupling intervals (assumes you have set dt for each model accordingly, and ocean dt is evenly divisble in land ice dt)
# (We may want to do more sophisticated time-keeping of the coupling if there is a need.)
numCoupleIntervals = 8

# Some parameters
rhoi = 900.0  # land ice model's ice density (kg m^-3)
rhow = 1000.0  # density of pure water (kg m^-3)
grav = 9.8101 # gravitational acceleration
secondsInYear = 3600.0 * 24.0 * 365.0  # number of seconds in a year

ocnThicknessMin = 50.0  # limit to how thin the ocean can be (m)

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

os.system('mkdir -p land_ice/output-files')
os.system('mkdir -p ocean/output-files')

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
        liceLowerSurface = liceSourceFile.variables['lowerSurface'][-1,:]
        liceBedTopography = liceSourceFile.variables['bedTopography'][-1,:]
        liceBasalTemp = liceSourceFile.variables['temperature'][-1,:,-1]
    except:
        print "Problem getting needed Source fields from the land ice files"
    try:
        ocnSurfTemp = ocnSourceFile.variables['temperature'][-1,:,0]
        ocnSurfDensity = ocnSourceFile.variables['density'][-1,:,0]
        ocnSurfSalinity = ocnSourceFile.variables['salinity'][-1,:,0]
        ocnSurfPressure = ocnSourceFile.variables['pressure'][-1,:,0]
        ocnBottomLayerDensity = ocnSourceFile.variables['density'][-1,:,-1]  # density of bottom layer
        ocnBottomLayerPressure = ocnSourceFile.variables['pressure'][-1,:,-1] # pressure at center of bottom layer
        ocnBottomLayerThickness = ocnSourceFile.variables['layerThickness'][-1,:,-1] # thickness of bottom layer
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

    # (these could be function calls)

    # --- Calculate the ocean surface pressure 
    # Where there is no ice, this will take a value of 0.0
    # Some areas will have ice that is thick enough to make a negative ocean thickness (grounded ice).  
    # In these locations instead assign the ocean surface pressure to leave some small amount of water in the ocean,
    # based on the parameter defined above called ocnThicknessMin.

    # compute pressure at bottom of ocean: take the pressure at the center of the bottom cell and add on pressure due to the 
    # lower half of the bottom cell.
    ocnBottomPressure = ocnBottomLayerPressure + ocnBottomLayerDensity * grav * (0.5 * ocnBottomLayerThickness)
    ocnSurfPressureLimit = (ocnBottomPressure - ocnBottomLayerDensity * grav * ocnThicknessMin)
    ocnSurfPressure = np.minimum(rhoi * grav * liceThickness, ocnSurfPressureLimit)  # Units: Pa


    # --- Calculate boundary layer fluxes
    ###calculateBdyLyrFluxes(liceBasalTemp, ocnSurfTemp, ocnSurfSalin,   liceBHF, ocnMassFlux, ocnHeatFlux, ocnSalinFlux) 

    # Create a mask of where the ice shelf is located
    isShelf = (liceLowerSurface > liceBedTopography)

    # Heat Flux: Eq. 3 from ISOMIP.  Tf should be in-situ freezing point of seawater at the local pressure.
    a = -.0573     # C/psu Salinity coefficient of freezing equation
    b = .0939      # C Salinity coefficient of freezing equation
    c = -7.53e-8   # C/Pa Salinity coefficient of freezing equation
    Tf = a * ocnSurfaceSalinity + b + c * ocnSurfacePressure     # Eq. 1 from Holland and Jenkins, Modeling Thermodynamic IceOcean Interactions at the Base of an Ice Shelf. J. Physical Oceanography, Aug 1999.
    heatFlux = ocnSurfDensity * waterSpecificHeat * exchangeVel * (ocnSurfTemp - Tf)  # upward positive, units: W/m2
    heatFlux[not(isShelf)] = 0.0  # Only allow a heat flux where you have an ice shelf
    ocnSurfTempFlux = heatFlux # TODO: What should this be exactly???
    
    # Mass Flux: Eq. 4, modified
    massFlux = heatFlux / latentHeatFusion  # positive=melting, units: kg/(m2 s); divide by appropriate density to get thickness rate of change for each model component.
    ocnSurfMassFlux = massFlux / rhow  # units: m of water per second.  Assuming pure freshwater density here - TODO is that right?
    liceBMB = -massFlux / rhoi * secondsInYear  # units: m of ice per yr


    # =======================================
    # 3. Write the new forcing data to the files
    # =======================================

    # Write new quantities to the destination files
    # TODO Does this need to be done explicitly?  I always forget how netCDF i/o works in python precisely.
    liceDestFile.close()
    ocnDestFile.close()

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
                subprocess.check_call('echo "0000-01-01_00:00:00" > restart_timestamp', shell=True, executable='/bin/bash') # Setup initial time stamp.
                subprocess.check_call('cp grid.nc input.nc', shell=True, executable='/bin/bash') # Copy the initial input file to input.nc
                subprocess.check_call('sed -i.SEDBACKUP "s/config_write_output_on_startup.*/config_write_output_on_startup = .true./" namelist.input' , shell=True, executable='/bin/bash')  # write the initial time
            else:
                # Get the old stop time so we can set it to be the new start time.
                p = os.popen('cat restart_timestamp')
                endtime = p.read().rstrip().replace(" ", "").replace(":",".")
                restart_filename = "restart.%s.nc"%endtime
                subprocess.check_call('cp %s input.nc'%restart_filename, shell=True, executable='/bin/bash')  # Copy the restart file to input.nc
                subprocess.check_call('sed -i.SEDBACKUP "s/config_write_output_on_startup.*/config_write_output_on_startup = .false./" namelist.input' , shell=True, executable='/bin/bash')  # write the initial time


            print 'Starting ocean model.'
            subprocess.check_call(runOCNcmd, shell=True, executable='/bin/bash')
            print 'Ocean model completed.'
            if t == 0:
                output_filename = "output.0000-01-01_00.00.00.nc"
            else:
                p = os.popen('cat restart_timestamp')
                endtime = p.read().rstrip().replace(" ", "").replace(":",".")
                output_filename = "output.%s.nc"%endtime

            print "File: %s"%output_filename
            subprocess.check_call('mv %s output-files/.'%output_filename, shell=True, executable='/bin/bash')
            print "Moved"

            subprocess.check_call('rm -f output.*.nc', shell=True, executable='/bin/bash')
            print "Deleted"
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
print 'Successful completion.'



