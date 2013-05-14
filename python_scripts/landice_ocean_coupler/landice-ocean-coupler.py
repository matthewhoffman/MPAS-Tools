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
#  ocean model is setup in BASEDIR/ocean
#  land ice model is setup in BASEDIR/land_ice
#  Each has a subdirectory called 'output-files' where all the output is accumulated.
#  Each directory has the appropriate executable, input .nc file, and graph.info files located there.
#  namelist.input.land_ice and namelist.input.ocean should be located in BASEDIR.  THe script will copy them to the run directories.


# =======================================

BASEDIR='/Users/mhoffman/documents/mpas-git/MPAS-Tools/python_scripts/landice_ocean_coupler/dome-island-testcase'

# Set if you are using data or running a model
runLICE = True
runOCN = False

# Setup needed executables
runLICEcmd = 'mpirun -np 4 land_ice_model.exe'
runOCNcmd = 'mpirun -np 4 ocean_model'

# Setup needed file names
liceInput = 'land_ice_grid.nc'
liceOutput = 'output.nc'
ocnInput = 'input.nc'
ocnOutput = 'output.nc'

# Define number of coupling intervals (assumes you have set dt for each model accordingly, and ocean dt is evenly divisble in land ice dt)
# (We may want to do more sophisticated time-keeping of the coupling if there is a need.)
numCoupleIntervals = 3

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
os.chdir(BASEDIR)  # actually move to the appropriate directory location and do all operations there and in its subdirs.
# clean up previous runs
os.system('rm -rf land_ice/output.*')
os.system('rm -rf land_ice/restart.*')
os.system('rm -rf land_ice/log.*')
os.system('rm -rf land_ice/output-files')  # completely remove this dir before starting
os.system('rm -rf ocean/output.*')
os.system('rm -rf ocean/restart.*')
os.system('rm -rf ocean/output-files')  # completely remove this dir before starting

os.system('mkdir -p land_ice/output-files')
os.system('mkdir -p ocean/output-files')

# copy in the namelist templates to the run dirs.
os.system('cp namelist.input.land_ice land_ice/namelist.input')
os.system('cp namelist.input.ocean ocean/namelist.input')

# =======================================
# =======================================
# 0. At time 0, run both models once without timestepping to fill in the necessary diagnostic variables
#    (land ice needs the lower surface, ocean model needs the density and pressure fields.)
#    We also want a restart file for restarting the models during the coupling loop.
#    This initial run could be an entire coupling interval, but it could just be really short, like 1 day.
# =======================================
# =======================================

print 'Performing initial run for time 0 to get diagnostic output fields and a restart file.'
t=0

# === LAND ICE MODEL ===
try:
    if runLICE:
        os.chdir('land_ice')
        # Set the namelist file properly for a restart or no and set the start and stop times
        print 'Setting up land ice namelist file'
        subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_do_restart.*/   config_do_restart = .false./" namelist.input' , shell=True, executable='/bin/bash')  # no restart
        print '  set no restart'
        subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_write_output_on_startup.*/   config_write_output_on_startup = .true./" namelist.input' , shell=True, executable='/bin/bash')  # write the initial time
        print '  set write initial time'
        subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_start_time.*/   config_start_time = \'0000-01-01_00:00:00\'/" namelist.input' , shell=True, executable='/bin/bash')  # set the start time
        print '  set start time'
        print 'Starting land ice model.'
        subprocess.check_call(runLICEcmd, shell=True, executable='/bin/bash')
        print 'Land ice model execution completed.'
        # Don't need to copy output file
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

        #####################  DOUG, FILL OUT THIS SECTION ##########################


        # Set the namelist file properly for a restart or no
        print 'Setting up ocean namelist file'
        print 'Starting ocean model.'
        subprocess.check_call(runOCNcmd, shell=True, executable='/bin/bash')
        print 'Ocean model completed.'
        # Don't need to copy output file if the namelist is setup properly
        ### Copy output file so we don't lose it when the model restarts - could use Python's shutil module.  TODO Also may want to change format of the string conversion of t
        ##subprocess.check_call('cp ' + ocnOutput + ' output.' + str(t) + '.nc', shell=True, executable='/bin/bash')
        os.chdir('..')
except:
    sys.exit('Ocean model failed!')


for t in np.arange(numCoupleIntervals)+1:
    print '============================================================'
    print 'Starting coupling number ', t

    # =======================================
    # =======================================
    # 1. Get the appropriate files and variables needed by the coupler
    # =======================================
    # =======================================

    # Get the appropriate source and destination files for the input/output of the coupler
    if runLICE:
        
        # Use values from output files as Source at other times.
        try:
          liceSourceFile = NetCDFFile('./land_ice/' + liceOutput, 'r')
        except:
          sys.exit('Problem opening liceSourceFile: '+liceOutput)

        # At other times we want to write to the restart files as the Destination
        # Need to choose the appropriate restart file - there may be smarter ways to do this, but an easy way is to take the most recent
        # os.listdir does not have a sort ability, so use a shell command (unpythonic but works)
        try:
          p = os.popen('ls -1t ./land_ice/restart.*.nc | head -n 1')
          liceRestart = p.read().rstrip()  # get stdout
          liceDestFile = NetCDFFile(liceRestart,'r+')
          print 'For land ice destination file, using:' + liceRestart
        except:
          sys.exit('Problem opening liceDestFile: ' + liceRestart)
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
        liceBasalTemp = liceSourceFile.variables['temperature'][-1,:,-1]
    except:
        print "Problem getting needed Source fields from the land ice file:" + liceOutput
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
        liceBedTopography = liceDestFile.variables['bedTopography'][-1,:]  # The bed topo is really a Source field (not a Destination field), but it is not included in output files, only in restart files!
        liceBMB = liceDestFile.variables['marineBasalMassBalTimeSeries'][:,0]  # There should only be on time level
        #liceBHF = liceDestFile.variables['basalHeatFluxTimeSeries'][:,0]  # NEED TO RESTRICT THIS TO CELLS WITH ICE
    except:
        print "Problem getting needed Destination fields from the land ice files" + liceRestart
    try:
        ocnSurfPressure = ocnDestFile.variables['seaSurfacePressure'][0,:]
        ocnSurfTempFlux = ocnDestFile.variables['surfaceTemperatureFlux'][0,:]
        ocnSurfMassFlux = ocnDestFile.variables['surfaceMassFlux'][0,:]
    except:
        print "Problem getting needed Destination fields from the ocean files"

    # Setup dummy values for the coupler fields if both models are not being run
    if runOCN == False:
        print "Setting up dummy ocean variables since ocean model is not being run."
        fieldshape = liceThickness.shape 
        ocnSurfTemp = np.zeros( fieldshape )
        ocnSurfDensity = np.zeros( fieldshape )
        ocnSurfSalinity = np.zeros( fieldshape )
        ocnSurfPressure = np.zeros( fieldshape )
        ocnBottomLayerDensity = np.zeros( fieldshape )
        ocnBottomLayerPressure = np.zeros( fieldshape )
        ocnBottomLayerThickness = np.zeros( fieldshape )

        ocnSurfPressure = np.zeros( fieldshape )
        ocnSurfTempFlux = np.zeros( fieldshape )
        ocnSurfMassFlux = np.zeros( fieldshape )

    if runLICE == False:
        print "Setting up dummy land ice variables since land ice model is not being run."
        fieldshape = ocnSurfTemp.shape
        liceThickness = np.zeros( fieldshape )
        liceLowerSurface = np.zeros( fieldshape )
        liceBedTopography = np.zeros( fieldshape )
        liceBasalTemp = np.zeros( fieldshape )

        liceBMB = np.zeros( fieldshape )

    # =======================================
    # =======================================
    # 2. Do the coupler calculations
    # =======================================
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
    Tf = a * ocnSurfSalinity + b + c * ocnSurfPressure     # Eq. 1 from Holland and Jenkins, Modeling Thermodynamic IceOcean Interactions at the Base of an Ice Shelf. J. Physical Oceanography, Aug 1999.
    heatFlux = ocnSurfDensity * waterSpecificHeat * exchangeVel * (ocnSurfTemp - Tf)  # upward positive, units: W/m2
    heatFlux[np.logical_not(isShelf)] = 0.0  # Only allow a heat flux where you have an ice shelf
    ocnSurfTempFlux = -exchangeVel * (ocnSurfTemp - Tf)  # downward positive, units: deg C * m/s
    ocnSurfTempFlux[np.logical_not(isShelf)] = 0.0 # Only allow a temp. flux where you have an ice shelf
    
    # Mass Flux: Eq. 4, modified
    massFlux = heatFlux / latentHeatFusion  # positive=melting, units: kg/(m2 s); divide by appropriate density to get thickness rate of change for each model component.
    ocnSurfMassFlux = massFlux / rhow  # units: m of water per second.  Assuming pure freshwater density here - TODO is that right?
    liceBMB = -massFlux / rhoi * secondsInYear  # units: m of ice per yr

    # =======================================
    # =======================================
    # 3. Write the new forcing data to the files
    # =======================================
    # =======================================

    # Write new quantities to the destination files
    # TODO Does this need to be done explicitly?  I always forget how netCDF i/o works in python precisely.
    if runLICE == True:
       liceDestFile.variables['marineBasalMassBalTimeSeries'][:,0] = liceBMB
       liceDestFile.close()
    if runOCN == True:
        ocnDestFile.variables['seaSurfacePressure'][0,:] = ocnSurfPressure
        ocnDestFile.variables['surfaceTemperatureFlux'][0,:] = ocnSurfTempFlux
        ocnDestFile.variables['surfaceMassFlux'][0,:] = ocnSurfMassFlux
        ocnDestFile.close()

    # =======================================
    # =======================================
    # 4. Run the models
    # =======================================
    # =======================================
    # TODO Do we want to do this in parallel?  Probably no reason to do so.


    # === LAND ICE MODEL ===
    try:
        if runLICE:
            os.chdir('land_ice')
            # Set the namelist file properly for a restart or no and set the start and stop times
            print 'Setting up land ice namelist file'
            subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_do_restart.*/   config_do_restart = .true./" namelist.input' , shell=True, executable='/bin/bash') # do a restart
            print '  set restart'
            subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_write_output_on_startup.*/   config_write_output_on_startup = .true./" namelist.input' , shell=True, executable='/bin/bash')  # don't write the initial time
            print '  disabled initial write'
            # Get the old stop time so we can set it to be the new start time. - it seems more reliable to search for the most recent restart file rather than the most recent output file.
            p = os.popen('ls -1t restart*nc | head -n 1')
            lastrestartfilename=p.read().rstrip()
            p = os.popen('ncks -v xtime -H -s "%c" ' + lastrestartfilename + ' | tail -r -c 65')
            oldstoptime=p.read().rstrip()
            print '  For land ice, getting the restart time from file: ' + lastrestartfilename + ' which gives a time of ' + oldstoptime
            subprocess.check_call('sed -i.SEDBACKUP "s/^.*config_start_time.*/   config_start_time = ' + oldstoptime + '/" namelist.input ' , shell=True, executable='/bin/bash')  # set the start time to the old stop time
            print '  set start time'

            print 'Starting land ice model.'
            subprocess.check_call(runLICEcmd, shell=True, executable='/bin/bash')
            print 'Land ice model completed.'
            # Don't need to copy output file if the namelist is setup properly
            ### Copy output file so we don't lose it when the model restarts - could use Python's shutil module.  TODO Also may want to change format of the string conversion of t
            subprocess.check_call('cp ' + liceOutput + ' ./output-files/output.' + '{0:04d}'.format(t+1) + '.nc', shell=True, executable='/bin/bash')  # the number is the ending time level (if there are more than one)
            subprocess.check_call('cp log.0000.out ./output-files/log.0000.' + '{0:04d}'.format(t+1) + '.out', shell=True, executable='/bin/bash')  # the number is the ending time level (if there are more than one)
            subprocess.check_call('cp log.0000.err ./output-files/log.0000.' + '{0:04d}'.format(t+1) + '.err', shell=True, executable='/bin/bash')  # the number is the ending time level (if there are more than one)

            os.chdir('..')
        else:
            print 'Land ice model is not selected for execution'
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
        else:
            print 'Ocean model is not selected for execution'
    except:
        sys.exit('Ocean model failed!')


# =======================================

print 'Coupled model run complete!'


# =======================================
# =======================================
# 5. Post-processing
# =======================================
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



