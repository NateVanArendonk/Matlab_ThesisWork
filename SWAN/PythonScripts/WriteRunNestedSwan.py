# -*- coding: utf-8 -*-
"""
Created on Mon May 14 15:39:32 2018

@author: NVanArendonk
"""

import subprocess
import os

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)



speed = range(5,35,5)
direction = range(0,360,10)
#tide = [1]
tide = [-2,-1,0,1,2,3,4,5.5]

forcing_list = []

## Create a set of lists inside of a list for each forcing set
for ii in range(0,len(tide)):
    for mm in range(0,len(speed)):
        for jj in range(0,len(direction)):
            forcing_list.append([speed[mm],direction[jj],tide[ii]])


os.chdir('E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\Model_Runs\NV_LUT')


for ii in range(0,len(forcing_list)):
    current_forcings = forcing_list[ii]
    unique_name = "SpatialPS_s%d_d%d_t%d.mat" % (current_forcings[0]*10,current_forcings[1],current_forcings[2]*10)
    with cd("RES2"):
        if os.path.isfile(unique_name): # If the output already exists, skip it
            print('File {0:s} exists, moving on to next file'.format(unique_name))
        else: # If the output doesn't exist, run the damn thing
            with cd("../"):
                swn_file = "PS.swn"
    
                current_forcings = forcing_list[ii] # Grab tide and wind data for this model run
                    # Note forcing goes [SPD, DIR, WL]
                
                print 'Running SWAN Model: Tide: %.1f, Speed: %.1f, Direction: %.1f' % (current_forcings[2],current_forcings[0],current_forcings[1])
               
    
                f = open(swn_file, 'w') #open file in write mode
                    ########################## HEADER #########################################
                header = """$ South Central PS Model

$******************************************************************************

PROJECT 'PG' '01'

MODE STAT

SET MAXERR=2 NAUT\n"""

                f.writelines(header) #write header
                f.close()
                f = open(swn_file, 'a')
                f.writelines("SET LEVEL %.1f\n" % current_forcings[2])# write tide level
                f.writelines("COORD SPHERICAL CCM\n")



    ########################## Model Setup ####################################
                model_setup = """\n$************************************* MODEL SETUP******************************
CGRID CURVI 1555 524 EXC 99.0 CIRCLE 90 0.03 2.0
READ  COOR  1. '../../INP_sph/pug8_sph_swn.grd' IDLA=4 NHEDF=0 NHEDVEC=1 FREE

INP  BOTTOM CURVI 0. 0. 1555 524     EXC -999
READ BOTTOM -1. '../../INP_sph/PugetSound2.BOT' IDLA=3 NHEDF=0 FREE

$ Wind input (no variation)"""
                boundary_setup = "\n$Wave Boundary Conditions from Salish Sea"
                boundary_forcing_file = "BND_PG_s%d_d%d_t%d.SA2" % (current_forcings[0]*10,current_forcings[1],current_forcings[2]*10)
                boundary_forcing = "RES1/" + boundary_forcing_file
                f.writelines(model_setup)
                f.writelines("\nWIND %.1f %.1f DRAG WU\n" % (current_forcings[0], current_forcings[1]))
                f.writelines(boundary_setup)
                f.writelines("\nBOUND NEST '%s' OPEN\n" % boundary_forcing)
    ########################## PHYSICS ########################################
                physics = """\n$****************************************** PHYSICS ****************************
GEN3GEN3 KOMEN
WCAP KOMEN      
QUAD iquad=2 lambda=0.250000 Cnl4=3.00000e+07
LIMITER ursell=10.0000 qb=1.00000
FRIC JONSWAP cfjon=0.038
BREA CON alpha=1.0 gamma=0.73
TRIAD itriad=1 trfac=0.8 cutfr=2.5
    """
                f.writelines(physics)

    ######################## NUMERICAL PARAMETERS #############################
                numerics = """\n$************************************ NUMERICAL PARAMETERS *********************
PROP BSBT
NUM STOPC dabs=0.00  drel=0.01  curvat=0.001  npnts=99.  STAT mxitst=80  &
  alfa=0.001
    """
                f.writelines(numerics)

    ####################### OUTPUT LOCS/RAYS
                out_locs = """\n$************************** OUTPUT LOCATIONS/RAYS*******************************
QUANT HS TMM10 TM01 TM02 DIR DSPR FMIN 0.03 FMAX 2

POINTS 'BND' FILE '../../INP_sph/bnd_jdaf.loc'
POINTS 'R1' FILE  '../../INP_sph/pg_ray1.loc'
POINTS 'R2' FILE  '../../INP_sph/pg_ray2.loc'
POINTS 'R3' FILE  '../../INP_sph/pg_ray3.loc'
    """

                f.writelines(out_locs)

    ######################### OUTPUT ##########################################
                out_name = "SpatialPS_s%d_d%d_t%d" % (current_forcings[0]*10,current_forcings[1],current_forcings[2]*10)
                output = """\n$ ************************************ OUTPUT  *********************************    
BLOCK  'COMPGRID' NOHEADER 'RES2/%s.mat' LAYOUT 3 &
        XP YP DEP BOTLEV WATLEV HS HSWELL RTP TMM10 TM01 TM02 DIR DSPR   &
        TPS DHSIGN DRTM01 WIND VEL WLENGTH DISBOT DISSURF DISWCAP GENE   & 
        REDI REDQ REDT PROPA

COMPUTE

STOP
    """ % out_name
                f.writelines(output)
                f.close()

    # Run SWAN Model
#                cwd = os.getcwd()
#                 print 'Running SWAN Model: Tide: %.1f, Speed: %.1f, Direction: %.1f' % (current_forcings[2],current_forcings[0],current_forcings[1])
                subprocess.check_call('swanrun PS',shell=True)
