# -*- coding: utf-8 -*-
"""
Created on Mon May 14 15:39:32 2018

@author: NVanArendonk
"""

import subprocess
import shutil
import os


speed = [5,15]
direction = range(0,350,10)
tide = [-2,5.5]

forcing_list = []

## Create a set of lists inside of a list for each forcing set
for ii in range(0,len(tide)):
    for jj in range(0,len(direction)):
        forcing_list.append([30,direction[jj],tide[ii]])
        
os.chdir('E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\Model_Runs\NV_LUT')

for ii in range(0,len(forcing_list)):
    swn_file = "JDF.swn"
        
    current_forcings = forcing_list[ii] # Grab tide and wind data for this model run 
    # Note forcing goes [SPD, DIR, WL]
    
    f = open(swn_file, 'w') #open file in write mode
    ########################## HEADER #########################################   
    header = """$ JDF and Strait Georgia Model
    
$******************************************************************************
        
PROJECT 'SaS' '01'
        
MODE STAT
        
SET MAXERR=2 NAUT\n"""
        
    f.writelines(header) #write header
    f.close()
    f = open(swn_file, 'a') 
    f.writelines("SET LEVEL %.1f\n" % current_forcings[2])# write tide level 
    f.writelines("COORD SPHERICAL CCM\n")
    
    
    
    ########################## Model Setup #################################### 
    model_setup = """\n$************************************* MODEL SETUP******************************
CGRID CURVI 659 1195 EXC 99.0 CIRCLE 90 0.03 2.0
READ  COOR  1. '../../INP_sph2/jdf_sph_swn.grd' IDLA=4 NHEDF=0 NHEDVEC=1 FREE
    
INP  BOTTOM CURVI 0. 0. 659 1195     EXC -999
READ BOTTOM -1. '../../INP_sph2/jdaf_new.BOT' IDLA=3 NHEDF=0 FREE
    
$ Wind input (no variation)"""
    f.writelines(model_setup)
    f.writelines("\nWIND %.1f %.1f\n" % (current_forcings[0], current_forcings[1])) 
    
    ########################## PHYSICS ########################################
    physics = """\n$****************************************** PHYSICS ****************************
GEN3 KOMEN
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
NUM STOPC dabs=0.00  drel=0.01  curvat=0.001  npnts=99. STAT mxitst=50  &
    alfa=0.001 
    """
    f.writelines(numerics)
    
    ####################### OUTPUT LOCS/RAYS
    out_locs = """\n$************************** OUTPUT LOCATIONS/RAYS*******************************
QUANT HS TMM10 TM01 TM02 DIR DSPR FMIN 0.03 FMAX .5

POINTS 'BND' FILE  '../../INP_sph/bnd_pg.loc'
POINTS 'OBS' FILE  '../../INP_sph/meas_jdf.loc'
POINTS 'R1' FILE  '../../INP_sph/jdaf_ray1.loc'
POINTS 'R2' FILE  '../../INP_sph/jdaf_ray2.loc'
    """
    
    f.writelines(out_locs)
    
    ######################### OUTPUT ##########################################
    out_name = "SpatialJDF_s%d_d%d_t%d" % (current_forcings[0]*10,current_forcings[1],current_forcings[2]*10)
    out_name2 = "_s%d_d%d_t%d" % (current_forcings[0]*10,current_forcings[1],current_forcings[2]*10)
    output = """\n$ ************************************ OUTPUT  *********************************
BLOCK  'COMPGRID' NOHEADER 'RES1/%s.mat' LAYOUT 3 &
    XP YP DEP BOTLEV WATLEV HS HSWELL RTP TMM10 TM01 TM02 DIR DSPR   &
    TPS DHSIGN DRTM01 WIND VEL WLENGTH DISBOT DISSURF DISWCAP GENE   & 
    REDI REDQ REDT PROPA
        
SPEC   'BND' SPEC2D ABS 'RES1/BND_PG%s.SA2'

COMPUTE 

STOP
    """ % (out_name,out_name2)
    f.writelines(output)
    f.close()
    
    # Run SWAN Model
    cwd = os.getcwd()
#    print(cwd)
    print 'Running SWAN Model: Tide: %.1f, Speed: %.1f, Direction: %.1f' % (current_forcings[2],current_forcings[0],current_forcings[1])
    subprocess.check_call('swanrun JDF',shell=True)    
    
    
       
    