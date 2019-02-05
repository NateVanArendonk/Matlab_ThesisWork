# -*- coding: utf-8 -*-
"""
Created on Mon May 14 16:22:31 2018

@author: NVanArendonk
"""


speed = 30
direction = range(0,350,10)
tide = [-2,5.5]

forcing_list = []

## Create a set of lists inside of a list for each forcing set
for ii in range(0,len(tide)):
    for jj in range(0,len(direction)):
        forcing_list.append([30,direction[jj],tide[ii]])
        


swn_file = "E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\Model_Runs\NV_LUT\JDF.swn"
    
current_forcings = forcing_list[ii] # Grab tide and wind data for this model run 
# Note forcing goes [DIR, SPD, WL]

f = open(swn_file, 'w') #open file in write mode
########################## HEADER #########################################   
header = """$ JDF and Strait Georgia Model

$******************************************************************************
    
PROJECT 'SaS' '01'
    
MODE STAT
    
SET MAXERR=2 NAUT
"""
    
f.writelines(header) #write header
f.close()

f = open(swn_file, 'a') 
f.writelines("Set LEvel %.1f" % current_forcings[2])
f.close()

