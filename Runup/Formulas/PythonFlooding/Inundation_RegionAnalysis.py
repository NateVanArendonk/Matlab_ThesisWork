######################################################################
#
# 'Inundation_RegionAnalysis.py'
#  ----------------------------
#
# Peter Horne, April 2013
#
# Finds inundated areas based on user specified water level values
# Uses the Spatial Analyst tool 'Region Group' to find connected areas
#
# Inputs:
# ------
#  - a geodatabase to place output
#  - a DEM
#  - a list of scenario names (string)
#  - a list of water level heights associated with each scenario (numeric)
#
#
# Outputs:
# -------
#  - a raster dataset with all areas in the dem that are lower than the scenario water height
#      -> cells that are lower than water height have value = 1 in output
#      -> cells that are higher than water height are nodata in output
#      -> output name = '[scenario_name]_all'
#
#  - a raster dataset with all connected regions
#      -> result of Spatial Analyst 'Region Group' function
#      -> output name = '[scenario_name]_regions'
#
#  - a raster with inundated areas
#      -> excludes low areas that are not connected to the main water body
#      -> values are either 1 or nodata
#      -> output name = '[scenario_name]_connected'
#
###########################################################################


import arcpy, os
from arcpy import env
from arcpy.sa import *

arcpy.CheckOutExtension("Spatial")


# ----- User Input ----- #

# Geodatabase to place output files - must include .gdb extension
# If the geodatabase does not already exist, it will be created
workspace = "E:\Projects\Deltas\Snohomish\Analyses\InundationAnalyses\RegionAnalysis\Snohomish_RegionAnalysis.gdb"

# DEM with elevations to query
dem = r'E:\GIS_DATA\Washington\RiverDeltas\Snohomish_Data\Grids\Mosaics\sno_mos3'


# Specify the names of the scenarios
# 'Scenario_Names' must be a list of strings that are in the same order as 'Scenario_Values'
# For example: Scenario_Names = ['MHHW_2010', 'MHHW_2050_Mean', 'HOW_2010', 'MHHW_2050_Max', 'HOW_2050_Mean', 'HOW_2050_Max', 'MHHW_2100_Mean', 'HOW_2100_Mean', 'MHHW_2100_Max', 'HOW_2100_Max']
# The final inundation grid will be Scenario_Name plus '_connected'

Scenario_Names = ['MHHW_2010', 'MHHW_2050_Mean', 'MHHW_2050_Max', 'MHHW_2100_Mean', 'HOW_2010', 'HOW_2050_Mean', 'MHHW_2100_Max', 'HOW_2050_Max', 'HOW_2100_Mean', 'HOW_2100_Max']
                  

# Specify the Values for each scenario
# 'Scenario_Values' must be a list of numbers inthe same order as 'Scenario_Names'
# For example: Scenario_Values = [2.75, 3.04, 3.05, 3.21, 3.34, 3.51, 3.6, 3.9, 4.05, 4.35]

Scenario_Values = [2.76, 3.04, 3.19, 3.61, 3.36, 3.64, 4.00, 3.79, 4.21, 4.60]


# If you want to run a series of scenarios with incrementing heights, you can 
# set up the 'Scenario_Names' and 'Scenario_Values' using a loop like the one below:

## Sea Level Values will increment by 0.1m
## Names will be numeric 1 - 38

##    i = -0.7
##
##    Scenario_Names = []
##    Scenario_Values = []
##
##    while i <= 4.4:
##        Scenario_Values.append(i)
##        Scenario_Names.append(str(len(SLR)))
##        i += 0.1


# ----- End user Input ----- #



# -- Functions -- #

# Select cells that meet criteria
def waterSelect(datumHeight, inGrid, outGrid):

    print '  Finding Areas Lower Than Scenario Height . . . '
    outGrid_tmp = Reclassify(inGrid, "Value", RemapRange([[-9999, datumHeight, 1]]), "NODATA")
    outGrid_tmp.save(outGrid)
    

# Run region analysis
def waterRegions(inGrid, outGrid):

    print '  Assessing Connectivity . . . '
    
    outRgnGrp = RegionGroup(inGrid, "EIGHT", "", "", "")
    outRgnGrp.save(outGrid)


# Find the largest region
# We assume that the connected water area is the largest region
def selectRegion(inGrid):

    rows = arcpy.SearchCursor(inGrid)

    region = []
    count = []

    for r in rows:
        region.append(r.value)
        count.append(r.count)

    index = count.index(max(count))

    return region[index]


#Isolates the region in a separate grid
def isolateRegion(region, inGrid, outGrid):

    print '  Saving Final Inundation Grid . . .'
    outGrid_tmp = Reclassify(inGrid, "Value", RemapRange([[region, region, 1]]), "NODATA")
    outGrid_tmp.save(outGrid)
    


def createDepthGrid(inGrid, outGrid, dem, datumHeight):

    print '  Creating Depth Grid . . .'
    inExpression = "Con(\"" + inGrid + "\"==1," + str(datumHeight) + "-\"" + dem + "\")"  #-\"" + DEM + "\")"
    arcpy.gp.RasterCalculator_sa(inExpression, outGrid)   
    


#Main Sequence of Functions
def main(datumHeight, dem, scenario):


    # 1. Select areas lower than datum height
    outGrid = scenario + '_all'
    waterSelect(datumHeight, dem, outGrid)


    # 2. Run Region analysis to find connected areas

    #rename inGrid and outGrid
    inGrid = outGrid
    outGrid = scenario + '_regions'

    waterRegions(inGrid, outGrid)


    # 3. Select Largest region as the connected water area

    inGrid = outGrid
    region = selectRegion(inGrid)


    # 4. Isolate the largest region in separate grid - this is the final inundate area
    outGrid = scenario + '_connected'
    isolateRegion(region, inGrid, outGrid)


    # 5. Create depth Grids From DEM
    
    inGrid = outGrid
    outGrid = scenario + '_connected_depth'
    createDepthGrid(inGrid, outGrid, dem, datumHeight)


if __name__ == "__main__":

    #Check to see if geodatabase exists
    #If not, then create it
    if not arcpy.Exists(workspace):
        print ' -- Creating GDB: ' + os.path.basename(workspace) + ' -- \n'
        arcpy.CreateFileGDB_management(os.path.dirname(workspace), os.path.basename(workspace))

    #Set working environment
    env.workspace = workspace

    #Run Analysis
    for i in range(len(Scenario_Values)):
        
        print 'Scenario: ' + Scenario_Names[i] + ' = ' + str(Scenario_Values[i])
        main(Scenario_Values[i], dem, Scenario_Names[i])
    

    
    



