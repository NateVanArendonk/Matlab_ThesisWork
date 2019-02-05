import arcpy
from arcpy import env
from arcpy.sa import *

arcpy.CheckOutExtension("Spatial")

workingDir = r'E:\Projects\CR_WebFiles\SLR_Scenarios'

inGDB = 'Skagit_Inundation.gdb'
outGDB = 'Skagit_Inundation_Reclass.gdb'

inGDB = 'Snohomish_Inundation.gdb'
outGDB = 'Snohomish_Inundation_Reclass.gdb'

#Define remap
remap = [[0,1,1],[1,2,2],[2,3,3],[3,4,4],[4,9999,5]]


#Create Geodatabase if one does not already exist
if arcpy.Exists(workingDir + '\\' + outGDB):
    print ' -- GDB \'' + outGDB + '\' Exists -- '
else:
    print ' -- Creating GDB: \'' + outGDB + ' -- \''
    arcpy.CreateFileGDB_management(workingDir, outGDB)



#Get grids in input directory

env.workspace = workingDir + '\\' + inGDB

grids = arcpy.ListDatasets()

for grid in grids:

    inGrid = grid
    outGrid = outGDB + '\\' + inGrid.split('\\')[-1] + '_RC'

    print 'Reclassifying: ' + inGrid.split('\\')[-1]

    #Reclassify Depth Grid
    outGrid_tmp = Reclassify(inGrid, "Value", RemapRange(remap), "NODATA")
    outGrid_tmp.save(outGrid)


