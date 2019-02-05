
# Inundation_GridProcessing.py
#
# Creates a geodatabase with raster datasets
# of the Ascii grids created by the connectivity script
#
#########################################################

import arcpy
from arcpy import env
from arcpy import sa
arcpy.CheckOutExtension("spatial")




#-- USER INPUT --#

#Working Directory - Directory with the Ascii grid output from connectivity code
workingDir = r'E:\Projects\Deltas\Skagit\Analyses\InundationAnalyses\Connectivity\DEM_Mos14\output'


#Name of GeoDatabase to be Created
gdbName = 'Skagit_Inundation_Mos14'


#Prefix for outfiles from Connectivity Script - same as designated in connectivity code
prefix = 'sk_'
suffix = ''  # can also have a suffix, but usually blank


#Tide Scenario List

#These lists are the same as the scenario and SLR lists in the connectivity code
#Used to calculate depths for each scenario

datumList = ['MHHW_2010', 'MHHW_2050_Mean', 'HOW_2010', 'MHHW_2050_Max', 'HOW_2050_Mean', 'HOW_2050_Max', 'MHHW_2100_Mean', 'HOW_2100_Mean', 'MHHW_2100_Max', 'HOW_2100_Max']
datumHeight = [2.75, 3.04, 3.05, 3.21, 3.34, 3.51, 3.6, 3.9, 4.05, 4.35] #Tide Heights - in same order as tideList


#DEM For Creating Inundation Depth Grids
#This is the same as the DEM used in the connectivity code - just point to the folder of the grid rather than the .adf file (since this script will run arcpy)
DEM = r"E:\Projects\Deltas\Skagit\Analyses\InundationAnalyses\Connectivity\DEM_Mos14\Grids\dem_mos14_cl"


#Mask Directory - if you want to clip the water grids by masks, point to the directory that holds the masks
#maskDir = r"F:\DeltasAnalyses\Nisqually\GIS\NisquallyProjectAreas"

#Mask List -- areas to be clipped out of grid for analysis - masks assumed to be shapefiles
#maskList = ['NNWR', 'Phase1', 'Phase2', 'Pilot']


# -- END USER INPUT -- #





def asciiGrid_to_Raster(workingDir, gdbName, prefix, datumList):

    print ' \n -- Converting Ascii Grids -- \n'

    for i in range(len(datumList)):

        inFile = workingDir + '\\' + prefix + datumList[i] + suffix + '.txt'
        print 'Converting ' + datumList[i] + ' Grid: ' + inFile
        outGrid = workingDir + '\\' + gdbName + '.gdb\\' + prefix + datumList[i]
        if arcpy.Exists(outGrid):
            print '  - File \'' + inFile + '\' Already Exists'
        else:
            arcpy.ASCIIToRaster_conversion(inFile, outGrid, "INTEGER")

    print ' \n -- End Grid Conversion -- \n'



def createDepthGrids(workingDir, gdbName, prefix, datumList, datumHeight, DEM):

    
    for i in range(len(datumList)):

        outDir = workingDir + '\\' + gdbName + '.gdb\\'

        inGrid = outDir + prefix + datumList[i]
        outGrid = outDir + prefix + datumList[i] + '_Depth'

        print 'Creating Depth Grid: ' + prefix + datumList[i] + '_Depth'

        if arcpy.Exists(outGrid):
            print ' - ' + gdbName + '\\' + datumList[i] + '_Depth' + ' Already Exists'
        else:
            inExpression = "Con(\"" + inGrid + "\"==1," + str(datumHeight[i]) + "-\"" + DEM + "\")"  #-\"" + DEM + "\")"
            print inExpression
            arcpy.gp.RasterCalculator_sa(inExpression, outGrid)
            print " Created"

    print ' \n -- End Depth Grids -- \n'


def clipWaterGrids(workingDir, gdbName, prefix, maskDir, maskList, datumList):

    print ' \n -- Clipping Grids to Study Areas -- \n'

    for i in range(len(maskList)):

        mask = maskDir + '\\' + maskList[i] + '.shp'
        outDir = workingDir + '\\' + gdbName + '.gdb\\'


        print 'Clipping to mask: ' + mask

        for j in range(len(datumList)):

            waterGrid = workingDir + '\\' + gdbName + '.gdb\\' + prefix + datumList[j]

            outWaterGrid = outDir + '\\' + maskList[i] + '_' + datumList[j] 


            if arcpy.Exists(outWaterGrid):
                print ' - ' + maskList[i] + '_' + datumList[j] + ' Already Exists'
            else:
                arcpy.Clip_management(waterGrid, "#", outWaterGrid, mask, "#", "ClippingGeometry")
                print   ' - ' + gdbName + '.gdb\\' + prefix + datumList[j] + ' clipped to: ' + maskList[i] + '_' + datumList[j]
                print ' - Rebuilding Attribute Table'
                arcpy.BuildRasterAttributeTable_management(outWaterGrid, "Overwrite")
               

                

    print ' \n -- End Grid Clipping -- \n'




def getCellArea(inGrid):
    
    #Get cellsize and calculate cell area (assumes square cells)
    cellSize = arcpy.GetRasterProperties_management(inGrid, 'CELLSIZEX')
    cellSize = round(float(cellSize[0]))
    cellArea = pow(int(cellSize), 2)

    return cellArea




def calcWaterArea(waterGrid, cellArea):

    Rows = arcpy.SearchCursor(waterGrid)

    for row in Rows:
        val = row.getValue('Value')

        count = 0
        area = 0
        
        if val == 1:

            count = row.getValue('Count')
            area = count * cellArea
            break

    #Area returned in units of input grid
    return count, area




def calcAreaForWaterGrids(workingDir, gdbName, prefix, maskDir, maskList, datumList, datumHeight):

    print ' -- Calculating Inundation Areas -- \n'


    outFile = workingDir + '\\' + workingDir.split('\\')[-1] + '_' + "waterAreas.csv"  #csv File to print results

    #Open outFile for writing and print header
    fout = open(outFile, 'w')

    print >> fout, "Tide_Datum, Datum_Height_m, Mask, Area_sqM, Area_sqKM, Cell_Count, Cell_Area_sqM, GridFile"

        

    for j in range(len(maskList)):

        print "--- Processing " + maskList[j] + " ---"


        for i in range(len(datumList)):

            waterGrid = workingDir + '\\' + gdbName + '.gdb\\' + maskList[j] + '_' + datumList[i]

            print " - Processing Grid: " + gdbName + '.gdb\\' + maskList[j] + '_' + datumList[i]
                
            cellArea = getCellArea(waterGrid)
            count, area = calcWaterArea(waterGrid, cellArea)
            area_km = float(area/1000000.0)

            fout.write(datumList[i] + ',')
            fout.write(str(datumHeight[i]) + ',')
            fout.write(maskList[j] + ',')
            fout.write(str(area) + ',')
            fout.write(str(area_km) + ',')
            fout.write(str(count) + ',')
            fout.write(str(cellArea) + ',')
            fout.write(gdbName + '.gdb\\' + maskList[j] + '_' + datumList[i] + "\n")

            print '  - Water Cells: ' + str(count)
            print '  - Cell Area: ' + str(cellArea)
            print '  - Water Area (sqM): ' + str(area)
            print '  - Water Area (sqKM): ' + str(area_km)

    fout.close()








if __name__ == "__main__":

    #Set Workin Directory
    env.workspace = workingDir
    
    #Create Geodatabase
    if arcpy.Exists(workingDir + '\\' + gdbName + '.gdb'):
        print ' -- GDB \'' + gdbName + '\' Exists -- '
    else:
        print ' -- Creating GDB: \'' + gdbName + '\''
        arcpy.CreateFileGDB_management(workingDir, gdbName)


    #Create Grids From Ascii Grids
    #Place Grids into Geodatabase
    asciiGrid_to_Raster(workingDir, gdbName, prefix, datumList)
    

    #Create Inundation Depth Grids
    createDepthGrids(workingDir, gdbName, prefix, datumList, datumHeight, DEM)
    

    #Uncomment the following 2 functions if you want to clip the grids by maks and calculate inundation areas within masks.

    #Clip Grids to Study Areas by Masks
    #clipWaterGrids(workingDir, gdbName, prefix, maskDir, maskList, datumList)


    #Calculate Inundation Area
    #calcAreaForWaterGrids(workingDir, gdbName, prefix, maskDir, maskList, datumList, datumHeight)


    

    
