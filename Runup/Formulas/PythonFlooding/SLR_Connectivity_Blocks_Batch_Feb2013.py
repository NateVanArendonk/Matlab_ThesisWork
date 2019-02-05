#-------------------------------------------------------------------------------#
#
# SLR_Connectivity_Blocks_Batch.py
#
# Evaluates hydrologic connectivity from a known water source
#---------------------------------------------------------------------------------------------------------
#
#  ** USER INPUTS ARE AT THE BOTTOM OF THE SCRIPT IN THE  if __name__ == '__main__':  SECTION **
#
#
#  Input Grids -- ESRI Grid Format
#  -------------------------------
#  - DEM - an elevation model for the study area
#  - Water - a binary grid identifying source of water (i.e. the ocean, sound, etc.)
#          - water cells are 1, land cells are 0
#   ** NOTE: Both input grids must have the same cell size, llcorner, number of rows, and number of colums
# 
#
#  Output Files -- ESRI ASCII Grid Format
#  --------------------------------------
#  - User Defined Name - binary grid with water cells = 1 and dry cells = 0
#
#
#----------------------------------------------------------------------------------------------------------


import numpy, sys, os, time, arcpy
from gdalconst import *
from numpy import *

try:
    import gdal
except importError:
    from osgeo import gdal


def openFiles(DEM, Water):

    print "Opening Data Files . . . \n"

    #DEM GRID
    inDEM = gdal.Open(DEM)

    ## Get Cell Info for DEM Grid
    inDEM_rows = inDEM.RasterYSize
    inDEM_cols = inDEM.RasterXSize
    inDEM_xMin = inDEM.GetGeoTransform()[0]
    inDEM_yMax = inDEM.GetGeoTransform()[3]
    inDEM_cellSize = inDEM.GetGeoTransform()[1]
    inDEM_yMin = inDEM_yMax - (inDEM_rows*inDEM_cellSize)

    ## WATER GRID
    inWater = gdal.Open(Water)

    inWater_rows = inWater.RasterYSize
    inWater_cols = inWater.RasterXSize
    inWater_xMin = inWater.GetGeoTransform()[0]
    inWater_yMax = inWater.GetGeoTransform()[3]
    inWater_cellSize = inWater.GetGeoTransform()[1]
    inWater_yMin = inWater_yMax - (inWater_rows*inWater_cellSize)

    ##COMPARE GRIDS . . .  NEED TO ADD THIS ##

    return inDEM, inWater, inDEM_rows, inDEM_cols, inDEM_xMin, inDEM_yMin, inDEM_cellSize



def SLR_Connectivity(dm, wtr, SLR):

    # Find Water Cells
    #wtrCells = where(wtr == 1,1,0)
    #dmCells = where(dm>SLR-2, 1, 0)
    #wtrCells2 = wtrCells + dmCells
    #slrCells = where(wtrCells2 == 2)

                    
    slrCells = where((wtr == 1)) # & (dm > -11))

    change = 0 #Will count the number of water cells added

    if any(slrCells):
        rows = slrCells[0]
        cols = slrCells[1]
        print "        Number of Water Cells in this Block: " + str(len(rows))
    else:
        print "        No Water Cells in this Block"
        return wtr, 0


    #Find lowCells that are not already water
    lowCells = where(dm< SLR, 1, 0)

    if not any(lowCells):
        print "        No Cells Below Elevation Threshold in this Block"
        return wtr, 0
    else:
        lowCells = where(wtr==1, 0, 1) + lowCells

        lowNum = len(where(lowCells == 2)[0])

        if any(lowCells):
            print "        Number of Low and Dry Cells in this Block: " + str(lowNum)
        else:
            print "        No Cells Below Elevation Threshold and Dry in this Block"
            return wtr, 0
    

    #If there are fewer dry cells, then use those to find boundary water cells
    if lowNum < len(rows):

        print "        Using Low/Dry Cells to find water boundary"

        dryRows = where(lowCells == 2)[0]
        dryCols = where(lowCells == 2)[1]

        newRows = array([], int32)
        newCols = array([], int32)

        for cell in range(len(dryRows)):
        
            i = dryRows[cell]
            j = dryCols[cell]

            if i == 0:
                kStart = i
            else:
                kStart = i - 1

            if j == 0:
                lStart = j
            else:
                lStart = j - 1


            for k in range(kStart, i+2):
                for l in range(lStart, j+2):

                    try:
                        if (wtr[k, l] == 1):
                    
                            newRows = append(newRows, k)
                            newCols = append(newCols, l)

                    except:
                        pass

        rows = newRows
        cols = newCols

    # - End Dry Cell Check - #


    #Iterate through rows and cols arrays to evaluate cells
    iteration = 1

    while rows.any():

        newRows = array([], int32)
        newCols = array([], int32)
           
        for cell in range(len(rows)):

            i = rows[cell]
            j = cols[cell]
                

            #The following if statements prevent negative indexing
            if i == 0:
                kStart = i
            else:
                kStart = i - 1

            if j == 0:
                lStart = j
            else:
                lStart = j - 1

            #Determine if surrounding cell meets criteria.
            #If so, change that cell to 1 and add to cell list
            for k in range(kStart, i+2):
                for l in range(lStart, j+2):

                    try:
                        if (wtr[k, l] == 0):
                            if (dm[k, l] < SLR):

                                wtr[k, l] = 1
                    
                                newRows = append(newRows, k)
                                newCols = append(newCols, l)

                                change += 1

                    except:
                        pass

        iteration+=1
        
        rows = newRows
        cols = newCols
           
    return wtr, change




#Determine how many blocks there are to process in the grid based on user defined block size
def blockNums(inArray, xBlockSize, yBlockSize):

    rows = len(inArray)
    cols = len(inArray[0])


    if rows%yBlockSize == 0:
        rowBlocks = int(rows/yBlockSize)
    else:
        rowBlocks = int(round(rows/yBlockSize)) + 1

    if cols%xBlockSize == 0:
        colBlocks = int(cols/xBlockSize)
    else:
        colBlocks = int(round(cols/xBlockSize)) + 1

    totalBlocks = rowBlocks * colBlocks

    return totalBlocks


#Iterate through blocks of dem grid
def slrBlocks(inDEM, inWater, rows, cols, xBlockSize, yBlockSize, SLR):

    #Create Water Output Array
    outWater = inWater.ReadAsArray()

    #Calculate Total Number of Blocks to Process
    totalBlocks = blockNums(outWater, xBlockSize, yBlockSize)

    print 'Total Blocks to Process:  '+ str(totalBlocks)
    

    #Read Input Raster By Block

    change = 1
    Iteration = 1

    while change == 1:

        print '\nIteration: ' + str(Iteration) + '\n'

        #Re-Set change to zero
        change = 0

        block = 1
        
        for i in range(0, rows, yBlockSize):  #Iterate through row blocks
           
            #Determine if block is last in the series
            if i + yBlockSize < rows:
                numRows = yBlockSize
            else:
                numRows = rows - i


            #Make Blocks Overlap
            if i > 0:
                rowStart = i - 5
                numRows = numRows + 5
            else:
                rowStart = i


            for j in range(0, cols, xBlockSize):  #For each row block, iterate through column blocks


                if block % 1 == 0:
                    print "    Processing Block: " + str(block)

                block += 1
                
                #Determine if block is last in the series
                if j + xBlockSize < cols:
                    numCols = xBlockSize
                else:
                    numCols = cols - j

                #Make Blocks Overlap
                if j > 0:
                    colStart = j - 5
                    numCols = numCols + 5

                else:
                    colStart = j


                ##Read data by block

                dm = array([], float32)
                wtrIn = array([], int32)
                wtrOut = array([], int32)
                
                #Get Blocks of data
                dm = inDEM.ReadAsArray(colStart, rowStart, numCols, numRows)
                wtrIn = outWater[rowStart:rowStart+numRows, colStart:colStart+numCols]

                ##Send data blocks to connectivity function and save results to outWater Array
                wtrOut, wtrChange = SLR_Connectivity(dm, wtrIn, SLR)

                #If there is a change in the wtr array, then set 'change' to 1
                # While loop with continue until there are no more changes in any of the blocks
                
                if wtrChange > 0:
                    change = 1
                    print '     ' + str(wtrChange) + ' Water Cells Added'

                #if any(logical_and(abs(wtrOut-wtrIn), 1)):  #If there is no change, then wtrIn - wtrOut will = 0
                #    change = 1                              #Otherise, wtrIn- wtrOut will have at least one value of 1
                #    print change

                #THIS WAS THE ORIGINAL CODE FOR CHANGE SWITCH:

                # -!- NEED TO REHINK THIS PHRASE -- MAYBE SHOULD == TRUE?? -!- #
                #if (wtrIn - wtrOut).all() == False:
                #   change = 1

                outWater[rowStart:rowStart+numRows, colStart:colStart+numCols] = wtrOut
                
        Iteration += 1
        
    return outWater



def writeAsciiGrid(inArray, outFile, rows, cols, xMin, yMin, cellSize):

    print "\nWriting Ascii Grid: '" + outFile + "'\n"

    #Open outfile for writing
    fout = open(outFile, 'w')

    #Print Header Info

    txt = "ncols         " + str(cols) + "\n"
    fout.writelines(txt)
    txt = "nrows         " + str(rows) + "\n"
    fout.writelines(txt)
    txt = "xllcorner     " + str(xMin) + "\n"
    fout.writelines(txt)
    txt = "yllcorner     " + str(yMin) + "\n"
    fout.writelines(txt)
    txt = "cellsize      " + str(cellSize) + "\n"
    fout.writelines(txt)
    txt = "NODATA_value  -9999\n"
    fout.writelines(txt)


    #Print Array to File
    for i in range(len(inArray)):
        inArray[i].tofile(fout, sep = " ", format = '%s')
        fout.write('\n')

    fout.close()
    


def SLR_Sim(DEM, Water, outFile, xBlockSize, yBlockSize, SLR, Nodata = -3.4028235e+38):

    ##Open Input Files - open the raster files for reading with gdal - get raster properties
    inDEM, inWater, rows, cols, xMin, yMin, cellSize  = openFiles(DEM, Water)

    #Analyze Input Rasters by Block
    outWater = slrBlocks(inDEM, inWater, rows, cols, xBlockSize, yBlockSize, SLR)

    #Write to outFile
    writeAsciiGrid(outWater, outFile, rows, cols, xMin, yMin, cellSize)

    

if __name__ == '__main__':


    startTime = time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime())
    print "Start:  " +  startTime + "\n"
    sys.stdout.flush()

    #--------------------------------------------------------------------------------#

    ### USER DEFINITIONS ###

    ##Infiles

    #Identify the DEM - the required file in the Arc Grid folder and labeled 'w001001.adf'
    DEM = r"E:\Projects\Deltas\Skagit\Analyses\InundationAnalyses\Connectivity\DEM_Mos14\Grids\dem_mos14_cl\w001001.adf"

    #Identify the grid that holds water cells -  Water = 1, Dry = 0
    Water = r'E:\\Projects\\Deltas\\Skagit\\Analyses\\InundationAnalyses\\Connectivity\\DEM_Mos14\\scratch\\sk_1\\w001001.adf'
 
    #Identify the directory for printing ASCII Grids
    outDir = r"E:\Projects\Deltas\Skagit\Analyses\InundationAnalyses\Connectivity\DEM_Mos14\output"

    #Identify the directory to use for scratch workspace
    scratchDir = r"E:\Projects\Deltas\Skagit\Analyses\InundationAnalyses\Connectivity\DEM_Mos14\scratch"

    #Prefix for Ascii Grid outfiles - outfiles will be named: prefix + SLR_Name + '.txt'
    prefix = "\\sk_"

    #Set Sea Level Value for Analysis
    #Names must be a list of strings in the same order as the SLR values
    #SLR is the NAVD88 height that corresponds to the scenario name - a list of floating point values
    # NOTE: The SLR values must be in ascending order

    SLR_Names = ['MHHW_2010', 'MHHW_2050_Mean', 'HOW_2010', 'MHHW_2050_Max', 'HOW_2050_Mean', 'HOW_2050_Max', 'MHHW_2100_Mean', 'HOW_2100_Mean', 'MHHW_2100_Max', 'HOW_2100_Max']
    SLR = [2.75, 3.04, 3.05, 3.21, 3.34, 3.51, 3.6, 3.9, 4.05, 4.35]


    #Can use the following code to generate SLR values that increment by constant value
    
    #Sea Level Values will increment by 0.1m
    #Names will be numeric 1 - 38
##    i = -0.7
##    SLR = []
##    SLR_Names = []
##    while i <= 4.4:
##        SLR.append(i)
##        SLR_Names.append(str(len(SLR)))
##        i += 0.1


    #Set NoData Value - this is standard for ESRI grids
    Nodata = -3.4028235e+38


    #Define Size of Data Blocks - The script processes the grid in blocks - they are generally too large to keep in memory.
    #Should not really need to change the block size unless you want to just play around with it.  At some point the script will crash if the blocks are too large.
    xBlockSize = 2000 #473
    yBlockSize = 2000 #496

    

##   --  DO NOT CHANGE CODE BELOW  --   ##
##########################################

    if not arcpy.Exists(scratchDir):
        os.mkdir(scratchDir)
        

    for i in range(2,len(SLR)):
        
        print "Evaluating Connectivity for " + SLR_Names[i] + ". . . \n"
        outFile = outDir + prefix + SLR_Names[i] + ".txt"


        SLR_Sim(DEM, Water, outFile, xBlockSize, yBlockSize, SLR[i], Nodata = -3.4028235e+38)

        #Save outFile as grid in scratch directory
##        outGrid = scratchDir + os.sep + "sk_" + str(i)  #SLR_Names[i]
##        print "Saving Raster . . . \n"
##        arcpy.ASCIIToRaster_conversion (outFile, outGrid, "Integer")
##
##        Water = outGrid + os.sep + "w001001.adf"

    



    #Print Time Message to Screen
    endTime = time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime())

    print "\n"
    print "Start Time: " + startTime
    print "End Time: " + endTime + "\n"
    sys.stdout.flush()



