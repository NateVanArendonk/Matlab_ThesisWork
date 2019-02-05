import arcpy, os
from arcpy import env


#Deltas = ['Nooksack', 'Skagit', 'Stillaguamish', 'Snohomish', 'Nisqually', 'Skokomish']
Deltas = ['Snohomish']


#Directory that holds output GDB from connectivity script
##Nooksack_Dir = r'E:\Projects\Deltas\Nooksack\Analyses\InundationAnalyses\Connectivity\Nook_1\Nooksack_Inundation.gdb'
##Skagit_Dir = r'E:\Projects\Deltas\Skagit\Analyses\InundationAnalyses\Connectivity\DEM_Mos14\output\Skagit_Inundation_Mos14.gdb'
##Stillaguamish_Dir = r'E:\Projects\Deltas\Stillaguamish\Analyses\InundationAnalyses\Connectivity\Still_1\Still_Inundation.gdb'
Snohomish_Dir = r' E:\Projects\Deltas\Snohomish\Analyses\InundationAnalyses\Connectivity\output\Snohomish_Inundation.gdb'
##Nisqually_Dir = r'E:\Projects\Deltas\Nisqually\Analyses\InundationAnalyses\Connectivity\SLR\Nisqually_Inundation.gdb'
##Skokomish_Dir = r'E:\Projects\Deltas\Skokomish\Analyses\InundationAnalyses\Connectivity\Skok_1\Skokomish_Inundation.gdb'

#inDir = [Nooksack_Dir, Skagit_Dir, Stillaguamish_Dir, Snohomish_Dir, Nisqually_Dir, Skokomish_Dir]
inDir = [Snohomish_Dir]
outDir = os.getcwd() # -- outDir will be location of script


#Cellsize for grid resampling
cellsize = '5'

#Projection Files for Reprojecting from UTM to Web Mercator
webMerc_prjFile = r"C:\ArcGIS\Desktop10.0\Coordinate Systems\Projected Coordinate Systems\World\WGS 1984 Web Mercator.prj"
utm_prjFile = r"C:\ArcGIS\Desktop10.0\Coordinate Systems\Projected Coordinate Systems\UTM\NAD 1983\NAD 1983 UTM Zone 10N.prj"

#Create gdb

for i in range(len(Deltas)):



    outGDB = Deltas[i] + '_Inundation.gdb'

    #Create Geodatabase 
    if arcpy.Exists(outGDB):
        print ' -- GDB \'' + outGDB + '\' Exists -- '
    else:
        print ' -- Creating GDB: \'' + outGDB + '\' -- '
        arcpy.CreateFileGDB_management(outDir, outGDB)



#Resample and Reproject Inundation Grids to Web Mercator and place into GDB

for i in range(len(Deltas)):

##    if Deltas[i] != 'Skagit':
##        continue

    print '\n - ' + Deltas[i] + ' - \n'

    inGDB = inDir[i]
    env.workspace= inGDB

    grids = arcpy.ListDatasets('*_depth')

    for grid in grids:

        print 'Reprojecting: ' + grid + ' -- Resampling to ' + cellsize + ' meters --'
        
        inGrid = grid
        outGrid = outDir + '\\' + Deltas[i] + '_Inundation.gdb' + '\\' + grid + '_Web'

        if arcpy.Exists(outGrid):
            print ' -- \'' + outGrid + '\' Already  Exists -- '
        else:
            arcpy.ProjectRaster_management(inGrid, outGrid, webMerc_prjFile ,"NEAREST", cellsize, "NAD_1983_To_WGS_1984_5;WGS_1984_Major_Auxiliary_Sphere_To_WGS_1984", "", utm_prjFile) # params2 )

            print ' -- \'' + outGrid + '\' Created' + ' -- '



        
        
