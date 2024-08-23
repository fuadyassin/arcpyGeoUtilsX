'''
Create instance of the small project gis processing python module
'''
import sys
sys.path.append(r"H:\Basin Operations\F_HYD\HYD\STAFF FOLDERS\FYAS\GIS_scripts")
from Small_Projects_GIS_processing import SmallProjectsGISprocessing
# Set the file path
file_path = r"H:\Basin Operations\F_HYD\HYD\STAFF FOLDERS\FYAS\SmallHydrologyProjects\11-ASSINIBOINE_RIVER_BASIN\H11-1_WSS_SW09-02-25-W2_PeterYoung-WarkenRanchingLtd\ProjectGISinput.txt"
# create instance for the class SmallProjectsGISprocessing
extractor = SmallProjectsGISprocessing(file_path)
# Check the extracted data
metadata = extractor.get_data()
print metadata

'''
# Optional use this section only if you want to clip raster skip if you arlready have
'''
# extracting DEM using provided input for min max of x and y (optional)
clipped_raster_path = extractor.extrdem()

'''
# Use default name clipped_raster.tif or provide name s
'''
# using raster process neccassary hydrological analysis such as flow direction, flow accumulation
custom_clipped_raster = "clipped_raster.tif"
flow_direction_raster, flow_accumulation_raster = extractor.prepareRaster(clipped_raster_name=custom_clipped_raster)
# with default clliped raster file
#flow_direction_raster, flow_accumulation_raster = extractor.prepareRaster()

'''
# Optional use this section only if you want to force the stream shape file
'''
custom_stream_file = "custom_stream.shp"
results = extractor.prepareRasterWithStream(stream_shapefile_name=custom_stream_file)
# With the default stream shapefile name
results = extractor.prepareRasterWithStream()

'''
# You can run and re-run this with updated location on existing flow accumulation
'''
# Define the parameters for the createPourPointAndDelineateWatershed method
pour_point_x = 478335.000   
pour_point_y = 5439520.000  
snap_distance = 10          
output_name = "watershed1"  
# Call the createPourPointAndDelineateWatershed method
extractor.createPourPointAndDelineateWatershed(pour_point_x, pour_point_y, snap_distance, output_name)

'''
# appending your current project, EDA, GDA to your main one that contains all project
'''
# Specify the names of the current project shapefiles
Projcurrent = "watershed1_snapped_pour_point.shp"
EDAcurrent = "watershed1_EDA.shp"
GDAcurrent = "watershed1_GDA.shp"
extractor.appendProjectShapefiles(Projcurrent, EDAcurrent, GDAcurrent)

'''
# Calculate 50% and 70% water availability from the runoff index
'''
# GDA 50 and 70 water availablility
inZoneData = "watershed1_GDA.shp"
zoneField = "Id"
outTable = "InflowRunoff_GDA1.dbf"
extractor.runoffIndex(inZoneData, zoneField, outTable)

# EDA 50 and 70
inZoneData = "watershed1_EDA.shp"
zoneField = "Id"
outTable = "InflowRunoff_EDA1.dbf"
extractor.runoffIndex(inZoneData, zoneField, outTable)


