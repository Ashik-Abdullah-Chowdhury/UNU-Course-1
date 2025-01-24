// SAR SENTINEL - 1 DATA FILTERING 
 //import study area shapefile
var ctg = table.filter(ee.Filter.eq("ADM1_EN","Chittagong"));
var upzla = table.filter(ee.Filter.eq("ADM3_EN","Teknaf"));
Map.centerObject(upzla, 10);
Map.addLayer(ctg);
// Import Sentinel-1 SAR GRD
var image = ee.ImageCollection('COPERNICUS/S1_GRD')
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) 
// Filter to get images with VV polarization “VV” stands for vertical transmit, vertical recieved.
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
// Filter to get images with VH polarization “VH” stands for vertical transmit, horizontal recieved.
.filter(ee.Filter.eq('instrumentMode', 'IW'))  // Filter to get images collected in interferometric wide swath mode.
.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')) // Also filter based on the orbit: descending or ascending mode
//.select('VV', 'VH', 'angle')
.filterMetadata('resolution_meters', 'equals' , 10)
.filterBounds(upzla) 
.filterDate('2022-08-01', '2022-08-31');
//.first();
//Print Collection and add it to the map
print(image);
print(image.size());
Map.addLayer(image, {min: -25, max: 0}, 'Search Sentine-1', true);
//mosaic Images
var mosaicS1 = image.mosaic();
print(mosaicS1);
Map.addLayer(mosaicS1.clip(upzla), {min: -25, max: 0}, 'Mosaic Sentinel 1', true);
//Select Bands from mosaic Collection and and calculate their difference for display
//(it will revert the colors)
var VV  = (mosaicS1.select('VV')).rename('VV');
var VH  = (mosaicS1.select('VH')).rename('VH');
var diff = ((VH).divide(VV)).rename('diffHV'); //calculate difference for different and better display (VH-VV), and for default (VV-VH)
var S1 = VH.addBands(VV).addBands(diff); //add bands
var S1_Viz = {min: -25, max: 5}; //visual parameters
//add it to the map by clipping it
Map.addLayer(S1.clip(upzla), S1_Viz, 'Clipped S1 Teknaf', true);