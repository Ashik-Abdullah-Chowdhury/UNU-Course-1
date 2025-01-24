// SAR SENTINEL - 1 DATA FILTERING 
 //import study area shapefile
var ctg = table.filter(ee.Filter.eq("ADM1_EN","Chittagong"));
var upzla = table.filter(ee.Filter.eq("ADM3_EN","Teknaf")).geometry();
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


//***********************************************************************
//CALCULATING FLOOD
function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10; // Bits 10 and 11 are dense clouds and cirrus, respectively.
  var cirrusBitMask = 1 << 11; // Bits 10 and 11 are dense clouds and cirrus, respectively.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0) // Both flags should be set to zero, indicating clear conditions.
               .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}
//Import Sentinel 2 Harmonized Collection
var dataset = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterDate('2022-06-01', '2022-08-31')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',20))
                  .map(maskS2clouds); //apply cloud mask function
//Visual Parameters
var visualization = {
  min: 0.0,
  max: 0.5,
  bands: ['B12', 'B8', 'B4'],
};
 //Add it to the map by creating composite (mean)
Map.addLayer(dataset.mean().clip(upzla), visualization, 'RGB');
//Calculate NDWI
var s2 = dataset.mean().clip(upzla);
var ndwi = s2.normalizedDifference(['B3', 'B8']);
// var ndsiViz = {min: 0, max: 1, palette: ['44c9f1', '1637f1']};
// Map.addLayer(ndsi.clip(aoi), ndsiViz, 'NDSI', false);
var ndwiViz = {min: 0.5, max: 1, palette: ['1637f1', '44c9f1']};
//Add SRTM
var srtm = ee.Image("CGIAR/SRTM90_V4").toInt();
var ndwiMasked = ndwi.updateMask(ndwi.gte(0)).updateMask(srtm.lte(500));
Map.addLayer(ndwiMasked.clip(upzla), ndwiViz, 'Permanent Water');
//Calculate Flood Water
var post_VH = S1.select('VH');
var flood = post_VH.lte(-21).and(srtm.lt(1700));
var flood_mask = flood.updateMask(flood.eq(1)).updateMask(ndwi.lt(0));
var flood_viz = {palette:['Red']};
Map.addLayer(flood_mask.clip(upzla),flood_viz,'Flood_Inundation');

//CALCULATE FLOOD AREA
//Area Calculation 
//var flood_area = ee.Image.pixelArea().divide(10000).addBands(flood_mask).reduceRegion({
  //reducer: ee.Reducer.sum(),
  //geometry: Sind,
  //maxPixels: 1e10,
  //scale: 10
//});
//print('Flood area (ha)', flood_area);
//Area Calculation
var flood_area = flood_mask.multiply(ee.Image.pixelArea().divide(1e6).round());
var area = flood_area.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: upzla,
  scale: 10,
  maxPixels: 1e10
}); //results will be in square meters regardless of the projection of the input image.
print(area);
//****************************************************************************************
//RASTER TO VECTOR CONVERSION
var flood_vectors = flood_mask.reduceToVectors({
  geometry: upzla,
  crs: flood_mask.projection(),
  scale:500,
  geometryType:'polygon',
  eightConnected: false,
  labelProperty: 'zone',
});
Map.addLayer(flood_vectors, {}, 'Flood - Raster to Vectors');
//********************************************************************************
Export.image.toDrive({
    image: flood_mask,
    description: 'July Masked Image',
    folder: 'earthengine',
    fileNamePrefix: 'July Masked Image',
    region: upzla,
    scale: 10,
    maxPixels: 1e9
});