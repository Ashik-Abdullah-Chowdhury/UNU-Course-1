var AOI=ee.Geometry.Polygon([[91.51530777374853,22.886944339801406],
                            [91.26262222687353,22.648883140435398],
                            [91.60319839874853,22.47133475658129],
                            [91.72404800812353,22.653952589145867],
                            [91.51530777374853,22.886944339801406]])
print(AOI) 
Map.centerObject(AOI);

var s1=ee.ImageCollection("COPERNICUS/S1_GRD")
        .filterDate('2022-08-01', '2022-08-31')
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        .filterBounds(AOI)
        .select(['VV'])
        .sort('system:time_start',false)
        .first();
print(s1)
Map.addLayer(s1,{bands: ['VV'], min: -25, max: 0}, '2022 first image');

var histogram=ui.Chart.image.histogram({
  image:s1,
  region:AOI,
  scale:10,
  // maxBuckets: 
  minBucketWidth:0.4 ,
  // maxRaw: 
  maxPixels:1e19
})
print(histogram)

var collectionVV = ee.ImageCollection("COPERNICUS/S1_GRD").filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')); 
print('number of VV images', collectionVV.size());
print('VV image information', collectionVV.first());
