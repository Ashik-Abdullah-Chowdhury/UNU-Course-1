var AOI=ee.Geometry.Polygon([[91.51530777374853,22.886944339801406],
                            [91.26262222687353,22.648883140435398],
                            [91.60319839874853,22.47133475658129],
                            [91.72404800812353,22.653952589145867],
                            [91.51530777374853,22.886944339801406]])
print(AOI) 
Map.centerObject(AOI);

var S1_IW = ee.ImageCollection('COPERNICUS/S1_GRD')
                         .filter(ee.Filter.eq('instrumentMode', 'IW'))
                         .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
                         .filterBounds(AOI);
print(S1_IW)

var S1_IW_first=ee.Image(S1_IW.sort('system:time_start',false).first());
print('S1 IW latest acquired image',S1_IW_first);
Map.addLayer(S1_IW_first, {'bands': 'VV', min: -15,max: 0}, 'Sentinel-1 IW VV', false);
Map.addLayer(S1_IW_first, {'bands': 'VH', min: -25,max: -5}, 'Sentinel-1 IW VH', false);
Map.addLayer(S1_IW_first, {'bands': 'VV,VH,VV', min: [-15, -25, -15],max: [0, -5, 0]}, 'Sentinel-1 IW VV,VH,VV false color',false);

// Seasonal composite
// Create a 3 band stack by selecting from different periods (months)
// this is useful for visualising temporal changes
var s1_VV = S1_IW.select(['VV']);
var VV1 = ee.Image(s1_VV .filterDate('2022-06-01', '2022-06-30').median());
var VV2 = ee.Image(s1_VV .filterDate('2022-07-01', '2022-07-31').median());
var VV3 = ee.Image(s1_VV .filterDate('2022-08-01', '2022-08-31').median());
//Add to map
Map.addLayer(VV1.addBands(VV2).addBands(VV3), {min: -12, max: -7}, 'Season composite', false);