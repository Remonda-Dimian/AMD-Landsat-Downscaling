// ============================================
// AMD Part 1: SLC-off Gap-Filling (RGB)
// Eqs. (1)-(2). Landsat-5 TM → Landsat-7 ETM+.
// ============================================

var CROSS_CAL = {
  'B1': {a: 0.943, b: 0.012}, // Blue
  'B2': {a: 0.956, b: 0.008}, // Green
  'B3': {a: 0.974, b: 0.006}  // Red
};

function crossCalibrateL5toL7(l5Image, coeffs) {
  var bands = ['B1', 'B2', 'B3'];
  var calibrated = bands.map(function(b) {
    var c = coeffs[b];
    return l5Image.select(b).multiply(c.a).add(c.b).rename(b);
  });
  return ee.Image.cat(calibrated)
    .copyProperties(l5Image, ['system:time_start', 'system:index']);
}

function getInvalidPixelMask(l7Image) {
  var qa = l7Image.select('QA_PIXEL');
  return qa.bitwiseAnd(1).eq(0); // 1=valid, 0=invalid/fill
}

function fillSLCGaps(l7Image, l5Collection, windowDays) {
  var t = l7Image.date();
  var start = t.advance(-windowDays, 'day');
  var end = t.advance(windowDays, 'day');
  
  var l5 = l5Collection
    .filterDate(start, end)
    .filterBounds(l7Image.geometry())
    .sort('system:time_start')
    .first();
  
  // Reproject L5 to L7 grid before calibration
  l5 = ee.Image(l5).reproject({crs: l7Image.projection(), scale: 30});
  
  var l5Cal = crossCalibrateL5toL7(l5, CROSS_CAL);
  var validMask = getInvalidPixelMask(l7Image);
  
  var bands = ['B1', 'B2', 'B3'];
  var filled = bands.map(function(b) {
    return l7Image.select(b).updateMask(validMask).unmask(l5Cal.select(b)).rename(b);
  });
  
  return ee.Image.cat(filled)
    .copyProperties(l7Image, l7Image.propertyNames())
    .set('SLC_FILLED', 1);
}

// Usage
/*
var aoi = ee.Geometry.Point([79.9, 37.1]);
var l5Col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(aoi)
  .filter(ee.Filter.lt('CLOUD_COVER', 20))
  .map(function(img) {
    return img.select(['SR_B1','SR_B2','SR_B3','QA_PIXEL'])
              .rename(['B1','B2','B3','QA_PIXEL']);
  });

var l7Col = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterBounds(aoi)
  .filter(ee.Filter.lt('CLOUD_COVER', 20))
  .map(function(img) {
    return img.select(['SR_B1','SR_B2','SR_B3','QA_PIXEL'])
              .rename(['B1','B2','B3','QA_PIXEL']);
  });

var l7Filled = l7Col.map(function(img) { return fillSLCGaps(img, l5Col, 32); });
*/
