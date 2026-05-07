// ============================================
// AMD Part 2: Adaptive Multi-Pixel Downscaling (RGB)
// Eqs. (3)-(6). 30m MS + 15m PAN → 10m RGB.
// PAN is used for texture-adaptive weighting only (not pan-sharpening fusion).
// ============================================

var AMD_PARAMS = {
  targetScale: 10,
  sourceScale: 30,
  sigmaP: 30,
  cvThreshold: 0.05,
  wtMin: 0.8,
  wtMax: 2.0,
  gamma: 1.0,
  kernelHomo: 11,
  kernelHetero: 3
};

function computePANTexture(pan15m) {
  var k3 = ee.Kernel.square(1, 'pixels');
  var mean = pan15m.reduceNeighborhood(ee.Reducer.mean(), k3);
  var sqMean = pan15m.pow(2).reduceNeighborhood(ee.Reducer.mean(), k3);
  var variance = sqMean.subtract(mean.pow(2));
  
  var meanVariance = variance.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: pan15m.geometry(),
    scale: 15, maxPixels: 1e13, bestEffort: true
  }).get('variance');
  
  return variance.divide(ee.Number(meanVariance)).rename('pan_var_norm');
}

function resamplePANvarTo30m(panVar15m) {
  return panVar15m.reduceResolution({reducer: ee.Reducer.mean(), bestEffort: true, maxPixels: 1024})
    .reproject({crs: panVar15m.projection(), scale: 30});
}

function computeCV_multispectral(ms30m) {
  var k3 = ee.Kernel.square(1, 'pixels');
  var mean = ms30m.reduceNeighborhood(ee.Reducer.mean(), k3);
  var sqMean = ms30m.pow(2).reduceNeighborhood(ee.Reducer.mean(), k3);
  var stdDev = sqMean.subtract(mean.pow(2)).sqrt();
  
  var stdNorm = stdDev.pow(2).reduce(ee.Reducer.sum()).sqrt();
  var meanNorm = mean.pow(2).reduce(ee.Reducer.sum()).sqrt();
  
  return stdNorm.divide(meanNorm).rename('CV');
}

function computeTextureWeight(panVar30m, params) {
  return panVar30m.multiply(params.gamma).add(1)
    .clamp(params.wtMin, params.wtMax).rename('wt');
}

function amdKernel(msImage, panImage, bandName, params, aoi) {
  var band = msImage.select(bandName);
  var proj30m = band.projection();
  var proj10m = proj30m.atScale(params.targetScale);
  
  var panVar15 = computePANTexture(panImage);
  var panVar30 = resamplePANvarTo30m(panVar15);
  var wt30 = computeTextureWeight(panVar30, params);
  
  var cv = computeCV_multispectral(msImage);
  var isHomogeneous = cv.lt(params.cvThreshold);
  
  var coords10m = ee.Image.pixelCoordinates(proj10m).rename(['x', 'y']);
  var coords30m = ee.Image.pixelCoordinates(proj30m).rename(['x', 'y']);
  
  function buildAMD(kernelSize) {
    var half = Math.floor(kernelSize / 2);
    var msBands = [], wtBands = [], wpBands = [], names = [];
    
    for (var i = -half; i <= half; i++) {
      for (var j = -half; j <= half; j++) {
        var name = 'n_' + i + '_' + j;
        names.push(name);
        
        var disp = ee.Image.constant([j * params.sourceScale, i * params.sourceScale]);
        
        var msD = band.displace({displacement: disp, projection: proj30m})
          .resample('nearest').reproject({crs: proj10m, scale: params.targetScale}).rename(name);
        msBands.push(msD);
        
        var wtD = wt30.displace({displacement: disp, projection: proj30m})
          .resample('nearest').reproject({crs: proj10m, scale: params.targetScale}).rename(name);
        wtBands.push(wtD);
        
        var c30n = coords30m.displace({displacement: disp, projection: proj30m})
          .resample('nearest').reproject({crs: proj10m, scale: params.targetScale}).rename(['x','y']);
        
        var dx = coords10m.select('x').subtract(c30n.select('x'));
        var dy = coords10m.select('y').subtract(c30n.select('y'));
        var distSq = dx.pow(2).add(dy.pow(2));
        
        var wp = distSq.multiply(-1).divide(2 * params.sigmaP * params.sigmaP).exp().rename(name);
        wpBands.push(wp);
      }
    }
    
    var msNbhd = ee.Image.cat(msBands).rename(names);
    var wtNbhd = ee.Image.cat(wtBands).rename(names);
    var wpNbhd = ee.Image.cat(wpBands).rename(names);
    
    var center10 = band.resample('nearest').reproject({crs: proj10m, scale: params.targetScale});
    var centerRep = center10.rename(names);
    
    // Eq. (4)
    var ws = msNbhd.subtract(centerRep).pow(2).add(1).pow(-1);
    
    // Eq. (3)
    var weights = ws.multiply(wpNbhd).multiply(wtNbhd);
    var num = msNbhd.multiply(weights).reduce(ee.Reducer.sum());
    var den = weights.reduce(ee.Reducer.sum());
    
    return num.divide(den).rename(bandName);
  }
  
  var amdHomo = buildAMD(params.kernelHomo);
  var amdHetero = buildAMD(params.kernelHetero);
  
  var cv10 = isHomogeneous.resample('nearest').reproject({crs: proj10m, scale: params.targetScale});
  
  return amdHomo.updateMask(cv10)
    .unmask(amdHetero.updateMask(cv10.not()))
    .unmask(band.resample('nearest').reproject({crs: proj10m, scale: params.targetScale}), false);
}

function runAMD(msImage, panImage, aoi) {
  var bands = ['B2', 'B3', 'B4'];
  var out = bands.map(function(b) {
    return amdKernel(msImage, panImage, b, AMD_PARAMS, aoi);
  });
  
  return ee.Image.cat(out)
    .copyProperties(msImage, ['system:time_start'])
    .set('AMD_PROCESSED', 1);
}

// Usage
/*
var aoi = ee.Geometry.Point([79.9, 37.1]).buffer(5000);
var l8 = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_145034_20200818');
var ms = l8.select(['SR_B2','SR_B3','SR_B4']).rename(['B2','B3','B4']);
var pan = l8.select(['SR_B8']).rename('pan');

var amdResult = runAMD(ms, pan, aoi);

Export.image.toAsset({
  image: amdResult, description: 'AMD_10m_RGB',
  assetId: 'users/your_username/AMD_10m_RGB', region: aoi, scale: 10, maxPixels: 1e13
});
*/
