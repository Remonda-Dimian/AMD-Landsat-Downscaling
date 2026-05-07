// -----------------------------------------
// Landsat-7 SLC-off Gap Filling using Landsat-5
// -----------------------------------------

// Load Landsat-7 and Landsat-5 images
var l7 = ee.Image('LANDSAT/LE07/C02/T1_L2/LE07_XXXXX');
var l5 = ee.Image('LANDSAT/LT05/C02/T1_L2/LT05_XXXXX');

// Select harmonized reflectance bands
var bands = ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'];

l7 = l7.select(bands);
l5 = l5.select(bands);

// Radiometric harmonization (example coefficients)
var slope = 1.02;
var intercept = -0.001;

var l5_adj = l5.multiply(slope).add(intercept);

// Identify valid L7 pixels
var l7_mask = l7.mask();

// Fill SLC-off gaps
var l7_filled = l7.unmask(l5_adj);

// Preserve original mask structure
l7_filled = l7_filled.updateMask(l7_mask.or(l5_adj.mask()));

// Visualization
Map.addLayer(l7_filled, {min:0, max:0.3}, 'Gap-filled L7');
