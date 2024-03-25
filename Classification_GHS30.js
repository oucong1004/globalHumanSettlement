//----------functions---------//
function ObjectBasedFilter(img) {
  var objectId = img.connectedComponents({
    connectedness: ee.Kernel.plus(1),
    maxSize: 256
  });
  var objectSize = objectId.select('labels').connectedPixelCount({
    maxSize: 1024, eightConnected: true
  });
  var objectArea = objectSize.multiply(ee.Image.pixelArea()).divide(1e6);
  var objectFiltered = objectId.updateMask(objectArea.lt(5));//lt 5km2
  return objectFiltered;
}

//----------Classification of GHS30---------//
//Load the global fishnet in 10 degree
var fishCol = ee.FeatureCollection("FeatureCollection ID of your stored fishnet");
var fishList = fishCol.toList(fishCol.size()); 
// Load ISA_KDE_CA_MP dataset
var year = 1985; //e.g 1985/1990/1995/2000/2005/2010/2015/2020
for (var i=0; i<<=fishList.size(); i++){
  var ISA_KDE_CA_MP = ee.Image('Your stored ISA_KDE_CA_MP path'+i);
  if (i===0){
    var ISA_KDE_CA_MP_imgs = ee.ImageCollection([ISA_KDE_CA_MP]);
  }
  else{
    ISA_KDE_CA_MP_imgs = ISA_KDE_CA_MP_imgs.merge(ISA_KDE_CA_MP);
  }
}
var ISA_KDE_CA_MP_mosaic = ISA_KDE_CA_MP_imgs.mosaic().rename('b0');
// Load GHS-POP dataset
var GHS_POP_year = ee.Image("JRC/GHSL/P2023A/GHS_POP/"+year); 
//Object-based patch cluster
var patchid = ISA_KDE_CA_MP_mosaic.connectedComponents(ee.Kernel.plus(1), 1024); 
patchid = patchid.select('labels').unmask(0, true).add(patchid.select('b0'));
//patch size
var patchSizeFiltered = ObjectBasedFilter(ISA_KDE_CA_MP_mosaic).select('b0');
patchSizeFiltered = ISA_KDE_CA_MP_mosaic.where(patchSizeFiltered.eq(1), 2); //size>=5--1, size<5--2
//patch POP density
var PopWithID = GHS_POP_year.addBands(patchid.select('labels'));
var patchPOPmean = PopWithID.reduceConnectedComponents({
  reducer: ee.Reducer.mean(),
  labelBand: 'labels',
  maxSize: 1024
}).rename('patchPOPmean'); 
//classify
var patchClass = patchSizeFiltered
                .where(patchSizeFiltered.eq(1).and(patchPOPmean.gte(15)), 1) //High population density urban settlement (HDU)
                .where(patchSizeFiltered.eq(1).and(patchPOPmean.lt(15)), 2)  //Low population density urban settlement (LDU)
                .where(patchSizeFiltered.eq(2).and(patchPOPmean.gte(3)), 3)  //High population density rural settlement (HDR)
                .where(patchSizeFiltered.eq(2).and(patchPOPmean.lt(3)), 4)   //Low population density rural settlement (HDR)
                .rename('patchClass');
//Batch output
for (var id=0; id<=fishList.size(); id++){
  var selFish =  ee.Feature(fishList.get(id));
  var roiBound = selFish.geometry();
  var patchClass_fish = patchClass.clip(roiBound).unmask(0, true);
  var exportname = 'patchClass_' + id;
  Export.image.toAsset({
    image: patchClass_fish.updateMask(patchClass_fish), 
    description: exportname, 
    assetId: 'Your stored path' + exportname, 
    region: roiBound, 
    scale: 30,
    maxPixels:1e13
  });
}

