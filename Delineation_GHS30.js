//----------functions---------//
// CA_neighborhood
exports.CA_neighborhood = function(windowSize, ISA_old){
  var selKernel = ee.Kernel.square(windowSize); 
  var nscale = ISA_old.projection().nominalScale(); 
  var ISA_Nei = ISA_old.neighborhoodToBands(selKernel);
  var bandName = ISA_Nei.bandNames();
  var idList = ee.List.sequence(0, bandName.length().subtract(1)); 
  var ISA_Nei_Col = ee.ImageCollection(idList.map(function(id){
    return ISA_Nei.select([id]).rename(['b0']).float();})); 
    return ISA_Nei_Col.mean().reproject({crs: 'EPSG:4326', scale: nscale});
}

//Erode
exports.erode = function(img, distance){
  var d = (img.not().unmask(1)
       .fastDistanceTransform(30).sqrt()
       .multiply(ee.Image.pixelArea().sqrt()));
  return img.updateMask(d.gt(distance));
}

//Dilate
exports.dilate = function(img, distance){
  var d = (img.fastDistanceTransform(30).sqrt()
       .multiply(ee.Image.pixelArea().sqrt()));
  return d.lt(distance);
} 

//----------kernel density estimation---------//
//Load the GISD30 dataset
var GISD_8520 = ee.Image("Image ID of your stored GISD30 dataset");
//For year 1985--1, 1990--2, 1995--3, 2000--4, 2005--5, 2010--6, 2015--7, 2020--8
var ISA_year = GISD_8520.gte(1).and(GISD_8520.lte(6)).unmask(0,true).float();
//Load the global fishnet in 10 degree
var fishCol = ee.FeatureCollection("FeatureCollection ID of your stored fishnet");
var fishList = fishCol.toList(fishCol.size()); 
//Batch output
for (var id=0; id<=fishList.size(); id++){
  var selFish =  ee.Feature(fishList.get(id));
  var selFishID = selFish.get('Id');
  var roiBound = selFish.geometry();
  var ISA_fish = ISA.clip(roiBound);
  var selISA = ISA.reproject({crs: GISD_8520.projection()}).reduceResolution({
    reducer:ee.Reducer.mean(),
    maxPixels: 65535
  }); 
  var ISA_KDE = selISA.multiply(100).round().uint8(); 
  var exportname = 'ISA_KDE_' + id;
  Export.image.toAsset({
    image: ISA_KDE, 
    description: exportname, 
    assetId: 'Your stored path' + exportname, 
    region: roiBound, 
    scale: 100,
    maxPixels:1e13
  });
}

//----------initial boundaries delineation---------//
//Load ISA_KDE dataset
for (var i=0; i<<=fishList.size(); i++){
  var ISA_KDE = ee.Image('Your stored ISA_KDE path'+i);
  if (i===0){
    var ISA_KDE_imgs = ee.ImageCollection([ISA_KDE]);
  }
  else{
    ISA_KDE_imgs = ISA_KDE_imgs.merge(ISA_KDE);
  }
}
var ISA_KDE_mosaic = ISA_KDE_imgs.mosaic();
//Batch output
for (var id=0; id<=fishList.size(); id++){
  var selFish =  ee.Feature(fishList.get(id));
  var roiBound = selFish.geometry().buffer(1000).bounds();
  var ISA_fish = ISA_year.clip(roiBound).unmask(0, true);
  var ISA_KDE_fish = ISA_KDE_mosaic.clip(roiBound).unmask(0, true); 
  var ISA_CA = CA_neighborhood(3, ISA_fish).gt(0.3);
  var ISA_KDE_CA = ISA_CA.gt(0).add(ISA_KDE_fish.gt(80)).gte(1);
  var exportname = 'ISA_KDE_CA_' + id;
  Export.image.toAsset({
    image: ISA_KDE_CA.updateMask(ISA_KDE_CA), 
    description: exportname, 
    assetId: 'Your stored path' + exportname, 
    region: roiBound, 
    scale: 30,
    maxPixels:1e13
  });
}

//----------morphological processing---------//
//Load ISA_KDE_CA dataset
for (var i=0; i<<=fishList.size(); i++){
  var ISA_KDE_CA = ee.Image('Your stored ISA_KDE_CA path'+i);
  if (i===0){
    var ISA_KDE_CA_imgs = ee.ImageCollection([ISA_KDE_CA]);
  }
  else{
    ISA_KDE_CA_imgs = ISA_KDE_CA_imgs.merge(ISA_KDE_CA);
  }
}
var ISA_KDE_CA_mosaic = ISA_KDE_CA_imgs.mosaic();
//Batch output
for (var id=0; id<=fishList.size(); id++){
  var selFish =  ee.Feature(fishList.get(id));
  var roiBound = selFish.geometry().buffer(1000).bounds();
  var ISA_KDE_CA_fish = ISA_KDE_CA_mosaic.clip(roiBound).unmask(0, true);
  var ISA_KDE_CA_MP = erode(ISA_KDE_CA_fish, 100).clip(roiBound).unmask(0,true);
  ISA_KDE_CA_MP = dilate(ISA_KDE_CA_MP, 100).clip(roiBound).unmask(0,true);
  var exportname = 'ISA_KDE_CA_MP' + id;
  Export.image.toAsset({
    image: ISA_KDE_CA_MP.updateMask(ISA_KDE_CA_MP), 
    description: exportname, 
    assetId: 'Your stored path' + exportname, 
    region: roiBound, 
    scale: 30,
    maxPixels:1e13
  });
}
  
