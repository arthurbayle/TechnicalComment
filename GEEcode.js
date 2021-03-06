/* 
 * Code for the manuscript "From white to green: Multidecadal trends of
 * snow cover and vegetation productivity in the European Alps" by
 * Rumpf et al., submitted December 2021 to Science
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

// simon.gascoin@cesbio.cnes.fr
// files added in shared repo 
var projectPath='users/sgascoin/greening/' 

var L4 = ee.ImageCollection("LANDSAT/LT04/C01/T1_SR"),
    L5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR"),
    L7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR"),
    geometry = 
    /* color: #d63000 */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[8.086392982858342, 46.60915355637467],
          [8.086392982858342, 46.536757352783646],
          [8.270671471017522, 46.536757352783646],
          [8.270671471017522, 46.60915355637467]]], null, false),
    treeCover = ee.ImageCollection("GLCF/GLS_TCC"),
    srtm = ee.Image("USGS/SRTMGL1_003"),
    L8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR"),
    glacier = ee.FeatureCollection("GLIMS/current"),
    alpsArea = ee.FeatureCollection(projectPath+"Alpine_Convention_Perimeter_2018_v2"),
    iceObservation = ee.Image(projectPath+"IceObservation");

var startSection=2018; // start

var stepSize=36; // full export

var firstNoIceObservation=ee.Image.constant(0).where(iceObservation.select('fNIO').mask(),iceObservation.select('fNIO'));

var quantile=75; // NDVI quantile
var thresholdNDSI=0.4; // clear enough :)


//start all Landsat collections
// shift teh band names of L8 to match the other once
var mergedCollection=L4.merge(L5).merge(L7).merge(L8.map(function(image){return image.select(['B2','B3','B4','B5','B6','B7','B10','pixel_qa'])
      .rename(['B1','B2','B3','B4','B5','B7','B6','pixel_qa'])})).select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7','pixel_qa']).filterDate(startSection+"-01-01",(startSection+stepSize+1)+"01-01");

//remove pixel with negative reflectance
mergedCollection=mergedCollection.map(function(image){return image.updateMask(image.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']).gt(0).toArray().arrayReduce(ee.Reducer.min(),[0]).arrayGet(0));});

// remove cloud, shadow ...
mergedCollection=mergedCollection.map(function(image){return image.updateMask(image.select('pixel_qa').bitwiseAnd(32+8+4).not())});

//remove glacier
mergedCollection=mergedCollection.map(function(image){return image.updateMask(firstNoIceObservation.lt(image.date().get('year')));});

//compute and add NDVI and NDSI
mergedCollection=mergedCollection.map(function(image){
  return image.addBands(image.normalizedDifference(['B4','B3']).rename('NDVI'))
    .addBands(image.normalizedDifference(['B2','B5']).rename('NDSI'));
});

mergedCollection=mergedCollection.select(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'NDVI', 'NDSI']);

// remove genative NDVI, and threshold NDSI
mergedCollection=mergedCollection.map(function(image){
  return image.addBands(image.select('NDVI').max(0).rename('corrcetedNDVI'))
              .addBands(image.select('NDSI').gt(thresholdNDSI).rename('Snow'));
});

// filter to the area of interest
mergedCollection=mergedCollection.filterBounds(alpsArea);
var forestCover10=treeCover.mosaic().select(0).lt(10);

//remove forst and threshold altitude 
mergedCollection=mergedCollection.map(function(image){
  return image.updateMask(forestCover10).updateMask(srtm.gt(1700));
});

//filter from june to september
mergedCollection=mergedCollection.filter(ee.Filter.calendarRange(6, 9, 'month'));


// apply sesnor corrections to harmonize to L7
mergedCollection=mergedCollection.map(function(im){
  var imL8=im;
  imL8=imL8.addBands(imL8.select('NDVI').multiply(0.9166515729479716).add(0.014362264022815808),null,true);
  imL8=imL8.addBands(imL8.select('NDSI').multiply(0.9284655117221521).add(-0.010708073607665443),null,true);
  return ee.Algorithms.If({
    condition:ee.String(im.get('SATELLITE')).equals('LANDSAT_8'),
    trueCase:imL8,
    falseCase:im,
  }) 
}).map(function(im){
  var imL5=im;
  imL5=imL5.addBands(imL5.select('NDVI').multiply(0.9848094305727885).add(0.017977388487747152),null,true);
  imL5=imL5.addBands(imL5.select('NDSI').multiply(0.8652263528673204).add(-0.07197670457494443),null,true);
  return ee.Algorithms.If({
    condition:ee.String(im.get('SATELLITE')).equals('LANDSAT_5'),
    trueCase:imL5,
    falseCase:im,
  })
}).map(function(im){
  var imL5=im;
  imL5=imL5.addBands(imL5.select('NDVI').multiply(0.9848094305727885).add(0.017977388487747152),null,true);
  imL5=imL5.addBands(imL5.select('NDSI').multiply(0.8652263528673204).add(-0.07197670457494443),null,true);
  return ee.Algorithms.If({
    condition:ee.String(im.get('SATELLITE')).equals('LANDSAT_4'),
    trueCase:imL5,
    falseCase:im,
  })
}).map(function(im){return im.toFloat();});


// Create a time filter to define a match as overlapping timestamps.
var timeFilter = ee.Filter.maxDifference({
    difference: 1000 * 3600 * 24 * 2 * 31,  //2 months each side
    leftField: 'system:time_start',
    rightField: 'system:time_start'
  });

// Define the join.
var saveAllJoin = ee.Join.saveAll({
  matchesKey: 'ImageOfTheYear',
  ordering: 'system:time_start',
  ascending: true
});

// create year event
var yearsList=ee.FeatureCollection(ee.List.sequence(startSection, +startSection+stepSize-1, 1)
  .map(function(y){return ee.Feature(null,
      ee.Dictionary.fromLists(['system:time_start'],[ee.Date.fromYMD(y,8,1).millis()]))}));

// Apply the join.
var YearlyMergedCollection = saveAllJoin.apply(yearsList, mergedCollection, timeFilter);

// aply computation fro each yearly set, proportion of snow, NDVI quantile ...
YearlyMergedCollection=ee.ImageCollection(YearlyMergedCollection.map(function(yearFeature){
  var imCol=ee.ImageCollection(ee.List(yearFeature.get('ImageOfTheYear')));
  var im=imCol.reduce(ee.Reducer.percentile([quantile])).rename(imCol.first().bandNames()).set('system:time_start', yearFeature.get('system:time_start'));
  im=im.addBands(imCol.reduce(ee.Reducer.count()).select('NDVI_count').rename('imageCount'))
  im=im.addBands(imCol.select('Snow').mean().rename('SnowProb'))
  im=im.addBands(im.select('SnowProb').gt(0.999).rename('PermanantSnow'));
  return im;
}));


//set the band to export
var bandsToExport=["corrcetedNDVI","imageCount","SnowProb","PermanantSnow"];

YearlyMergedCollection=YearlyMergedCollection.select(bandsToExport);

//full export
/*
YearlyMergedCollection.size().evaluate(function(size){
  for(var year=0; year<size; year++){
    var im=ee.Image(YearlyMergedCollection.toList(1000).get(year));
    var y=ee.Date(ee.Number(im.get('system:time_start'))).get('year').getInfo()
      for( var bandNameIndex=0; bandNameIndex< bandsToExport.length; bandNameIndex++){
        Export.image.toCloudStorage({
            image: im.toFloat().select(bandsToExport[bandNameIndex]),
            description: 'alps_'+bandsToExport[bandNameIndex]+'_'+y,
            bucket: 'alps-image-ee/'+bandsToExport[bandNameIndex],
            fileNamePrefix: 'alps_'+bandsToExport[bandNameIndex]+'_'+y,
            scale: 30,
            maxPixels:1e10,
            region: alpsArea.geometry().bounds(),
            fileDimensions:256*256,
          });
      }
  }
});

// for R
var im=YearlyMergedCollection.first();
print(YearlyMergedCollection)
var y=ee.Date(ee.Number(im.get('system:time_start'))).get('year').getInfo()
  for( var bandNameIndex=0; bandNameIndex< bandsToExport.length; bandNameIndex++){
    Export.image.toCloudStorage({
        image: im.toFloat().select(bandsToExport[bandNameIndex]),
        description: 'alps_'+bandsToExport[bandNameIndex]+'_'+y,
        bucket: 'alps-image-ee',
        fileNamePrefix: bandsToExport[bandNameIndex]+'/alps_'+bandsToExport[bandNameIndex]+'_'+y,
        scale: 30,
        maxPixels:1e10,
        region: alpsArea.geometry().bounds(),
        fileDimensions:256*256,
      });
  }
*/

// simon.gascoin@cesbio.cnes.fr
// export loop  
for(var y=2020; y<2022; y++){ 
var exportYear = ee.String(ee.Number(y))
var im = YearlyMergedCollection
          .filterDate(exportYear.cat('-08-01'), exportYear.cat('-08-02'))
          .first()
          print(im)
Export.image.toDrive({
        image: im.select('imageCount'),
        description: 'alps_imageCount_'+exportYear.getInfo(),
        folder:'GEE/greening/',
        scale: 30,
        maxPixels:1e10,
        region: alpsArea.geometry().bounds()
      });
}
