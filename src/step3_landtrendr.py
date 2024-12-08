import ee
import pandas as pd
import geemap
from datetime import date
import os 

# Initialize Earth Engine
ee.Initialize()

# Read in fire-insect co-occurrence data
occurrence = pd.read_csv('/home/goldma34/fire_insect_co-occurence/data/outputs/qc/qc_fire_insect_co-occurrence_v1.csv')
fire_names = occurrence['Fire_ID'].unique()

# Read in fires from GEE
fires = ee.FeatureCollection("users/jandrewgoldman/qc-fire-perims-shield-2").filter(ee.Filter.inList('Fire_ID', fire_names))
# Get required columns and add them together
fires = fires.select('Fire_ID', 'Fire_Year')


########### functions
def getYearStr(year):
  return(ee.String('yr_').cat(ee.Algorithms.String(year).slice(0,4)))

def getYearNumber(y):
  y = y.format('YYYY')
  y = y.getInfo()
  y = ee.Number.parse(y)
  return(y)

def getYear(ls, number):
  y = yoiStr[number]
  return(y)

def calcBS(img, ft1):
  ft1 = ee.Feature(ft1)
  dnbr = img.expression( "(b('preNBR') - b('postNBR')) * 1000").rename('dnbr').toFloat()
  ring  = ft1.buffer(180).difference(ft1);
  offset = ee.Image.constant(ee.Number(dnbr.select('dnbr').reduceRegion(**{'reducer': ee.Reducer.mean(),'geometry': ring.geometry(),'scale': 30,'maxPixels': 1e13}).get('dnbr')))
  offset = offset.rename('offset').toFloat()
  dnbr = dnbr.addBands(offset)
  dnbr = dnbr.addBands(img)
  dnbr_offset = dnbr.expression("b('dnbr') - b('offset')") \
          .rename('dnbr_w_offset').toFloat()
  dnbr_offset = dnbr_offset.addBands(img).addBands(dnbr.select('dnbr'))
  rbr = dnbr_offset.expression("b('dnbr') / (b('preNBR') + 1.001)").rename('rbr').toFloat().addBands(dnbr_offset) 
  rbr_offset = rbr.expression("b('dnbr_w_offset') / (b('preNBR') + 1.001)").rename('rbr_w_offset').toFloat().addBands(rbr)
  rbr_offset = rbr_offset.set('fireID' , ft1.get('Fire_ID'),'fireName' , ft1.get('Fire_Name'), 'fireYear' ,  ft1.get('Fire_Year')) 
  return(rbr_offset)



######### Get Fires
fids = fires.aggregate_array('Fire_ID').getInfo()



# Define function to run LandTrendr for a list of fire IDs
def run_landtrendr_for_fires(fire_ids, fires, folder="chapter3_lt_nbr_v2", scale=30):
    """
    Run LandTrendr for a list of fire IDs and export the results to Google Drive.

    Parameters:
    fire_ids (list): List of fire IDs to process.
    fires (ee.FeatureCollection): FeatureCollection of fires.
    start_year (int): Start year for LandTrendr analysis.
    end_year (int): End year for LandTrendr analysis.
    folder (str): Google Drive folder to save the exported images. Default is "chapter3_lt_nbr_v2".
    scale (int): Scale in meters for the exported images. Default is 30.
    """
    for fire_id in fire_ids:
        # Filter the fire by ID
        fire = fires.filterMetadata('Fire_ID', 'equals', fire_id).first()
        fire_feature = ee.Feature(fire)
        fire_name = fire_feature.get("Fire_ID").getInfo()
        fire_year = fire_feature.get("Fire_Year").getInfo()

         # turn fire year into a string and get prefire year
        fire_year = str(fire_year)
        prefire_year = str(int(fire_year) - 1)
        postfire_year = str(int(fire_year) + 11)

        # Define LandTrendr parameters
        lt_params = {
            "start_date": date(int(prefire_year), 4,1),
            "end_date": date(int(postfire_year), 11,1),
            "index": 'NBR',
            "ftv_list": ['NBR'],
            "mask_labels": ['cloud', 'shadow', 'snow', 'water'],
            "area_of_interest": fire_feature.geometry(),
            "run_params": {
                "maxSegments": 6,
                "spikeThreshold": 0.9,
                "vertexCountOvershoot": 3,
                "preventOneYearRecovery": True,
                "recoveryThreshold": 0.25,
                "pvalThreshold": 0.05,
                "bestModelProportion": 0.75,
                "minObservationsNeeded": 2,
            },
        }

        # Initialize LandTrendr
        lt = geemap.landtrendr.LandTrendr(**lt_params)

        # Run LandTrendr
        lt_result = lt.run()

        ftv_nbr = lt_result.data.select('ftv_nbr_fit').clip(lt_params['area_of_interest'])

        # creat list of years to get imagery
        years = ee.List.sequence(1984, 2021)
        
        # make it a string
        yearsStr = years.map(getYearStr)
        
        ### get RBR image
        # create sequence for each fire to filter image to years we want
        start = ee.Date(fire_year.advance(-1, 'year'))
        start = getYearNumber(start)
        end = ee.Date(fire_year.advance(11, 'year'))
        end = getYearNumber(end)
        yoi = ee.List.sequence(start, end)
        yoiStr = yoi.map(getYearStr)
        yoiStr = yoiStr.getInfo()
    
        
        #nbr
        ftv_nbr2 = ftv_nbr.arrayFlatten([yearsStr])
        nbr = ftv_nbr2.select(yoiStr)

        ## get nbr
        preNBR = nbr.select(getYear(yoiStr, 0)).rename('preNBR').toFloat()
        yof = nbr.select(getYear(yoiStr, 1)).rename('yof')
        postNBR = nbr.select(getYear(yoiStr, 2)).rename('postNBR').toFloat()
        nbr2 = nbr.select(getYear(yoiStr, 3)).rename('nbr_2yr')
        nbr3 = nbr.select(getYear(yoiStr, 4)).rename('nbr_3yr')
        nbr4 = nbr.select(getYear(yoiStr, 5)).rename('nbr_4yr')
        nbr5 = nbr.select(getYear(yoiStr, 6)).rename('nbr_5yr')
        nbr6 = nbr.select(getYear(yoiStr, 7)).rename('nbr_6yr')
        nbr7 = nbr.select(getYear(yoiStr, 8)).rename('nbr_7yr')
        nbr8 = nbr.select(getYear(yoiStr, 9)).rename('nbr_8yr')
        nbr9 = nbr.select(getYear(yoiStr, 10)).rename('nbr_9yr')
        nbr10 = nbr.select(getYear(yoiStr,11)).rename('nbr_10yr')

        
        # add prefire and post fire
        indices = preNBR.addBands(postNBR)
        
        # calculate burn severity metrics
        burnIndices = calcBS(indices, fire_feature)
        
        # add nbr bands 
        nbrBands = nbr2.addBands(nbr3)
        nbrBands = nbrBands.addBands(nbr4)
        nbrBands = nbrBands.addBands(nbr5)
        nbrBands = nbrBands.addBands(nbr6)
        nbrBands = nbrBands.addBands(nbr7)
        nbrBands = nbrBands.addBands(nbr8)
        nbrBands = nbrBands.addBands(nbr9)
        nbrBands = nbrBands.addBands(nbr10)
        nbrBands = nbrBands.addBands(yof)
        
        # add nbr bands to burn severity
        #burnIndices = burnIndices.addBands(nbrBands)
        
        burnIndices = burnIndices.select('preNBR', 'postNBR', 'rbr', 'rbr_offset')


        # Export the result to Google Drive
        export_nbr_name = f'clipped_image_{fire_name}'
        geemap.ee_export_image_to_drive(
            image=lt_result,
            description=export_nbr_name,
            folder=folder,
            region=fire_feature.geometry().bounds().getInfo()['coordinates'],
            scale=scale
        )

# Example usage
run_landtrendr_for_fires(fids, fires)