import geopandas as gpd
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
import ee
import geemap
import os
from ltgee import LandTrendr
from datetime import date
import ltgee

# Intialize Earth Engine
ee.Initialize()


# read in fire-insect co-occurrence data
occurrence = pd.read_csv('/home/goldma34/fire_insect_co-occurence/data/outputs/qc/qc_fire_insect_co-occurrence_v1.csv')
fire_names = occurrence['Fire_ID'].unique()

# Read in fires from GEE
fires = ee.FeatureCollection("users/jandrewgoldman/qc-fire-perims-shield-2").filter(ee.Filter.inList('Fire_ID', fire_names))
# get required columns and add them together
fires = fires.select('Fire_ID', 'Fire_Year')

############### START OF LANDTRENDR AND BURN SEVERITY ANALYSIS ###############

# get Fires
fids = fires.aggregate_array('Fire_ID')
 # turn it to python list
fids = fids.getInfo()

print(fids)
fids1 = str(fids)


#test_fire = fires.filter(ee.Filter.eq('Fire_ID', 'QC_533_1996'))
#ft = test_fire.first()
# loop through list of fire ids

for j in fids:
  Name = j
  fiya= fires.filterMetadata('Fire_ID', 'equals', Name).first();
  ft = ee.Feature(fiya)
  fName = ft.get("Fire_ID")
  fire = ft


  fireBounds = ft.geometry().bounds()
  fyr = ft.get('Fire_Year')
  fyr = fyr.getInfo()
  fyr = str(fyr)
  fireYear = ee.Date(fyr)
  preFireYear = fireYear.advance(-1, 'year')
  
  
  lt_params = {
      "start_date": date(1984, 4,1),
      "end_date": date(2021, 11,1),
      "index": 'NBR',
      "ftv_list": ['TCB', 'TCG', 'TCW', 'TCA', 'NBR'],
      "mask_labels": ['cloud', 'shadow', 'snow', 'water'],
      "area_of_interest": fireBounds 
      ,"run_params": {
              "maxSegments": 6,
              "spikeThreshold": 0.9,
              "vertexCountOvershoot":  3,
              "preventOneYearRecovery":  True,
              "recoveryThreshold":  0.25,
              "pvalThreshold":  .05,
              "bestModelProportion":  0.75,
              "minObservationsNeeded": 2,
          }}

  lt = LandTrendr(debug=True, **lt_params)
  lt.run()
  
  #lt.data.getInfo()
  #lt.sr_collection.getInfo()
 
  #lt.lt_collection.getInfo()
   
 # lt.clear_pixel_count_collection
   
  ftv_nbr = lt.data.select('ftv_nbr_fit').clip(lt_params['area_of_interest'])

  # creat list of years to get imagery
  years = ee.List.sequence(1984, 2021)
  
  # make it a string
  yearsStr = years.map(getYearStr)
  
   ### get RBR image
  # create sequence for each fire to filter image to years we want
  start = ee.Date(fireYear.advance(-1, 'year'))
  start = getYearNumber(start)
  end = ee.Date(fireYear.advance(11, 'year'))
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
  burnIndices = calcBS(indices, ft)
  
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
  
  burnIndices = burnIndices.select('preNBR', 'postNBR', 'rbr', 'rbr_w_offset')
  
  export_bi_name = Name + "_bi"
  export_nbr_name = Name + "_nbr"
  
  # export  burn Indices imagery
  geemap.ee_export_image_to_drive(burnIndices, description=export_bi_name , folder="chapter3_lt_bi_v2", region=fireBounds, scale=30)

  # export nbr recovery imagery
  geemap.ee_export_image_to_drive(nbrBands, description=export_nbr_name, folder="chapter3_lt_nbr_v2", region=fireBounds, scale=30)