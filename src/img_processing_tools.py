import ee


def check_image_info(img):
  imgInf = img.getInfo()
  print(imgInf)
  
  # Returns vegetation indices for LS8 / change to CO2_T2_L2
def ls8_Indices(lsImage):
  nbr = lsImage.normalizedDifference(['SR_B5', 'SR_B7']).toFloat()
  qa = lsImage.select(['QA_PIXEL'])
  return nbr.addBands([qa]).select([0,1], ['nbr', 'QA_PIXEL']).copyProperties(lsImage, ['system:time_start'])

  
#// Returns vegetation indices for LS4, LS5 and LS7 . change to CO2_T2_L2
def ls4_7_Indices(lsImage):
  nbr = lsImage.normalizedDifference(['SR_B4', 'SR_B7']).toFloat()
  qa = lsImage.select(['QA_PIXEL'])
  return nbr.addBands([qa]).select([0,1], ['nbr', 'QA_PIXEL']).copyProperties(lsImage, ['system:time_start'])
  
 
# Mask Landsat surface reflectance images
# Creates a mask for clear pixels 
def lsCfmask(lsImg):
  qa = lsImg.select("QA_PIXEL")
  cloudBitMask = 1 << 3
  cirrusBitMask = 1 << 4
  waterBitMask = 1 << 7
  snowBitMask = 1 << 5
  mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0)).And(
    qa.bitwiseAnd(waterBitMask).eq(0)).And(
    qa.bitwiseAnd(snowBitMask).eq(0))
  return lsImg.updateMask(mask) \
      .select("nbr") \
      .copyProperties(lsImg, ["system:time_start"])
   
   
#  check to see if there are pixels: this prints to console
def check_pixel_count(img, geom):
  res = img.reduceRegion(**{
  'reducer': ee.Reducer.count(),\
  'geometry':geom,\
  'scale': 30})
  print(res.getInfo())

# check value
def check_pixel_mean(img, geom):
  res = img.reduceRegion(**{
  'reducer': ee.Reducer.mean(),\
  'geometry':geom,\
  'scale': 30})
  print(res.getInfo())

# check to see if there are pixels: this returns a number object
def return_pixel_count(img, geom, time):
  res = img.reduceRegion(**{
  'reducer': ee.Reducer.count(),\
  'geometry':geom,\
  'scale': 30})
  if time == "preNBR":
    res = res.getNumber('preNBR')
  elif time == "postNBR":
    res = res.getNumber('postNBR')
  return(res.getInfo())

# check dates - converts ee date to y-m-d and prints in console
def check_dates(date):
  date = date.format('Y-M-d')
  print(date.getInfo())


# check if image contains key
def get_keys(img, geom):
  res = img.reduceRegion(**{
  'reducer': ee.Reducer.count(),\
  'geometry':geom,\
  'scale': 30})
  res = res.keys()
  return(res.getInfo())


def water_mask(img):
  gsw = ee.Image("JRC/GSW1_1/GlobalSurfaceWater")
  water = gsw.select('max_extent')
  mask = water.eq(0)
  maskedImg = img.updateMask(mask)
  return(maskedImg)

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

# for each img calculate nbr
def nbrPerYear(img):
  nbr = img.normalizedDifference(['nir', 'swir']).rename('nbr').toFloat()
  return nbr.copyProperties(img,['system:time_start'])


# def function, for each year, create new img for that year wth two bands
def mergeByYear(img1, img2, key):
  yr = getYear(yoiStr, key)
  nir2 = img1.select(yr).rename('nir')
  swir2 = img2.select(yr).rename('swir')
  merged = nir2.addBands(swir2)
  nbr= nbrPerYear(merged)
  if (key == 0):
    img = merged.addBands(nbr).select('nbr').rename('preNBR')
  elif (key == 2):
    img = merged.addBands(nbr).select('nbr').rename('postNBR')
  elif (key == 3):
    img = merged.addBands(nbr).select('nbr').rename('nbr_2yr')
  elif (key == 4):
    img = merged.addBands(nbr).select('nbr').rename('nbr_3yr')
  elif (key == 5):
    img = merged.addBands(nbr).select('nbr').rename('nbr_4yr')
  elif (key == 6):
    img = merged.addBands(nbr).select('nbr').rename('nbr_5yr')
  elif (key == 7):
    img = merged.addBands(nbr).select('nbr').rename('nbr_6yr')
  elif (key == 8):
    img = merged.addBands(nbr).select('nbr').rename('nbr_7yr')
  elif (key == 9):
    img = merged.addBands(nbr).select('nbr').rename('nbr_8yr')
  elif (key == 10):
    img = merged.addBands(nbr).select('nbr').rename('nbr_9yr')
  elif (key == 11):
    img = merged.addBands(nbr).select('nbr').rename('nbr_10yr')
  return(img)


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


# for each image in image collection map through the collection, if the year of the image is > 