import ee
import pandas as pd
import geemap
from datetime import date
import os 
from ltgee import LandTrendr, LandsatComposite, LtCollection
from datetime import date

# Initialize Earth Engine
ee.Initialize()

# Define Helpter Functions

def get_year_range(prefire_year, postfire_year):
    """
    Create a list of year strings in the format 'yrYYYY' for a given range of years.
    prefire_year: str - The prefire year.
    postfire_year: str - The postfire year.
    """
    years = []

    # Iterate over the range of years and append the formatted strings to the list
    for i in range(int(prefire_year), int(postfire_year) + 1):
        years.append('yr' + str(i))
    
    return years
    

def rename_bands(image):
    """
    Renames each band in the image by stripping "yr" from the name.

    Parameters:
    image (ee.Image): The input image with bands to be renamed.

    Returns:
    ee.Image: The image with renamed bands.
    """
    # Get the list of band names
    band_names = image.bandNames()

    # Create a list of new band names by stripping "yr" from each name
    new_band_names = band_names.map(lambda name: ee.String(name).replace('yr', '', 'g'))

    # Rename the bands in the image
    renamed_image = image.rename(new_band_names)

    return renamed_image

    


def calcBS(img, ft1):
    """
    Calculate burn severity indices for a given image and feature.
    img: ee.Image - Image containing NBR bands.
    ft1: ee.Feature - Feature containing the fire perimeter.
    """
    ft1 = ee.Feature(ft1)
    
    # Select the bands by their indices
    preNBR = img.select([0]).rename('preNBR')
    postNBR = img.select([2]).rename('postNBR')
    
    # Calculate dNBR
    dnbr = preNBR.subtract(postNBR).multiply(1000).rename('dnbr').toFloat()
    
    # Create a ring buffer around the feature
    ring = ft1.buffer(180).difference(ft1)
    
    # Calculate the offset
    offset = ee.Image.constant(ee.Number(dnbr.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=ring.geometry(),
        scale=30,
        maxPixels=1e13
    ).get('dnbr'))).rename('offset').toFloat()
    
    # Add the offset and original bands to the dNBR image
    dnbr = dnbr.addBands(offset).addBands(preNBR).addBands(postNBR)
    
    # Calculate dNBR with offset
    dnbr_offset = dnbr.expression("b('dnbr') - b('offset')").rename('dnbr_w_offset').toFloat()
    dnbr_offset = dnbr_offset.addBands(dnbr)
    
    # Calculate RBR
    rbr = dnbr_offset.expression("b('dnbr') / (b('preNBR') + 1.001)").rename('rbr').toFloat().addBands(dnbr_offset)
    
    # Calculate RBR with offset
    rbr_offset = rbr.expression("b('dnbr_w_offset') / (b('preNBR') + 1.001)").rename('rbr_w_offset').toFloat().addBands(rbr)
    
    # Set properties
    rbr_offset = rbr_offset.set('fireID', ft1.get('Fire_ID'), 'fireName', ft1.get('Fire_Name'), 'fireYear', ft1.get('Fire_Year'))
    
    return rbr_offset


# Define Main Processing Functions

def export_landtrendr_image(feature):
    """
    Export LandTrendr imagery for a given feature.

    input:
    feature: ee.Feature - The feature to process.

    Subsequent variables:
    fire_id: str - The fire ID.
    geometry: ee.Geometry - The geometry of the feature.
    fire_year: str - The year of the fire.
    prefire_year: str - The year before the fire.
    postfire_year: str - The year after the fire.
    composite_params: dict - Dictionary of parameters for LandsatComposite.
    lt_collection_params: dict - Dictionary of parameters for LtCollection.
    lt_params: dict - Dictionary of parameters for LandTrendr.
    lt: LandTrendr - The LandTrendr object.
    nbr_stack: ee.Image - The NBR stack image.
    clean_nbr_stack: ee.Image - The cleaned NBR stack image.
    results: ee.Image - The burn severity indices image.
    export_bi_name: str - The name of the burn indices image to export.
    export_nbr_name: str - The name of the NBR recovery image to export.

    output:
    writes out band to drive

    
    
    """
    # get the feature's ID, geometry and year
    fire_id = feature.get('Fire_ID')
    geometry = feature.geometry().bounds()
    fire_year = feature.get('FireYear')

    # turn fire year into a string and get prefire year
    fire_year = str(fire_year)
    prefire_year = str(int(fire_year) - 1)
    postfire_year = str(int(fire_year) + 11)

    # Example arguments to pass to the LandTrendr class. See docstring of LandTrendr for more information'
    composite_params = {
    "start_date": date(int(prefire_year), 4,1),
    "end_date": date(int(postfire_year), 11,1),
    "area_of_interest": geometry,
    "mask_labels": ['cloud', 'shadow', 'snow', 'water'],
    "debug": True
}
    lt_collection_params = {
        "sr_collection": LandsatComposite(**composite_params),
        # "sr_collection": composite_params, # - you may also just pass in your own collection or the params directly
        "index": 'NBR',
        "ftv_list": ['NBR']}
    lt_params = {
        "lt_collection": LtCollection(**lt_collection_params),
        # "lt_collection": lt_collection_params, # - you may also just pass in your own collection or the params directly
        "run_params": {
                "maxSegments": 6,
                "spikeThreshold": 0.9,
                "vertexCountOvershoot":  3,
                "preventOneYearRecovery":  True,
                "recoveryThreshold":  0.25,
                "pvalThreshold":  .05,
                "bestModelProportion":  0.75,
                "minObservationsNeeded": 2,
            }
    }

    lt = LandTrendr(**lt_params)
    lt.run()
    
    ftv_nbr = lt.data.select('ftv_nbr_fit').clip(lt_params['area_of_interest'])

    # get the year range
    years = get_year_range(prefire_year, postfire_year)

    # convert 1-D fitted nbr array to an image with bands representing years
    nbr_stack = ftv_nbr.arrayFlatten([years])

    # clean the names of the nbr stack bands
    clean_nbr_stack = rename_bands(nbr_stack)

    # calculate burn severity indices
    results = calcBS(clean_nbr_stack, feature)

    # make the file names to export each image
    export_bi_name = fire_id + "_bi"
    export_nbr_name = fire_id + "_nbr"
    
    # export  burn Indices imagery
    geemap.ee_export_image_to_drive(results, description=export_bi_name , folder="chapter3_lt_bi_v2", region=geometry, scale=30)

    # export nbr recovery imagery
    geemap.ee_export_image_to_drive(clean_nbr_stack, description=export_nbr_name, folder="chapter3_lt_nbr_v2", region=geometry, scale=30)
    

# Main Function
def main():
    occurrence = pd.read_csv('/home/goldma34/fire_insect_co-occurence/data/outputs/qc/qc_fire_insect_co-occurrence_v1.csv')
    fire_names = occurrence['Fire_ID'].unique()

    fires = ee.FeatureCollection("users/jandrewgoldman/qc-fire-perims-shield-2").filter(ee.Filter.inList('Fire_ID', fire_names))
    fires = fires.select('Fire_ID', 'Fire_Year')

    fire_ids = fires.aggregate_array('Fire_ID').getInfo()

    for fire_id in fire_ids:
        print(fire_id)
        fire = fires.filterMetadata('Fire_ID', 'equals', fire_id).first()
        results = export_landtrendr_image(fire)
        print(f"Exported LandTrendr imagery for fire {fire_id}")

if __name__ == "__main__":
    main()