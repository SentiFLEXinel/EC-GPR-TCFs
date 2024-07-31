import ee
import json
import datetime
import importlib
import sys
import requests
import zipfile
import os
from os import listdir
from os.path import isfile,join
import numpy as np
import netCDF4
from netCDF4 import Dataset
from Models import RECO_model_14VarsNoPANoLCTGFGUY as model_RECO
import ee
PROJECTID = 'ee-pablosrmz'

#Initializing the GEE python api
ee.Initialize()

## To be defined by the user
StartDate = ee.Date('2004-01-01')

EndDate = ee.Date('2024-01-01')

timeStep = 8
 
Iterations = round(EndDate.difference(StartDate, 'day').getInfo()/timeStep)

assetPath='projects/ee-user/assets/collection/'

fc = ee.FeatureCollection("projects/ee-user/assets/AreaOfInterest")

#region = fc.geometry()

OutputFileName = 'GPP_AOI'

variables_GREEN = {'GPP':['GPP', limitUpQualityFlagF, limitDownQualityFlag]}

##Until this line custom settings

#Defining ancilliary functions.
def sequence_GREEN(variable):
	sequence_GREEN = []
	model = globals()['model_' + variable]
	for i in range(0, model.XTrain_dim_GREEN):
		sequence_GREEN.append(str(i))
	return sequence_GREEN

def getInputDates(i):
	fecha_inicio = StartDate.advance(ee.Number(i*timeStep), 'day')
	fecha_fin = fecha_inicio.advance(timeStep, 'day')
	fecha_str = datetime.datetime.utcfromtimestamp(fecha_inicio.getInfo()['value']/1000.0).strftime('%Y%m%d')
	return {'fecha_inicio':fecha_inicio, 'fecha_fin':fecha_fin, 'fecha_str':fecha_str}

def addTimeProp(image):
	image_timed = image.set('system:time_start', (image.getString('system:index').slice(-8)))
	return image_timed
	
def maskS3badPixels(image):
  qa = ee.Image(image.select('quality_flags'));
  coastLine = 1 << 30;
  inLandWater = 1 << 29;
  bright = 1 << 27;
  invalid = 1 << 25;
  Oa12Sat = 1 << 9;
  mask = qa.bitwiseAnd(coastLine).eq(0).And(qa.bitwiseAnd(inLandWater).eq(0)).And(qa.bitwiseAnd(bright).eq(0))
  
  return image.updateMask(mask);

def addVariables(image):
  date = ee.Date(image.get("system:time_start"));
  years = date.difference(ee.Date('1970-01-01'),'days');
  return image.addBands(ee.Image(years).rename('t').float());

def dewPointVaporPressure(Tdew):
  return ee.Image(0.6108).multiply(ee.Image(17.27).multiply(Tdew.subtract(ee.Image(273.15))).divide((Tdew).add(ee.Image(-35.85))).exp());

#We define EC-GPR-TCFs input variables and compose the input image in GEE for estimating TCFs through GPR

def compose_image(fecha_inicio, fecha_fin, fecha_inicio_str, n, S3TOAGPR):
  
  LAIMOD = ee.Image(ee.ImageCollection("MODIS/061/MCD15A3H")
	.filterDate(fecha_inicio, fecha_fin)
	.select('Lai')
	.median()
	.divide(10)
	)

  ERA = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")

  PAR = ERA.filterDate(fecha_inicio, fecha_fin).select('surface_solar_radiation_downwards_hourly').mean().divide(3600)

  TMean = ERA.filterDate(fecha_inicio, fecha_fin).select('temperature_2m').mean().multiply(1.8).multiply(0.7)

  TMeanC = ERA.filterDate(getInputDates(0)['fecha_inicio'], getInputDates(0)['fecha_fin']).select('temperature_2m').mean().subtract(273.3).multiply(0.7)

  dewPointT = ERA5.filterDate(fecha_inicio, fecha_fin).select('dewpoint_temperature_2m').mean()

  TS1 = ERA.filterDate(fecha_inicio, fecha_fin).select('soil_temperature_level_1').mean().subtract(273.3)
  TS2 = ERA.filterDate(fecha_inicio, fecha_fin).select('soil_temperature_level_2').mean().subtract(273.3)
  TS3 = ERA.filterDate(fecha_inicio, fecha_fin).select('soil_temperature_level_3').mean().subtract(273.3)

  SWC1 = ERA.filterDate(fecha_inicio, fecha_fin).select('volumetric_soil_water_layer_1').mean().multiply(100)
  SWC2 = ERA.filterDate(fecha_inicio, fecha_fin).select('volumetric_soil_water_layer_2').mean().multiply(100)
  SWC3 = ERA.filterDate(fecha_inicio, fecha_fin).select('volumetric_soil_water_layer_3').mean().multiply(100)

  LE = ERA.filterDate(fecha_inicio, fecha_fin).select('surface_latent_heat_flux').mean().divide(86400).multiply(-1)
  
  HE = ERA.filterDate(fecha_inicio, fecha_fin).select('surface_sensible_heat_flux').mean().divide(86400).multiply(-1)

  P = ERA.filterDate(fecha_inicio, fecha_fin).select('total_precipitation').mean()

  vpair = dewPointVaporPressure(dewPointT)
        
  vpsat = dewPointVaporPressure(TMean)
        
  VPD = vpair.subtract(vpsat)

  WSu = ERA.filterDate(fecha_inicio, fecha_fin).select('u_component_of_wind_10m').first()#.resample('bilinear').reproject(crs=LAIMOD.select('Lai').projection().crs(), scale=500)
  WSv = ERA.filterDate(fecha_inicio, fecha_fin).select('v_component_of_wind_10m').first()#.resample('bilinear').reproject(crs=LAIMOD.select('Lai').projection().crs(), scale=5000)

  WS = ((WSu.multiply(WSu)).add(WSv.multiply(WSv))).sqrt()

  return ee.Image(LAIMOD.addBands(PAR).addBands(VPD).addBands(TS1).addBands(TS2).addBands(TS3).addBands(SWC1).addBands(SWC2).addBands(SWC3).addBands(HE).addBands(LE).addBands(P).addBands(TMeanC).addBands(WS))

#We apply the GPR function

def calculate_GREEN(fecha_inicio, fecha_fin, fecha_inicio_str, variable, limitUp, limitDown, n, S3TOAGPR): 
  
  model = globals()['model_' + variable]
  
  image = compose_image(fecha_inicio, fecha_fin, fecha_inicio_str, n, S3TOAGPR).clipToCollection(config.fc);

  im_norm_ell2D_hypell = image.subtract(model.mx_GREEN).divide(model.sx_GREEN).multiply(model.hyp_ell_GREEN).toArray().toArray(1); 
  im_norm_ell2D = image.subtract(model.mx_GREEN).divide(model.sx_GREEN).toArray().toArray(1); 
  PtTPt  = im_norm_ell2D_hypell.matrixTranspose().matrixMultiply(im_norm_ell2D).arrayProject([0]).multiply(-0.5); #OK

  PtTDX  = ee.Image(model.X_train_GREEN).matrixMultiply(im_norm_ell2D_hypell).arrayProject([0]).arrayFlatten([sequence_GREEN(variable)]);
  arg1   = PtTPt.exp().multiply(model.hyp_sig0_GREEN);
  k_star = PtTDX.subtract(model.XDX_pre_calc_GREEN.multiply(0.5)).exp().toArray();
  mean_pred = k_star.arrayDotProduct(model.alpha_coefficients_GREEN.toArray()).multiply(arg1);
  mean_pred = mean_pred.toArray(1).arrayProject([0]).arrayFlatten([[variable + '_GREEN']]);
  mean_pred = mean_pred.add(model.mean_model_GREEN);#.updateMask(mean_pred.gt(0));
  filterDown = mean_pred.gt(limitDown)
 
  filterUp = mean_pred.lt(limitUp)
  
  quality_flags = (filterDown.multiply(filterUp)).Not().toArray().arrayFlatten([[variable + '_QUALITY_FLAG']])

  k_star_uncert = PtTDX.subtract(model.XDX_pre_calc_GREEN.multiply(0.5)).exp().multiply(arg1).toArray();
  Vvector = ee.Image(model.LMatrixInverse).matrixMultiply(k_star_uncert.toArray(0).toArray(1)).arrayProject([0])
  Variance = ee.Image(model.hyp_sig_GREEN).toArray().subtract(Vvector.arrayDotProduct(Vvector)).sqrt()
  Variance = Variance.toArray(1).arrayProject([0]).arrayFlatten([[variable + '_UNCERTAINTY_GREEN']])

  image= image.addBands(mean_pred);
  image = image.addBands(Variance);
  image = image.addBands(quality_flags);

  return image.select(variable + '_GREEN', variable + '_UNCERTAINTY_GREEN', variable + '_QUALITY_FLAG');

# We iterate to create a collection over a temporal range of interest

def maploop():

	for i in range(0,Iterations):
    
		imageHolder = ee.Image();

		for variable_GREEN in variables_GREEN:
			params = variables_GREEN[variable_GREEN]
			variable = params[0]
			limitUp = params[1]
			limitDown = params[2]
			imagen = calculate_GREEN(getInputDates(i)['fecha_inicio'], getInputDates(i)['fecha_fin'], getInputDates(i)['fecha_str'], variable, limitUp, limitDown, i, S3TOAGPR)
			imageHolder = imageHolder.addBands(imagen) 

		image_export = imageHolder.select(variable + '_GREEN', variable + '_UNCERTAINTY_GREEN')
		exportar = ee.batch.Export.image.toAsset(
		assetId=assetPath + OutputFileName + getInputDates(i)['fecha_str'],
		image=image_export,
		maxPixels=1e10,
		crs='EPSG:4326',
		description=OutputFileName + getInputDates(i)['fecha_str'],
		scale=5000,
		region=fc.geometry()
		)
		exportar.start()
		exportar.status()  
  

