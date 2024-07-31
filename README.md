# EC-GPR-TCFs
Workflow for upscaling terrestrial carbon fluxes (TCFs) from a GPR model trained with LAI and climate variables at eddy-covariance towers.

Author: Pablo Reyes-Muñoz

Code: Pablo Reyes-Muñoz

Workflow for upscaling TCFs from MCD15A3H and ERA-5 through Google Earth Engine, from the paper "Tower-to-global upscaling of terrestrial carbon fluxes driven by Copernicus ERA-5 and MODIS-LAI data".

Please, download the source code available in this site (TCFs_upscaling / main).

The code is formed by a set of Python functions

The compose_image function format the predictors input for the GPR function below

Calculate_Green function is the core of the GPR algorithm implemented in GEE

    <\br>
The map_loop function iterates over the defined temporal windows

    <\br>

<p style="text-align:center;"> <img src="https://github.com/psreyes/EC-GPR-TCFs/blob/main/TCFs_Global.png"></p> 

