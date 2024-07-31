# EC-GPR-TCFs
Workflow for upscaling terrestrial carbon fluxes (TCFs) from a GPR model trained with LAI and climate variables at eddy-covariance towers.

Author: Pablo Reyes-Muñoz

Code: Pablo Reyes-Muñoz

Workflow for upscaling TCFs from MCD15A3H and ERA-5 through Google Earth Engine, from the paper "Tower-to-global upscaling of terrestrial carbon fluxes driven by Copernicus ERA-5 and MODIS-LAI data".

<ol style='list-style-type:disc'> 

<li> Please, download the source code available in this site (TCFs_upscaling / main). </li>

 </br>

The code is formed by a set of Python functions

 </br>

<li> The compose_image function format the predictors input for the GPR function below </li>

 </br>

<li> Calculate_Green function is the core of the GPR algorithm implemented in GEE </li>

</br>
<li></li> The map_loop function iterates over the defined temporal windows </li>

</br>

<p style="text-align:center;"> <img src="https://github.com/psreyes/EC-GPR-TCFs/blob/main/TCFs_Global.png"></p> 

For <b> training and exporting customized models from ARTMO to GEE </b>, please follow the guidelines in https://github.com/msalinero/ARTMOtoGEE.git

A workflow in Python to produce <b> time series mapping </b> over a region of interest can be found in <a href="https://colab.research.google.com/github/daviddkovacs/Global-EVT-maps/blob/main/Main%20Python%20script.ipynb"> this link</a>

