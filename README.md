##
Alexa Villaume

####PURPOSE:

#### FILE RUNDOWN:
* **createModelGrid.py** -
* **userpar.inc/lambda_grid.dat** - Internal DUSTY files to control the wavelength grid. Used the BaSeL wavelength grid as input for DUSTY. 
* **CGrid.txt/OGrid.txt** - The grids of circumstellar dust models to be included in FSPS. The first line is the number of wavelength points in the
  grid. The second line are the wavelengths in microns. After that the lines alternate between the temperature (K) and tau value (in log) of the model
  on the following line.
