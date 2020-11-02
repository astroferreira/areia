# areia
Artificial Redshift Effects for IA

This Python code is losely based on the IDL FERENGI code (https://www2.mpia-hd.mpg.de/FERENGI/).

Currently it has only a single band mode, so bandpass-shifting and k-corrections are yet to be implemmented. Also, there is no PSF reconstruction, as the method used by FERENGI does not work well with realistic psfs (the reconstruction is done in fourier space with a wiener filter). 

Checkout the IPython Notebook for an usage example. A feature to use it from the command line to be implemmented soon.
