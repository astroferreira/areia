import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io import fits
from astropy import units as u
from photutils import CircularAperture
from photutils import aperture_photometry

gds = pd.read_pickle('/data/captain/BGS/gs_complete.pk')

field125 = fits.open('/data/captain/BGS/hlsp_candels_hst_wfc3_gs-tot_f125w_v1.0_drz.fits', 
                  memmap=True)

field160 = fits.open('/data/captain/BGS/hlsp_candels_hst_wfc3_gs-tot_f160w_v1.0_drz.fits', 
                  memmap=True)

wcs160 = WCS(field160[0].header)
wcs125 = WCS(field125[0].header)
center_coord = SkyCoord(ra=53.12272274*u.degree, dec=-27.805064*u.degree)
ra_min = (center_coord.ra - 120*u.arcsec).to(u.degree).value
ra_max = (center_coord.ra + 120*u.arcsec).to(u.degree).value
dec_min = (center_coord.dec - 120*u.arcsec).to(u.degree).value
dec_max = (center_coord.dec + 120*u.arcsec).to(u.degree).value
subset = np.where((gds.DEC >= dec_min) & (gds.DEC <= dec_max) & (gds.RA >= ra_min) & (gds.RA <= ra_max) & (gds.FLAGS == 0))
RAs = gds.RA.values*u.degree
DECs = gds.DEC.values*u.degree

source_coords = SkyCoord(ra=RAs, dec=DECs)


def find_bg_section(bg, img):

    bg = random_flips(bg)

    a, b = bg.shape
    c, d = img.shape

    if((c <= a) & (d <= b)):
        x = np.random.randint(0, a-c)
        y = np.random.randint(0, b-d)
    else:
        #return np.zeros_like(img) 
        raise(IndexError('Galaxy Image larger than BG'))
    
    bg_section = bg[x:(x+c), y:(y+d)]

    return bg_section


def find_bg_section_from_goodss(size, bg_coord=None):

    offset = size

    if(bg_coord is None):
        reiterate = True
        num_iters = 0
        while(reiterate):    

            x_c = np.random.randint(offset, field160[0].data.shape[0]-offset)
            y_c = np.random.randint(offset, field160[0].data.shape[0]-offset)

            x_i = int(x_c-offset/2)
            x_f = int(x_c+offset/2)
            y_i = int(y_c-offset/2)
            y_f = int(y_c+offset/2)

            bg_coord = SkyCoord(ra=wcs160.wcs_pix2world(x_c, y_c, 0)[0]*u.degree, dec=wcs160.wcs_pix2world(x_c, y_c, 0)[1]*u.degree)

            try:
                bg160 = Cutout2D(field160[0].data, position=bg_coord, size=offset*u.pixel, wcs=wcs160)
            except:
                continue

            distances = bg_coord.separation(source_coords[subset]).to(u.arcsec)/0.06
            within_reach = np.where(distances.value < offset/2)
            too_close = np.where(distances[within_reach].value <= offset/2.5)[0]

            if(distances.min().value > offset/1.5):
                continue
    
            aperture = CircularAperture((offset/2, offset/2), r=offset/3)
            LT = aperture_photometry(bg160.data, aperture)['aperture_sum'][0]
            light_outside = bg160.data.sum() - LT
        
            if(LT == 0 or LT > 20 or LT > light_outside):
                continue
                
            if(len(too_close)==0):
                reiterate = False

            num_iters += 1

        bg125 = Cutout2D(field160[0].data, position=bg_coord, size=offset*u.pixel, wcs=wcs160)
        bg = bg125.data
    else:
        bg160 = Cutout2D(field160[0].data, position=bg_coord, size=offset*u.pixel, wcs=wcs160)
        bg = bg160.data

    return bg, bg_coord

    
def find_bg_coords(bg, img):

    bg = random_flips(bg)

    a, b = bg.shape
    c, d = img.shape

    if((c <= a) & (d <= b)):
        x = np.random.randint(0, a-c)
        y = np.random.randint(0, b-d)
    else:
        raise(IndexError('Galaxy Image larger than BG'))

    return x, c, y, d