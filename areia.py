__title__ = 'AREIA: Artificial Redshift Effects for IA'
__author__ = 'Leonardo Ferreira & Clar-Brid Tohill'
__version__ = '0.0.2'

import sys
import argparse
import numpy as np
import glob


from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM  
from astropy.convolution import convolve
from astropy import constants as const
from astropy import units as u

from photutils import detect_sources, detect_threshold

from .galclean import central_segmentation_map, measure_background

from scipy.ndimage import zoom
from scipy.interpolate import interpn

from matplotlib import pyplot as plt


def zoom_contents(image, scale, image_axes=[0, 1], method='linear', conserve_flux=True, fill_value=0):
    """
       Change the scale with fixed image size. Implemmentation by Aiden Rolfe (aidenrolfe) and
       Steven Bamford (bamford) from https://github.com/aidenrolfe/ARG.

       Suggestion by Bamford in Issue https://github.com/astroferreira/areia/issues/1
    """
    # Resize contents of image relative to fixed image size using interpolation
    in_coords = [np.arange(s) - (s - 1) / 2 for s in image.shape]
    out_coords = np.array(np.meshgrid(*in_coords, indexing='ij'))
    # match shape of input scale to allow broadcasting
    scale = np.atleast_1d(scale)
    scale = scale.reshape(scale.shape + (1,) * (image.ndim - scale.ndim))
    out_coords[image_axes] /= scale
    out_coords = np.transpose(out_coords)
    output = interpn(in_coords, image, out_coords,
                     method=method, bounds_error=False, fill_value=fill_value)
    output = output.T
    if conserve_flux:
        output /= scale**2

    return output



class Config(object):
    '''
        Tracks all the flags used by Artificial Redshift.
        If none is provided, it loads the default defined below.
    '''

    h = 0.7
    cosmo = FlatLambdaCDM(H0=100 * h, Om0=0.3, Tcmb0=2.725)
    add_background = True
    rebinning = True
    convolve_with_psf = True
    make_cutout = True
    dimming = True
    shot_noise = True
    size_correction = True
    evo = True
    evo_alpha = -0.13
    output_size = 101
    bg_centered = True
    fixed_grid = True

class ObservationFrame(object):
    '''
        Class that represents one observation frame with a given
        instrument setup.
    '''

    def __init__(self, redshift, pixelscale, exptime):
        self.pixelscale = pixelscale
        self.redshift = redshift
        self.exptime = exptime

class ArtificialRedshift(object):
    '''
        This handles all transformations and effects selected
        in the Config class to be applied to the input data,
        from initial_frame to target_frame. It keeps track
        of the transformation in the input image, it is possible 
        to retrieve partial results between each step, ideal for
        debugging.
    '''

    def __init__(self, image, psf, background, initial_frame, target_frame, MAG, bg_position=None, 
                 config=None):

        self.image = image
        self.psf = psf
        self.background = background
        self.initial_frame = initial_frame
        self.target_frame = target_frame
        self.MAG = MAG
        self.bg_position = bg_position

        if config is None:
            self.config = Config()
        else:
            self.config = config

        self.cosmo = self.config.cosmo

        self.final = self.image

        self.cutout_source() 
        self.geometric_rebinning() 
        self.apply_dimming()
        self.evolution_correction()
        self.convolve_psf()
        self.apply_shot_noise()
        self.add_background()
        self.crop_for_network(size=self.config.output_size)
       #self.to_candelizer()

    @classmethod
    def fromrawdata(cls, image,
                         psf,
                         background,
                         initial_redshift,
                         target_redshift, 
                         initial_pixelscale, 
                         target_pixelscale,
                         obs_exptime, 
                         target_exptime):

        current_frame = ObservationFrame(initial_redshift, initial_pixelscale, obs_exptime)
        target_frame = ObservationFrame(target_redshift, target_pixelscale, target_exptime)
    
        return cls(image, psf, background, initial_frame, target_frame)

    def cutout_source(self):

        if self.config.make_cutout:
            segmentation = central_segmentation_map(self.image)
            self.cutout = self.image.copy()
            self.masked = self.image.copy()
            self.masked[segmentation == True] = 0
            self.cutout[segmentation == False] = 0
            self.final = self.cutout.copy()

    def geometric_rebinning(self):

        def _size_correction(redshift):
            return (1 + redshift)**(-0.97)

        self.scale_factor = 1
        self.size_correction_factor = 1
        
        self.flux = self.final.sum()
        self.rebinned = self.final / self.flux

        if self.config.rebinning:           
            initial_distance = self.cosmo.luminosity_distance(self.initial_frame.redshift).value   
            target_distance = self.cosmo.luminosity_distance(self.target_frame.redshift).value   
            self.scale_factor = (initial_distance * (1 + self.target_frame.redshift)**2 * self.initial_frame.pixelscale) / (target_distance * (1 + self.initial_frame.redshift)**2 * self.target_frame.pixelscale)
        
        if self.config.size_correction:
            self.size_correction_factor = _size_correction(self.target_frame.redshift)
        
        if self.config.fixed_grid:
            self.rebinned = zoom_contents(self.rebinned, self.scale_factor * self.size_correction_factor, conserve_flux=False)
        else:
            self.rebinned = zoom(self.rebinned, self.scale_factor * self.size_correction_factor, order=0, prefilter=True)
            
        self.rebinned /= self.rebinned.sum()
        self.rebinned *= self.flux
        self.final = self.rebinned.copy()

    def apply_dimming(self):

        if self.config.dimming:
            self.dimming_factor = (self.cosmo.luminosity_distance(self.initial_frame.redshift) / self.cosmo.luminosity_distance(self.target_frame.redshift))**2
            self.dimming_factor = self.dimming_factor.value
            self.dimmed = self.final * self.dimming_factor
            self.final = self.dimmed.copy()


    def flat_evolution_correction(self):
        
        if self.config.evo:
            self.evo_factor = 10**(-0.4 * self.config.evo_alpha * (self.target_frame.redshift))
            self.with_evolution = self.final * self.evo_factor
            self.final = self.with_evolution.copy()

    def evolution_correction(self):
        
        if self.config.evo:
            self.mag_correction = (self.MAG - self.MAG*(1 + self.target_frame.redshift)**(self.config.evo_alpha))
            self.evo_factor = 10**(-0.4 * self.mag_correction)
            self.with_evolution = self.final * self.evo_factor
            self.final = self.with_evolution.copy()

    def convolve_psf(self):
        
        if self.config.convolve_with_psf:
            
            original_flux = self.final.sum()

            self.psf /= self.psf.sum()
            self.convolved = convolve(self.final, self.psf)   
            self.final = self.convolved.copy()

        
    def apply_shot_noise(self):         
        
        if self.config.shot_noise:
            self.shot_noise = np.sqrt(abs(self.final * self.target_frame.exptime)) * np.random.randn(self.final.shape[0], self.final.shape[1]) / self.target_frame.exptime         
            self.with_shot_noise = self.final + self.shot_noise
            self.final = self.with_shot_noise.copy()


    def _find_bg_section(self, bg, img, CENTER=False):
    
        a, b = bg.shape
        c, d = img.shape 
        c = int(c)
        d = int(d)
        if CENTER:
            xi = int(a/2 - c/2)
            xf = int(b/2 + d/2)
            bg_section = bg[xi:xf, xi:xf] 
        else:
            if self.bg_position is None:
                if((c <= a) & (d <= b)):
                    x = np.random.randint(0, a-c)
                    y = np.random.randint(0, b-d)
                else:
                    raise(IndexError('Galaxy Image larger than BG'))
            else:
                x, y = self.bg_position

            bg_section = bg[x:(x+c), y:(y+d)]

        return bg_section

    def add_background(self):



        if self.config.add_background:
            if self.background is None:
                background_shape = self.final.shape
                if self.scale_factor > 1:
                    background_shape = self.final.shape
                
                mean, median, std = measure_background(self.final, 2, np.zeros((background_shape)))
                self.background = np.random.normal(0, std, size=background_shape)
            else:
                self.background = self._find_bg_section(self.background, self.final, CENTER=self.config.bg_centered) if self.scale_factor <= 1 else self._find_bg_section(self.background, self.final, CENTER=self.config.bg_centered)
                
            source_shape = self.final.shape
            
            offset = 1
            if source_shape[0] % 2 == 0:
                offset = 0

            offset_min = int(self.background.shape[0]/2) - int(np.floor(source_shape[0]/2)) 
            offset_max = int(self.background.shape[0]/2) + int(np.floor(source_shape[0]/2)) + offset


            self.with_background = self.background.copy()
            self.with_background[offset_min:offset_max, offset_min:offset_max] += self.final
            
            self.final = self.with_background.copy()

    def crop_for_network(self, size):

        if self.final.shape[0] > size:

            half_size = int(size/2)    
            xc = int(self.final.shape[0]/2)
            yc = int(self.final.shape[1]/2)

            self.final_crop = self.final[yc-half_size:yc+half_size, xc-half_size:xc+half_size]
        else:
            self.final_crop = self.final


    def writeto(self, filepath, data, overwrite=False):

        hdr = fits.Header()

        if self.config.rebinning:
            hdr['REBIN'] = self.scale_factor
        
        if self.config.dimming:
            hdr['DIM_FACT'] = self.dimming_factor

        if self.config.evo:
            hdr['EVO_FACT'] = self.evo_factor
            hdr['EV_ALPHA'] = self.config.evo_alpha

        hdr['HISTORY'] = 'Image simulated with AREIA: Artificial Redshift Effects for IA'
        hdr['HISTORY'] = 'Source Extracted with Galclean'
        hdr['HISTORY'] = f'From z = {self.initial_frame.redshift} to z = {self.target_frame.redshift}'
        hdr['HISTORY'] = 'Any issues forward it to leonardo.ferreira@nottingham.ac.uk'
        hdr['HISTORY'] = 'Or seek help at https://github.com/astroferreira/areia'

        fits.writeto(filename=filepath, data=data, header=hdr, overwrite=overwrite)
        












    

    

