'''
Module for visually inspecting and interacting with UVOT images using DS9, including selection of DS9 regions for aperture photometry.
'''

import numpy as np
import os.path as path
import pyds9 as ds9
import pyregion

from astropy.coordinates import SkyCoord
from astropy import units as u


class SourceImageViewer(ds9.DS9):
    '''Inherits from a `pyds9`_ object. 
    Initialize the class with ``filepath`` input. The radii of source and background are set to suggested UVOT values.
    These are hardcoded at the moment. If they require to be changed, `aperture correction in UVOTSOURCE`_ will be necessary,
    which is not currently included in these wrappers.

    Attributes:
        filepath (str): path to image to display in DS9
        source_coords (`astropy.coordinates.SkyCoord`_): Coordinates of the source
        bkg_coords (`astropy.coordinates.SkyCoord`_): Coordinate of background region center

    .. _pyds9:
        https://github.com/ericmandel/pyds9
    .. _aperture correction in UVOTSOURCE:
        https://heasarc.nasa.gov/lheasoft/ftools/headas/uvotsource.html
    .. _astropy.coordinates.SkyCoord:
        http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord
    '''

    def __init__(self,filepath=''):
        #this will launch a DS9 window if there isn't one already or connect to an existing one
        super(SourceImageViewer,self).__init__()
        self.filepath = filepath
        self.source_coords = None
        self.bkg_coords = None

    def open_fits(self):
        '''Tell DS9 to open ``filepath``
        ''' 
        opencmd = "file %s" %self.filepath
        self.set(opencmd)

    def display_source_region(self,radius = 5.):
        '''Tell DS9 to display a region corresponding to source coordinates.

        args:
            radius (float): radius of the source region
        ''' 
        reg = 'fk5; circle(%s,%s,%s")' %(self.source_coords.ra.value,self.source_coords.dec.value,radius)
        self.set('regions', reg);

    def display_bkg_region(self,radius = 20.):
        '''Tell DS9 to display a region corresponding to background center coordinates.

        args:
            radius (float): radius of the background region
        ''' 

        bkgreg = 'fk5; circle(%s,%s,%s")' %(self.bkg_coords.ra.value,self.bkg_coords.dec.value,radius)
        self.set('regions', bkgreg)

    def remove_regions(self):
        '''Tell DS9 to delete all regions currently open in a frame.
        '''
        self.set('regions delete all')

    def load_regions(self):
        '''Tell DS9 to load regions based on the provided ``filepath``. The `filepath`` gets parsed to 
        construct the name of the region files based on the observation ID and band. If the region files
        are not found in the parent path of ``filepath``, the default coordinates will get used for 
        the source (likely from SIMBAD) and for the background.
        '''

        #parse filepath to get obs ID and band
        dirpath,filename = path.split(self.filepath)
        base,extn = filename.split('_')
        obs = base[2:-3]
        band = base[-2:]

        #constructs the source and background region file paths
        regfile = path.join(dirpath,'detect_%s_%s.reg' %(obs,band))
        bkgregfile = path.join(dirpath,'back_%s_%s.reg' %(obs,band))

        #trying to load the background region. If the attempt fails, whatever coordinates are
        #currently stored in the object bkg_ra and bkg_dec attributes will be displayed.
        try:
            bkgregion = pyregion.open(bkgregfile)
            bkg_ra = bkgregion[0].coord_list[0]
            bkg_dec = bkgregion[0].coord_list[1]
            self.bkg_coords = SkyCoord('%s %s' %(bkg_ra, bkg_dec),unit=(u.deg, u.deg))
        except IOError:
            print 'Background region file missing. Will try default coordinates.'
        
        self.display_bkg_region()

        #trying to load the source region. If the attempt fails, whatever coordinates are
        #currently stored in the object source_ra and source_dec attributes will be displayed.
        try:
            region = pyregion.open(regfile)
            ra = region[0].coord_list[0]
            dec = region[0].coord_list[1]

            self.source_coords = SkyCoord('%s %s' %(ra,dec),unit=(u.deg, u.deg))
        except IOError:
            print 'Region file missing. Only showing default coordinates, if any.'

        self.display_source_region()

    def center_on_source(self):
        '''Tell DS9 to center the frame on the source location.
        '''
        centercmd = 'pan to %s %s wcs fk5' %(self.source_coords.ra.value,self.source_coords.dec.value)
        self.set(centercmd)

    def zoom_in(self,zoom = 5):
        '''Tell DS9 to zoom in.
        '''
        self.set('zoom to fit')
        self.set('zoom %s' %zoom)

    def prettify(self):
        '''Tell DS9 to scale the image to log and minmax and use the heat color map.
        '''
        self.set('scale mode minmax')
        self.set('scale log')
        self.set('cmap heat')

    def format_frame(self):
        '''Calls ``zoom_in`` and ``prettify`` methods.
        '''
        self.zoom_in()
        self.prettify()


    def setup_frame(self):
        '''Prepares a frame by opening a file, calling ``format_frame`` method, loading regions, and centering on the source position.
        '''
        self.open_fits()
        self.format_frame()
        self.load_regions()
        self.center_on_source()

    def get_user_coords(self):
        '''Tell DS9 to go into an interactive mode and get the coordinates of the position the user clicks.
        '''
        return self.get('iexam coordinate wcs fk5')

    def prime_bkg(self):
        '''Primes the user to select a position for the center of the background region.
        An image file will be displayed and formatted. The frame will be centered on the source position and display a source region.
        DS9 will then go into interactive mode and prompt the user to select a position for the background region center on the image.
        Once the user clicks on the image, a background region will be drawn at the clicked position and the user's happiness level 
        will be assessed.

        Returns:
            list: background region center RA and Dec coordinates
        '''

        self.open_fits()
        self.format_frame()
        self.center_on_source()

        #open loop to keep prompting the user for background selection if they're not happy after an initial attempt
        while True:
            self.display_source_region()
            print 'Use the DS9 window to select a location for the background region.'
            coords = self.get_user_coords()
            bkg_ra, bkg_dec = coords.split()
            self.bkg_coords = SkyCoord('%s %s' %(bkg_ra,bkg_dec),unit=(u.deg, u.deg))
            self.display_bkg_region()

            response = raw_input('Happy with the background selection (y / anything_else)?\n') 
            if response == 'y':
                print 'Great! That is all we need. Bye.'
                break
            else:
                print 'Fine, try again...'
                self.remove_regions()

        self.set('quit')
        return self.bkg_coords 

    def prime_source(self):
        '''Primes the user to select a position for a source. 
        An image file will be displayed and formatted. DS9 will then go into interactive mode and prompt the 
        user to click on the location for the source of interest. Once the user clicks on the image, a source region
        will be drawn at the clicked position and the user's happiness level will be assessed.

        Returns:
            list: source region center RA and Dec coordinates
        '''

        self.open_fits()
        self.format_frame()

        #open loop to keep prompting the user for source selection if they're not happy after an initial attempt
        while True:

            print 'Use the DS9 window to select a source region.'
            coords = self.get_user_coords()
            ra, dec = coords.split()
            self.source_coords = SkyCoord('%s %s' %(ra,dec),unit=(u.deg, u.deg))
            self.display_source_region()

            response = raw_input('Happy with source location selection (y/n)?\n') 
            if response == 'y':
                print 'Great! That is all we need. Bye.'
                break
            else:
                print 'Fine, try again...'
                self.remove_regions()

        return self.source_coords 
