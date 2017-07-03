'''
This is the main module for an automated Swift-UVOT analysis.

The process consists of three primary steps:
    * Run ``UVOTDETECT`` to get a proper centroid for the source for each observation, and generate region files for source and background regions,
    * Visually check images for potential issues/artifacts (and manually reject observations),
    * Run ``UVOTSOURCE`` to obtain photometry and parse data for light curve and SED use.
'''
import os.path
import sys
import glob
import numpy as np
import argparse

from astropy.coordinates import SkyCoord
from astropy import units as u

def correct_extinction(val, filtr, source_coords, EBminV = None, mag = False):

    '''Function to correct for Galactic extinction using values from IRSA

    
    Note:
        Central wavelengths for UVOT filters are taken from `Poole et al. (2008)`_.               
        :math:`R_{\lambda}` values are derived using the `York Extinction Solver`_.

    
    Args:
        val (float): flux or mag requiring extinction correction (flux is assumed, unless ``mag = True``)
        filtr (str): UVOT filter name (vv, uu, bb, w1, m2, w2)
        source_coords (`astropy.coordinates.SkyCoord`_): Source position to be used for querying the amount of extinction
    
    Returns:
        float: Extinction-corrected flux (erg/s/cm2) or magnitude (mag)

    .. _astropy.coordinates.SkyCoord:
        http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord
    .. _Poole et al. (2008):
        http://adsabs.harvard.edu/abs/2008MNRAS.383..627P
    .. _York Extinction Solver: 
        http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/YorkExtinctionSolver/coefficients.cgi        
    '''

    central_wav = {'uu':3465.,'w1':2600.,'m2':2246.,'w2':1928.,'bb':4392.,'vv':5468.}
    R_lambda = {'uu':4.89172,'w1':6.55663,'m2':9.15389,'w2':8.10997,'bb':4.00555,'vv':2.99692}


    #query for the E(B-V) value, unless user specifies one
    if not EBminV:
        from astroquery.irsa_dust import IrsaDust
        extTable = IrsaDust.get_extinction_table(source_coords)
        EBminV = np.median(extTable['A_SandF']/extTable['A_over_E_B_V_SandF'])
    
    #calculate extinction magnitude
    A_lambda = R_lambda[filtr]*EBminV

    if mag:
        return val - R_lambda[filtr]*EBminV
    else:
        return val*central_wav[filtr]*10**(R_lambda[filtr]*EBminV/2.5)


class uvot_runner(object):
    '''Wrapper class for the wrappers! Class to coordinate running the different stages of the UVOT analysis.

    Attributes:
        filepath (str): path of image file used for analysis
        source_ra (float): right ascension coordinate of source
        source_dec (float): declination coordinate of source
        bkg_ra (float): right ascension coordinate of background region center
        bkg_dec (float): declination coordinate of background region center
        ebminv (float): E(B-V) value for the source.
    '''
    def __init__(self):
        self.filepaths = None
        self.source_coords = None
        self.bkg_coords = None
        #self.source_ra = None
        #self.source_dec = None
        #self.bkg_ra = None
        #3self.bkg_dec = None
        self.ebminv = None

    def query_for_source_coords(self, source_name):
    
        import astroquery
        from astroquery.ned import Ned
        try:
            coords = Ned.query_object(source_name)
            self.source_coords = SkyCoord('%s %s' %(coords['RA(deg)'][0], coords['DEC(deg)'][0]),unit=(u.deg, u.deg))
            #self.source_ra, self.source_dec = coords['RA(deg)'][0], coords['DEC(deg)'][0]
        except astroquery.exceptions.RemoteServiceError:
            from astroquery.simbad import Simbad
            coords = Simbad.query_object(source_name)
            self.source_coords = SkyCoord('%s %s' %(coords['RA'][0], coords['DEC'][0]),unit=(u.hourangle, u.deg))
            #self.source_ra, self.source_dec = coords['RA'][0], coords['DEC'][0]
            if not coords:
                raise astroquery.exceptions.RemoteServiceError('Both NED and SIMBAD queries failed for the provided source name. Maybe try supplying coordinates or selecting the source position interactively.')
    
    def uvot_primer(self, prime_source = False):
        '''Makes user identify a location for the background region (and source location) using a DS9 window.

        Args:
            prime_source (bool): If true, in addition to prompting user to interactively select a center for background region, will also ask user to select location for a source.
        '''
        
        import pyds9 as ds9
        from check_images import SourceImageViewer

        #launching a ds9 window

        iv = SourceImageViewer()

        # picking first optical filter image for user background selection,
        # or just the first image, if there aren't optical ones.
        try:
            iv.filepath = next(filepath for filepath in self.filepaths if 'ubb' in filepath or 'uvv' in filepath)
        except:
            iv.filepath = self.filepaths[0]

        if prime_source:
            #self.source_ra, self.source_dec = iv.prime_source()
            self.source_coords = iv.prime_source()

        #iv.source_ra,iv.source_dec = self.source_ra, self.source_dec
        iv.source_coords = self.source_coords

        #self.bkg_ra, self.bkg_dec = iv.prime_bkg()
        #bkg_ra, bkg_dec = iv.prime_bkg()
        self.bkg_coords = iv.prime_bkg() 


    def uvot_detecter(self):
        '''Wrapper for running ``UVOTDETECT`` on all requested observations and extracting position information.
        '''

        from source_position import PositionExtractor

        pos = PositionExtractor()

        pos.source_coords = self.source_coords
        pos.bkg_coords = self.bkg_coords
        #pos.source_ra = self.source_ra
        #pos.source_dec = self.source_dec

        for filepath in self.filepaths:
            print 'working on %s...' %filepath 
            pos.filepath = filepath
            out, err = pos.run_uvotdetect()
            pos.getNearestSource()
            pos.cleanup()


    def uvot_checker(self):
        '''Wrapper for viewing images from all requested observations to verify the source and background regions and image quality.
        '''

        from check_images import SourceImageViewer

        #launching a ds9 window

        for filepath in self.filepaths:

            iv = SourceImageViewer()
        
            iv.source_coords = self.source_coords
            iv.bkg_coords = self.bkg_coords

            iv.filepath = filepath
            #iv.source_ra = self.source_ra 
            #iv.source_dec = self.source_dec
            #iv.bkg_ra = self.bkg_ra 
            #iv.bkg_dec = self.bkg_dec
            iv.setup_frame()
            x = raw_input('Viewing %s. Hit Enter to continue.' %filepath)


    def uvot_measurer(self,measure=True):
        '''Wrapper for running ``UVOTDETECT`` and extracting photometry information from all requested observations.

        Args:
            measure (bool): If True, will run ``UVOTDETECT`` otherwise will simply parse existing output files.

        Returns:
            `astropy.table.Table`_: Astropy Table with photomoetry information grouped by observation filter

        .. _astropy.table.Table: 
            http://docs.astropy.org/en/stable/api/astropy.table.Table.html
        '''
        
        from astropy.table import Table
        from uvot_photometry import MeasureSource
    
        #creating Astropy table for storing the photometry information
        ptab = Table(names=('filter','MJD','Mag','MagErr','FluxDensity','FluxDensityErr','FluxDensityJy','FluxDensityJyErr','FluxExtCorr','FluxExtCorrErr'),
                                dtype=('S2','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
    
        #defining units for each column
        ptab['MJD'].unit = u.d
        ptab['Mag'].unit = u.mag
        ptab['MagErr'].unit = u.mag
        ptab['FluxDensity'].unit = u.erg/u.cm/u.cm/u.second/u.AA
        ptab['FluxDensityErr'].unit = u.erg/u.cm/u.cm/u.second/u.AA
        ptab['FluxDensityJy'].unit = u.Jy
        ptab['FluxDensityJyErr'].unit = u.Jy
        ptab['FluxExtCorr'].unit = u.erg/u.cm/u.cm/u.second 
        ptab['FluxExtCorrErr'].unit = u.erg/u.cm/u.cm/u.second
    
        #run through all image files to perform and/or extract photometry
        for filepath in self.filepaths:
            measurer = MeasureSource(filepath)
            filtr = measurer.band
            if measure:
                tmp = measurer.run_uvotsource()
             
            objphot = measurer.get_observation_data()
    
            #getting extinction-corrected fluxes
            flux = correct_extinction(objphot[-4],filtr,self.source_coords,EBminV = self.ebminv)
            fluxerr = correct_extinction(objphot[-3],filtr,self.source_coords,EBminV = self.ebminv)
            objphot.append(flux)
            objphot.append(fluxerr)
    
            ptab.add_row(objphot)
        
        return ptab.group_by('filter')
        

def main():

    parser = argparse.ArgumentParser(description='Quick verification of image quality and regions selection for UVOT analysis.')
    parser.add_argument('-p',required=True, help='Parent directory of UVOT data structures (directories with names like 00030901252).')
    parser.add_argument('-obs',default=None,help='If specified, will only look at data for a single observations: e.g. 00030901252')
    parser.add_argument('-f',default=None, help='File with list of (full or relative) directory paths of UVOT data structures. Use this option if UVOT data structure directory names do not start with 000 for some reason...')
    parser.add_argument('-o',default='photometry.fits',help='Name of output fits file for storing resulting photometry.')
    parser.add_argument('-c',default=None,help='Specify coordinates for the source of interest.')
    parser.add_argument('-s',default=None,help='Specify the name for the source of interest. Ignored if -c is specified.')
    parser.add_argument('-ebv',default=None,help='User-specified E(B-V) value. Will query IRSA database, if not specified.')
    parser.add_argument('--check',action='store_true',default=False,help='View all images with the source marked to check observation quality.')
    parser.add_argument('--detect',action='store_true',default=False,help='Run UVOTDETECT on all images.')
    parser.add_argument('--measure',action='store_true',default=False,help='Run UVOTSOURCE on all images to get photometry.')
    parser.add_argument('--extract_only',action='store_true',default=False,help='Extract photometry information from existing UVOTSOURCE output files, without rerunning UVOTSOURCE.')
    args = parser.parse_args()

    #use all observations in the provided directory unless a specific obs is requested
    if args.f:
        filedirs = np.genfromtxt(args.f,dtype=str)
        tmppaths = [glob.glob('%s/uvot/image/*sk.img.gz'%d) for d in filedirs ]
        filepaths = [path for paths in tmppaths for path in paths]
        print filepaths
    elif not args.obs:
        filepaths = glob.glob('%s/000*/uvot/image/*sk.img.gz' %args.p)#setting up filepaths
    else:
        obspath = os.path.join(args.p,str(args.obs))
        if os.path.exists(obspath):
            filepaths = glob.glob('%s/uvot/image/*sk.img.gz'%obspath)
            print 'Will use images from observation %s only' %args.obs
        else:
            raise NameError('%s does not exist.' %obspath)

    #raise an error if an operation is not selected.
    if not ( args.measure or args.check or  args.detect or  args.extract_only ):
        raise parser.error('Nothing to do. Use --detect, --check, --measure, or --extract_only flags for desired operations or -h for help.')


    #run all the wrappers, though unless --detect --check or --measure are given, nothing will happen!
    runner = uvot_runner()
    runner.filepaths = filepaths


    #use coordinates either provided by the user, coordinates looked up based on the source name provided by the user, or coordinates selected interactively by the user.
    prime_source = False
    if args.c:
        #need to get the parsing right
        runner.source_coords = SkyCoord(args.c)
        #runner.source_ra, runner.source_dec = args.c
    elif args.s:
        #query NED for the source position using provided source name and set source_ra and source_dec
        runner.query_for_source_coords(args.s)
    else:
        prime_source = True
        print 'Source not specified by either -c or -s option. Will make user select the source location interactively.'
        print 'Quit (ctrl-C) and specifiy source name with -s or source coordinates with -c options if the source is known.'

    #get user to identify background region (and source region if coordinates or source name are not supplied)

    if args.detect:
        runner.uvot_primer(prime_source=prime_source)
        runner.uvot_detecter()

    if args.check:
        print 'User selected coordinates will only be used if region files are missing.'
        runner.uvot_primer(prime_source=prime_source)
        runner.uvot_checker()

    #can set up options for multiple formats here, but will probably default to fits
    if args.measure or args.extract_only:
        photometry = runner.uvot_measurer(measure = not args.extract_only)
        #photometry.write('OJ287_uvot_photometry_wExtCorr.dat',format='ascii.commented_header')
        photometry.write(args.o)


if __name__ == '__main__':
    main()

