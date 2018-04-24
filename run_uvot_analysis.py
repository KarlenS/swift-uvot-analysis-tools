#!/local/gammasoft/anaconda2/bin/python

'''
This is the main module for an automated Swift-UVOT analysis.

The process consists of three primary steps:
    * Run ``UVOTDETECT`` to get a proper centroid for the source for each observation, and generate region files for source and background regions,
    * Visually check images for potential issues/artifacts (and manually reject observations),
    * Run ``UVOTSOURCE`` to obtain photometry and parse data for light curve and SED use.

    TODO:
        Add automatic data downloading (possibly with astroquery.heasarc)
'''
import os.path
import sys
import glob
import numpy as np
import argparse
from tqdm import tqdm

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time


class uvot_runner(object):
    '''Wrapper class for the wrappers! Class to coordinate running the different stages of the UVOT analysis.

    Attributes:
        filepath (str): path of image file used for analysis
        source_coords (`astropy.coordinates.SkyCoord`_): coordinates of the source
        bkg_coords (`astropy.coordinates.SkyCoord`_): coordinates of the background region center
        ebminv (float): E(B-V) value for the source.

    .. _astropy.coordinates.SkyCoord:
        http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord
    '''
    def __init__(self):
        self.filepaths = None
        self.source_coords = None
        self.bkg_coords = None
        self.ebminv = None

    def query_for_source_coords(self, source_name):
        '''Queries NED or SIMBAD databases for coordinates using a source name.

        Args:
            source_name (str): source name to use for querying for coordinates 
        '''
    
        import astroquery
        from astroquery.ned import Ned

        try:
            coords = Ned.query_object(source_name)
            self.source_coords = SkyCoord('%s %s' %(coords['RA(deg)'][0], coords['DEC(deg)'][0]),unit=(u.deg, u.deg))
        except astroquery.exceptions.RemoteServiceError:
            from astroquery.simbad import Simbad
            coords = Simbad.query_object(source_name)
            self.source_coords = SkyCoord('%s %s' %(coords['RA'][0], coords['DEC'][0]),unit=(u.hourangle, u.deg))
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
            self.source_coords = iv.prime_source()

        iv.source_coords = self.source_coords

        self.bkg_coords = iv.prime_bkg() 


    def uvot_detecter(self,default_fs = True):
        '''Wrapper for running ``UVOTDETECT`` on all requested observations and extracting position information.
        '''

        from source_position import PositionExtractor

        pos = PositionExtractor(default_fs = default_fs)

        pos.source_coords = self.source_coords
        pos.bkg_coords = self.bkg_coords

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
            iv.setup_frame()
            x = raw_input('Viewing %s. Hit Enter to continue.' %filepath)


    def uvot_measurer(self,measure=True,default_fs = True):
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
        for filepath in tqdm(self.filepaths):
            measurer = MeasureSource(filepath,default_fs = default_fs)
            measurer.source_coords = self.source_coords

            if measure:
                tmp = measurer.run_uvotsource()
             
            objphot = measurer.get_observation_data()

            ptab.add_row(objphot)
        
        return ptab.group_by('filter')
        

    def uvot_sed_maker(self,startdate,enddate,outbase):
        from uvot_sed import MakeSED

        #preparing files
        sed = MakeSED(self.filepaths)
        sortedpaths = sed.sortFilelistByFilter()
        sed.sortedpaths = sortedpaths
        self.filepaths = sed.combineFits(startdate,enddate)

        self.uvot_primer()
        self.uvot_detecter(default_fs = False)
        photometry = self.uvot_measurer(default_fs = False)
        photometry.write('%s_%s-%s_sed.fits'%(outbase,startdate.mjd,enddate.mjd),overwrite=True)


def main():

    parser = argparse.ArgumentParser(description='Quick verification of image quality and regions selection for UVOT analysis.')
    parser.add_argument('-p',required=True, help='Parent directory of UVOT data structures (directories with names like 00030901252).')
    parser.add_argument('-obs',default=None,help='If specified, will only look at data for a single observations: e.g. 00030901252')
    parser.add_argument('-f',default=None, help='File with list of (full or relative) directory paths of UVOT data structures. Use this option if UVOT data structure directory names do not start with 000 for some reason...')
    parser.add_argument('-o',default='photometry.fits',help='Name of output fits file for storing resulting photometry.')
    parser.add_argument('-c',default=None,help='Specify coordinates for the source of interest.')
    parser.add_argument('-s',default=None,help='Specify the name for the source of interest. Ignored if -c is specified.')
    parser.add_argument('-ebv',default=None,help='User-specified E(B-V) value. Will query IRSA database, if not specified.')
    parser.add_argument('-date_range',default='1990-01-01,2020-01-01',help='Comma-separated beginning and end of date(time) range to average images for SED construction. Format is somewhat flexible, but try \'YYYY-MM-DD HH:MM,YYYY-MM-DD HH:MM\'')
    parser.add_argument('--check',action='store_true',default=False,help='View all images with the source marked to check observation quality.')
    parser.add_argument('--detect',action='store_true',default=False,help='Run UVOTDETECT on all images.')
    parser.add_argument('--measure',action='store_true',default=False,help='Run UVOTSOURCE on all images to get photometry.')
    parser.add_argument('--extract_only',action='store_true',default=False,help='Extract photometry information from existing UVOTSOURCE output files, without rerunning UVOTSOURCE.')
    parser.add_argument('--sed',action='store_true',default=False,help='Extract source SED by combining all observations/images located in the path (specified with -p) and falling within the date range (specified by -date_range)')
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
    if not ( args.measure or args.check or  args.detect or  args.extract_only or args.sed ):
        raise parser.error('Nothing to do. Use --detect, --check, --measure, or --extract_only flags for desired operations or -h for help.')

    if args.sed and args.date_range == '1990-01-01,2020-01-01':
        print 'Option -date_range not specified. Will use all provided observations for SED.'


    #run all the wrappers, though unless --detect --check or --measure are given, nothing will happen!
    runner = uvot_runner()
    runner.filepaths = filepaths


    #use coordinates either provided by the user, coordinates looked up based on the source name provided by the user, or coordinates selected interactively by the user.
    prime_source = False
    if args.c:
        #need to get the parsing right
        runner.source_coords = SkyCoord(args.c)
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
        photometry.write(args.o,overwrite=True)

    if args.sed:
        start,end = args.date_range.split(',')
        if args.s:
            outbase = args.s
        else:
            outbase = 'Source'

        try:
            runner.uvot_sed_maker(Time(start.lstrip()),Time(end.lstrip()),outbase)
        except ValueError:
            if start.lstrip().isdigit() and start.lstrip().isdigit():
                runner.uvot_sed_maker(Time(np.float(start.lstrip()),format='mjd'),Time(np.float(end.lstrip()),format='mjd'),outbase)
            else:
                raise ValueError('Invalid date_range format. Try providing dates in MJD or ISO formats.')




if __name__ == '__main__':
    main()

