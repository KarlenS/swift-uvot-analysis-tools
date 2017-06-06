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
import numpy
import argparse

from astroquery.simbad import Simbad
from astroquery.ned import Ned


def correct_extinction(val, filtr, source_coords, EBminV = None, mag = False):

    '''Function to correct for Galactic extinction using values from IRSA

    
    .. note:: Central wavelengths for UVOT filters are taken from `Poole et al. (2008)`_. 
              
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
        EBminV = np.median(extTable['A_SandF']/x['A_over_E_B_V_SandF'])
    
    A_lambda = R_lambda[filtr]*EBminV

    if mag:
        return val - R_lambda[filtr]*EBminV
    else:
        return val*central_wav[filtr]*10**(R_lambda[filtr]*EBminV/2.5)

def class uvot_runner(object):

    def __init__(self):
        self.filepaths = None
        self.source_ra = None
        self.source_dec = None
        self.ebminv = None

    def uvot_detecter(self):

        from source_position import PositionExtractor

        pos = PositionExtractor()

        pos.source_ra = self.source_ra
        pos.source_dec = self.source_dec

        for filepath in self.filepaths:
            print 'working on %s...' %filepath 
            pos.filepath = filepath
            detect_out = pos.run_uvotdetect()
            nearest = pos.getNearestSource()
            pos.cleanup()


    def uvot_checker(self):

        import pyds9 as ds9
        from check_images import SourceImageViewer

        #launching a ds9 window
        d = ds9.DS9()

        for filepath in self.filepaths:

            iv = SourceImageViewer()

            iv.filepath = filepath
            iv.source_ra = self.source_ra 
            iv.source_dec = self.source_dec
            iv.bkg_ra = '8:54:48.772'
            iv.bkg_dec = '+20:05:32.576'
            iv.setup_frame(d)
            x = raw_input('Viewing %s. Hit Enter to continue.' %filepath)


    def uvot_measurer(self,measure=True):
    
        from astropy.table import Table
        from astropy import units as u
        from uvot_photometry import MeasureSource
    
    
        ptab = Table(names=('filter','MJD','Mag','MagErr','FluxDensity','FluxDensityErr','FluxDensityJy','FluxDensityJyErr','FluxExtCorr','FluxExtcorrErr'),
                                dtype=('S2','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
    
        ptab['MJD'].unit = u.d
        ptab['Mag'].unit = u.mag
        ptab['MagErr'].unit = u.mag
        ptab['FluxDensity'].unit = u.erg/u.cm/u.cm/u.second/u.AA
        ptab['FluxDensityErr'].unit = u.erg/u.cm/u.cm/u.second/u.AA
        ptab['FluxDensityJy'].unit = u.Jy
        ptab['FluxDensityJyErr'].unit = u.Jy
        ptab['Flux'].unit = u.erg/u.cm/u.cm/u.second 
        ptab['FluxErr'].unit = u.erg/u.cm/u.cm/u.second
    
        for filepath in self.filepaths:
            measurer = MeasureSource(filepath)
            filtr = measurer.band
            if measure:
                tmp = measurer.run_uvotsource()
             
            objphot = measurer.get_observation_data()
    
            #getting extinction-corrected fluxes
            flux = correct_extinction(objphot[-4],filtr,EBminV = self.ebminv)
            fluxerr = correct_extinction(objphot[-3],filtr,EBminV = self.ebminv)
            objphot.append(flux)
            objphot.append(fluxerr)
    
            ptab.add_row(objphot)
    
        #flux_extcorr = correct_exctinction(ptab['FluxDensity'],ptab['filter'],ebminv)
        #fluxerr_extcorr = correct_exctinction(ptab['FluxDensityErr'],ptab['filter'],ebminv)
        
    
        return ptab.group_by('filter')
        

def main():

    parser = argparse.ArgumentParser(description='Quick verification of image quality and regions selection for UVOT analysis.')
    parser.add_argument('-p',required=True, help='Parent directory of UVOT data structures (directories with names like 00030901252).')
    parser.add_argument('-obs',default=None,help='If specified, will only look at data for a single observations: e.g. 00030901252')
    parser.add_argument('-f', help='File with list of directory names of UVOT data structures. Use this option if directory names do not start with 000 for some reason...')
    parser.add_argument('-o',default='photometry.fits',help='Name of output fits file for storing resulting photometry.')
    parser.add_argument('-c',default=None,help='Specify coordinates for the source of interest.')
    parser.add_argument('-s',default=None,help='Specify the name for the source of interest. Ignored if -c is also specified.')
    parser.add_argument('-ebv',default=None,help='User-specified E(B-V) value. Will query IRSA database, if not specified.')
    parser.add_argument('--check',action='store_true',default=False,help='View all images with the source marked to check observation quality.')
    parser.add_argument('--detect',action='store_true',default=False,help='Run uvotdetect on all images.')
    parser.add_argument('--measure',action='store_true',default=False,help='Run uvotsource on all images to get photometry.')
    parser.add_argument('--print_only',action='store_true',default=False,help='Only plot light curves.')
    args = parser.parse_args()

    if not args.obs:
        filepaths = glob.iglob('%s/000*/uvot/image/*sk.img.gz' %args.p)#setting up filepaths
    else:
        obspath = os.path.join(args.p,str(args.obs))
        if os.path.exists(obspath):
            filepaths = glob.iglob('%s/uvot/image/*sk.img.gz'%obspath)
        else:
            raise NameError('%s does not exist.' %obspath)

    if args.c:
        coords = args.c
    elif args.s:
        coords = ned.query_object(args.s) 
    else:
        raise parser.error('Please specify a source name with -s option or coordinates with -c option.')

    runner = uvot_runner()
    runner.filepaths = filepaths
    runner.source_coords = coords

    if args.detect:
        runner.uvot_detecter()

    if args.check:
        runner.uvot_checker()


    if args.measure:
        photometry = runner.uvot_measurer(measure = not args.print_only)
        photometry.write('OJ287_uvot_photometry_wExtCorr.dat',format='ascii.commented_header')
        photometry.write(args.o)

if __name__ == '__main__':
    main()

