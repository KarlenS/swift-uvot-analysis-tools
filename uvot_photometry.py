#!/local/gammasoft/anaconda2/bin/python

'''
Use ``UVOTSOURCE`` to get aperture photometry from UVOT images.
'''
import os.path as path
import numpy as np
import subprocess

from astropy.time import Time
from astropy.io import fits

class MeasureSource(object):
    '''Class to measure photometry on UVOT data using the ``UVOTSOURCE`` tool and to parse the results.

    Attributes:
        filepath (str): path of fits image to input into ``UVOTSOURCE`` 
        source_coords (`astropy.coordinates.SkyCoord`_): coordinates of the source
        dirpath (str): directory containing ``filepath``
        obs (str): observation ID parsed from ``filepath``
        band (str): observation band parsed from ``filepath``

    .. _astropy.coordinates.SkyCoord:
        http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord
    '''

    def __init__(self,filepath='',default_fs = True):

        self.filepath = filepath
        self.source_coords = None
        
        #parse filename, directory name, obs ID, and band from filepath
        dirpath,filename = path.split(self.filepath)
        self.dirpath = dirpath
        if default_fs:
            base,extn = filename.split('_')
            self.obs = base[2:-3]
            self.band = base[-2:]
        else:
            self.obs = filename[8:23]
            self.band = filename[5:7]

    def correct_extinction(self, val, filtr, EBminV = None, mag = False):
    
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
            extTable = IrsaDust.get_extinction_table(self.source_coords)
            EBminV = np.median(extTable['A_SandF']/extTable['A_over_E_B_V_SandF'])
    
        #calculate extinction magnitude
        A_lambda = R_lambda[filtr]*EBminV
    
        if mag:
            return val - R_lambda[filtr]*EBminV
        else:
            return val*central_wav[filtr]*10**(R_lambda[filtr]*EBminV/2.5)
    

    def run_uvotsource(self):
        '''Wrapper to run the ``UVOTSOURCE`` tool for extracting photometry. The input/output file information
        is taken from the class attributes, which are obtained by parsing the provided ``filepath`` argument.

        Warning:
            Running this will override previously generated output files with the same names.

        Returns:
            list: 2 elements corresponding to the output and error messages from ``UVOTSOURCE``
        '''
        
        srcregfile = path.join(self.dirpath,'detect_%s_%s.reg' %(self.obs,self.band))
        bkgregfile = path.join(self.dirpath,'back_%s_%s.reg' %(self.obs,self.band))
        uvotsourcefile = path.join(self.dirpath,'uvotsource_%s_%s.fits' %(self.obs,self.band))

        tmp = subprocess.Popen(['uvotsource','image=%s' %self.filepath,'srcreg=%s' %srcregfile, 'bkgreg=%s' %bkgregfile, 
                                'outfile=%s' %uvotsourcefile,'chatter=0','sigma=3','clobber=YES'], stdout=subprocess.PIPE)
        
        return tmp.communicate()

    def get_observation_time(self):
        '''Uses the header of the image corresponding to ``filepath`` to get the start and end time of the observation.

        Returns:
            float: the middle of the observation time in ISOt format
        '''

        mainfits = fits.open(self.filepath)
        obs_start = Time(mainfits[0].header['DATE-OBS'],format='isot')
        obs_end = Time(mainfits[0].header['DATE-END'],format='isot')


        try:
            if mainfits[1].header['ASPCORR'] == 'DIRECT':
                aspflag = True
            else:
                aspflag = False
        except KeyError:
            print 'Assuming apsflag is fine for %s' %self.filepath
            aspflag = True

        mainfits.close()

        return obs_start + (obs_end - obs_start)/2,aspflag 

    def get_observation_data(self):
        '''Parse the output of ``UVOTSOURCE`` to extract essential photometry information.

        Returns:
            list: observing band, observing time (MJD), magnitude and 1$\sigma$ error in UVOT system, flux and 1$\sigma$ error in erg/s/cm$^{2}$/$\AA$, flux and 1$\sigma$ error in milliJanskies
        '''

        uvotsourcefile = path.join(self.dirpath,'uvotsource_%s_%s.fits' %(self.obs,self.band))
        try:
            data = fits.getdata(uvotsourcefile)
        except IOError:
            print '''%s not found. Make sure to run uvotsource beforehand or this will get annoying.\n
                   Will try running uvotsource now and will skip this observation if it fails.''' %uvotsourcefile
            tmp = self.run_uvotsource()
            try:
                data = fits.getdata(uvotsourcefile)
                print 'Running uvotsource worked.'
            except IOError:
                print 'Failed again... Skipping %s' %self.filepath
                return None


        obstime,corrflg = self.get_observation_time()
        mag = data['MAG'][0]
        magerr = data['MAG_ERR'][0]
        flux = data['FLUX_AA'][0]
        fluxerr = data['FLUX_AA_ERR'][0]
        fluxj = data['FLUX_HZ'][0]
        fluxjerr = data['FLUX_HZ_ERR'][0]
        fluxextcorr = self.correct_extinction(fluxj,self.band)
        fluxextcorrerr = self.correct_extinction(fluxjerr,self.band)
        return [self.band,obstime.mjd,mag,magerr,flux,fluxerr,fluxj,fluxjerr,fluxextcorr,fluxextcorrerr]


