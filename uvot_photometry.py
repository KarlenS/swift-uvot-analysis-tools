'''
Use ``UVOTSOURCE`` to get aperture photometry from UVOT images.
'''
import os.path as path
import subprocess

from astropy.time import Time
from astropy.io import fits

class MeasureSource(object):
    '''Class to measure photometry on UVOT data using the ``UVOTSOURCE`` tool and to parse the results.

    Attributes:
        filepath (str): path of fits image to input into ``UVOTSOURCE`` 
        dirpath (str): directory containing ``filepath``
        obs (str): observation ID parsed from ``filepath``
        band (str): observation band parsed from ``filepath``
    '''

    def __init__(self,filepath=''):

        self.filepath = filepath
        
        #parse filename, directory name, obs ID, and band from filepath
        dirpath,filename = path.split(self.filepath)
        base,extn = filename.split('_')
        self.dirpath = dirpath
        self.obs = base[2:-3]
        self.band = base[-2:]

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
        mainfits.close()

        return obs_start + (obs_end - obs_start)/2

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


        obstime = self.get_observation_time()
        mag = data['MAG'][0]
        magerr = data['MAG_ERR'][0]
        flux = data['FLUX_AA'][0]
        fluxerr = data['FLUX_AA_ERR'][0]
        fluxj = data['FLUX_HZ'][0]
        fluxjerr = data['FLUX_HZ_ERR'][0]
        return [self.band,obstime.mjd,mag,magerr,flux,fluxerr,fluxj,fluxjerr]
