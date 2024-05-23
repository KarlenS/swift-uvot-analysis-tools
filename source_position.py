#!/Users/karlen/anaconda2/envs/astroconda/bin/python

'''
Module for determining the source position using ``UVOTDETECT`` and for generating source and background region files used in the photometry.
'''
import os
import os.path as path
import glob
import subprocess
import numpy as np

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord


class PositionExtractor(object):
    '''Class to handle position-related requests for the UVOT analysis chain, including running ``UVOTDETECT`` and parsing its results.

    Attributes:
        source_coords (`astropy.coordinates.SkyCoord`_): coordinate of the source
        bkg_coords (`astropy.coordinates.SkyCoord`_): coordinate of the background region center
        filepath (str): path of image file used for analysis
        detect (str): path of ``UVOTDETECT`` output file
        regfile (str): path of source region file
        bkgregfile (str): path of background region file

    .. _astropy.coordinates.SkyCoord:
        http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord

    '''

    def __init__(self,default_fs=True):
        self.source_coords = None
        self.bkg_coords = None
        self.filepath = None
        self.detect = None
        self.regfile = None
        self.bkgregfile = None
        self.default_fs = default_fs

    def cleanup(self):
        '''Cleans up ``UVOTDETECT`` output files after necessary information has been extracted from them.
        '''
        try:
            os.remove(self.detect)
        except OSError:
            print('File %s does not exist. Cannot delete what does not exist.' %self.detect)
            pass

    def createRegionFiles(self):
        '''Creates source and background region files, using ``source_ra``, ``source_dec``, ``bkg_ra``, and ``bkg_dec``  attributes. Filenames are generated using the ``filepath`` attribute.
        '''
        
        dirpath,filename = path.split(self.filepath)

        # write source region file
        region = 'fk5;circle(%s,%s,5")' %(self.source_coords.ra.value,self.source_coords.dec.value)
        regfile = open(self.regfile,'w')
        regfile.write(region)
        regfile.close()

        # write background region file if it doesn't already exist 
        #if not path.isfile(self.bkgregfile):
        bkgregion = 'fk5;circle(%s,%s,20")' %(self.bkg_coords.ra.value,self.bkg_coords.dec.value)
        bkgregfile = open(self.bkgregfile,'w')
        bkgregfile.write(bkgregion)
        bkgregfile.close()
        #else:
        #    print 'Background region file already exists. Not recreating.'

  
    def run_uvotimsum(self, inputFile, outputFile):
        '''Wrapper for running ``UVOTIMSUM`` for combining multiple extensions in a fits file.

        args:
            inputFile (str): name of the input file for ``UVOTIMSUM``
            outputFile (str): name of the output file for ``UVOTIMSUM``
        '''

        # using  subprocess to run uvotimsum
        tmp = subprocess.Popen(["uvotimsum",inputFile,outputFile,'chatter=0','clobber=yes','exclude=NONE'], stdout=subprocess.PIPE)
        tmp = tmp.communicate()

    def getNearestSource(self):
        '''Uses the results of ``UVOTDETECT`` to get the refined source position per image. If this is unsuccessful, the default coordinates provided to the class are used. Once this is sorted out, the function orders the region file creation.
        '''

        try:
            data = fits.getdata(self.detect)
        except:
            print('File %s does not exist. %s is likely a multiple extension file. Attempting to combine extensions and rerunning.' %(self.detect,self.filepath))
            try:
                # running uvotimsum to combine multiple extensions
                # to keep things relatively clean, changing filepath to combined while remembering orignal name
                comb = '%s.comb' %self.filepath
                orig = self.filepath 
                self.run_uvotimsum(self.filepath,comb)
                self.filepath = comb
                # trying to run uvotdetect again
                detect_out = self.run_uvotdetect(exp=None)
                data = fits.getdata(self.detect)

                # once we have the data from combined, changing back filepath to original
                self.filepath = orig
                print('Combining and rerunning seems to have worked!')

            except:
                print('Attempt to combine extensions failed. Defaulting to NED/SIMBAD position.')
                self.createRegionFiles()

        if np.size(data) > 0:
            # finds the nearest source in UVOTDETECT to default coordinates, and updates coordinate attributes
            source_ind = np.argmin(np.abs(data['RA']-self.source_coords.ra.value) + np.abs(data['DEC']-self.source_coords.dec.value))
            c_detect = SkyCoord(data[source_ind]['RA'],data[source_ind]['DEC'],unit='deg')

            # only accept the nearest source if it is within about 2 arcseconds of the default position and order creation of region files
            if c_detect.separation(self.source_coords).value < 6.E-4:
                self.source_coords = SkyCoord('%s %s' %(data[source_ind]['RA'],data[source_ind]['DEC']),unit=(u.deg, u.deg))
                self.createRegionFiles()
            else:
                print('Nearest source too far away (more than 2 arcsecs)! Defaulting to SIMBAD coords.')
                self.createRegionFiles()

        else:
            self.createRegionFiles()
            print('No sources were detected in %s! Defaulting to SIMBAD coords.' %self.filepath)

    def run_uvotdetect(self,exp='exp.fits'):
        '''Wrapper to run ``UVOTDETECT``. The ``filepath`` attribute is used for output filename generation. 

        args:
            exp (str): name of exposure file, assumed to be in the same directory as the primary image.
        '''
        dirpath,filename = path.split(self.filepath)
        

        #this is terrible - streamline it more
        expfile = ''
        if self.default_fs:
            base,extn = filename.split('_')
            obs = base[2:-3]
            band = base[-2:]
            expfile = path.join(dirpath,'%s_ex.img.gz' %base)
        else:
            obs = filename[8:23]
            band = filename[5:7]

        # optional usage of exposure files
        if not path.isfile(expfile) or not exp:
            expfile='NONE'

        outfile = path.join(dirpath,'detect_%s_%s.fits' %(obs,band))
        regfile = path.join(dirpath,'detect_%s_%s.reg' %(obs,band))
        bkgregfile = path.join(dirpath,'back_%s_%s.reg' %(obs,band))

        #using subprocess to run UVOTDETECT
        tmp = subprocess.Popen(['uvotdetect','infile=%s' %self.filepath,'outfile=%s' %outfile,'chatter=0','plotsrc=NO',
                                'expfile=%s'%expfile,'threshold=3','clobber=YES'], stdout=subprocess.PIPE)

        # sets filename attributes for region files and UVOTDETECT output filename
        self.detect = outfile
        self.regfile = regfile
        self.bkgregfile = bkgregfile

        return tmp.communicate()
