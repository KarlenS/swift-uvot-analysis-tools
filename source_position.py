import os
import os.path as path
import numpy as np
import glob
import subprocess
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord


class PositionExtractor(object):

    def __init__(self):
        #hardcoded
        self.source_ra = '8:54:48.867'
        self.source_dec = '+20:06:30.97'
	    self.bkg_ra =  '8:54:48.772'
	    self.bkg_dec = '+20:05:32.576'
        self.filepath = ''
        self.detect = ''
	    self.regfile = ''
	    self.bkgregfile = ''

    def cleanup(self):
        try:
            os.remove(self.detect)
        except OSError:
            print 'File %s does not exist. Cannot delete what does not exist.' %self.detect
            pass

    def createRegionFiles(self,c):
        dirpath,filename = path.split(self.filepath)
        region = 'fk5;circle(%s,%s,5")' %(c.ra.value,c.dec.value)
        regfile = open(self.regfile,'w')
        regfile.write(region)
        regfile.close()

        if not path.isfile(self.bkgregfile):
                bkgregion = 'fk5;circle(%s,%s,20")' %(self.bkg_ra,self.bkg_dec)
                bkgregfile = open(self.bkgregfile,'w')
                bkgregfile.write(bkgregion)
                bkgregfile.close()
        else:
            print 'Background region file already exists. Not recreating.'

  
    def run_uvotimsum(self,inputFile, outputFile):
    	tmp = subprocess.Popen(["uvotimsum",inputFile,outputFile,'chatter=0','clobber=no','exclude=NONE'], stdout=subprocess.PIPE)
  	    tmp = tmp.communicate()

    def getNearestSource(self):
        c = SkyCoord('%s %s' %(self.source_ra,self.source_dec),unit=(u.hourangle, u.deg))
        try:
            data = fits.getdata(self.detect)
        except:
            print 'File %s does not exist. %s is likely a multiple extension file. Attempting to combine extensions and rerunning.' %(self.detect,self.filepath)
	    try:
		#to keep things relatively clean, changing filepath to combined while remembering orignal name
	        comb = '%s.comb' %self.filepath
		    orig = self.filepath 
	        self.run_uvotimsum(self.filepath,comb)
		    self.filepath = comb
	        detect_out = self.run_uvotdetect(exp=None)
	        data = fits.getdata(self.detect)
		    #once we have the data from combined, changing back filepath to original
		    self.filepath = orig
	        print 'Combining and rerunning seems to have worked!'
	    except:
            print 'Attempt to combine extensions failed. Defaulting to SIMBAD position.'
            self.createRegionFiles(c)
            return None

        if np.size(data) > 0:
            source_ind = np.argmin(np.abs(data['RA']-c.ra.value) + np.abs(data['DEC']-c.dec.value))
            c_detect = SkyCoord(data[source_ind]['RA'],data[source_ind]['DEC'],unit='deg')

            if c_detect.separation(c).value < 6.E-4:
                self.createRegionFiles(c_detect)
            else:
                print 'Nearest source too far away (more than 2 arcsecs)! Defaulting to SIMBAD coords.'
                self.createRegionFiles(c)

            return c_detect
        else:
            self.createRegionFiles(c)
            print 'No sources were detected in %s! Defaulting to SIMBAD coords.' %self.filepath
            return None

    def run_uvotdetect(self,exp='exp.fits'):
        dirpath,filename = path.split(self.filepath)
        base,extn = filename.split('_')
        obs = base[2:-3]
        band = base[-2:]

        expfile = path.join(dirpath,'%s_ex.img.gz' %base)
        if not path.isfile(expfile) or not exp:
            expfile='NONE'

        outfile = path.join(dirpath,'detect_%s_%s.fits' %(obs,band))
        regfile = path.join(dirpath,'detect_%s_%s.reg' %(obs,band))
	    bkgregfile = path.join(dirpath,'back_%s_%s.reg' %(obs,band))

        tmp = subprocess.Popen(['uvotdetect','infile=%s' %self.filepath,'outfile=%s' %outfile,'chatter=0','plotsrc=NO',
                                'expfile=%s'%expfile,'threshold=3','clobber=YES'], stdout=subprocess.PIPE)
        self.detect = outfile
        self.regfile = regfile
        self.bkgregfile = bkgregfile

        return tmp.communicate()
