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
        self.source_ra = '8:54:48.867'
        self.source_dec = '+20:06:30.97'
        self.filepath = ''
        self.detect = ''

    def cleanup(self):
        try:
            os.remove(self.detect)
        except OSError:
            print 'File %s does not exist. Cannot delete what does not exist.' %self.detect
            pass

    def createRegionFile(self,c):
        dirpath,filename = path.split(self.filepath)
        region = 'fk5;circle(%s,%s,5")' %(c.ra.value,c.dec.value)
        regfile = open(self.regfile,'w')
        regfile.write(region)
        regfile.close()


    def getNearestSource(self):
        c = SkyCoord('%s %s' %(self.source_ra,self.source_dec),unit=(u.hourangle, u.deg))
        try:
            data = fits.getdata(self.detect)
        except:
            print 'File %s does not exist. %s is likely a multiple extension file. Defaulting to SIMBAD position.' %(self.detect,self.filepath)
            self.createRegionFile(c)
            return None

        if np.size(data) > 0:
            source_ind = np.argmin(np.abs(data['RA']-c.ra.value) + np.abs(data['DEC']-c.dec.value))
            c_detect = SkyCoord(data[source_ind]['RA'],data[source_ind]['DEC'],unit='deg')

            if c_detect.separation(c) < 3.E-4:
                self.createRegionFile(c_detect)
            else:
                print 'Nearest source too far away! Defaulting to SIMBAD coords.'
                self.createRegionFile(c)

            return [c_detect,c_detect.separation(c)]
        else:
            self.createRegionFile(c)
            print 'No sources were detected in %s! Defaulting to SIMBAD coords.' %self.filepath
            return None

    def run_uvotdetect(self):
        dirpath,filename = path.split(self.filepath)
        base,extn = filename.split('_')
        obs = base[2:-3]
        band = base[-2:]

        expfile = path.join(dirpath,'%s_ex.img.gz' %base)
        if not path.isfile(expfile):
            expfile='NONE'

        outfile = path.join(dirpath,'detect_%s_%s.fits' %(obs,band))
        regfile = path.join(dirpath,'detect_%s_%s.reg' %(obs,band))

        tmp = subprocess.Popen(['uvotdetect','infile=%s' %self.filepath,'outfile=%s' %outfile,'chatter=0','plotsrc=NO',
                                'expfile=%s'%expfile,'threshold=3','clobber=YES'], stdout=subprocess.PIPE)
        self.detect = outfile
        self.regfile = regfile

        return tmp.communicate()
