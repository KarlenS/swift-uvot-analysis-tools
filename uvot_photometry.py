import os.path as path
import subprocess

from astropy.io import fits

class MeasureSource(object):

    def __init__(self,filepath=''):
        self.filepath = filepath
        
        dirpath,filename = path.split(self.filepath)
        base,extn = filename.split('_')
        self.dirpath = dirpath
        self.obs = base[2:-3]
        self.band = base[-2:]

    def run_uvotsource(self):

	srcregfile = path.join(self.dirpath,'detect_%s_%s.reg' %(self.obs,self.band))
	bkgregfile = path.join(self.dirpath,'back_%s_%s.reg' %(obs,band))
        uvotsourcefile = path.join(self.dirpath,'uvotsource_%s_%s.reg' %(self.obs,self.band))
        
        tmp = subprocess.Popen(['uvotsource','image=%s' %self.filepath,'srcreg=%s' %srcregfile, 'bkgreg=%s' %bkgregfile,
				'outfile=%s' %uvotsourcefile,'chatter=0','sigma=3','clobber=YES'], stdout=subprocess.PIPE)

        return tmp.communicate()

    def get_observation_data(self):

        uvotsourcefile = path.join(self.dirpath,'uvotsource_%s_%s.reg' %(self.obs,self.band))
        try:
            data = fits.getdata(uvotsourcefile)
        except IOError:
            print '%s not found. Make sure to run uvotsource beforehand or this will get annoying.\n
                   Will try running uvotsource now and will skip this observation if it fails.' %uvosourcefile
            tmp = self.run_uvotsource()
            try:
                data = fits.getdata(uvotsourcefile)
                print 'Running uvotsource worked.'
            except IOError:
                print 'Failed again... Skipping %s' %self.filepath
                return None


