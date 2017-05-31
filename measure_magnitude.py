import os.path as path
import subprocess

class MeasureSource(object):

    def __init__(self,filepath=''):
        self.source_ra = None
        self.source_dec = None
        self.bkg_ra = None # 8:54:48.772
        self.bkg_dec = None # +20:05:32.576
        self.filepath = filepath


    def run_uvotsource(self):

        dirpath,filename = path.split(self.filepath)
        base,extn = filename.split('_')
        obs = base[2:-3]
        band = base[-2:]
	
	srcregfile = path.join(dirpath,'detect_%s_%s.reg' %(obs,band))
	self.makeBackgroundRegionFile(bkgregfile)
        
        tmp = subprocess.Popen(['uvotsource','infile=%s' %self.filepath,'srcreg=%s' %srcregfile, 'bkgreg=%s' %bkgregfile,
				'chatter=0','threshold=3','clobber=YES'], stdout=subprocess.PIPE)
	tmp = tmp.communicate()

