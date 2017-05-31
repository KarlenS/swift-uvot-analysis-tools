import numpy as np
import os.path as path
import pyregion

class SourceImageViewer(object):

    def __init__(self,filepath=''):
        self.filepath = filepath
        self.source_ra = 0 
        self.source_dec = 0
        self.bkg_ra = 0 
        self.bkg_dec = 0 

    def open_fits(self,d):
        opencmd = "file %s" %self.filepath
        d.set(opencmd)

    def load_source_region(self,d,radius = 5):
        
        dirpath,filename = path.split(self.filepath)
        base,extn = filename.split('_')
        obs = base[2:-3]
        band = base[-2:]

        #reg_default = 'fk5; circle(%s,%s,%s")' %(self.source_ra,self.source_dec,radius)
        #d.set('regions', reg_default);

        regfile = path.join(dirpath,'detect_%s_%s.reg' %(obs,band))
        bkgregfile = path.join(dirpath,'back_%s_%s.reg' %(obs,band))

        try:
            bkgregion = pyregion.open(bkgregfile)
            bkg_ra = bkgregion[0].coord_list[0]
            bkg_dec = bkgregion[0].coord_list[1]
            self.bkg_ra = bkg_ra
            self.bkg_dec = bkg_dec
        except IOError:
            print 'Background region file missing. Will try default coordinates.'
        
        bkgreg = 'fk5; circle(%s,%s,%s")' %(self.bkg_ra,self.bkg_dec,20)
        d.set('regions', bkgreg)


        try:
            region = pyregion.open(regfile)
            ra = region[0].coord_list[0]
            dec = region[0].coord_list[1]
            self.source_ra = ra
            self.source_dec = dec
        except IOError:
            print 'Region file missing. Only showing default coordinates, if any.'


        reg = 'fk5; circle(%s,%s,%s")' %(self.source_ra,self.source_dec,radius)
        d.set('regions', reg);

    def center_on_source(self,d):
        centercmd = 'pan to %s %s wcs fk5' %(self.source_ra,self.source_dec)
        d.set(centercmd)

    def zoom_in(self,d,zoom = 5):
        d.set('zoom to fit')
        d.set('zoom %s' %zoom)

    def prettify(self,d):
        d.set('scale mode minmax')
        d.set('scale log')
        d.set('cmap heat')

    def setup_frame(self,d):
        self.open_fits(d)
        self.load_source_region(d)
        self.zoom_in(d)
        self.center_on_source(d)
        self.prettify(d)