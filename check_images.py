import numpy as np
import os.path as path
import pyregion

class SourceImageViewer(object):

    def __init__(self,filepath=''):
        self.filepath = filepath
        self.source_ra = None
        self.source_dec = None

    def open_fits(self,d):
        opencmd = "file %s" %self.filepath
        d.set(opencmd)

    def load_source_region(self,d,radius = 5):
        
        dirpath,filename = path.split(self.filepath)
        base,extn = filename.split('_')
        obs = base[2:-3]
        band = base[-2:]

        reg_default = 'fk5; circle(%s,%s,%s")' %(self.source_ra,self.source_dec,radius)
        d.set('regions', reg_default);

        regfile = path.join(dirpath,'detect_%s_%s.reg' %(obs,band))

        try:
            region = pyregion.open(regfile)
            ra = region[0].coord_list[0]
            dec = region[0].coord_list[1]
            self.source_ra = ra
            self.source_dec = dec
        except IOError:
            print 'Region file missing. Only showing default coordinates'


        reg = 'fk5; circle(%s,%s,%s")' %(self.source_ra,self.source_dec,radius)
        d.set('regions', reg);

    def center_on_source(self,d):
        centercmd = 'pan to %s %s wcs fk5' %(self.source_ra,self.source_dec)
        d.set(centercmd)

    def zoom_in(self,d,zoom = 4):
        d.set('zoom to fit')
        d.set('zoom %s' %zoom)

    def setup_frame(self,d):
        self.open_fits(d)
        self.load_source_region(d)
        self.zoom_in(d)
        self.center_on_source(d)
