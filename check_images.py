import numpy as np
import os.path as path
import pyregion

class SourceImageViewer(object):

    def __init__(self,filepath=''):
        self.filepath = filepath
        self.source_ra = None 
        self.source_dec = None
        self.bkg_ra = None
        self.bkg_dec = None

    def open_fits(self,d):
        opencmd = "file %s" %self.filepath
        d.set(opencmd)

    def display_source_region(self,d,radius = 5):

        reg = 'fk5; circle(%s,%s,%s")' %(self.source_ra,self.source_dec,radius)
        d.set('regions', reg);

    def display_bkg_region(self,d,radius = 20):

        bkgreg = 'fk5; circle(%s,%s,%s")' %(self.bkg_ra,self.bkg_dec,radius)
        d.set('regions', bkgreg)

    def remove_regions(self,d):
        d.set('regions delete all')

    def load_regions(self,d):
        
        dirpath,filename = path.split(self.filepath)
        base,extn = filename.split('_')
        obs = base[2:-3]
        band = base[-2:]

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
        
        self.display_bkg_region(d)

        try:
            region = pyregion.open(regfile)
            ra = region[0].coord_list[0]
            dec = region[0].coord_list[1]
            self.source_ra = ra
            self.source_dec = dec
        except IOError:
            print 'Region file missing. Only showing default coordinates, if any.'

        self.display_source_region(d)

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

    def format_frame(self,d):
        self.zoom_in(d)
        self.prettify(d)


    def setup_frame(self,d):

        self.open_fits(d)
        self.format_frame(d)
        self.load_regions(d)
        self.center_on_source(d)

    def get_user_coords(self,d):
        return d.get('iexam coordinate wcs fk5')

    def prime_bkg(self,d):

        self.open_fits(d)
        self.format_frame(d)
        self.center_on_source(d)

        while True:

            self.display_source_region(d)
            print 'Use the DS9 window to select a background region.'
            coords = self.get_user_coords(d)
            self.bkg_ra, self.bkg_dec = coords.split()
            self.display_bkg_region(d)

            response = raw_input('Keep background selection (y/n)?\n') 
            if response == 'y':
                print 'Great! That is all we need. Bye.'
                break
            else:
                print 'Fine, try again...'
                self.remove_regions(d)

        d.set('quit')

    def prime_source(self,d):

        self.open_fits(d)
        self.format_frame(d)

        while True:

            print 'Use the DS9 window to select a source region.'
            coords = self.get_user_coords(d)
            self.source_ra, self.source_dec = coords.split()
            self.display_source_region(d)

            response = raw_input('Keep source selection (y/n)?\n') 
            if response == 'y':
                print 'Great! That is all we need. Bye.'
                break
            else:
                print 'Fine, try again...'
                self.remove_regions(d)

        d.set('quit')
