import glob
import argparse
import pyds9 as ds9

class SourceImageViewer(object):

    def __init__(self,filepath=''):
        self.filepath = filepath
        self.source_ra = '8:54:48.867'
        self.source_dec = '+20:06:30.97'

    def open_fits(self,d):
        opencmd = "file %s" %self.filepath
        d.set(opencmd)

    def load_source_region(self,d,radius = 5):
        regcmd = 'fk5; circle(%s,%s,%s")' %(self.source_ra,self.source_dec,radius)
        d.set('regions', regcmd);

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

def main():

    parser = argparse.ArgumentParser(description='Quick verification of image quality and regions selection for UVOT analysis.')
    parser.add_argument('-p',required=True, help='Parent directory of UVOT data structures (directories with names like 00030901252).')
    parser.add_argument('-f', help='File with list of directory names of UVOT data structures. Use this option if directory names do not start with 000 for some reason...')
    parser.add_argument('--pause_at_obs',action='store_true',default=False)
    args = parser.parse_args()

    filepaths = glob.glob(args.p+'/000*/uvot/image/*sk.img.gz')#setting up filepaths

    #launching a ds9 window
    d = ds9.DS9()

    iv = SourceImageViewer()
    iv.source_ra = '8:54:48.867'
    iv.source_dec = '+20:06:30.97'

    for filepath in filepaths:
        iv.filepath = filepath
        iv.setup_frame(d)
        x = raw_input('')

if __name__ == '__main__':
    main()
