import os.path
import sys
import glob
import argparse 

def uvot_detecter(filepaths):
    from source_position import PositionExtractor
    pos = PositionExtractor()
    pos.source_ra = '8:54:48.867'
    pos.source_dec = '+20:06:30.97'
    for filepath in filepaths:
        print 'working on %s...' %filepath 
        pos.filepath = filepath
        detect_out = pos.run_uvotdetect()
        nearest = pos.getNearestSource()
        pos.cleanup()


def uvot_checker(filepaths):

    import pyds9 as ds9
    from check_images import SourceImageViewer
    #launching a ds9 window
    d = ds9.DS9()

    for filepath in filepaths:
        iv = SourceImageViewer()
        iv.filepath = filepath
        iv.source_ra = '8:54:48.867'
        iv.source_dec = '+20:06:30.97'
        iv.bkg_ra = '8:54:48.772'
        iv.bkg_dec = '+20:05:32.576'
        iv.setup_frame(d)
        x = raw_input('Viewing %s. Hit Enter to continue.' %filepath)


def uvot_measurer(filepaths,measure=True):
    #can try replacing this whole thing with a function on a single file and map iterator

    from astropy.table import Table
    from uvot_photometry import MeasureSource

    ptab = Table(names=('filter','MJD','Mag','MagErr'),dtype=('S2','f8','f8','f8'))

    for filepath in filepaths:
        measurer = MeasureSource(filepath)
        filtr = measurer.band
        if measure:
            tmp = measurer.run_uvotsource()
         
        objphot = measurer.get_observation_data()
        ptab.add_row(objphot)


    return ptab
        

def main():

    parser = argparse.ArgumentParser(description='Quick verification of image quality and regions selection for UVOT analysis.')
    parser.add_argument('-p',required=True, help='Parent directory of UVOT data structures (directories with names like 00030901252).')
    parser.add_argument('-obs',default=None,help='If specified, will only look at data for a single observations: e.g. 00030901252')
    parser.add_argument('-f', help='File with list of directory names of UVOT data structures. Use this option if directory names do not start with 000 for some reason...')
    parser.add_argument('-o',default='photometry.fits',help='Name of output fits file for storing resulting photometry.')
    parser.add_argument('--check',action='store_true',default=False,help='View all images with the source marked to check observation quality.')
    parser.add_argument('--detect',action='store_true',default=False,help='Run uvotdetect on all images.')
    parser.add_argument('--measure',action='store_true',default=False,help='Run uvotsource on all images to get photometry.')
    parser.add_argument('--print_only',action='store_true',default=False,help='Only plot light curves.')
    args = parser.parse_args()

    if not args.obs:
        filepaths = glob.iglob('%s/000*/uvot/image/*sk.img.gz' %args.p)#setting up filepaths
    else:
        obspath = os.path.join(args.p,str(args.obs))
        if os.path.exists(obspath):
            filepaths = glob.iglob('%s/uvot/image/*sk.img.gz'%obspath)
        else:
            raise NameError('%s does not exist.' %obspath)


    if args.check:
        uvot_checker(filepaths)

    if args.detect:
        uvot_detecter(filepaths)

    if args.measure:
        photometry = uvot_measurer(filepaths,measure = not args.print_only)
        photometry.write(args.o)

if __name__ == '__main__':
    main()

