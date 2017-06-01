import os.path
import sys
import glob
import argparse 


def correct_extinction(val, filtr, EBminV, mag = False):

    '''assumes val expects flux densities in erg/s/cm2/Angstrom units, unless mag is True
    
    returns:
        erg/s/cm2
    '''

    #from Poole et al. 2008 (MNRAS 383, 627-645)
    central_wav = {'uu':3465.,'w1':2600.,'m2':2246.,'w2':1928.,'bb':4392.,'vv':5468.}
    # derived from http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/YorkExtinctionSolver/coefficients.cgi
    R_lambda = {'uu':4.89172,'w1':6.55663,'m2':9.15389,'w2':8.10997,'bb':4.00555,'vv':2.99692}
    
    A_lambda = R_lambda[filtr]*EBminV

    if mag:
        return val - R_lambda[filtr]*EBminV
    else:
        return val*central_wav[filtr]*10**(R_lambda[filtr]*EBminV/2.5)

def uvot_detecter(filepaths):
    from source_position import PositionExtractor
    pos = PositionExtractor()
    #hardcoded
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
        #hardcoded
        iv.source_ra = '8:54:48.867'
        iv.source_dec = '+20:06:30.97'
        iv.bkg_ra = '8:54:48.772'
        iv.bkg_dec = '+20:05:32.576'
        iv.setup_frame(d)
        x = raw_input('Viewing %s. Hit Enter to continue.' %filepath)


def uvot_measurer(filepaths,measure=True):
    #can try replacing this whole thing with a function on a single file and map iterator

    from astropy.table import Table
    from astropy import units as u
    from uvot_photometry import MeasureSource


    ptab = Table(names=('filter','MJD','Mag','MagErr','FluxDensity','FluxDensityErr','FluxDensityJy','FluxDensityJyErr','FluxExtCorr','FluxExtcorrErr'),
                            dtype=('S2','f8','f8','f8','f8','f8','f8','f8','f8','f8'))

    ptab['MJD'].unit = u.d
    ptab['Mag'].unit = u.mag
    ptab['MagErr'].unit = u.mag
    ptab['FluxDensity'].unit = u.erg/u.cm/u.cm/u.second/u.AA
    ptab['FluxDensityErr'].unit = u.erg/u.cm/u.cm/u.second/u.AA
    ptab['FluxDensityJy'].unit = u.Jy
    ptab['FluxDensityJyErr'].unit = u.Jy
    ptab['Flux'].unit = u.erg/u.cm/u.cm/u.second 
    ptab['FluxErr'].unit = u.erg/u.cm/u.cm/u.second

    #hardcoded
    ebminv = 0.077/3.1

    for filepath in filepaths:
        measurer = MeasureSource(filepath)
        filtr = measurer.band
        if measure:
            tmp = measurer.run_uvotsource()
         
        objphot = measurer.get_observation_data()

        #getting extinction-corrected fluxes
        flux = correct_extinction(objphot[-4],filtr,ebminv)
        fluxerr = correct_extinction(objphot[-3],filtr,ebminv)
        objphot.append(flux)
        objphot.append(fluxerr)

        ptab.add_row(objphot)

    #flux_extcorr = correct_exctinction(ptab['FluxDensity'],ptab['filter'],ebminv)
    #fluxerr_extcorr = correct_exctinction(ptab['FluxDensityErr'],ptab['filter'],ebminv)
    

    return ptab.group_by('filter')
        

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
        photometry.write('OJ287_uvot_photometry_wExtCorr.dat',format='ascii.commented_header')
        photometry.write(args.o)

if __name__ == '__main__':
    main()

