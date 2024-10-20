#!/Users/karlen/anaconda2/envs/astroconda/bin/python

'''
Quickly plots up the UVOT light curve in each available band
'''
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from future import standard_library
standard_library.install_aliases()
from builtins import zip
from builtins import *
from builtins import object
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt

import seaborn as sns

from astropy.io import fits


class UVOTLightCurve(object):

    def __init__(self,filename):
        self.filename = filename

    def readFile(self):

        try:
            return fits.getdata(self.filename)
        except:
            raise IOError('Only doing this with fits files. Dump your photometry into a fits file and try again')

    def plotLightCurves(self, dat, multipanel=True):

        filters = np.unique(dat['filter'])

        if not multipanel:
            fig, ax = plt.subplots(1)
            plt.ylim([0.1E-14,2E-14])

            for f in filters:
                ax.errorbar(dat['mjd'][dat['filter']==f],1E14*dat['FluxDensity'][dat['filter']==f],yerr=1E14*dat['FluxDensityErr'][dat['filter']==f],label='f',fmt='.')
            ax.legend()

        else:
            fig, axes = plt.subplots(np.size(filters),sharex=True,sharey=True)
            plt.subplots_adjust(hspace=0)
            plt.ylim([0.1,1.8])
            cols = cycle('bgrcmk')
            for ax,f in zip(axes,filters):
                #ax.errorbar(dat['mjd'][dat['filter']==f],dat['Mag'][dat['filter']==f],yerr=dat['MagErr'][dat['filter']==f],fmt='.')
                ax.errorbar(dat['mjd'][dat['filter']==f],1E14*dat['FluxDensity'][dat['filter']==f],yerr=1E14*dat['FluxDensityErr'][dat['filter']==f],label=f,color=next(cols),fmt='.')
                ax.legend()
                
        fig.text(0.04, 0.5,r'Flux Density [ 10$^{-14}$ erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ]',va='center',rotation='vertical')
        plt.xlabel('MJD')
        #plt.savefig("lightcurve.png")
        plt.show() 

def main():

    sns.set_style('ticks')

    filename = 'photometry.fits'
    lc = UVOTLightCurve(filename=filename)
    lcdata = lc.readFile()
    lc.plotLightCurves(lcdata)

if __name__ == '__main__':
    main()
