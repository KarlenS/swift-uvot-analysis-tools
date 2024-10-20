#!/Users/karlen/anaconda2/envs/astroconda/bin/python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals
from future import standard_library
standard_library.install_aliases()
from builtins import zip
from builtins import *
from past.utils import old_div
from astropy.io import fits
import matplotlib.pyplot as plt
import argparse
import numpy as np

def readData(filename):
    return fits.getdata(filename)

def plotSED(dat,axs,color='black',label=None):


    central_wav = {'uu':3465.,'w1':2600.,'m2':2246.,'w2':1928.,'bb':4392.,'vv':5468.}

    c = 2.9979E8*1E10 #speed of light in angstroms

    axs.errorbar(None,None,fmt='o',color=color,label=label,alpha=0.7)

    for keys,vals in list(central_wav.items()):
        ind = dat['filter'] == keys
        flux = dat[ind]['FluxDensity']
        ferr = dat[ind]['FluxDensityErr']
        fluxcorr = dat[ind]['FluxExtCorr']
        ferrcorr = dat[ind]['FluxExtCorrErr']
        if True in ind:
            axs.errorbar(old_div(c,vals),fluxcorr,yerr=ferrcorr, fmt='o',color=color,alpha=0.7)

    axs.set_ylabel(r'Flux (ExtCorr) [ erg cm$^{-2}$ s$^{-1}$ ]')
    axs.set_xlabel(r'Frequency [ Hz ]')

    axs.legend(ncol=4)

def main():

    parser = argparse.ArgumentParser(description='Quick plotter for UVOT SEDs.')
    parser.add_argument('-f', default=None, help='Name of fits file to plot.')
    parser.add_argument('-l', default=None, help='List of fits files to plot.')
    args = parser.parse_args()

    if not args.f and not args.l:
        raise parser.error('Either provide a fits file or a file with a list of fits files to plot')

    if args.f:
        dat = readData(args.f)
        fig, axs = plt.subplots(1,sharex=True,figsize=(10, 6))
        plotSED(dat,axs)

    if args.l:
        files = np.genfromtxt(args.l,dtype=str)
        colors = plt.cm.rainbow(np.linspace(0,1,np.size(files)))
        fig, axs = plt.subplots(1,sharex=True,figsize=(10, 6))
        for f,c in zip(files,colors):
            dat = readData(f)
            plotSED(dat,axs,color=c,label=f[6:21])

    plt.show()

if __name__ == '__main__':
    main()
