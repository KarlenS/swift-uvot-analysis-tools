#!/local/gammasoft/anaconda2/bin/python
from astropy.io import fits
import matplotlib.pyplot as plt
import argparse
import numpy as np

def readData(filename):
    return fits.getdata(filename)

def plotSED(dat,axs,color='black',label=None):


    central_wav = {'uu':3465.,'w1':2600.,'m2':2246.,'w2':1928.,'bb':4392.,'vv':5468.}

    c = 2.9979E8*1E10 #speed of light in angstroms
    ebminv = 0.085
    R_lambda = {'uu':4.89172,'w1':6.55663,'m2':9.15389,'w2':8.10997,'bb':4.00555,'vv':2.99692}

    (ax, ax_corr) = axs

    ax.errorbar(None,None,fmt='o',color=color,label=label,alpha=0.7)
    ax_corr.errorbar(None,None,fmt='o',color=color,label=label,alpha=0.7)

    for keys,vals in central_wav.items():
        ind = dat['filter'] == keys
        flux = dat[ind]['FluxDensity']
        ferr = dat[ind]['FluxDensityErr']
        fluxcorr = flux*vals*10**(0.4*R_lambda[keys]*ebminv)
        ferrcorr = ferr*vals*10**(0.4*R_lambda[keys]*ebminv)
        if True in ind:
            ax.errorbar(c/vals,flux*vals,yerr=ferr, fmt='o',color=color,alpha=0.7)
            #print '%.5e %.5e %.5e' %(c/vals, flux[0]*vals, ferr[0]*vals, fluxcorr[0], ferrcorr[0])
            ax_corr.errorbar(c/vals,fluxcorr,yerr=ferrcorr, fmt='o',color=color,alpha=0.7)

    ax.set_ylabel(r'Flux [ erg cm$^{-2}$ s$^{-1}$ ]')
    ax_corr.set_ylabel(r'Flux (ExtCorr) [ erg cm$^{-2}$ s$^{-1}$ ]')
    ax_corr.set_xlabel(r'Frequency [ Hz ]')

    #ax.set_ylim([0.8E-11,7.2E-11])
    #ax_corr.set_ylim([1E-11,9E-11])


    ax.legend(ncol=4)
    ax_corr.legend(ncol=4)

def main():

    parser = argparse.ArgumentParser(description='Quick plotter for UVOT SEDs.')
    parser.add_argument('-f', default=None, help='Name of fits file to plot.')
    parser.add_argument('-l', default=None, help='List of fits files to plot.')
    args = parser.parse_args()

    if not args.f and not args.l:
        raise parser.error('Either provide a fits file or a file with a list of fits files to plot')

    if args.f:
        dat = readData(args.f)
        fig, axs = plt.subplots(2,sharex=True,figsize=(10, 6))
        plotSED(dat,axs)
    
    if args.l:
        files = np.genfromtxt(args.l,dtype=str)
        colors = plt.cm.rainbow(np.linspace(0,1,np.size(files)))
        fig, axs = plt.subplots(2,sharex=True,figsize=(10, 6))
        for f,c in zip(files,colors):
            dat = readData(f)
            plotSED(dat,axs,color=c,label=f[6:21])

    #plt.savefig('3C66A_seds.png',bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main()
