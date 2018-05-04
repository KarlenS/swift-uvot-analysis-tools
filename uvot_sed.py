#!/Users/karlen/anaconda2/envs/astroconda/bin/python
import os
import re
import subprocess
import gzip
import numpy as np
from shutil import copyfile,copyfileobj
from source_position import PositionExtractor
from uvot_photometry import MeasureSource
from astropy.time import Time
from tqdm import tqdm


class MakeSED(object):

    def __init__(self,filelist):
        self.filelist = filelist
        self.filters = ['vv','bb','uu','w1','w2','m2']
        self.sortedpaths = None

    def sortFilelistByFilter(self):
        
        fdict = dict.fromkeys(self.filters,np.array([],dtype=str))
        filelistarr = np.array(self.filelist)

        for filtr in fdict.keys():
            indx = np.where([re.search(filtr, path) for path in filelistarr])[0]
        
            if np.size(indx) > 0:
                fdict[filtr] = np.append(fdict[filtr],filelistarr[indx])

        return fdict 
        

    def combineFits(self,startdate, enddate):

        sumfile = 'summed_image.fits'
        sumfilegz = 'summed_image.gz'
        filepaths = []
        for filtr,paths in self.sortedpaths.items():
            firstfile = True
            extfile = '%s_%s-%s_%s' %(filtr,startdate.mjd,enddate.mjd,sumfile)
            combfile = 'comb_%s' % extfile
            nfiles = 0
            print 'Working on %s filter' % (filtr)

            for f in tqdm(paths):
                measurer = MeasureSource(f)
                obstime,aspflag = measurer.get_observation_time()
                if (obstime-startdate).value > 0 and (enddate-obstime).value > 0:
                    #print '%s is within the time range.' %f
                    nfiles += 1
                    if firstfile:
                        copyfile(f,'%s_%s-%s_%s' %(filtr,startdate.mjd,enddate.mjd,sumfilegz))
                        with gzip.open('%s_%s-%s_%s' %(filtr,startdate.mjd,enddate.mjd,sumfilegz),'rb') as f_in:
                            with open(extfile,'wb') as f_out:
                                copyfileobj(f_in,f_out)
                        os.remove('%s_%s-%s_%s' %(filtr,startdate.mjd,enddate.mjd,sumfilegz))

                        firstfile = False
                    elif aspflag:
                        self.runFappend(f,extfile)
                    else:
                        print 'FUUUUCK THIS FILE: %s' %f
                        continue

            if nfiles == 0:
                print 'Filter %s had no files to combine' %filtr
                continue
            else:

                sp = PositionExtractor() #make a utils file so don't have to use a whole class for a simple function
                sp.run_uvotimsum(extfile,combfile)
                filepaths.append(combfile)
        
        return filepaths



    def runFappend(self,inputFile,mergedFile):
        #print 'appending %s to %s' %(inputFile,mergedFile)
        tmp = subprocess.Popen(["fappend",inputFile,mergedFile], stdout=subprocess.PIPE)
        tmp = tmp.communicate()
