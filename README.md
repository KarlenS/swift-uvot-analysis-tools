# Swift UVOT analysis wrappers

Python helper wrappers for HEASoft UVOT analysis tools.

```
- Runs uvotdetect and generates region files for all observations
- Rudimentary data quality check by displaying all images in ds9
- Runs uvotsource on all observations
 -- currently set up for producing observation-by-observation light curve data (if multiple exposures exist per observation, they are summed)
```
## Required Software

[XPA](http://ds9.si.edu/site/XPA.html) , [DS9](http://ds9.si.edu/site/Download.html), [HEASoft (at least for Swift)](https://heasarc.gsfc.nasa.gov/lheasoft/download.html), [CALDB](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_install.html)[(with data for at least for Swift)](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/swift/)

Python packages:
```bash
$> pip -r install requirements.txt
```

## Generating photometry and light curve data

To see all available options:
```bash
$> python run_uvot_analysis.py --help
```

### Step1: Run uvotdetect and generate region files for all observations

```bash
$> python run_uvot_analysis.py -p /path/to/Swift/data --detect
```
The directory */path/to/Swift/data/* should contain UVOT data strucutres from all observations (typically directories titles e.g., 00034934001)

To run with a single observation:

```bash
$> python run_uvot_analysis.py -p /path/to/Swift/data -obs 00034934001 --detect
```

### Step2: Visually check observation from each image (this will require DS9 and XPA)
```bash
$> python run_uvot_analysis.py -p /path/to/Swift/data --check
```
To view images for a single observations, run:
```bash
$> python run_uvot_analysis.py -p /path/to/Swift/data -obs 00034934001 --check
```

### Step3: Run uvotsource and extract photometry
The following will run uvotsource on all observations and store photometry data in *photometry.fits* (in fits format)
```bash
$> python run_uvot_analysis.py -p /path/to/Swift/data --measure
```

To only parse photometry data and store it in *MySource_photometry.cvs* (in csv format; same formats as supported by [astropy.table](http://docs.astropy.org/en/stable/table/io.html)).

***NOTE*: Do this if *uvotsource* has already been run with --measure**
```bash
$> python run_uvot_analysis.py -p /path/to/Swift/data --measure -o MySource_phometry.csv --print_only
```

As in previous steps, this works for a single observation as well.
```bash
$> python run_uvot_analysis.py -p /path/to/Swift/data --measure -obs 00034934001 -o MySource_photometry.fits
```
