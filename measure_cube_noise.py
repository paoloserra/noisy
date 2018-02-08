#!/usr/bin/env python

import sys,os
import noisy



###############################
### READ COMMAND LINE INPUT ###
###############################

# Print help message and exit
if 'help' in sys.argv or '-h' in sys.argv:
    print ' Run as follows'
    print 'python measure_cube_noise.py <cube1.fits> [cube2.fits [cube3.fits ... [cubeN.ms] ... ]]'
    print '       [-plot <plot name with extension>] [-title <plot title>] [-random]'
    sys.exit()
else: arg=sys.argv

# make plots (default: False)
if '-plot' in arg:
    plotName=arg[arg.index('-plot')+1]
    del(arg[arg.index('-plot'):arg.index('-plot')+2])
else: plotName=None

# plots title (default: None)
if '-title' in arg:
    title=arg[arg.index('-title')+1]
    del(arg[arg.index('-title'):arg.index('-title')+2])
else: title=None

# randomize channel order
if '-random' in arg:
    randomizeChanOrder=True
    del(arg[arg.index('-random')])
else: randomizeChanOrder=False

# Get input files from command line
FITS=arg[1:]
checkfiles=[os.path.exists(jj) for jj in FITS]
if not len(FITS):
    print ' CATASTROPHE!'
    print ' No input .FITS file provided'
    print ' Aborting ...'
    sys.exit()
elif False in checkfiles:
    print ' CATASTROPHE!'
    for jj in range(len(FITS)):
        if not checkfiles[jj]: print ' Input file {0:s} does not exist'.format(FITS[jj])
    print ' Aborting ...'
    sys.exit()
else:
    print ''
    print '--- Will work on the following .FITS files ---'
    for jj in range(len(FITS)): print '',FITS[jj]



#######################
### GET NATURAL RMS ###
#######################

noisy.MeasureCubeNoise(FITS,plotName,randomizeChanOrder,title=title)
