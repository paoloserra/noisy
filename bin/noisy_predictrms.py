#!/usr/bin/env python

# Import modules
import noisy as noisypkg # critical imports and settings, should be first

import numpy as np
import pyrap.tables as tables
import sys,os
from noisy.helpers import noisy

###############################
### READ COMMAND LINE INPUT ###
###############################
# Print help message and exit
print 'Noisy: Predict Natural Noise'
print '(C) Paolo Siera, 2017'
print ''
print 'Version: %s' % noisypkg.__version__
print 'Package installed at: %s' % noisypkg.__path__[0]
print ''

if 'help' in sys.argv or '-h' in sys.argv or '--help' in sys.argv:
    print 'Usage:'
    print 'predict_natural_rms.py <file1.ms> [file2.ms [file3.ms ... [fileN.ms] ... ]] [-tsyseff <Tsys/eff (K) OR file>] '
    print '       [-diam <antenna diameter (m)>] [-field <field name>] [-plot <plot name with extension>]'
    print ''
    print ' If you give Tsys/eff as a file it should have two columns: frequency (Hz) and Tsys/eff (K)'
    sys.exit()
else: arg=sys.argv

# Get telescope parameters (default or from command line) and calculate derived quantities
# Tsys/eff in K (default: 22 K)
if '-tsyseff' in arg:
    tsyseff=arg[arg.index('-tsyseff')+1]
    del(arg[arg.index('-tsyseff'):arg.index('-tsyseff')+2])
else: tsyseff='22' # leave as string
# antenna diameter in m (default: 13.5 m)
if '-diam' in arg:
    diam=np.float(arg[arg.index('-diam')+1])
    del(arg[arg.index('-diam'):arg.index('-diam')+2])
else: diam=13.5
# field name (default: None = all fields)
if '-field' in arg:
    selectFieldName=arg[arg.index('-field')+1]
    del(arg[arg.index('-field'):arg.index('-field')+2])
else: selectFieldName=None
# make plots? (default: False)
if '-plot' in arg:
    plotName=arg[arg.index('-plot')+1]
    del(arg[arg.index('-plot'):arg.index('-plot')+2])
else: plotName=None

# Get input files from command line
MS=arg[1:]
checkfiles=[os.path.exists(jj) for jj in MS]
if not len(MS):
    print ' CATASTROPHE!'
    print ' No input .MS file provided'
    print ' Aborting ...'
    sys.exit()
elif False in checkfiles:
    print ' CATASTROPHE!'
    for jj in range(len(MS)):
        if not checkfiles[jj]: print ' Input file {0:s} does not exist'.format(MS[jj])
    print ' Aborting ...'
    sys.exit()
else:
    print ''
    print '--- Will work on the following .MS files ---'
    for jj in range(len(MS)): print '',MS[jj]



#######################
### GET NATURAL RMS ###
#######################

noisy.PredictNoise(MS,tsyseff,diam,plotName,selectFieldName)
