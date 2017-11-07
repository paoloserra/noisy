#!/usr/bin/env python

# Import modules
import numpy as np
import pyrap.tables as tables
import sys,os
import noisy

# To be done
# - frequency dependent Tsys
# - MS files with 4 polarisation products (currently uses all available pols)


###############################
### READ COMMAND LINE INPUT ###
###############################

# Print help message and exit
if 'help' in sys.argv or '-h' in sys.argv:
    print ' Run as follows'
    print 'python predict_natural_rms.py <file1.ms> [file2.ms [file3.ms ... [fileN.ms] ... ]] [-tsys <tsys/K>] '
    print '       [-eff <antenna efficiency>] [-diam <antenna diameter/m>] [-field <field name>] [-plot]'
    sys.exit()
else: arg=sys.argv

# Get telescope parameters (default or from command line) and calculate derived quantities
# tsys in K (default: 22 K)
if '-tsys' in arg:
    tsys=np.float(arg[arg.index('-tsys')+1])
    del(arg[arg.index('-tsys'):arg.index('-tsys')+2])
else: tsys=22
# antenna efficiency (default: 1)
if '-eff' in arg:
    eff=np.float(arg[arg.index('-eff')+1])
    del(arg[arg.index('-eff'):arg.index('-eff')+2])
else: eff=1.
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
    makePlot=True
    del(arg[arg.index('-plot')])
else: makePlot=False

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

noisy.combineMS(MS,tsys,eff,diam,makePlot,selectFieldName)
