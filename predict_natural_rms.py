#!/usr/bin/env python

# Import modules
import numpy as np
import pyrap.tables as tables
import sys,os

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




########################
### DEFINE FUNCTIONS ###
########################

# Get single-MS flags, intervals, channel widths, channel frequencies and calculate natural rms (ignoring flags)
def processMS(ms,kB,tsys,eff,Aant,selectFieldName):
    print ''
    print '--- Working on file {0:s} ---'.format(ms)
    t=tables.table(ms)
    fieldIDs=t.getcol('FIELD_ID')
    ant1=t.getcol('ANTENNA1')
    ant2=t.getcol('ANTENNA2')
    fieldNames=tables.table(ms+'/FIELD').getcol('NAME')
    spw=tables.table(ms+'/SPECTRAL_WINDOW')
    channelWidths=spw.getcol('CHAN_WIDTH')
    channelFreqs=spw.getcol('CHAN_FREQ')

    if selectFieldName:
        try:
            selectFieldID=fieldNames.index(selectFieldName)
        except ValueError:
            print ' CATASTROPHE!'
            print ' Cannot find the field you want to process, {0:s}'.format(selectFieldName)
            print ' Available fields are',fieldNames
            print ' Aborting ...'
            sys.exit()
        print 'Successfully selected Field with name {0:s} (Field ID = {1:d})'.format(selectFieldName,selectFieldID)
        selection=fieldIDs==selectFieldID
    else:
        print 'Will process all available fields:',fieldNames
        selection=fieldIDs>=fieldIDs.min()

    autoCorr=ant1==ant2
    if autoCorr.sum(): print 'Successfully selected crosscorrelations only'
    else: print 'Found crosscorrelations only'
    selection*=ant1!=ant2
    nrAnt=np.unique(np.concatenate((ant1,ant2))).shape[0]
    nrBaseline=nrAnt*(nrAnt-1)/2
    print 'Number of antennas  = {0:d}'.format(nrAnt)
    print 'Number of baselines = {0:d}'.format(nrBaseline)
    print 'Frequency coverage  = {0:.5e} Hz - {1:.5e} Hz'.format(channelFreqs.min(),channelFreqs.max())
    if np.unique(channelWidths).shape[0]==1: print 'Channel width = {0:.5e} Hz'.format(np.unique(channelWidths)[0])
    else: print 'The channel width takes the following values:',channelWidths,'Hz'

    print 'Loading flags and intervals ...'
    flag=t.getcol('FLAG')[selection]          # flagged data have flag = True
    interval=t.getcol('INTERVAL')[selection]
    t.close()

    print 'The *flag* array has shape (Nr_integrations, Nr_channels, Nr_polarisations) =',flag.shape
    print 'The *interval* array has shape (Nr_integrations) =',interval.shape
    print 'The *channel* width array has shape (-, Nr_channels) =',channelWidths.shape

    rms=np.sqrt(2)*kB*tsys/eff/Aant/np.sqrt(channelWidths*interval.sum()*flag.shape[2])
    if len(rms.shape)==2 and rms.shape[0]==1: rms=rms[0]

    print 'Total Integration on selected field(s) = {0:.2f} h ({1:d} polarisations)'.format(interval.sum()/nrBaseline/3600,flag.shape[2])
    print 'The Stokes I theoretical natural rms ignoring flags is in the range:    *** ({0:.3e} - {1:.3e}) Jy ***'.format(np.nanmin(rms),np.nanmax(rms))

    return flag,interval,channelWidths,channelFreqs,rms

# Combine all input .MS files and predict natural rms
def combineMS(MS,tsys,eff,diam,makePlot,selectFieldName):

    # Derive quantities
    kB=1380.6                    # Boltzmann constant (Jy m^2 / K)
    Aant=np.pi*(diam/2)**2       # collecting area of 1 antenna (m^2)
    SEFD=2*kB*tsys/eff/Aant      # system equivalent flux density (Jy)

    # Print assumptions
    print ''
    print '--- Assumptions ---'
    print '  Tsys/efficiency      = {0:.1f} K (frequency independent)'.format(tsys/eff)
    print '  Dish diameter        = {0:.1f} m'.format(diam)
    print '    and therefore SEFD = {0:.1f} Jy'.format(SEFD)

    # Read MS files to get the flags and calculate single-MS natural rms values (ignoring flags)

    # Start with first file ...
    flag0,interval0,channelWidths0,channelFreqs0,rms0=processMS(MS[0],kB,tsys,eff,Aant,selectFieldName)
    rmsAll=[rms0]

    # ... and do the same for all other MS's appending to the flag array, checking that the channelisation is the same
    for ii in range(1,len(MS)):
        flagi,intervali,channelWidthsi,channelFreqsi,rmsi=processMS(MS[ii],kB,tsys,eff,Aant,selectFieldName)

        if channelWidths0.shape!=channelWidthsi.shape or (channelWidths0!=channelWidthsi).sum() or (channelFreqs0!=channelFreqsi).sum():
            print ''
            print ' CATASTROPHE!'
            print ' The input .MS file {1:s} has different channelization than the first input .MS file {2:s}'.format(ii,MS[ii],MS[0])
            print ' Cannot combine files to estimate their joint theoretical noise'
            print ' Aborting ...'
            sys.exit()
        else:
            flag0=np.concatenate((flag0,flagi),axis=0)
            interval0=np.concatenate((interval0,intervali),axis=0)
            rmsAll.append(rmsi)

    # Concatenate files
    print ''
    print '--- All input tables concatenated ---'
    print 'The concatenated *flag* array has shape (Nr_integrations, Nr_channels, Nr_polarisations) =',flag0.shape
    print 'The concatenated *interval* array has shape (Nr_integrations) =',interval0.shape
    print 'The concatenated *channel* width array has shape (-, Nr_channels) =',channelWidths0.shape
    
    print ''
    print '--- Result ---'
    print 'Calculating natural rms of all .MS files combined ...'

    # Reshape arrays
    interval0.resize((interval0.shape[0],1,1))
    channelWidths0.resize((channelWidths0.shape[1]))
    channelFreqs0.resize((channelFreqs0.shape[1]))

    # Calculate theoretical natural rms
    rmsAll=np.array(rmsAll)
    rmsAll=1./np.sqrt( (1./rmsAll**2).sum(axis=0) )
    unflaggedIntegration=(interval0*(1-flag0.astype(int))).sum(axis=(0,2)) # total integration per channel adding up all UNFLAGGED integrations and polarisations (sec)
    unflaggedIntegration[unflaggedIntegration==0]=np.nan
    rmsUnflagged=np.sqrt(2)*kB*tsys/eff/Aant/np.sqrt(channelWidths0*unflaggedIntegration)

    print 'The Stokes I theoretical natural rms ignoring flags is in the range:    *** ({0:.3e} - {1:.3e}) Jy ***'.format(np.nanmin(rmsAll),np.nanmax(rmsAll))
    print 'The Stokes I theoretical natural rms applying flags is in the range:    *** ({0:.3e} - {1:.3e}) Jy ***'.format(np.nanmin(rmsUnflagged),np.nanmax(rmsUnflagged))

    # Plot
    if makePlot==True:
        print ''
        print '--- Plot ---'
        print 'Busy plotting ...'
        import matplotlib.pyplot as plt
        plt.subplot(111)
        plt.plot(channelFreqs0,rmsAll*1e+3,'k:')
        plt.plot(channelFreqs0,rmsUnflagged*1e+3,'r-')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('rms (mJy)')
        plt.legend(('all data','unflagged data'),loc='lower right')
        plt.savefig('rms.png')
        print 'rms.png saved in working directory'



################
### RUN CODE ###
################

combineMS(MS,tsys,eff,diam,makePlot,selectFieldName)
