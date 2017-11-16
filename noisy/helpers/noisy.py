# This module defines functions useful for various scripts in the noisy repository

# To be done
# - MS files with 4 polarisation products (currently uses all available pols)


######################
### IMPORT MODULES ###
######################



import numpy as np
import pyrap.tables as tables
import sys,os



########################
### DEFINE FUNCTIONS ###
########################


### Get Tsys/eff (possibly from file)
def GetTsyseff(tsyseff):
    if os.path.exists(tsyseff):
        print ' ( Tsys/eff from file {0:s} )'.format(tsyseff)
        tsyseffFile,tsyseff=tsyseff,np.loadtxt(tsyseff)
    else:
        try: tsyseff=float(tsyseff)
        except ValueError:
            print ''
            print ' CATASTROPHE!'
            print ' You set Tsys/eff = {0:s}'.format(tsyseff)
            print ' This is either a file that cannot be found or a value that cannot be converted to float'
            print ' Correct any mistakes and try again'
            print ' Aborting ...'
            sys.exit()
        tsyseffFile=None
    return tsyseffFile,tsyseff

### Interpolate input Tsys/eff table to observed frequencies
def InterpolateTsyseff(tsyseff,chans):
    print 'Interpolating Tsys/eff table to observed frequencies ...'
    return np.interp(np.ravel(chans),tsyseff[:,0],tsyseff[:,1])

### Get single-MS flags, intervals, channel widths, channel frequencies and calculate natural rms (ignoring flags)
def ProcessSingleMS(ms,kB,tsyseff,tsyseffFile,Aant,selectFieldName):
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

    print 'Total Integration on selected field(s) = {0:.2f} h ({1:d} polarisations)'.format(interval.sum()/nrBaseline/3600,flag.shape[2])
    if tsyseffFile!=None:
        rms=np.sqrt(2)*kB*InterpolateTsyseff(tsyseff,channelFreqs)/Aant/np.sqrt(channelWidths*interval.sum()*flag.shape[2])
    else:
        rms=np.sqrt(2)*kB*tsyseff/Aant/np.sqrt(channelWidths*interval.sum()*flag.shape[2])
    if len(rms.shape)==2 and rms.shape[0]==1: rms=rms[0]

    print 'The Stokes I theoretical natural rms ignoring flags is in the range:    *** ({0:.3e} - {1:.3e}) Jy ***'.format(np.nanmin(rms),np.nanmax(rms))

    return flag,interval,channelWidths,channelFreqs,rms



### Predict natural rms for an arbitrary number of MS files (both ignoring and applying flags)
def PredictNoise(MS,tsyseff,diam,plotName,selectFieldName):

    # Get Tsys/eff either from table (col1 = frequency, col2 = Tsys/eff) or as a float values (frequency independent Tsys/eff value)
    tsyseffFile,tsyseff=GetTsyseff(tsyseff)

    # Derive quantities
    kB=1380.6                                   # Boltzmann constant (Jy m^2 / K)
    Aant=np.pi*(diam/2)**2                      # collecting area of 1 antenna (m^2)
    if tsyseffFile==None:
        SEFD=2*kB*tsyseff/Aant                  # frequency independent system equivalent flux density (Jy)
    else:
        SEFD=2*kB*np.median(tsyseff[:,1])/Aant  # median system equivalent flux density (Jy)

    # Print assumptions
    print ''
    print '--- Assumptions ---'
    if tsyseffFile==None:
        print '  Tsys/efficiency      = {0:.1f} K (frequency independent)'.format(tsyseff)
    else:
        print '  Tsys/efficiency      = ({0:.1f} - {1:.1f}) K (range over input table {2:s})'.format(tsyseff[:,1].min(),tsyseff[:,1].max(),tsyseffFile)
        print '                         (frequency range = ({0:.3e} - {1:.3e}) Hz'.format(tsyseff[:,0].min(),tsyseff[:,0].max())
    print     '  Dish diameter        = {0:.1f} m'.format(diam)
    print     '    and therefore SEFD = {0:.1f} Jy'.format(SEFD)

    # Read MS files to get the flags and calculate single-MS natural rms values (ignoring flags)

    # Start with first file ...
    flag0,interval0,channelWidths0,channelFreqs0,rms0=ProcessSingleMS(MS[0],kB,tsyseff,tsyseffFile,Aant,selectFieldName)
    rmsAll=[rms0]

    # ... and do the same for all other MS's appending to the flag array, checking that the channelisation is the same
    for ii in range(1,len(MS)):
        flagi,intervali,channelWidthsi,channelFreqsi,rmsi=ProcessSingleMS(MS[ii],kB,tsyseff,tsyseffFile,Aant,selectFieldName)

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

    # Message concatenated files
    print ''
    print '--- All input tables concatenated ---'
    print 'The concatenated *flag* array has shape (Nr_integrations, Nr_channels, Nr_polarisations) =',flag0.shape
    print 'The concatenated *interval* array has shape (Nr_integrations) =',interval0.shape
    print 'The concatenated *channel* width array has shape (-, Nr_channels) =',channelWidths0.shape

    print ''
    print '--- Result ---'

    # Reshape arrays
    print 'Reshaping arrays ...'
    interval0.resize((interval0.shape[0],1,1))
    channelWidths0.resize((channelWidths0.shape[1]))
    channelFreqs0.resize((channelFreqs0.shape[1]))

    # Interpolate Tsys
    if tsyseffFile!=None:
       tsyseff=InterpolateTsyseff(tsyseff,channelFreqs0)

    # Calculate theoretical natural rms
    print 'Calculating natural rms of all .MS files combined (with and without flags)...'
    rmsAll=np.array(rmsAll)
    rmsAll=1./np.sqrt( (1./rmsAll**2).sum(axis=0) )
    unflaggedIntegration=(interval0*(1-flag0.astype(int))).sum(axis=(0,2)) # total integration per channel adding up all UNFLAGGED integrations and polarisations (sec)
    unflaggedIntegration[unflaggedIntegration==0]=np.nan
    rmsUnflagged=np.sqrt(2)*kB*tsyseff/Aant/np.sqrt(channelWidths0*unflaggedIntegration)

    print 'The Stokes I theoretical natural rms ignoring flags is in the range:    *** ({0:.3e} - {1:.3e}) Jy ***'.format(np.nanmin(rmsAll),np.nanmax(rmsAll))
    print 'The Stokes I theoretical natural rms applying flags is in the range:    *** ({0:.3e} - {1:.3e}) Jy ***'.format(np.nanmin(rmsUnflagged),np.nanmax(rmsUnflagged))

    # Plot
    if plotName!=None:
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
        plt.savefig(plotName)
        print 'Plot {0:s} saved (full path)'.format(plotName)
