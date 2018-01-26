#!/usr/bin/env python

# This module defines functions useful for various scripts in the noisy repository

# To be done
# - MS files with 4 polarisation products (currently uses all available pols)



######################
### IMPORT MODULES ###
######################



import numpy as np
import pyrap.tables as tables
import sys,os
try: from astropy.io import fits
except ImportError:
    print ''
    print ' WARNING:'
    print ' Could not find astropy python module; certain functions may crash'
    print ''
import math as mt
import scipy as sp
from scipy import stats
from scipy import optimize
import scipy.ndimage as nd
from distutils.version import StrictVersion, LooseVersion

## check numpy and scipy version numbers for the nanmedian function import
if LooseVersion(np.__version__)>=LooseVersion('1.9.0'): from numpy import nanmedian
elif LooseVersion(sp.__version__)<LooseVersion('0.15.0'): from scipy.stats import nanmedian
else: from scipy import nanmedian


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
    else: print 'The channel width takes the following values:',np.unique(channelWidths),'Hz'

    print 'Loading flags and intervals ...'
    flag=t.getcol('FLAG')[selection]          # flagged data have flag = True
    interval=t.getcol('INTERVAL')[selection]
    if np.unique(interval).shape[0]==1: print 'Interval = {0:.5e} sec'.format(np.unique(interval)[0])
    else: print 'The channel width takes the following values:',interval,'sec'
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

    print 'The Stokes I theoretical natural rms ignoring flags has median and range:    *** {0:.3e} Jy/b, ({1:.3e} - {2:.3e}) Jy/b ***'.format(np.nanmedian(rms),np.nanmin(rms),np.nanmax(rms))

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
    print     '    and therefore SEFD = {0:.1f} Jy'.format(SEFD),
    if tsyseffFile==None: print ''
    else: print '(median)'

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

    print 'The Stokes I theoretical natural rms ignoring flags has median and range:    *** {0:.3e} Jy/b, ({1:.3e} - {2:.3e}) Jy/b ***'.format(np.nanmedian(rmsAll),np.nanmin(rmsAll),np.nanmax(rmsAll))
    print 'The Stokes I theoretical natural rms applying flags has median and range:    *** {0:.3e} Jy/b, ({1:.3e} - {2:.3e}) Jy/b ***'.format(np.nanmedian(rmsUnflagged),np.nanmin(rmsUnflagged),np.nanmax(rmsUnflagged))

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



### Measure noise of single cube
def GetRMS(cube,rmsMode='negative',zoomx=1,zoomy=1,zoomz=1,nrbins=10000,verbose=0,min_hist_peak=0.05,sample=1):
    sh=cube.shape

    if len(sh) == 2:
        # add an extra dimension to make it a 3d cube
        cube = np.array([cube])
    sh=cube.shape

    x0,x1=int(mt.ceil((1-1./zoomx)*sh[2]/2)),int(mt.floor((1+1./zoomx)*sh[2]/2))+1
    y0,y1=int(mt.ceil((1-1./zoomy)*sh[1]/2)),int(mt.floor((1+1./zoomy)*sh[1]/2))+1
    z0,z1=int(mt.ceil((1-1./zoomz)*sh[0]/2)),int(mt.floor((1+1./zoomz)*sh[0]/2))+1
    if verbose: print (' Estimating rms on subcube (x,y,z zoom = %.0f,%.0f,%.0f) ...' % (zoomx, zoomy, zoomz))
    if verbose: print (' Estimating rms on subcube sampling every %i voxels ...' % (sample))
    if verbose: print (' ... Subcube shape is ' + str(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample].shape) + ' ...')

    if rmsMode=='negative':
        nrbins=max(100,int(mt.ceil(float(np.array(cube.shape).prod())/1e+5))) # overwrites nrbins value!!!
        cubemin=np.nanmin(cube)
        bins=np.arange(cubemin,abs(cubemin)/nrbins-1e-12,abs(cubemin)/nrbins)
        fluxval=(bins[:-1]+bins[1:])/2
        rmshisto=np.histogram(cube[z0:z1:sample,y0:y1:sample,x0:x1:sample][~np.isnan(cube[z0:z1:sample,y0:y1:sample,x0:x1:sample])],bins=bins)[0]

        nrsummedbins=0
        while rmshisto[-nrsummedbins-1:].sum()<min_hist_peak*rmshisto.sum():
            nrsummedbins+=1
        if nrsummedbins:
            if verbose: print ('  ... adjusting bin size to get a fraction of voxels in central bin >= ' + str(min_hist_peak))
            nrbins/=(nrsummedbins+1)
            bins=np.arange(cubemin,abs(cubemin)/nrbins-1e-12,abs(cubemin)/nrbins)
            fluxval=(bins[:-1]+bins[1:])/2
            rmshisto=np.histogram(cube[z0:z1:sample,y0:y1:sample,x0:x1:sample][~np.isnan(cube[z0:z1:sample,y0:y1:sample,x0:x1:sample])],bins=bins)[0]

        rms=abs(optimize.curve_fit(GaussianNoise,fluxval,rmshisto,p0=[rmshisto.max(),-fluxval[rmshisto<rmshisto.max()/2].max()*2/2.355])[0][1])

    elif rmsMode=='mad':
        rms=nanmedian(abs(cube[z0:z1:sample,y0:y1:sample,x0:x1:sample]-nanmedian(cube[z0:z1:sample,y0:y1:sample,x0:x1:sample],axis=None)),axis=None)/0.6745

    elif rmsMode=='std':
        rms=np.nanstd(cube[z0:z1:sample,y0:y1:sample,x0:x1:sample],axis=None,dtype=np.float64)

    if verbose: print('  ... %s rms = %.2e (data units)' % (rmsMode, rms))
    return rms



### Measure noise of N cubes

def MeasureCubeNoise(FITS,plotName,randomizeChanOrder):

    for ii in range(0,len(FITS)):

        print ''
        print '--- Working on file {0:s} ---'.format(FITS[ii])

        # Open cube
        f=fits.open(FITS[ii])
        cube=f[0].data
        head=f[0].header
        f.close()
        cube=np.squeeze(cube)

        # Flag all-zero channels
        print(' Flagging all-zero channels ...')
        zeroChans=(cube==0).prod(axis=(1,2))
        for cc in range(zeroChans.shape[0]):
            if zeroChans[cc]: cube[cc]=np.nan
        print(' ... {0:d} channels flagged'.format(zeroChans.sum()))

        # Measure global rms
        globalRms=GetRMS(cube,rmsMode='mad',verbose=1)

        # Measure rms per channel and take median
        print(' Measuring median and range of single-channel {0:s} rms values (excluding all-nan channels) ...'.format('mad'))
        singleChanNoise=[]
        for cc in cube:
            if not np.isnan(cc).prod(): singleChanNoise.append(GetRMS(cc,rmsMode='mad'))

        print('  median: {0:.2e} (data units)'.format(nanmedian(np.array(singleChanNoise))))
        print('  range : ({0:.2e}-{1:.2e}) (data units)'.format(np.nanmin(singleChanNoise),np.nanmax(singleChanNoise)))

        # Sum channels and measure noise
        SumChanNoise(cube,plotName,randomizeChanOrder)


### Measure noise as a function of number of averaged channels

def SumChanNoise(cube,plotName,randomizeChanOrder):
    print ' Measuring noise increase as a function of number of summed channels ...'

    # Take first unmasked channel and measure its rms
    c0=0
    while np.isnan(cube[c0]).prod(): c0+=1
    avCube=cube[c0]
    xx,yy,zz,aa=[1,],[GetRMS(avCube,rmsMode='mad')],[GetRMS(avCube,rmsMode='mad')],[0]
    safelist=[c0,]

    # Add unmasked channels one by one and measure the noise
    #print('  {0:7d} {1:.3e} {2:.3e}  *  {3:10.3e}'.format(xx[-1],yy[-1],yy[0]*np.sqrt(xx[-1]),A))
    for cc in range(c0+1,cube.shape[0]):
        if randomizeChanOrder:
            cc=np.random.randint(c0+1,cube.shape[0])
            while cc in safelist: cc=np.random.randint(c0+1,cube.shape[0])
            safelist.append(cc)
        if not np.isnan(cube[cc]).prod():
            A=np.corrcoef(np.ravel(avCube),y=np.ravel(cube[cc]))[0,1]
            avCube+=cube[cc]
            xx.append(xx[-1]+1)
            yy.append(GetRMS(avCube,rmsMode='mad'))
            zz.append(GetRMS(cube[cc],rmsMode='mad'))
            aa.append(A)
            #print('  {0:7d} {1:.3e} {2:.3e}  *  {3:10.3e}'.format(xx[-1],yy[-1],yy[0]*np.sqrt(xx[-1]),A))

    # Plot
    if plotName!=None:
        print ''
        print '--- Plot ---'
        print 'Busy plotting ...'
        import matplotlib.pyplot as plt
        xx,yy=np.array(xx),np.array(yy)
        plt.subplot(211)
        plt.loglog(xx,yy,'ko')
        plt.loglog(xx,zz,'r.')
        plt.loglog(xx,yy[0]*np.sqrt(xx),'r-')
        plt.xlabel('Number summed channels')
        plt.ylabel('rms (data units)')
        plt.subplot(212)
        plt.semilogx(xx,aa,'ko')
        plt.xlabel('Number summed channels')
        plt.ylabel('corr coeff')
        plt.savefig(plotName)
        print 'Plot {0:s} saved (full path)'.format(plotName)

