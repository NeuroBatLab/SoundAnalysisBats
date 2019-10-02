# Import math, plotting and sound libraries

import numpy as np
import matplotlib.pyplot as plt
from soundsig.sound import BioSound 
from soundsig.sound import WavFile
import os

def BioSound1(dir_path, plotMe=True, normalize = False, spec_sample_rate=1000,freq_spacing=50, min_freq=100, max_freq=40000, cutoff_freq=75, amp_sample_rate=1000, f_high=60000, maxFund=4000, minFund=300, lowFc=500, highFc=10000, minSaliency=0.5, debugFig=0, minFormantFreq=500, maxFormantBW=500, method='Stack', N = 50):
    # Go to the folder that has the wav files
    os.chdir(dir_path)
    # This will be the output directory
    if not os.path.exists('h5files'):
        os.makedirs('h5files')
   


    # Find all the wave files 
    isound = 0   
    for fname in os.listdir('.'):
        if fname.endswith('.wav'):
            isound += 1;
            
            if isound == N+1:
                break
            
            # Read the sound file
            #print ('Processing sound %d:%s\n' % (isound, fname))*
            soundIn = WavFile(file_name=fname) 
            filename, file_extension = os.path.splitext(fname)
            
            
            # Normalize if wanted
            if normalize :
                maxAmp = np.abs(soundIn.data).max() 
            else :
                maxAmp = 1.0

            # Create BioSound Object
            myBioSound = BioSound(soundWave=soundIn.data.astype(float)/maxAmp, fs=float(soundIn.sample_rate), emitter = filename)
                 
            # Calculate the spectrogram and the rms
            myBioSound.spectroCalc(spec_sample_rate=spec_sample_rate, freq_spacing = freq_spacing, min_freq=min_freq, max_freq=max_freq)
            myBioSound.rms = myBioSound.sound.std() 
           
            # Calculate amplitude enveloppe
            myBioSound.ampenv(cutoff_freq = cutoff_freq, amp_sample_rate = amp_sample_rate)
            
            # Calculate the power spectrum
            myBioSound.spectrum(f_high=f_high)
        
            # Calculate fundamental and related values.  These are the default values.
            # For the estimation of the fundamental, four methods are available: 
            # 'AC' - Peak of the auto-correlation function
            # 'ACA' - Peak of envelope of auto-correlation function 
            # 'Cep' - First peak in cepstrum 
            # 'Stack' - Fitting of harmonic stacks (default - works well for zebra finches)
        
            myBioSound.fundest(maxFund = maxFund, minFund = minFund, lowFc = lowFc, highFc = highFc, minSaliency = minSaliency, debugFig = debugFig, minFormantFreq = minFormantFreq, maxFormantBW = maxFormantBW,  method=method)
          
            # Calculate the MPS
            myBioSound.mpsCalc(window=0.1, Norm = True)
    
                    
            if plotMe and isound <= 2:
               # print('                Bird %s    Call Type %s' % (myBioSound.emitter, myBioSound.type))
                myBioSound.plot(DBNOISE=50, f_low=250, f_high=f_high)  
    

            # Save the results
            fh5name = 'h5files/%s.h5' % (filename)
            myBioSound.saveh5(fh5name)
