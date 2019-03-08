'''
Created on May 14, 2018
CAMERON EDIT
Edited on March 6th, 2019
@author: xuwang and Cameron Amos
'''
import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt, freqz
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import argparse
import csv
#------------------------------------------------------------------------
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data, method = "gust")
    return y
#------------------------------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("-s", "--srcPath", required=True,
    help="source folder with canopy rate file")
#==============
# ap.add_argument("-p", "--pltFile", required=True,
#     help="plot layout file")
#==============
# ap.add_argument("-t", "--targetPath", required=True,
#     help="output path")
#==============
args = ap.parse_args()
sourceFile = args.srcPath+"\\CanopyCoverRate.csv"
# outputPath = args.targetPath
# Number of ranges 2018AM3 28, 2018IPSR 51, 2017AYN 56
rn = 28
#------------------------------------------------------------------------
df = pd.read_csv(sourceFile, usecols=[1,2], dtype=float)
ampCanopy = df['Canopy_Rate_Index'].values.tolist() # -1 to 1
allCanopy = df['Canopy_Rate'].values.tolist() # percentage
imgFileNames = df
normAmp = [float(i)/max(allCanopy) for i in allCanopy]
#------------------------------------------------------------------------
Fs = 29.99 # Fs, sampling rate, frame rate, Frequency
Ts = 1.0/Fs # sampling interval
n = len(ampCanopy) # length of the signal(canopy cover)
T = n/Fs # total sampling period
t = np.arange(0, T, Ts) # time vector
# if len(t) != n:
#     t=np.arange(0,T+Ts,Ts)
#------------------------------------------------------------------------
# FFT (Fast Fourier Transformation)
k = np.arange(n)
frq = k/T
frq = frq[range(int(n/2))]
Y = np.fft.fft(ampCanopy)/n
Y = Y[range(int(n/2))]
#------------------------------------------------------------------------
# Design Low-pass filter
order = 8
fs = Fs       # sample rate, Hz
cutoff = 1.0  # desired cutoff frequency of the filter, Hz, 1.6 1-AM3 2.4-IPSR
b, a = butter_lowpass(cutoff, fs, order)
#------------------------------------------------------------------------
# Filter the signal by applying the low-pass filter
y1 = butter_lowpass_filter(ampCanopy, cutoff, fs, order)
y2 = butter_lowpass_filter(normAmp, cutoff, fs, order)
#------------------------------------------------------------------------
# Find local peaks, y1
peaks_y1, _ = find_peaks(y1)
# peaks_y_sd = peaks_y.sort(reverse=True)
y1_peaks_sd = sorted(y1[peaks_y1], reverse=True)
if len(y1_peaks_sd)>=rn:
    height_peaks_th = y1_peaks_sd[rn-1]
    peaks_y1_rn, _ = find_peaks(y1, height=height_peaks_th)
    
else:
    peaks_y1_rn, _ = find_peaks(y1)
t1_peaks = peaks_y1_rn/Fs
#------------------------------------------------------------------------
# Find local peaks, y2
peaks_y2, _ = find_peaks(y2)
# peaks_y_sd = peaks_y.sort(reverse=True)
y2_peaks_sd = sorted(y2[peaks_y2], reverse=True)
if len(y2_peaks_sd)>=rn:
    height_peaks_th = y2_peaks_sd[rn-1]
    peaks_y2_rn, _ = find_peaks(y2, height=height_peaks_th)
else:
    peaks_y2_rn, _ = find_peaks(y2)
t2_peaks = peaks_y2_rn/Fs
#------------------------------------------------------------------------
#------------------------------------------------------------------------
fig, ax = plt.subplots(4,1)
# Plot original signal
print(len(t))
print(len(ampCanopy))
if (len(t)<=len(ampCanopy)):
    tx = len(t)
else:
    tx = len(ampCanopy)
ax[0].plot(t[0:tx-1],ampCanopy[0:tx-1],'r-',linewidth=2)
ax[0].plot(t[0:tx-1],normAmp[0:tx-1],'b',linewidth=1)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('C. Cover Indicator')
ax[0].set_ylim(-2,2)
#------------------------------------------------------------------------
# Plot frequency domain
ax[1].plot(frq,abs(Y),'r')
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('|Y(freq)|')
#------------------------------------------------------------------------
# Plot frequency response of the low-pass filter
w, h = freqz(b, a, worN=2000)
plt.subplot(4, 1, 3)
plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
plt.plot(cutoff, 0.5*np.sqrt(2), 'ko')
plt.axvline(cutoff, color='k')

plt.xlim(0, 0.5*fs)
plt.title("Lowpass Filter Frequency Response")
plt.xlabel('Frequency [Hz]')
plt.ylabel('Attenuation')
plt.grid()
#------------------------------------------------------------------------
# Plot filtered signal
plt.subplot(4, 1, 4)
# plt.plot(t[0:tx-1], y1[0:tx-1], 'r-', linewidth=2, label='Filtered Index')
# plt.plot(t[0:tx-1], y2[0:tx-1], 'b', linewidth=1, label='Filtered Rate')
plt.plot(t[0:tx-1], y1[0:tx-1], 'r-', linewidth=2)
plt.plot(t[0:tx-1], y2[0:tx-1], 'b', linewidth=1)
plt.plot(t1_peaks, y1[peaks_y1_rn], 'gx')
plt.plot(t2_peaks, y2[peaks_y2_rn], 'kx')

plt.xlabel('Time [sec]')
plt.ylabel('CCover Indicator')
plt.grid()
plt.legend()
#------------------------------------------------------------------------
plt.subplots_adjust(hspace=1)
plt.savefig(args.srcPath+"\\Fig_process.png")
plt.show()
plt.close()
#========================================================================
# sourceFile = args.srcFile
# outputPath = args.targetPath
st_range=sourceFile.find("_C0")+3
rangeNum = int(sourceFile[st_range:st_range+2])
print("Range Number: %d" % rangeNum)

df2 = pd.read_csv(sourceFile, usecols=[0], dtype=str)
imgFileNames = df2['Image_file'].values.tolist() # image file nams
finalFile = open(args.srcPath+"\\peaks.csv",'wt')

try:
    # Create final output file
    writer = csv.writer(finalFile, delimiter=',', lineterminator='\n')
    # Header row if needed
    writer.writerow(('Column','Image_file','Image_file_1','Image_file_2','Canopy_Rate_Index','Canopy_Rate'))

    
    # Construct peaks_merged_rn:
    #------------------------------------------------------------------------
    similar_threshold = 16
    # This is a constant that has been roughly assumed to be half the distance (in images) from one plot to the next.
    # If two peaks are only x images apart or less, assume the same microplot, and average together
    # If two peaks are greater than x images apart, assume to be of different microplots and remain separated
    #------------------------------------------------------------------------

    peaks_merged_rn = []
    for o in range(len(peaks_y1_rn)):
        averaged = False
        for t in range(len(peaks_y2_rn)):
            # If any number in y1 is within the similarity threshold to any number in y2, average and append
            if abs(peaks_y1_rn[o] - peaks_y2_rn[t]) <= similar_threshold:
                peaks_merged_rn.append(int(0.5*(peaks_y1_rn[o]+peaks_y2_rn[t]))) 
                averaged = True 
                break

        # If not one match between y1 and y2 was made, append the number in y1    
        if averaged == False:
            peaks_merged_rn.append(int((peaks_y1_rn[o]))) 
    
    for t in range(len(peaks_y2_rn)):
        similar = False
        for m in range(len(peaks_merged_rn)):
            # If the number in y2 has already been included within the threshold, do not add it
            if abs(peaks_y2_rn[t] - peaks_merged_rn[m]) <= similar_threshold:
                similar = True
    
        # If the number in y2 is unique to the merged imf list, add it
        if similar == False:
            peaks_merged_rn.append(int((peaks_y2_rn[t])))
          
    peaks_merged_rn.sort() #Sorts merged so that it is from least to greatest

    # If merged list does not have enough ranges,
    # and there is an obvious gap of at least 2.5x the amount of the similarity threshold,
    # then insert the average of the two numbers between them (starting at the largest gap)

    while len(peaks_merged_rn) < rn:
        differences = []
        for (x, y) in list(zip(peaks_merged_rn[1:], peaks_merged_rn[:-1])):
            differences.append(abs(x-y))
        max_gap = max(differences)
        max_gap_index = differences.index(max_gap)
        if max_gap >= (2.5*similar_threshold):
            peaks_merged_rn.insert(max_gap_index + 1, (int(0.5*(peaks_merged_rn[max_gap_index] + peaks_merged_rn[max_gap_index+1]))))
            print("Peak gap detected. Inserted averaged peak at index ", (max_gap_index + 1))
        else:
            break
    print("Peaks 1: ", peaks_y1_rn)
    print("Peaks 2: ", peaks_y2_rn)
    print("Peaks M: ", peaks_merged_rn)       
    #------------------------------------------------------------------------
    
    # Assign peaks_len to the length of the largest list
    if len(peaks_y1_rn) > len(peaks_y2_rn) and len(peaks_y1_rn) > len(peaks_merged_rn):
        peaks_len = len(peaks_y1_rn)
    elif len(peaks_y2_rn) > len(peaks_y1_rn) and len(peaks_y2_rn) > len(peaks_merged_rn):
        peaks_len = len(peaks_y2_rn)
    else:
        peaks_len = len(peaks_merged_rn)
    
    if (rangeNum % 2) == 0: # even ranges, go from north to south
        for o in range(peaks_len):
            col = o+1
            try:
                imf = imgFileNames[peaks_merged_rn[o]] 
            except IndexError:
                imf = "Null"
            
            try:
                imf1 = imgFileNames[peaks_y1_rn[o]]
            except IndexError:
                imf1 = "Null"
                
            try:
                imf2 = imgFileNames[peaks_y2_rn[o]]
            except IndexError:
                imf2 = "Null"
            
            try:
                plotRate = ampCanopy[peaks_y1_rn[o]]
            except IndexError:
                plotRate = "Null"
                
            try:
                plotRateAll = ampCanopy[peaks_y1_rn[o]]
            except IndexError:
                plotRateAll = "Null"

            #---------------------------------------------------
            writer.writerow((col, imf, imf1, imf2, plotRate, plotRateAll))
    
    else: # odd ranges, go from south to north, flip
        peaks_merged_rn.reverse() #Reverse the merged list, as it is always constructed least to greatest regardless of y1 and y2 order
        for o in range(peaks_len):
            col = o+1
            try:
                imf = imgFileNames[peaks_merged_rn[o]] 
                print (o,imf)
            except IndexError:
                imf = "Null"
                
            try:
                imf1 = imgFileNames[peaks_y1_rn[rn - 1 - o]]
            except IndexError:
                imf1 = "Null"
                
            try:
                imf2 = imgFileNames[peaks_y2_rn[rn - 1 - o]]
            except IndexError:
                imf2 = "Null"
            
            try:
                plotRate = ampCanopy[peaks_y1_rn[rn - 1 - o]]
            except IndexError:
                plotRate = "Null"
                
            try:
                plotRateAll = ampCanopy[peaks_y1_rn[rn - 1 - o]]
            except IndexError:
                plotRateAll = "Null"
            
            #---------------------------------------------------
            writer.writerow((col, imf, imf1, imf2, plotRate, plotRateAll))
finally:
    finalFile.close()     
