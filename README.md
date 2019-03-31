# Plot Identification Pipeline

This image processing pipeline is designed isolation of compact, accurate images of cultivars in a micro-plot field for use in high-throughtput deep-learning-assisted phenotyping. Raw data is collected using a high-resolution camera mounted on a low-astitude UAV guided by gps markers. 

## Features / Summary of Files
### 4 distinct python files are employed during the pipeline, listed in order of use

- dngConversion.py: After an mp4 file of a pass down the microplot has been converted to .dng files at a certain fps, this python script simply converts themn to .jpg. The old .dng files may be moved to a local file for readability if desired.

- boundaryTest.py: Processes each image (through a process detailed below) and rates them according to where the position of the alleyway is relative to a micro-plot. These findings are stored in CanopyCoverRate.csv.

- canopyBPFilter.py: Utilizes CanopyCoverRate.csv to select .jpg files that are directly over the microplot, and show as little amount of the alleyway between the microplots as possible. These results are ouputted into peaks.csv. __NOTE:__ The variable 'similar_threshold' is a variable that is used in the script to determine if to hits are recording the same or different microplot. This variable *MUST* be manually set, as flight speed, the fps used to create the .dng files, as well as microplot dimensions all factor into this length. To find the the value of the variable 'similar_threshold', simply find 2 images that are centered above adjacent microplots and take the number of images between them divided by 2. 

- plotMatch.py: Crops 6 images centered horizontally inside each image specified in peaks.csv for use in deep learning algorithms. These images are stored in a folder labeled 'cropped'.

The CanopyFilter interprets the peaks of the wavelength generated from how much green there is in the middle of the image. A peak therefore can be assumed to be the location of an image that is directly over a microplot. There are two formulas that generate two different waveforms, and the position data from these peaks are averaged if non-unique and inserted if unique into a merged peaks list.

If the quadcopter is offset and not directly over the microplot, a good amount of the field of view is taken up by the alleyway between columns either to the left or right. This can result in the cropped images of that side composing of dirt or having more edge variance than intended.

Occasionally, the camera feed (and therefore the images) will have a ‘skip’ in the otherwise smooth progression. These skips can fly past anywhere from half of a microplot to 2 full micro plots.

Occasionally, an entire column can go missing for an unknown reason. This could be due to either human error, when the quadcopter stops in the middle of the field for a recharge and is sent to two columns ahead instead of the next one. It could also simply be due to a column file missing.

Occasionally, the camera feed will either begin or end in the middle of the field. This is due to lag from the camera controls or by human error.

The two algorithms, about a third of the time, will misalign the image ‘positions’, and one algorithm will be offset compared to the other one. This is due to one algorithm not having a peak high enough to register as a microplot in it’s list. The plot lists are mostly one short than normal in this situation. This problem was fixed.

Occasionally the peaks.csv file will not work and fail to yield anything. The reason for this is unknown.

Occasionally datasets will be blurry and unusable. This is either due to human error or due to the camera not being reset before recording.
