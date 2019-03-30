# Plot-Identification-Pipeline

The CanopyFilter interprets the peaks of the wavelength generated from how much green there is in the middle of the image. A peak therefore can be assumed to be the location of an image that is directly over a microplot. There are two formulas that generate two different waveforms, and the position data from these peaks are averaged if non-unique and inserted if unique into a merged peaks list.

If the quadcopter is offset and not directly over the microplot, a good amount of the field of view is taken up by the alleyway between columns either to the left or right. This can result in the cropped images of that side composing of dirt or having more edge variance than intended.

Occasionally, the camera feed (and therefore the images) will have a ‘skip’ in the otherwise smooth progression. These skips can fly past anywhere from half of a microplot to 2 full micro plots.

Occasionally, an entire column can go missing for an unknown reason. This could be due to either human error, when the quadcopter stops in the middle of the field for a recharge and is sent to two columns ahead instead of the next one. It could also simply be due to a column file missing.

Occasionally, the camera feed will either begin or end in the middle of the field. This is due to lag from the camera controls or by human error.

The two algorithms, about a third of the time, will misalign the image ‘positions’, and one algorithm will be offset compared to the other one. This is due to one algorithm not having a peak high enough to register as a microplot in it’s list. The plot lists are mostly one short than normal in this situation. This problem was fixed.

Occasionally the peaks.csv file will not work and fail to yield anything. The reason for this is unknown.

Occasionally datasets will be blurry and unusable. This is either due to human error or due to the camera not being reset before recording.
