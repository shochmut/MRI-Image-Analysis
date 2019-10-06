# MRI-Image-Analysis
This algorithm implements a seed-growing based segmentation algorithm in order to segment out the nucleus pulposus from the rest of the intervertebral disc MRI images. The region growing algorithm is written by Dirk-Jan Kroon and can be found on MathWorks. The Image Volume Algorithm Condensed implements a series of pre-processing steps combined with the region grow and post processing steps to segment the images. The image slices are then written into a single .tif file as well as .raw that can be taken for further processing. A hyper elastic warping finite element analysis will be done on these final processed images to develop a stress strain analysis of the intervertebral discs using FEBio.
