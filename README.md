# track_analyzer
A MATLAB class for segmentation of image data (2D, 3D, cells, particles, etc). 

This was designed to analyze experiments from my PhD work on YAP distribution dynamics and nascent transcription dynamics (https://www.biorxiv.org/content/10.1101/539049v2). The main purposes of this class is to track cells, and information about cells over time such as protein distribution and nascent transcribing regions within the cell nucleus. Although these are specific cases, there are general utilities for analyzing miscroscopy data. 

Segmentation functions (segment_2D, segment_3D, segment_3D_gpu) have many parameters and possible steps used for segmenting miscroscopy images. Examples will be provided to make sense of all options and parameters. 

