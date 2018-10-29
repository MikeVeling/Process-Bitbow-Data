This folder contains all the data and the python code we used to analyze the fly brainbow data for the PNS. There are a few key files that I will describe below.

Input (folder): This folder contains all the input files in CSV format for analysis. Each of these input files contain coded brainbow data for each of the neurons for each hemisegment analyzed. The first line describes the columns. Date is the date the larva was imaged, Larva # denotes the larva that was imaged on that given day, Segment refers to the hemisgement that was analyzed, and Body side refers to the side of the body wall the hemisgement came from. Beyond that there are columns numbered 1-45. These correspond to the 45 neurons in the fly PNS. Each file put into this folder will be analyzed separately.

key.csv file: This file converts between the 45 separate neurons in the input files and the neuron groupings for distinguishable neurons. The paper used the “Peach” layer conversion nomenclature.
Active_hemisegments.csv file: This file allows you to limit the hemisemgments you would like to study. This allows you to test hypothesizes about development timing per hemisegment by limiting the data analysis to only particular hemisgements.

Process brainbow.py: code for running the analysis. This code was run in python 2.7.11 with the following dependency packages:

	Warnings, itertools, sys, json, statsmodels, operator, numpy, scipy, random, collections, csv, math, os, multiprocessing, and time.

After processing, the program builds the output folder will all the information we used in the publication. If you are interested in these outputs, please contact us for more information.
