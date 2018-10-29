# Process-Brainbow-Data
This git contains the script and sample input files for the Brainbow analysis associated with the Veling Et. Al. 2019 paper (link when published). The main analysis allows for an estimation of neuron lineage relationships based on Brainbow data. We used it to establish neuron relationships in drosophila melanogaster's peripheral nervous system (PNS). 

This estimation was based on data from a 5 color Brainbow system. This system stochastically labels neurons during development with one of 31 potential colors. This color is maintained through development to maturity. This allows us to identify related neurons based shared color expression. Neurons that express the same color likely developed from the same progenitor cell (early neuron precursor cells we labeled early in development). Observations of shared color frequencies allows us estimate how closely related sets of neurons are. These relationships can be used to estimate neuron developmental patterning as described in more details in the paper.

The following sections describe the input, execution, and outputs of this analysis. Feel free to contact [me](mailto:mike.veling@gmail.com) for further clarification if you would like.

## Input files for analysis

### Input (folder)
This folder contains all the input files in CSV format for analysis. Each of these input files contain coded brainbow data for each of the neurons for each hemisegment analyzed. The first line describes the columns.

<b>Date:</b> the date the larva was imaged
<br><b>Larva #:</b> denotes the larva ID that was imaged on that given day 
<br><b>Segment:</b> the hemisgement that was analyzed
<br><b>Body side:</b> the side of the body wall the hemisgement came from.

Beyond that there are columns numbered 1-45. These correspond to the 45 neurons in the fly PNS. Each file put into this folder will be analyzed separately.

### key.csv
This file converts between the 45 separate neurons in the input files and the neuron groupings for distinguishable neurons. The paper used the “Veling Neuron Groupings” layer conversion nomenclature.

### Active_hemisegments.csv 
This file allows you to limit the hemisemgments you would like to study. This allows you to test hypothesizes about development timing per hemisegment by limiting the data analysis to only particular hemisgements.

### Process brainbow.py
code for running the analysis. This code was run in python 2.7.11 with the following dependency packages:

itertools, sys, json, statsmodels, operator, numpy, scipy, random, collections, csv, math, os, multiprocessing, and time.

## Running the code
