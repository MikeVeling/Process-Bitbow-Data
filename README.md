# Process-Brainbow-Data
This git contains the script and sample input files for the Brainbow analysis associated with the Veling Et. Al. 2019 paper (link when published). The main analysis allows for an estimation of neuron lineage relationships based on Brainbow data. We used it to establish neuron relationships in <i>Drosophila melanogaster's</i> peripheral nervous system (PNS). 

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

### Process_brainbow.py
code for running the analysis. More details will be provided in the next section.

## Running the code
This code runs on python 2.7.11 with several dependency packages. These packages can be installed by running the following code in your terminal

```
python -m pip install json
python -m pip install statsmodels
python -m pip install numpy
python -m pip install scipy
python -m pip install multiprocessing
```

Once these packages are installed to your python environment, you should be able to simply run the code with no arguments assuming it is in the same folder as the provided <b>Input</b> folder as well as the provided sample input files. The <b>key.csv</b> and <b>Active_hemisegments.csv</b> also needs to be present in that same folder as is in the git.

```
python Process_brainbow.py
```

This code unfortunately runs with no arguments because variables are coded within the script (see first lines in the script). If this code becomes more popular, I could easily edit the code to accept arguments from the command line. If you would like to edit the running parameters, please see the code for further details on the editable variables or contact [me](mailto:mike.veling@gmail.com) for further clarification.

## Outputs

This program outputs an <b>output</b> folder and several subfolders. These subfolders contain detailed information on the analysis performed on subsets of the data. This is important for establishing the reproducibility of the data across replicates. Therefore, it is important to describe the repeating folder structure in some detail as to familiarize yourself with the most relevant outputs for your study.

<b>Output [active hemisegments]</b>
   <b>[input_file_name folders]</b>
      <b>all_data</b>
         <b>[see all data below]⋅</b>
      <b>[date folders]</b>
         <b>all_data</b>
            <b>[see all data below]⋅</b>
         <b>[larva folders]</b>
            <b>[see all data below]⋅</b>

<b>all_data</b>
   <b>color_number_stats.csv</b>
   <b>color_stats.csv</b>
   <b>group_color_stats.csv</b>
   <b>[neuron grouping layer folder]</b>
      <b>by_the_colors.csv</b>
      <b>FDR_matrix.csv</b>
      <b>pval_dic.json</b>
      <b>pval_matrix.csv</b>
      <b>Match_fraction_matrix.csv</b>
      <b>Match_precent_matrix.csv</b>


<b>Output [active hemisegments]</b>: This is the root folder for an analysis performed on a set of hemisegments. By default this is simply set to all hemisegments. To change this simply provide your hemisgements list in the <b>Active_hemisegments.csv</b>.
<b>[input_file_name folders]</b>: These folders are nexted within the 
