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

This program outputs an <b>output</b> folder and several subfolders. These subfolders contain detailed information on the analysis performed on subsets of the data, but it may simply be better to run the code with the default settings to generate an output to explore. This is important for establishing the reproducibility of the data across replicates. Therefore, it is important to describe the repeating folder structure in some detail as to familiarize yourself with the most relevant outputs for your study. Below is a summary of the data structure followed by a brief description of each of the folder and files.

### Data Structure
```
Output [active hemisegments]
   [input_file_name folders]
      all_data
         [see all data below]⋅
      [date folders]
         all_data
            [see all data below]⋅
         [larva folders]
            [see all data below]⋅

all_data
   color_number_stats.csv
   color_stats.csv
   group_color_stats.csv
   [neuron grouping layer folder]
      FDR_matrix.csv
      pval_dic.json
      pval_matrix.csv
      Match_fraction_matrix.csv
      Match_precent_matrix.csv
      by_the_colors.csv
```
### Output Folder description
<b>Output [active hemisegments]</b>: This is the root folder for an analysis performed on a set of hemisegments. By default, this is simply set to all hemisegments. To change this simply provide your hemisgements list in the <b>Active_hemisegments.csv</b>.

<b>[input_file_name folders]</b>: These folders are nested within the <b>Output [active hemisegments]</b> folder. They contain all analysis based on data provided within a particular <b>input_file_name</b> from the <b>input</b> folder. Within this folder, there is a subfolder called <b>all_data</b> containing the analysis performed on everything. It also contains subfolders based on dates.

<b>[date folders]</b>: These folders are nested within the <b>[input_file_name folders]</b> folder. They contain all analysis based on data provided from a particular date within a particular <b>input_file_name</b> from the <b>input</b> folder. Within this folder, there is a subfolder called <b>all_data</b> containing the analysis performed on everything. It also contains subfolders based on larva.

<b>[larva folders]</b>: These folders are nested within the <b>[date folders]</b> folder. They contain all analysis based on larva from a particular data provided from a particular date within a particular <b>input_file_name</b> from the <b>input</b> folder. Within this folder, there are several files that are all shared between this and the <b>all_data</b> folder. See that section for more details.

<b>all_data folder</b>: These folder contain an iteration of the analysis based on all the data from a particular data set. For example, if the <b>all_data folder</b> was nested within an <b>[input_file_name folders]</b> it would contain an analysis based on a particular input file from the <b>input</b> folder. Below are details about the files within this folder.

### all_data contents description

<b>color_number_stats.csv</b>: This file contains information about the number of different florescent proteins that were activated by all neurons within a given dataset defined by the folder this file is found within.

<b>color_stats.csv</b>:This file contains information about the identity of the florescent proteins that were activated by all neurons within a given dataset defined by the folder this file is found within.

<b>group_color_stats.csv</b>: This file contains information about particular color sets that were activated by all neurons within a given dataset defined by the folder this file is found within. Of note, this file can be used to calculate the probability χ<sub>c</sub> values for calculating the test statistic based on the colors neuron pairs share.

<b>[neuron grouping layer folder]</b> This folder contains information based on the neuron grouping definitions provided in the <b>key.csv</b> input file.

<b>FDR_matrix.csv</b> This matrix provides FDR values for the relatedness test performed between every pair of neurons. Each neuron appears as a row and a column. The intersection of the two neurons in the upper right corner indicates its FDR (see paper for methods)

<b>pval_dic.json</b> This file is a storage of a dictionary that was generated during the 100000 iterations of data randomization to establish the random distribution of matching colors based in your data. If this file is deleted Python will attempt to regenerate it, but it will cost a lot of CPU time.

<b>pval_matrix.csv</b>This matrix provides non corrected p-values for the relatedness test performed between every pair of neurons. Each neuron appears as a row and a column. The intersection of the two neurons in the upper right corner indicates its FDR (see paper for methods)

<b>Match_fraction_matrix.csv</b>This matrix provides the fraction of matches between every pair of neurons. Each neuron appears as a row and a column. The intersection of the two neurons in the upper right corner indicates its FDR (see paper for methods)

<b>Match_precent_matrix.csv</b>This matrix provides the percent of matches between every pair of neurons. Each neuron appears as a row and a column. The intersection of the two neurons in the upper right corner indicates its FDR (see paper for methods)

<b>by_the_colors.csv</b> <b>This is the main output of the analysis</b> This output contains all information about the relationship between pairs of neurons. Each line is a pair of neurons and the columns contain various data about that pair. As this file is so important, I will delineated the data within this file below.

### by_the_colors.csv column explanation.

<b>Neuron 1</b>: Name of the first neuron or neuron group based on the <b>key.csv</b> file

<b>Neuron 2</b>: Name of the second neuron or neuron group based on the <b>key.csv</b> file

<b>Neuron 1 group size</b>: Size of the first neuron or neuron group based on the <b>key.csv</b> file

<b>Neuron 2 group size</b>: Size of the second neuron or neuron group based on the <b>key.csv</b> file

<b>Simple Match Name</b>: Simple concatenated name of both neuron 1 and neuron 2 to come up with an identifier for the pair.

<b>Total number of times both neurons were observed</b>: This counts the total number of times both neuron 1 and neuron 2 were observed within the same hemisegment. It does not consider if either of them were colored just if they could be identified.

<b>Total number of times both neurons were observed and at least one had a color</b>: This counts the total number of times both neuron 1 and neuron 2 were observed within the same hemisegment and at least one of them was not blank (not "00000" as a color). We consider this to be our total number of match trials for a pair as they should have a chance to match at this point. 

<b>Total number of times both neurons had color</b>: This counts the total number of times both neuron 1 and neuron 2 had a non "00000" color.

<b>Total number of times both neurons had the same color (matches)</b>: This counts the total number of times two neurons shared the same non "00000" color code. This column does not consider the particular color the match is made in only if there is a match.

<b>Total number of times both neurons were observed and at least one had a color - Total number of times both neurons had the same color (non-matches)</b>: This is a calculation of the number of times two neurons had an opportunity to match but did not. We therefore consider this a non-match column. The title describes the equation (total-matches).

<b>Fractional match</b>: This is simply the fraction of times these two neurons matched.

<b>Relatedness test statistic</b>: This is the test statistic calculated for this pair of neurons based on the number of matches and identity of the colors that the matches were made in (see paper for more details)

<b>Relatedness test P-Value</b>: This is a conversion of the test statistic into a p-value based on the randomization of the data, recalculation of the Relatedness test statistic for this pair of neurons, and then comparison of the given test statistic to the background distribution.

<b>FDR</b>: This is an FDR correction for the P-Values based on the parameters set in the program. By default α=0.05 using the Benjamini-Hochberg FDR correction (see https://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html).

<b>Total number of times both neurons were observed and at least one had a color randomly</b>: This column counts the number of times this pair of neurons was observed where at least one of them had color in the randomized dataset.

<b>Total number of times both neurons had the same color randomly (matches)</b>: This column counts the number of times this pair of neurons had a matching color in the randomized dataset.

<b>Total number of times both neurons were observed and at least one had a color randomly - Total number of times both neurons had the same color randomly (non-matches)</b>: This column counts the number of times this pair of neurons had a matching non matching color in the randomized dataset. The title describes the equation (total-matches).

<b>[set of colors]</b>: Beyond the last column there is a set of color codes. These codes correspond to the colors in the analysis and the values of the cells indicates how many times this particular neuron pair had a match in that particular color.

## Concluding remarks

As these analyses are somewhat complex, I expect there may be some confusion about all the inputs and outputs. As such, I would be happy to work with you to make it understandable and or help you tailor the code to your needs. Feel free to contact [me](mailto:mike.veling@gmail.com) if you would like my help.
