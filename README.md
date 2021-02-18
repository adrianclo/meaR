# meaR

## Introduction
Multi-electrode arrays (MEAs) are becoming increasingly popular to investigate activity of neuronal assemblies upon genetic, pharmacological and/or optogenetic manipulation. Neuronal activity can be assessed on acute brain slices, neuronal uniform or mixed-cell cultures, and human pluripotent stem cells, to name a few. It is a means for high-throughput screening and these MEAs come in different sizes, resolution and presets from different companies. However, MEA data in general are very complex, and the output files generated from these MEAs are rarely immediately usable for statistical analysis. As such, a tabular form (<i>i.e.</i>, a dataframe or tibble) is more recommended which also included the required parameters for statistical analyses.

<hr>

## Purpose
For these reasons, this repository contains a collection of functions to extract, visualize and analyze MEA data, specifically from <a href = "https://www.multichannelsystems.com/products/vitro-mea-systems">MultiChannel Systems</a>. <br>
<b>Note, that the script is in continuous development</b> and earlier versions have been used for the following papers:

- <b>Lo AC</b>, Rajan N, Gastaldo D, Telley L, Hilal ML, Buzzi A, Simonato M, Achsel T, Bagni C (accepted). Absence of RNA-binding protein FXR2P prevents prolonged phase of kainate-induced seizures. <i>EMBO Rep</i>
- Dominguez-Iturza N, <b>Lo AC</b>, Shah D, Armendariz M, Vannelli A, Mercaldo V, Trusel M, Li KW, Gastaldo D, Santos AR, Callaerts-Vegh Z, D'Hooge R, Mameli M, Van der Linden A, Smit AB, Achsel T, Bagni C (2019). The autism- and schizophrenia-associated protein CYFIP1 regulates bilateral brain connectivity and behaviour. <i>Nat Commun</i> <b>10</b>: 3454. <a href = "https://www.nature.com/articles/s41467-019-11203-y">[link to article]</a>

<hr>

## Step-by-step guide
### MEA output file
Settings for spike detection threshold are done within the Multi Channel Experimenter Spike Analyzer itself. Data output files are in .txt format. These contain for each separate electrode (N=60) spike data on the occurrence (s), interspike interval (isi in ms) and frequency (in Hz). The electrode ID is the combination of its xy-coordinates (e.g. 12: column 1, row 2).

<img src = "img/mea-layout.png"></img>

### Meta file

### Import MEA data

### Apply filters
- sample filter: This filter excludes these recordings for technical reasons, e.g., contamination. The user has marked these samples to be excluded in the meta file in the column exclude. 
- channel filter: This filter excludes particular electrodes for technical reasons, e.g., too noisy. The user has marked these channels in the meta file in the column "channels2exclude". Multiple channels are separated by a comma.
- recording length filter: This filter takes into account the max recording time to take into account. The user indicates the max recording time for all samples.
- active channel filter: This filter only keeps channels that are active. We define active as having at least a spike rate of 0.01 Hz.

- custom filters: region, layer, etc
### Burst detection
Burst are detected based on user-defined parameters. As default, we have set that a burst is identified as contain at least 5 spikes, with each being max 50 ms apart from each other. Other parameter includes the max time frame these spikes have to be in.

### Spike features

### Burst features

### Pipeline
A template processing pipeline is provided in `mea_pipeline.R`.
