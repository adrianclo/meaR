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
Settings for spike detection threshold are done within the Multi Channel Experimenter Spike Analyzer itself. Data output files are in .txt format. These contain for each separate electrode (N=60) spike data on the occurrence (s), interspike interval (isi in ms) and frequency (in Hz). The electrode ID is the combination of its xy-coordinates.

MEA channel layout:

     21  31  41  51  61  71 <br>
 12  22  32  42  52  62  72  82 <br>
 13  23  33  43  53  63  73  83 <br>
 14  24  34  44  54  64  74  84 <br>
 15  25  35  45  55  65  75  85 <br>
 16  26  36  46  56  66  76  86 <br>
 17  27  37  47  57  67  77  87 <br>
     28  38  48  58  68  78  

### Meta file

### Import MEA data

### Apply filters
- region filter:
- layer filter:
- recording length filter:
- active channel filter:

### Burst detection
Burst are detected based on user-defined parameters. As default, we have set that a burst is identified as contain at least 5 spikes, with each being max 50 ms apart from each other. Other parameter includes the max time frame these spikes have to be in.

2
### Spike features

### Burst features

- Pipeline to analyze data are provided in `mea_pipeline.R`.
