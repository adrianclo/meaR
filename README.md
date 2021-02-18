# meaR

## Introduction
Multi-electrode arrays (MEAs) are becoming increasingly popular to investigate activity of neuronal assemblies upon genetic or pharmacological manipulation. These can be performed on acute brain slices, neuronal cultures or human pluripotent stem cells. It is a means for high-throughput screening. These MEAs come in different sizes and presets from different companies. However, MEA data are very complex, and the output files from these MEAs can often not immediately be used for analysis.

## Purpose
For the aforementioned reason, this repository contains a collection of functions to extract and analyze data from the multi-electrode setup (MEA) from <a href = "https://www.multichannelsystems.com/products/vitro-mea-systems">MultiChannel Systems</a>. However,

<b>Note</b>: The script is in continuous development and earlier versions have been used for the following papers:

- <b>Lo AC</b>, Rajan N, Gastaldo D, Telley L, Hilal ML, Buzzi A, Simonato M, Achsel T, Bagni C (accepted). Absence of RNA-binding protein FXR2P prevents prolonged phase of kainate-induced seizures. <i>EMBO Rep</i>

- Dominguez-Iturza N, <b>Lo AC</b>, Shah D, Armendariz M, Vannelli A, Mercaldo V, Trusel M, Li KW, Gastaldo D, Santos AR, Callaerts-Vegh Z, D'Hooge R, Mameli M, Van der Linden A, Smit AB, Achsel T, Bagni C (2019). The autism- and schizophrenia-associated protein CYFIP1 regulates bilateral brain connectivity and behaviour. <i>Nat Commun</i> <b>10</b>: 3454. <i>doi: 10.1038/s41467-019-11203-y</i> <a href = "https://www.nature.com/articles/s41467-019-11203-y">[link to article]</a>

## Guide
- Pipeline to analyze data are provided in `mea_pipeline.R`.
