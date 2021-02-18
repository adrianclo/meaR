# meaR

## Introduction
Multi-electrode arrays (MEAs) are becoming increasingly popular to investigate activity of neuronal assemblies upon genetic, pharmacological and/or optogenetic manipulation. Neuronal activity can be assessed on acute brain slices, neuronal uniform or mixed-cell cultures, and human pluripotent stem cells, to name a few. It is a means for high-throughput screening and these MEAs come in different sizes, resolution and presets from different companies. However, MEA data in general are very complex, and the output files generated from these MEAs are rarely immediately usable for statistical analysis. As such, a tabular form is more recommended which should also include the parameters required to conduct a statistical analysis.

## Purpose
For these reasons, this repository contains a collection of functions to extract, visualize and analyze MEA data, specifically from <a href = "https://www.multichannelsystems.com/products/vitro-mea-systems">MultiChannel Systems</a>. <br>
<b>Note</b>: The script is in continuous development and earlier versions have been used for the following papers:

- <b>Lo AC</b>, Rajan N, Gastaldo D, Telley L, Hilal ML, Buzzi A, Simonato M, Achsel T, Bagni C (accepted). Absence of RNA-binding protein FXR2P prevents prolonged phase of kainate-induced seizures. <i>EMBO Rep</i>

- Dominguez-Iturza N, <b>Lo AC</b>, Shah D, Armendariz M, Vannelli A, Mercaldo V, Trusel M, Li KW, Gastaldo D, Santos AR, Callaerts-Vegh Z, D'Hooge R, Mameli M, Van der Linden A, Smit AB, Achsel T, Bagni C (2019). The autism- and schizophrenia-associated protein CYFIP1 regulates bilateral brain connectivity and behaviour. <i>Nat Commun</i> <b>10</b>: 3454. <a href = "https://www.nature.com/articles/s41467-019-11203-y">[link to article]</a>

## Guide
- Pipeline to analyze data are provided in `mea_pipeline.R`.
