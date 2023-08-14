# Research analyses

This repository is meant to host analysis templates for the various research projects that I have worked on.

I will try to provide the necessary input data where I can, but some projects may be too sensitive/far from publishing that I cannot make the data freely available.

Additionally, as projects become published, I will post all the necessary analyses here.

## Model RNA mods XGBoost

I work on many RNA modification projects. RNA modidifications can be predicted from regular RNA-seq with algorithms like [HAMR](https://github.com/GregoryLab/HAMR) and [ModTect](https://github.com/ktan8/ModTect). 

ModTect detects substantially more modifications relative to HAMR, but HAMR was build to classify RNA modifications using a reverse transcriptase error profile modeled on Yeast tRNA modifications. I wanted to the scale of RNA modifications predicted by ModTect but have the RNA modification class that HAMR reports. This script uses XGBoost to predict ModTect modification classes from HAMR output. Briefly, the important bits of data to model are the reference nucleotide of the position, the distribution of observed nucleotides at that position, and the relative mutation ratio of the position. Input data is present in this repo. The analysis HTML can be accessed 
