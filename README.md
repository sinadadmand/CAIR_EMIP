[![HitCount](http://hits.dwyl.com/synaptic-proteolab/CAIR_EMIP.svg)](http://hits.dwyl.com/synaptic-proteolab/CAIR_EMIP)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-red.svg?style=flat)](https://github.com/synaptic-proteolab/CAIR_EMIP/pulls)
[![GitHub issues](https://img.shields.io/github/issues/synaptic-proteolab/CAIR_EMIP.svg)](https://github.com/synaptic-proteolab/CAIR_EMIP/issues)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


<img src="Logo/Proteolab logo.jpg" width=970/>



[Website](https://www.synaptic.one/)



# An extension of Shannon’s entropy agrees with the tree of life and justifies human disease occurrences #

## Description
This Project involves the python and R implementation of the above entitled project.

## Table of Contents
* [Links to literature](#Links-to-literature)
* [The CAIR algorithm](#CAIR-algorithm)  
* [CAIR comparison analyses](#CAIR-comparison-analyses)  
* [The EMIP algorithm](#EMIP-algorithm)  
* [EMIP comparison analyses](#EMIP-comparison-analyses)  
* [Usage](#Usage)  
* [Toy data](#Toy-data)  


## Links to literature 

* **runCAIR** (Calculating the CAIR of all proteomes)  
_Prepared manuscript._  
Python codes: [runCAIR](https://github.com/synaptic-proteolab/CAIR_EMIP/tree/master/ShannoProt/runCAIR).  
* **CAIRcomp** (Proteomes CAIR comparisons through the tree of life)  
_Prepared manuscript._  
R codes: [CAIR Comparisons](https://github.com/synaptic-proteolab/CAIR_EMIP/tree/master/R_codes/CAIR_comparisons).  
* **runEMIP** (Human disease category indicators -- EMIP, etc.)  
_Prepared manuscript_  
Python codes: [runEMIP](https://github.com/synaptic-proteolab/CAIR_EMIP/tree/master/ShannoProt/runEMIP).  
* **EMIPanalys** (Analysis of the relation between disease occurances and calculated(gathered) parameters in human)  
_Prepared manuscript_  
R codes: [EMIP analyses](https://github.com/synaptic-proteolab/CAIR_EMIP/tree/master/R_codes/EMIP_analyses).

## The CAIR algorithm
CAIRs of each given protein is calculated using the following function of the shanon entropy:  

<img src="img/CAIR_Eq.JPG" height="260"> 

## CAIR comparison analyses
Violin plots of each evolutionary node and their corresponding mean comparison tests of the brunner-munzel:  
<img src="Figures/Figure1.jpg" height="1600"> 

## The EMIP algorithm
Estimation of mutual information for a protein (EMIP) is calculated using the following function:  
<img src="img/EMIP_Eq.JPG" height="150"> 

## EMIP comparison analyses
xxxxxxxxxxxxx:  
<img src="Figures/xxxxxxxxx.jpg" height="300"> 


## Usage


## Toy data
xxxxxxxxxxxxxxxx.



The easiest way to run CAIR-EMIP Python codes is by typing these lines in cmd/terminal:

for windows:
for linux:
.  
.  
.  

```

git clone https://github.com/synaptic-proteolab/CAIR_EMIP.git
cd ShannoProt/ToyData/runCAIR
python runCAIR.py
python 3 runCAIR.py
```
```
git clone https://github.com/synaptic-proteolab/CAIR_EMIP.git
cd ShannoProt/ToyData/runEMIP
python Human_disease_category_indicators.py
python 3 Human_disease_category_indicators.py
```

.
.
.

.
```
Reference Paper
-------
Not published yet

License
-------
© Copyright 2020 Synaptic Proteolab. Licensed under the MIT License. See LICENSE file for more details.
