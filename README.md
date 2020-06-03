[![HitCount](http://hits.dwyl.com/synaptic-proteolab/CAIR_EMIP.svg)](http://hits.dwyl.com/synaptic-proteolab/CAIR_EMIP)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-red.svg?style=flat)](https://github.com/synaptic-proteolab/CAIR_EMIP/pulls)
[![GitHub issues](https://img.shields.io/github/issues/synaptic-proteolab/CAIR_EMIP.svg)](https://github.com/synaptic-proteolab/CAIR_EMIP/issues)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<img src="https://img.shields.io/badge/Study%20Status-Results%20Available-yellow.svg" alt="Study Status: Results Available"> 



<img src="Logo/Proteolab logo.jpg" width=1100/>



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
  * [Install requirements](#Install-requirements)  
  * [Run from terminal](#run-from-terminal)  
  * [Run from python](#run-from-python)  
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
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/synaptic-proteolab/CAIR_EMIP/blob/master/ShannoProt/runCAIR/runCAIR.ipynb)

CAIRs of each given protein is calculated using the following function of the shanon entropy:  

<img src="img/CAIR_Eq.JPG" height="150"> 

## CAIR comparison analyses
Violin plots of each evolutionary node and their corresponding mean comparison tests of the brunner-munzel:  
<img src="Figures/Figure1.jpg" height="1400"> 

## The EMIP algorithm
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/synaptic-proteolab/CAIR_EMIP/blob/master/ShannoProt/runEMIP/Human_disease_category_indicators.ipynb)

Estimation of mutual information for a protein (EMIP) is calculated using the following function:  
<img src="img/EMIP_Eq.JPG" height="250"> 

## EMIP comparison analyses
xxxxxxxxxxxxx:  
<img src="Figures/Figure3.jpg" height="650"> 


## Usage

### Install requirements
```bash
pip install requests biopython numpy pandas import_ipynb
```

#### Run from terminal
CAIR_EMIP can be executed directly from terminal:


#### Run from python
Fire up your python shell or ipython notebook. 


## Toy data

The example data files to run all the functions are available here. Of note, these are just small subsets of species and human dataset. We provided these ToyData so that the user can test the method.

As an example you can re-calculate the Proteome CAIRs for the 16 mentioned organisms in table 1 of article. The required sample data files is in the ToyData folder. The program will count all the residues in ToyData FASTA files, group them by organism IDs and calculate CAIRs for each complete protome.

However, if you plan to run the complete algorithm, you should download all the required files (~100GB for CAIR+ ~100MB for EMIP).

1. CAIR Project
### Prerequisites

All required packages to run the scripts can be installed from terminal/CMD using:
```sh
pip3 install requests biopython numpy pandas import_ipynb
```

### Run from:

#### Windows

1. Clone the repo
```
git clone https://github.com/synaptic-proteolab/CAIR_EMIP.git
```
2. reset the currnet directory to user directory
```
cd %HOMEDRIVE%%HOMEPATH%
```
3. Change running directory to the project environment
```
cd ShannoProt/ToyData/runCAIR
```
4. Run the code file
```py
python runCAIR.py
or
python3 runCAIR.py
```
#### Linux
1. Clone the repo
```sh
git clone https://github.com/synaptic-proteolab/CAIR_EMIP.git
``` 

2. reset the currnet directory to user directory
```sh
cd ~
```
3. Change running directory to the project environment
```sh
cd ShannoProt/ToyData/runCAIR
```
4. Run the code file
```sh
python runCAIR.py
or
python3 runCAIR.py
```
2. EMIP Project
### Prerequisites

All required packages to run the scripts can be installed from terminal/CMD using:
```sh
pip3 install requests numpy pandas import_ipynb
```

### Run from:

#### Windows

1. Clone the repo
```
git clone https://github.com/synaptic-proteolab/CAIR_EMIP.git
```
2. reset the currnet directory to user directory
```
cd %HOMEDRIVE%%HOMEPATH%
```
3. Change running directory to the project environment
```
cd ShannoProt/ToyData/runEMIP
```
4. Run the code file
```py
python Human_disease_category_indicators.py
or
python3 Human_disease_category_indicators.py
```
#### Linux
1. Clone the repo
```sh
git clone https://github.com/synaptic-proteolab/CAIR_EMIP.git
``` 

2. reset the currnet directory to user directory
```sh
cd ~
```
3. Change running directory to the project environment
```sh
cd ShannoProt/ToyData/runEMIP
```
4. Run the code file
```sh
python Human_disease_category_indicators.py
or
python3 Human_disease_category_indicators.py
```

Reference Paper
-------
Not published yet

License
-------
© Copyright 2020 Synaptic Proteolab. Licensed under the MIT License. See LICENSE file for more details.
