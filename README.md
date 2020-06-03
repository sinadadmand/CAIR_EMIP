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
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/synaptic-proteolab/CAIR_EMIP/blob/master/ShannoProt/runCAIR/runCAIR.ipynb)

Violin plots of each evolutionary node and their corresponding mean comparison tests of the brunner-munzel:  
<img src="Figures/Figure1.jpg" height="1400"> 

## The EMIP algorithm
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/synaptic-proteolab/CAIR_EMIP/blob/master/ShannoProt/runEMIP/Human_disease_category_indicators.ipynb)

Estimation of mutual information for a protein (EMIP) is calculated using the following function:  
<img src="img/EMIP_Eq.JPG" height="250"> 

## EMIP comparison analyses
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/synaptic-proteolab/CAIR_EMIP/blob/master/R_codes/)

xxxxxxxxxxxxx:  
<img src="Figures/Figure3.jpg" height="650"> 


## Usage

### Install requirements
```bash
pip install requests numpy pandas biopython import_ipynb
```

#### Run from terminal
PyPuma can be run directly from the terminal with the following options:
```

```
To run PyPuma on the included Toy example:
```
python run_puma.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -i ToyData/ToyMiRList.txt -o output_puma.txt
```
To run LIONESS on PUMA networks, use the flag -q (note that this can take a long time and use considerable computing resources):

```python
python run_puma.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -i ToyData/ToyMiRList.txt -o output_puma.txt -q output_lioness.txt
```
Finally, note that running PUMA without importing motif and expression data will estimate a co-expression network using Pearson correlation.

#### Run from python
Fire up your python shell or ipython notebook. 
Import the classes in the PyPuma library:
```python
from pypuma.puma import Puma
from pypuma.lioness import Lioness
```
Then run PUMA:
```python
puma_obj = Puma('ToyData/ToyExpressionData.txt', 'ToyData/ToyMotifData.txt', 'ToyData/ToyPPIData.txt','ToyData/ToyMiRList.txt')
```
Save the results:
```python
puma_obj.save_puma_results('Toy_Puma.pairs.txt')
```
Example of returning a network visualization of the top edges:

```python
puma_obj.top_network_plot(top=70, file='top_genes.png')
```
<!--
or
```python
from PyPuma.analyze_puma import AnalyzePuma
plot = AnalyzePuma(puma_obj)
plot.top_network_plot(top=100, file='top_100_genes.png')
```-->
Calculate indegrees for further analysis:
```python
indegree = puma_obj.return_puma_indegree()
```


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
