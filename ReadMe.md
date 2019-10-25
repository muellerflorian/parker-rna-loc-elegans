<a href="https://imjoy.io" target="_blank" ><img src="https://muellerflorian.github.io/parker-rna-loc-elegans/img/imjoy-logo-powered.png" width="140"></img>
</a>
![GitHub](https://img.shields.io/github/license/muellerflorian/parker-rna-loc-elegans)

# Analyze RNA localization in c. elegans
We provide an analyze workflow to analyze **RNA localization patterns** in c. elegans 
in single molecule FISH (smFISH images).

Analysis is performed in Matlab and in Python scripts running in ImJoy.

<img src="https://muellerflorian.github.io/parker-rna-loc-elegans/img/dist_transform.png" width="400px"></img>

- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Instructions for use](#instructions-for-use)
- [License](#license)

# Documentation
Here we provide only a brief overview. The officiel documentation with usage is at: 
<a href="https://muellerflorian.github.io/parker-rna-loc-elegans/" target="_blank">**https://muellerflorian.github.io/parker-rna-loc-elegans/.**</a>

# System Requirements
## Hardware requirements
`parker-rna-loc-elegans` package requires only a standard computer.

## Software requirements

### OS Requirements
This package has been tested on *macOS* on the following systems:
+ macOS: High Sierra (10.13.6) on a Mac Pro (late 2010)

### Matlab dependencies
RNA detection is performed with our prevously published Matlab package 
<a href="https://bitbucket.org/muellerflorian/fish_quant" target="_blank">**FISH-quant.**</a>

FISH-quant requires the following **toolboxes**:
* Optimization toolbox
* Statistics toolbox
* Image processing toolbox
* (Optional) Parallel processing toolbox

FISH-quant has been tested on **Matlab 2017b**.

### ImJoy
ImJoy plugin were tested on **ImJoy (v0.9.93)** running on **Google Chrome (Version 77)** with the **ImJoy plugin engine (0.8.22)**.

More information on ImJoy: [https://arxiv.org/abs/1905.13105](https://arxiv.org/abs/1905.13105)

### Python Dependencies
`parker-rna-loc-elegans` mainly depends on the Python scientific stack, a few smaller libraries 
are required to read FIJI region definition files. 

Complete list of dependencies can be found in 
[here](https://github.com/muellerflorian/parker-rna-loc-elegans/blob/master/setup.py)

# Installation guide

## FISH-quant
You can obtain the latest version of FISH-quant
<a href="https://bitbucket.org/muellerflorian/fish_quant" target="_blank">**here.**</a>

Installation instruction are provided. Installation time is rapid and requires only
to download the most recent version.

## ImJoy
The provided workflow requires the installation of the ImJoy plugin engine, which
has to be installed only once, which takes several minutes.

Installation links to all plugins are provided and installation time for each plugin
is in the range of several minutes for the first install.

# Instructions for use

## Demo data
We provide an example data-set, with detailed instructions for how
data has to be organised. Then we detail in step-by-step instructions how to
perform the analysis.

Demo data also contains the expected results, which are obtained
in a few minutes.

## New data
The provided instructions for the demo data are directly applicable for new data.

# License
[MIT License](https://github.com/muellerflorian/parker-rna-loc-elegans/blob/master/LICENSE)
