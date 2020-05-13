# The phylogenetic range of bacterial and viral pathogens of vertebrates

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3823900.svg)](https://doi.org/10.5281/zenodo.3823900)

**Associated paper**: The phylogenetic range of bacterial and viral pathogens of vertebrates, *Molecular Ecology* (2020)  
**doi**: [10.1111/mec.15463](https://doi.org/10.1111/mec.15463)  

**Authors**: Liam P. Shaw\*, Alethea Wang\*, David Dylus, Magda Meier, Grega Pogacnik, Christophe Dessimoz, Francois Balloux (\* co-first authors)
![alt text][header]

[header]: data/pathogen-host-range-header-image.png

This repository contains code to reproduce all figures and analyses in the associated paper analysing the host range of pathogen species. This is bundled together in the Rmarkdown notebook (Supplementary Text 1 of the manuscript). It also contains several data files. The image above is a schematic for the manual literature review used to construct the database of host-pathogen associations. For more information on the dataset collection methods and analysis, please see the paper. 

Code can be altered to change any of the figures.  Some explanatory text is provided, but not much. Any questions can be addressed via email (liam.philip.shaw at gmail dot com). 

Unless otherwise stated, data related to viruses is plotted in <span style="color:red">red</span> and bacteria in black. The code should take <10 minutes to run on a standard laptop. 

**Note:** At several points code is adapted/reused from the fantastic supplementary code repository made available by Olival et. al. (2017) for an analysis of viral pathogens:

* Olival et al. (2017) paper: *Host and viral traits predict zoonotic spillover from mammals*  https://doi.org/10.1038/nature22975 
* Code repository: https://github.com/ecohealthalliance/HP3

This mainly applies to several scripts in the `scripts' directory relating to the fitting and plotting of generalized additive models (GAMs) -- more information is included where appropriate. All code here is also made available under the MIT License. 
