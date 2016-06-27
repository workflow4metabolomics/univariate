# Univariate parametric and non-parametric hypothesis testing with correction for multiple testing

A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) project.

Status: [![Build Status](https://travis-ci.org/workflow4metabolomics/univariate.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/univariate).

## Description

**Version:** 2.1.1   
**Date:** 2015-09-30  
**Author:** Marie Tremblay-Franco (INRA, MetaToul, MetaboHUB, W4M Core Development Team) and Etienne A. Thevenot (CEA, LIST, MetaboHUB, W4M Core Development Team)    
**Email:** [marie.tremblay-franco(at)toulouse.inra.fr](mailto:marie.tremblay-franco@toulouse.inra.fr); [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:** Thevenot E.A., Roux A., Xu Y., Ezan E. and Junot C. (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *Journal of Proteome Research*, **14**:3322-3335. [doi:10.1021/acs.jproteome.5b00354](http://dx.doi.org/10.1021/acs.jproteome.5b00354)  
**Reference history:** [W4M00001a_sacurine-subset-statistics](http://galaxy.workflow4metabolomics.org/history/list_published), [W4M00004_mtbls1](http://galaxy.workflow4metabolomics.org/history/list_published)  
**Licence:** CeCILL  
**Funding:** Agence Nationale de la Recherche ([MetaboHUB](http://www.metabohub.fr/index.php?lang=en&Itemid=473) national infrastructure for metabolomics and fluxomics, ANR-11-INBS-0010 grant)

## Installation

 * Configuration file: **univariate_config.xml**  
 * Image file: 
  + **static/images/univariate_workflowPositionImage.png**   
 * Wrapper file: **univariate_wrapper.R**  
 * Script file: **univariate_script.R**  
 * R packages  
  + **batch** from CRAN: `install.packages("batch", dep=TRUE)`.
  + **PMCMR** from Bioconductor: `install.packages("PMCMR", dep=TRUE)`.

## Tests

The code in the wrapper can be tested by running the **tests/univariate_tests.R** in R  

## News

## CHANGES IN VERSION 2.1.1  

 * Internal handling of 'NA' p-values (e.g. when intensities are identical in all samples).
