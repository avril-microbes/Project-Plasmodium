# Project-Plasmodium
Repository for my EEB397 project to understand the implications of plastic gametocyte conversion rate on artemisinin treatment for malaria parasites.

## Currently working on:
Creating a series of functions for non-age structured single infection and co-infection (without drugs for now) models where conversion rate depends on specific "cues" rather than time. Modelling conversion rate as a plastic trait allows for more realistic infection dynamics. This cue-based approach prevents parasites from "anticipating" drug administration/co-infection and relaxes the assumption that cue correlates with time post-infection (especially useful when we consider that drug adherence rate might be low in certain areas). 

## Function folder
Contains various functions that allows models to be customized and built rapidly. 
### Model functions
Naming of these types of functions follow the convention:
1. chabaudi = Plasmodium chabaudi-based 
2. si = single infection
3. ci = coinfection
4. opt = optimization function to obtain optimal conversion rate strategy

## Journal folder
This folder contains various Rmarkdown files where I run various analysis to test the functions and get results.
