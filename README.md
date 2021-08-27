# Project-Plasmodium
Repository for my EEB397 project to understand the implications of plastic gametocyte conversion rate on artemisinin treatment for malaria parasites.

Malaria is one of the deadliest infectious disease in the world, leading to 409,000 deaths in 2019 alone. In recent decades, global malaria burden has drastically declined, in part, due to greater accessibility of artemisinin-based combination treatment (ACT).  However, this encouraging trend is threatened by the emergence of artemisinin-resistant Plasmodium falciparum strains in Southeast Asia.  Understanding the mechanism of antimalarial resistance, especially to the artemisinins, is pivotal to the ultimate goal of malaria elimination.

Antimalarial resistance can occur via classical routes (e.g. efflux pump) or non-classical routes. The latter could involve alterations to the parasite’s life-history  traits, such as conversion rate. Conversion rate represents the proportion of merozoite that commits to gametocyte production and is shown to vary in response to different external stimuli. Lower conversion rate could entail higher within-host parasite density, which is associated with increased resistance to pyrimethamine and artemisinin in cell cultures. This project aims to better understand the mechanism that underlies conversion rate-mediated antimalarial resistance, with specific focus on artemisinin and its derivatives. Mathematical modelling will be used to derive the optimal conversion rate reaction norm in response to artemisinin treatment. Interactions between drug dosage, competing strains, and host immune response will also be examined.

## Progress so far
1. Streamlined optimization strategy search and improved run time drastically by implementing parallel computing and vectorizing sub-functions.
    * Allows sensitivity analysis to be conducted using plastic conversion rate strategy rather than static strategy.
    * Configurable code allows each aspect of modelling process (DDE solver, parameter sets, degrees of freedom...) to be altered with minimum manipulation. 
2. Incorporated option to allow conversion rate to depend on states (etc. density of infected RBC) rather than time. 
3. Incorporated several modes of innate immunity:
    * Greischar's model for saturating immunity (targeted infected RBC removal)
    * Kochin's model for saturating immunity (targeted infected RBC removal + explicit modelling of effector cell population)
    * Kamiya's model for innate immunity (targeted infected RBC removal, indiscriminant removal of RBC)

## Next steps
1. Creating configurable models that simulate realistic within-host parasite dynamics and can rapidly derive the optimal conversion rate strategy.
    * Incorporate adaptive immunity into the model. 
    * Allowing cues to be based on derivative of states (e.g. change in infected RBC density).
    * Further optimizing sensitivity analysis run time. 
2. Testing the effects of modelling parameters (e.g. degrees of freedom, cue range) on optimal conversion rate. 


## Currently working on:
Creating a series of functions for non-age structured single infection and co-infection (without drugs for now) models where conversion rate depends on specific "cues" rather than time. Modelling conversion rate as a plastic trait allows for more realistic infection dynamics. This cue-based approach prevents parasites from "anticipating" drug administration/co-infection and relaxes the assumption that cue correlates with time post-infection (especially useful when we consider that drug adherence rate might be low in certain areas). 

The main focus for the next few days is to incorporate adaptive immune response into the model to simulate more realistic infection dynamics.
