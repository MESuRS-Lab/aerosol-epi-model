
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

## Accounting for the long-distance transmission route: an epidemiological model of airborne disease transmission in hospitals 

### Context 

Respiratory infections are a major public health issue in hospital settings. They can be transmitted from patient to patient, as well as between patients and healthcare workers, resulting in nosocomial or hospital-acquired infections. These infections contribute to patient aggravation and complexify their care management.

Pathogens can be transmitted in different ways:
- During close-proximity interactions at short distances;
- Through aerosols emitted during breathing, speech, or coughing which can remain suspended in the air and cause long-distance transmission.

A thorough understanding of the respective roles of these transmission modes is essential for developing effective strategies to prevent nosocomial respiratory infections. To date, no model simultaneously integrates these two modes of transmission, and the relative importance of each remains to be evaluated.

The aim of this project is to develop a mathematical model of respiratory pathogen transmission in a hospital ward, combining inter-individual transmission during close-proximity interactions and long-range transmission due to persistent aerosols.


## Instructions

### Dependencies
<details>
  <summary> R </summary>
  
  * cowplot
  * data.table
  * DescTools
  * doParallel
  * extrafont
  * extrafontdb
  * foreach
  * gganimate
  * ggrepel
  * gifski
  * ggnetwork
  * ggpattern
  * ggpubr
  * ggtext
  * gt
  * gtsummary
  * hms
  * igraph
  * performance
  * rcompanion
  * Rcpp
  * rstatix
  * tidyverse
  * viridis

</details>

<details>
  <summary> c++ </summary>

  * C++ 11

  </details>


### Nods-Cov-2 scripts (R/nodscov2/)
The project is divided into numbered R scripts available in the `R/nodscov2/` directory. The scripts can be used as follows:

<details>
  <summary> 1.compare_synthetic_observed_contacts.R : compare network metrics and structure in the observed Nods-Cov-2 data and in the synthetic temporal networks. </summary>
  This script loads Nods-Cov-2 local data (private), filters it to keep only intensive care units, then does an analysis of the distribution of individuals' categories, type of interactions etc. An analysis on healthcare worker cumulative time spent interacting with patients is also done.
</details>

<details>
  <summary> 2.trim-synthetic-contacts.R : trim interactions in synthetic networks that do not respect the presence of individuals in the ward. </summary>
</details>

<details>
  <summary> 3.compare-truncated-observed-contacts.R : comparison of observed and synthetic networks after truncation of interactions. </summary>
</details>

<details>
  <summary> 4.localization-reconstruction-parallel.R : reconstruct locations in the observed and synthetic temporal networks. </summary>
    This script loads Nods-Cov-2 local data (private) and reconstructs the locations of individuals using interactions, individual category and comportemental rules.
</details>

<details>
  <summary> 5.localization-reconstruction-figures.R : plot reconstructed locations for observed and synthetic contact networks </summary>
</details>

<details>
  <summary> 6.generate-gif.R : code to generate a movie of individual movements in ICU2 for over one day. </summary>
</details>

<details>
  <summary> 7.generate-input-data-model.R : create object that contains data and parameters that are used as inputs by the simulation workflow. </summary>
</details>

<details>
  <summary> 8.grid-search.R : analysis of simulations from the grid search approach. </summary>
</details>

<details>
  <summary> 9.R0-calculation.R : summarize R_0 values across transmission conditions. </summary>
</details>

<details>
  <summary> 10.baseline-scenario.R : plot results from the baseline scenario without intervention. </summary>
</details>

<details>
  <summary> 11.intervention-analysis.R : plot results from scenarios with control interventions. </summary>
</details>

<details>
  <summary> 12.sensitivity-analysis.R : plot results from the univariate sensitivity analysis. </summary>
</details>

<details>
  <summary> 13.influenza.R : plot results from the analysis using influenza virus-like parameter values. </summary>
</details>

# Generate simulations

Simulations were generated using a Nextflow pipeline. Files to launch the pipelines are in the `nextflow/grid_search`, `nextflow/interventions`, `nextflow_sensitivity` and `nextflow_influenza` directories and R scripts that generate unique simulations are available in the `bin` subfolders. 

Epidemic simulations can be launched using the Rcpp code `dev-sensibility-analysis.cpp` in the `cpp` folder. Importantly, `dev-model-nodscov2.cpp` is the same code but also returns the trajectories of aerosol concentration in ward rooms at each time step (very demanding in terms of memory). 

## License

Distributed under the GNU General Public License. See `COPYING` for more information.


## Contact

Project Link: [https://github.com/MESuRS-Lab/aerosol-epi-model](https://github.com/MESuRS-Lab/aerosol-epi-model/tree/main)


<!-- MARKDOWN LINKS & IMAGES -->
[contributors-shield]: https://img.shields.io/github/contributors/nollive/projet-pacri.svg?style=for-the-badge
[contributors-url]: https://github.com/nollive/projet-pacri/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/nollive/projet-pacri.svg?style=for-the-badge
[forks-url]: https://github.com/nollive/projet-pacri/network/members
[stars-shield]: https://img.shields.io/github/stars/nollive/projet-pacri.svg?style=for-the-badge
[stars-url]: https://github.com/nollive/projet-pacri/stargazers
[license-shield]: https://img.shields.io/github/license/nollive/projet-pacri.svg?style=for-the-badge
[license-url]: https://github.com/nollive/projet-pacri/blob/main/COPYING
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://fr.linkedin.com/in/olivier-gaufr%C3%A8s
