## Diluted fecal community transplant restores *Clostridioides difficile* colonization resistance to antibiotic perturbed murine communities

### Abstract

Fecal communities transplanted into individuals can eliminate recurrent *Clostridioides difficile* infection (CDI) with high efficacy. However, this treatment is only used once CDI becomes resistant to antibiotics or has recurred multiple times. We sought to investigate whether a fecal community transplant (FCT) pre-treatment could be used to prevent CDI altogether. We treated male C57BL/6 mice with either clindamycin, cefoperazone, or streptomycin, and then inoculated them with the microbial community from untreated mice before challenging with *C. difficile*. We measured colonization and sequenced the V4 region of the 16S rRNA gene to understand the dynamics of the murine fecal community in response to the FCT and *C. difficile* challenge. Clindamycin-treated mice became colonized with *C. difficile* but cleared it naturally and did not benefit from the FCT. Cefoperazone-treated mice became colonized by *C. difficile*, but the FCT enabled clearance of *C. difficile*. In streptomycin-treated mice, the FCT was able to prevent *C. difficile* from colonizing. Then we diluted the FCT and repeated the experiments. Cefoperazone-treated mice no longer cleared *C. difficile*. However, streptomycin-treated mice colonized with 1:10^2^ dilutions resisted *C. difficile* colonization. Streptomycin-treated mice that received a FCT diluted 1:10^3^, *C. difficile* colonized but later was cleared. In streptomycin-treated mice, inhibition of *C. difficile* was associated with increased relative abundance of a group of bacteria related to *Porphyromonadaceae* and *Lachnospiraceae*. These data demonstrate that *C. difficile* colonization resistance can be restored to a susceptible community with a FCT as long as it complements the missing populations.  
  
  
### Importance  
  
Antibiotic use, ubiquitous with the healthcare environment, is a major risk factor for *Clostridioides difficile* infection (CDI), the most common nosocomial infection. When *C. difficile* becomes resistant to antibiotics, a fecal microbiota transplant from a healthy individual can effectively restore the gut bacterial community and eliminate the infection. While this relationship between the gut bacteria and CDI is well established, there are no therapies to treat a perturbed gut community to prevent CDI. This study explored the potential of restoring colonization resistance to antibiotic-induced susceptible gut communities. We described the effect gut bacteria community variation has on the effectiveness of a fecal community transplant for inhibiting CDI. These data demonstrated that communities susceptible to CDI can be supplemented with fecal communities but the effectiveness depended on the structure of the community following the perturbation. Thus, a reduced bacterial community may be able to recover colonization resistance to patients treated with antibiotics.  
  
  
### Overview

	project
	|- README         # the top level description of content (this doc)
	|- CONTRIBUTING   # instructions for how to contribute to your project
	|- LICENSE        # the license for this project
	|
	|- submission/	  # files for manuscript submission
	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered
	| |- mothur/      # mothur processed data
	| +- process/     # cleaned data, will not be altered once created
	|
	|- code/          # any programmatic code
	| |- ml/          # R scripts to run machine learning model and process data
	| |- mothur/      # dir for local copy of mothur
	|
	|- results        # all output from workflows and analyses
	| |- tables/      # text version of tables to be rendered with kable in R
	| |- figures/     # graphs, likely designated for manuscript figures
	| +- pictures/    # diagrams, images, and other non-graph graphics
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
	| +- scratch/     # temporary files that can be safely deleted or lost
	|
	+- Makefile       # executable Makefile for this study, if applicable


### How to regenerate this repository

#### Dependencies and locations
* Gnu Make should be located in the user's PATH
* mothur (v1.XX.0) should be located in the user's PATH
* R (v. 3.X.X) should be located in the user's PATH
* etc.


#### Running analysis

```
git clone https://github.com/SchlossLab/Lesniak_Severity_mBio_2022.git
make submission/manuscript.pdf
```