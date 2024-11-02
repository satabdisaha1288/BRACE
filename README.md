Bayesian compositional regression with flexible microbiome feature aggregation and selection
# Authors
Satabdi Saha, Liangliang Zhang, Kim Ahn Do, Christine B. Peterson
# BRACE
The paper "Flexible aggregation of compositional predictors with shared effects" presents a novel method, Bayesian Regression with Agglomerated Compositional Effects (BRACE) for rare feature clustering and variable selection in high dimensional compositional regression models. 

Ongoing advances in microbiome profiling have allowed unprecedented insights into the molecular activities of microbial communities. This has fueled a strong scientific interest in understanding the critical role the microbiome plays in governing human health, by identifying microbial features associated with clinical outcomes of interest. Several aspects of microbiome data limit the applicability of existing variable selection approaches. In particular, microbiome data are high-dimensional, extremely sparse, and compositional. Importantly, many of the observed features, although categorized as different taxa, may play related functional roles. To address these challenges, BRACE leverages the data-adaptive clustering and variable selection properties of the spiked Dirichlet process to identify taxa that exhibit similar functional roles, enabling the identification of a sparse set of features with shared impacts on the outcome and facilitating dimension reduction and model interpretation.

# Required Functions
This directory has all the functions required for running scripts Simulation.R and RunBRACElet.R

# Other files
* Simulation.R simulates data for an example case study.
* RunBRACElet.R runs BRACElet for the simulated data example with compositional predictors and a continuous response.
* ProcessOutputs.R processes the outputs returned by RunBRACElet.R and calculates the performance metrics.


