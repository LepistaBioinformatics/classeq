## v0.2.0 (2023-08-20)

### Feat

- implements the output file persistence after predictions
- implements the annotation port to update phylogeny nodes on runtime
- implements the phylojson parsing of tree to allow annotations

### Fix

- fix calculation of the adherence test and start the creation of the app port

## v0.1.0 (2023-07-22)

### Feat

- implements the full adherence test with outgroup checking
- implements the full adherence test with outgroup checkint
- wip - implements the multi child adherence test calculation
- wip - start the implementation of the refactored prediction of clades based on ingroup clades elements only
- replace default left returns by the new mappederrors left return using call class methods
- replace native elements of clean arch from elements from chean-base package
- wip - start the implementation of the code to turn possible the manual check of the phylogenetic adherence tests
- resolve deprecations from typings and implements the results-set hierarchical tree building
- implements the hierarchical tree creation and viewing methods
- include linear trees persistence into loading tree step
- update single clade adherence test calculation use-case to allow train to ingroups and outgroups
- implements the priors parsing from json file
- wip - implements the inclusion of groups into the training step to be used during predictions
- remove estimate_global_kmer_specific_priors and calculate_single_kmer_likelihood from use-cases bacause it was not longer used
- update main train use-case function to include clade specific train and their results persistence to file
- implements a specific use case to calculate priors of all tree clades
- wip - implements the conditional probilties calculation for single clades
- include a custom format at the logging initialization
- move the likelihood calculation funciton to a dedicated use-case to allow reuse of then
- include a extended node type check to include outgroup checks
- include an aditional node type for outgroup nodes
- update tree to force monophily of the outgroup nodes
- simplify the data loading process centralizing it to the dtos initialization cycle and create a new clade wrapper
- implements the reverse index generator to be used during probabilities calculations
- implements all methods for kmer indices generation and create methods for the msa initial parsing
- implement the basis for the bayesian calculation of the kmers posterior probabilities of the groups
- implements the fasta sequences loading use-cases and update msa and tree wrapper dtos
- mapping use-case main folders to the project
- initialize kmer-generation and models train use-cases files
- implements the first use-case to load phylogeny into memory
- initial commit

### Fix

- fix error raised on use trees with outgroups as the root tree
- wip - upgrade elements of recursive prediction test
- fix inverted either type hints of all project
- remove import of estimate_global_kmer_specific_priors from the train_from_single_phylogeny use-case test
- upgrade the hash indices storage to use immutable tuple instead sets

### Refactor

- move the calculate-recursive-priors to a dedicated use-case file
- mark all internal sub-use-cases files with underscore to identify their internal use
- move dto tests to a dedicated directory inside the dtos folder
