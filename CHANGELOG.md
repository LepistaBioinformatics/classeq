## v0.13.2 (2023-12-28)

### Fix

- include the response type for the models endpoint

## v0.13.1 (2023-12-10)

### Fix

- replace name selection from models to use the model id instead

### Refactor

- move manuscript figures to a dedicated directory

## v0.13.0 (2023-12-10)

### Feat

- implements a simple server to predict sequences from the web

## v0.12.4 (2023-11-02)

### Fix

- remove unused kmer-inverse-indices from prediction use-cases

## v0.12.3 (2023-10-29)

### Fix

- decrease the matches-coverage allowed lower value to allow predict over genome sequences

## v0.12.2 (2023-10-22)

### Fix

- reduce the minimum clade size to allow low leaf sister evaluation

## v0.12.1 (2023-10-22)

### Fix

- fix clades id generation that create repeaded uuids for simbling clades

## v0.12.0 (2023-10-14)

### Feat

- include the options to user ignore the usage of outgroups during predicion

## v0.11.1 (2023-10-14)

### Fix

- include a search break to avoid classify sequences with low similarity to train set

## v0.11.0 (2023-10-13)

### Feat

- create a cli command to sanitize phylogenies

## v0.10.2 (2023-10-13)

### Fix

- fix app to use the new pyside6 and replace tree root to use ete3

## v0.10.1 (2023-10-13)

### Fix

- replace middlepoint reroot from biopython features to use ete3 functionality

### Refactor

- remove unused code from the main cli port

## v0.10.0 (2023-10-08)

### Feat

- upgrade indexing and predictiong to include outgroup from out of the phylogeny

## v0.9.0 (2023-09-03)

### Feat

- implements a simple command to generate kmers from command line directely

## v0.8.1 (2023-09-02)

### Fix

- fix strand both as default option

## v0.8.0 (2023-09-02)

### Feat

- improve filtration of the predictions results during prediction phase to include new arguments
- include the minimum query match size at the clades prediciton

### Fix

- include conclusive-ingroup classification status as a break rule into predictions

## v0.7.0 (2023-08-30)

### Feat

- configure default values for the output predictions
- include the full prediciton results at the output table and json files

## v0.6.2 (2023-08-29)

### Fix

- improve the tree sanitization and prediction use-cases to allow index and predict using trees with non polithomic outgroups

## v0.6.1 (2023-08-28)

### Fix

- fix use-case that extends hierarchical tree over the original size

## v0.6.0 (2023-08-28)

### Feat

- create the strand attribute of reference-sets to control kmer generation and turn the tree clades ids as predictable uuid3

## v0.5.0 (2023-08-27)

### Feat

- upgrade project to work with tar balls as input instead to use individual gz files

## v0.4.2 (2023-08-26)

### Fix

- remove unused commands from cli port

## v0.4.1 (2023-08-26)

### Fix

- fix ui

## v0.4.0 (2023-08-26)

### Feat

- upgrade enums to use ordered-tuple as tuple manager default and update prediction use-cases to better performance

### Fix

- remove unused code of joint probability calculate

### Refactor

- remove joint probability strategy from application

## v0.3.6 (2023-08-22)

### Fix

- fix url of images that compose user-interface

## v0.3.5 (2023-08-22)

### Fix

- fix error on manage directories using cwd and refatorate predict use cases to modularize sub-use-cases

## v0.3.4 (2023-08-21)

### Fix

- fix the search by uuids of clades

## v0.3.3 (2023-08-21)

### Fix

- include the center option to focus on the first element found into a single search

## v0.3.2 (2023-08-21)

### Fix

- fix the nodes search to allow subsearches and case insensitive searches
- include clear selection option on search nodes into treeview of app

## v0.3.1 (2023-08-21)

### Fix

- fix user interface issues to improve the user experience

## v0.3.0 (2023-08-20)

### Feat

- include prediction final status into the user output file

## v0.2.1 (2023-08-20)

### Refactor

- remove trash temp file

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
