# [New Zealand Parkinson's Epidemiology Research](http://nzbri.org/Labs/parkinsons/Epidemiology/)

Any queries please contact:

Daniel Myall <daniel.myall@nzbri.org>

New Zealand Brain Research Institute

The data used in this research can not be released as it contains identifiable information. However, the code fully documents the methodology used in the analysis. If you are interested in a fake dataset please let me know.

## Preprocessing of medication and diagnosis data

A lot of this code is specific to our particular data sources.

### National Minimal Dataset

[python/nmds.py](python/nmds.py)

Processes hospital admissions and mortality data to extract diagnoses

### Diagnoses

[python/diagnoses.py](python/diagnoses.py)

Combine diagnoses from multiple sources

### Prescription data

[python/pharmacdata.py](python/pharmacdata.py)

Reads raw prescription data, basic processing, saves into a database, export into format for drug-based classification.

### Drug-based classification algorithm

[python/process.py](python/process.py)

Classifies individuals as very probable, probable, possible, or unlikely.


## Statistical Code

Requires R http://www.r-project.org/ and Stan http://www.mc-stan.org/

### Common base statistical code:

[R/model-base.R](R/model-base.R)

### Parkinson's ethnic analysis:

[R/model-ethnicity.R](R/model-ethnicity.R)

Stan model: [R/pd_epi_model_ethnicity_v2.stan](R/pd_epi_model_ethnicity_v2.stan)

## Census data

Contained in directory [input](input/) loaded by [R/model-base.R](R/model-base.R)
