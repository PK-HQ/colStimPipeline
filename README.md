# colStimPipeline

MATLAB pipeline for designing and analyzing optogenetic stimulation of 
orientation columns in macaque V1, used to elicit feature-specific 
visual percepts in behaving animals.

## What this does
- Generates column-targeted light patterns from widefield orientation maps
- Coregisters stimulation bitmaps to cortical imaging space
- Analyzes behavioral (psychometric) and neural responses to optogenetic stimulation
- Fits psychometric models to quantify perceptual effects under stimulation vs. control

## Key scripts
- `mainColStimPipeline.m` — top-level pipeline entry point
- `getColumnarBitmapV5.m` — generates stimulation bitmaps targeting specific columns
- `analyzePsychometrics.m` — psychometric analysis of stimulation effects
- `analyzeRecruitment.m` — quantifies network recruitment beyond stimulated site

## Context
Core analysis code for Tan et al. (in prep), Seidemann Lab, UT Austin.
