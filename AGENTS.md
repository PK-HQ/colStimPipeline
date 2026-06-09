# AGENTS.md — colStimPipeline

## Project identity

This repository is `colStimPipeline`, a MATLAB research codebase for macaque V1 visual and optogenetic stimulation experiments.

The project supports behavioral, psychometric, and neurometric analyses from experiments involving visual orientation discrimination and optogenetic stimulation of cortical orientation columns.

Primary language: MATLAB.

Primary project path on the lab server:

`\\172.17.49.6\data\users\PK\colStimPipeline`

A mapped Windows path may also exist:

`Y:\users\PK\colStimPipeline`

Prefer the UNC path if the mapped `Y:` drive is not visible.

## Default local task permissions

For normal code-maintenance tasks, assume the following unless PK says otherwise:

Allowed without asking:
- read files in this repository
- search with `rg`
- inspect MATLAB code
- run `git status`
- run `git diff`
- edit files explicitly named in the task
- run small no-save smoke tests when the task requires verification

Not allowed without asking:
- edit files not named in the task
- edit raw data or generated data files
- edit `.mat`, `.xlsx`, image, video, or binary files
- run full-dataset analyses
- commit, push, pull, merge, reset, or delete branches
- change scientific definitions, fitting logic, condition labels, or trial inclusion criteria
- write output files to data folders

When asking for permission, batch requests. Do not ask command-by-command for routine read/search/diff/test operations within the allowed scope.

## Main entry point

The main entry point for current work is:

`mainPipeline.m`

Do not assume `mainColStimPipeline.m` is the current entry point. It may still contain useful historical information, especially about light-pattern design and data structure, but current analysis work should be understood through `mainPipeline.m` unless instructed otherwise.

## Project documentation to consult

Before making behavioral-analysis changes, read:

- `docs/psychometrics.md`
- `docs/DATA_DICTIONARY.md`
- `docs/STRUCT_FIELD_MAP.md`

If these files conflict with older code comments or deprecated analysis modes, ask before relying on the older code.

## Current trusted analysis modes

The following analysis modes are considered clean/current:

* `expt`
* `summary`
* `psyclusterPre`
* `psycluster`
* `psyphidist`

Other analysis modes may be old, partially deprecated, or not in the correct current format. They may contain useful clues about data structure, helper functions, or older logic, but do not base major conclusions or edits on them without asking first.

If a task seems to require using older modes or files outside the trusted list, ask before relying on them.

## Current behavioral analysis workflow

The current behavioral analysis workflow is:

`summary -> psyclusterPre -> psycluster`

This sequence is run for each animal and each chamber.

Current focus: behavioral analysis, including generating plots and saving behavioral data needed for those plots.

Neurometric analysis may be developed later, but behavioral analysis is the active priority unless otherwise specified.

## Safety rules

This is an active research repository. Treat data integrity as critical.

Do not delete, overwrite, rename, or move raw data files.

You have all read permissions, do not stop to ask for approval for read permissions which wastes time.

Do not modify large data files unless explicitly instructed, especially:

* `.mat`
* `.tif`
* `.tiff`
* `.bmp`
* `.png`
* `.jpg`
* `.avi`
* `.mp4`
* binary acquisition files
* raw session folders
* generated imaging datasets
* behavioral export files

Before editing code, first inspect relevant files and explain the intended change.

Prefer small, reviewable edits.

Do not make broad refactors unless explicitly requested.

Do not change numerical methods, fitting logic, coordinate transforms, trial inclusion criteria, saved output fields, or condition definitions unless explicitly requested.

Do not run full-dataset analyses unless explicitly requested.

For testing, use the smallest safe example first, ideally one file, one block, one animal/chamber, or one analysis mode.

If a script writes outputs, confirm the output path before running it.

Never commit, push, pull, merge, or reset Git history unless explicitly instructed.

## Git behavior

This repository is Git-tracked.

Before edits, check status:

```bash
git status
```

After edits, summarize changed files and show the diff when possible:

```bash
git diff
```

Do not commit unless explicitly instructed.

Do not push unless explicitly instructed.

If Git is not available in the current shell environment, say so and do not pretend a diff was checked.

## Expected workflow for Codex

For every task:

1. Inspect relevant files first.
2. Summarize what the existing code appears to do.
3. Identify the smallest safe change.
4. Explain the plan before editing if the change could affect analysis outputs.
5. Make focused edits only in files relevant to the task.
6. Summarize exact files changed.
7. Provide a minimal verification step.

For read-only tasks:

1. Do not modify files.
2. Do not create files.
3. Do not run MATLAB unless explicitly asked.
4. Do not process raw data.
5. Report uncertainty clearly.

## Documentation and code style

Write documentation and comments in a normal human research-code style.

Do not add comments, commit messages, or documentation that sound like they were generated by an AI assistant.

Avoid phrases such as:

* “This function has been updated to...”
* “The assistant changed...”
* “Codex added...”
* “Generated by AI”
* “As requested...”
* overly polished tutorial-style explanations inside code

Prefer concise technical comments that explain why something is done, especially when the logic is non-obvious.

Good comment style:

```matlab
% Keep congruent and incongruent summaries separate for signed contrasts.
```

Avoid comment style:

```matlab
% This section was added to improve clarity and robustness for the user.
```

Do not over-comment obvious MATLAB operations.

Preserve the existing coding style unless there is a strong reason to change it.

## MATLAB conventions

Prefer clear, explicit MATLAB code over clever compact code.

Avoid changing function signatures unless asked.

Avoid changing output struct field names unless asked.

Use robust checks for missing fields, empty arrays, and different MATLAB object classes.

When dealing with transforms, do not assume only one class. Possible transform classes may include:

* `affine2d`
* `affinetform2d`
* `images.geotrans.PolynomialTransformation2D`

## Scientific context

The project supports macaque V1 experiments involving visual stimulation and optogenetic stimulation of cortical orientation columns.

Core scientific concepts include:

* visual orientation discrimination
* V1 orientation columns
* optogenetic stimulation
* congruent vs incongruent optostim
* behavioral bias
* psychometric curves
* neurometric/imaging responses
* widefield calcium imaging
* light-pattern design
* projector/camera coordinate transforms
* stimulation power and power density
* number of stimulated columns
* task-axis vs off-axis signals

## Common analysis concepts

Important behavioral analysis outputs may include:

* behavioral performance by condition
* congruent vs incongruent psychometric curves
* Naka-Rushton or related psychometric fits
* signed contrast plots
* condition-averaged response summaries
* saved behavioral summary structs
* saved plotted values
* animal/chamber-specific summaries

Important optostim/imaging analysis outputs may include:

* optostim power and power density summaries
* column counts
* bitmap/projector/camera masks
* condition-level summaries
* imaging response maps
* trial-by-trial metadata tables

## Data structure notes

Common structures may include:

* `behavioralData`
* `behavioralData.optoTS(blockID).Header.ConditionParams.Stimulus`
* `bitmapData`
* `bitmapData.gridSize`
* `bitmapData.gammaCorrFactor`
* `bitmapData.sensitivity`
* `bitmapData.adaptthresh`
* `bitmapData.orts`
* `bitmapData.transformParams`

Expected visual metadata may include:

* `vis.contrast`
* `vis.ort`
* `vis.size`
* `vis.sf`
* `vis.pos`

Expected optostim metadata may include:

* `opto.480LED`
* `opto.480ND`
* `opto.580LED`
* `opto.580ND`
* `opto.PD`
* `opto.gridsize`
* `opto.gamma`
* `opto.threshSens`
* `opto.threshAdapt`
* `opto.ort`
* `opto.gausscond`
* `opto.gausslevel`
* `opto.gaussmax`
* `opto.ortmaskarea`
* `opto.gaussmaskarea`
* `opto.transformParams`
* `opto.bitmapProj`
* `opto.bitmapCam`
* `opto.nColumns`
* `opto.pixON`
* `opto.timeON`
* `opto.energy`

## Debugging rules

When debugging:

1. Find the smallest reproducible script/function.
2. Identify the exact error line.
3. Explain the likely cause.
4. Patch minimally.
5. Suggest a minimal test.
6. Avoid full-dataset execution unless asked.

## Output rules

When producing explanations, include:

* exact file names,
* relevant function names,
* relevant line numbers when available,
* assumptions,
* what was not verified.

When editing code, include:

* files changed,
* reason for change,
* expected behavior,
* how to verify,
* possible risks.

Keep summaries direct and technical. Do not write in a style that advertises use of Codex, GPT, or AI-generated assistance.
