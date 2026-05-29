# AGENTS.md — colStimPipeline

## Project identity

This repository is `colStimPipeline`, a MATLAB research codebase for designing optogenetic light patterns targeting cortical columns and analyzing behavior/neurophysiology from macaque V1 visual + optogenetic stimulation experiments.

The main project path on the lab server is:

`\\172.17.49.6\data\users\PK\colStimPipeline`

A mapped Windows path may also exist:

`Y:\users\PK\colStimPipeline`

Prefer the UNC path if the mapped `Y:` drive is not visible.

## Primary language

MATLAB.

Relevant file types:

- `.m`
- `.mlx`
- `.fig`
- `.mat`
- `.csv`
- `.xlsx`
- `.md`
- `.txt`

## Safety rules

This is an active research repository. Treat data integrity as critical.

Do not delete, overwrite, rename, or move raw data files.

Do not modify large data files unless explicitly instructed, especially:

- `.mat`
- `.tif`
- `.tiff`
- `.bmp`
- `.png`
- `.jpg`
- `.avi`
- `.mp4`
- binary acquisition files
- raw session folders
- generated imaging datasets
- behavioral export files

Before editing code, first inspect relevant files and explain the intended change.

Prefer small, reviewable edits.

Do not make broad refactors unless explicitly requested.

Do not change numerical methods, fitting logic, coordinate transforms, or trial inclusion criteria unless explicitly requested.

Do not run full-dataset analyses unless explicitly requested.

For testing, use the smallest safe example first, ideally one file, one block, or one session.

If a script writes outputs, confirm the output path before running it.

Never commit, push, pull, merge, or reset Git history unless explicitly instructed.

## Git behavior

This repository appears to be a Git repository on branch `main`.

If Git is available, check status before any edit:

```bash
git status