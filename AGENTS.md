# AGENTS.md — colStimPipeline

# Agent Execution Rules

- Never edit MATLAB source using PowerShell string replacement, regex replacement,
  here-strings, ReadAllLines/WriteAllLines function-region replacement, or equivalent
  whole-block textual surgery.

- Use apply_patch or the native editor patch mechanism only.

- If apply_patch cannot access the network repository, stop and provide the proposed
  patch to the user. Do not invent a PowerShell workaround.

- After one failed edit command, do not attempt a second editing mechanism without
  explicit user approval.

- A request for “surgical edits only” means no helper-region replacement spanning
  more than 20 lines.

- Static inspection may use at most three targeted file-reading commands before
  returning control to the user.
  
## Priority: fast interactive iteration

This repository is on a slow network/UNC filesystem and MATLAB startup is expensive. Optimize for short edit-test cycles and avoid unattended validation loops.

### MATLAB execution

1. Do not launch MATLAB, `matlab.exe -batch`, `checkcode`, or a full pipeline run unless the user explicitly requests an automated run.

2. Prefer having the user run code in their existing interactive MATLAB session and return the first error.

3. When an automated MATLAB run is explicitly requested:

   * Run at most one MATLAB batch process per response.
   * Use a maximum timeout of 120 seconds unless the user explicitly authorizes longer.
   * If it times out, terminate only the batch MATLAB process that the agent started.
   * Do not wait an additional grace period.
   * Do not launch a second MATLAB process while the first is still alive.
   * Report the timeout and return control to the user.

4. Never terminate all MATLAB processes. The user may have an important interactive MATLAB session open.

5. For reloading edited MATLAB functions, recommend:

```matlab
clear functions
```

Do not recommend `clear all` unless it is specifically required.

6. After a runtime failure, fix only the first demonstrated error. Do not repeatedly relaunch MATLAB to discover successive errors unless the user explicitly requests that loop.

### Network filesystem

1. Avoid repeated `Test-Path`, `Get-Content`, or recursive searches against the UNC repository.

2. Perform at most one filesystem preflight.

3. Prefer targeted searches with `rg`, `git grep`, or `git diff` over broad recursive PowerShell commands.

4. Do not retry an unavailable mapped drive after the UNC repository is confirmed accessible.

5. Group related file reads into one command where practical.

### Editing

1. Prefer direct patch tools or small, deterministic edits.

2. Do not use large PowerShell here-strings, giant regex replacements, or full-function textual replacement for MATLAB files.

3. Do not bypass PowerShell execution policy merely to apply a patch.

4. If a helper or function requires extensive replacement, create a temporary local file and replace it using a simple deterministic script, or make smaller edits.

5. After editing, run only fast static checks:

   * `git diff --check`
   * Search for stale field names.
   * Search for duplicate local function names.
   * Inspect the edited regions.
   * Verify MATLAB `function`/`end` structure.

6. Do not treat a successful static check as scientific validation.

### Existing pipeline behavior

1. Reuse existing loaders, condition-indexing functions, transformations, and filtering helpers.

2. Do not create a parallel custom loader or condition parser until the existing implementation has been traced and shown to be insufficient.

3. Preserve existing conventions, including all equivalent projector labels such as:

   * `O000` and `O045` for the 0° target family.
   * `O090` and `O135` for the 90° target family.
   * Existing no-opto/baseline conventions such as `L000` where applicable.

4. Prefer `getUsableTrials` and established condition mappings over reimplementing condition parsing.

5. Do not refactor unrelated code during an urgent analysis task.

### User communication

1. Before starting a command likely to take more than 60 seconds, state exactly what will run and why.

2. If a command stalls or times out, stop and report it immediately.

3. Do not spend multiple minutes fixing shell quoting, PowerShell formatting, or agent-environment problems. Return a command for the user to run manually instead.

4. Clearly distinguish:

   * Static checks completed.
   * MATLAB execution completed.
   * Scientific QC inspected.
   * Items not yet verified.

5. When the user has an end-of-day deadline, prioritize:

   * Working core analysis.
   * First QC figures.
   * Clear runtime errors.
   * Minimal changes.

Defer optional refactoring, trial-level extensions, and exhaustive validation.

## Highest-priority operating contract

This repository is an active MATLAB research codebase. Preserve scientific correctness while keeping routine maintenance fast.

Unless PK explicitly requests broader execution, the default behavior is:

**Inspect narrowly, make the smallest correct patch, self-review it, perform static checks, and stop.**

Do not turn a local request into a full pipeline validation, architectural redesign, environment-debugging session, or broad audit.

The user normally runs the final MATLAB analysis and visually checks generated outputs.

A higher reasoning-effort setting permits more careful reasoning. It does not by itself grant permission to run broader analyses, edit more files, regenerate PDFs, or troubleshoot the MATLAB environment.

User instructions in the current task override this file.

---

## Network-share preflight and failure behavior

The repository is stored on a mounted Windows network share.

Repository paths:

* Interactive mapped path: `Y:\users\PK\colStimPipeline`
* Canonical UNC path: `\\172.17.49.6\data\users\PK\colStimPipeline`

Before editing files or launching MATLAB, perform one lightweight access check:

1. Test whether the mapped repository path is readable.
2. If `Y:` is unavailable, test the UNC repository path once.
3. If either path is accessible, use the accessible path and continue.
4. If both paths are inaccessible, stop immediately and report that the server/network share appears disconnected.

Do not keep retrying unavailable paths.

### Distinguish failure types

If the error indicates a missing or disconnected share, including messages such as:

* path not found;
* drive not found;
* network path not found;
* network name no longer available;
* specified network name is no longer available;
* device unavailable;
* connection reset;
* too many open files after failed network operations;

then:

* stop all work immediately;
* do not edit files;
* do not launch MATLAB;
* do not try alternate shells or repeated reconnect commands;
* tell the user the exact path and error so they can restore the server connection.

If the repository is reachable but access is denied by the sandbox or permissions:

* request one narrowly scoped approval for read/write access only beneath:
  `Y:\users\PK\colStimPipeline`
  or its UNC equivalent;
* do not request write access to animal-data directories;
* if approval is denied, stop and report the exact blocked operation.

### Data-directory protection

Treat these locations as strictly read-only:

* `Y:\Chip`
* `Y:\Pepper`
* `\\172.17.49.6\data\Chip`
* `\\172.17.49.6\data\Pepper`

Never create, modify, overwrite, rename, move, or delete files beneath those directories.

### Process behavior

* Do not launch MATLAB until repository access has passed the preflight check.
* Do not launch multiple MATLAB processes.
* If network access disappears during a task, stop immediately rather than continuing with partial files or repeated retries.
* Never claim completion when the network share became inaccessible before validation.

## Task classes and time budgets

Classify the task silently before acting.

### QUICK task

Examples:

* labels
* colors
* titles
* line styles
* mean ticks
* annotations
* display formatting
* metadata printed on a figure
* small plotting bugs
* localized source-field corrections
* a clearly bounded bug in one function

Expected behavior:

* target completion: approximately 5–10 minutes
* inspect no more than the directly relevant implementation and immediate call sites
* normally touch no more than 1–4 files
* do not launch MATLAB
* do not run the pipeline
* do not regenerate PDFs
* use static verification only
* stop after the patch and concise report

If a QUICK task cannot be safely completed within this scope, stop and state the exact dependency or blocker. Do not expand it into a larger task without permission.

### STANDARD task

Examples:

* a contained new analysis
* a multi-file behavioral change
* adding a reusable export
* a localized data-processing change
* a targeted debugging task with a reproducible failure

Expected behavior:

* target completion: approximately 15–30 minutes
* focused repository inspection
* one targeted smoke test may be used when genuinely useful
* do not run a full animal/chamber pipeline unless explicitly requested
* do not troubleshoot unrelated environment problems
* stop if the task grows into a substantial refactor or full validation campaign

### DEEP task

Examples:

* a new multi-chamber analysis
* a substantial statistical analysis
* pipeline restructuring
* a difficult data-integrity audit
* a broad refactor
* end-to-end validation requested by PK

Expected behavior:

* may take approximately 30–60 minutes
* only enter this mode when the request is clearly substantial or PK explicitly asks for deep/full validation
* remain within the requested scientific scope
* do not retry indefinitely when execution fails

### Budget overrun rule

If the task exceeds its expected scope or time:

1. Stop running commands.
2. Preserve the working tree.
3. Report what is complete.
4. State the precise blocker.
5. Do not keep troubleshooting autonomously.

Do not spend most of a task on environment problems, broad searches, repeated validation attempts, or command retries.

---

## Project identity

This repository is `colStimPipeline`, a MATLAB research codebase for macaque V1 visual and optogenetic stimulation experiments.

Primary language: MATLAB.

Primary project path:

`\\172.17.49.6\data\users\PK\colStimPipeline`

Mapped path:

`Y:\users\PK\colStimPipeline`

Prefer the UNC path only when the mapped drive is unavailable. Do not repeatedly switch between UNC and mapped paths while debugging an unrelated task.

---

## Default task permissions

### Allowed without asking

* read relevant repository files
* use narrowly scoped `rg` searches
* inspect immediate callers and callees
* run `git status`
* run `git diff`
* run `git diff --check`
* run `git diff --stat`
* edit files explicitly named in the task
* edit directly necessary immediate call sites
* add one small helper file when clearly justified
* perform a small, no-save, deterministic smoke test when the task is STANDARD and the test is expected to finish quickly

### Not allowed without explicit instruction

* run `mainPipeline.m`
* run `summary`, `psyclusterPre`, `psycluster`, or another full analysis mode
* regenerate full PDFs
* run all chambers or animals
* launch long MATLAB batch jobs
* modify raw or generated data
* write output into data folders
* commit
* push
* pull
* merge
* reset
* rebase
* create or delete branches
* alter the active animal, chamber, model, analysis mode, or loop settings
* change scientific definitions
* change fitting logic
* change statistical tests
* change cluster definitions or membership
* change condition labels
* change trial inclusion rules
* change contrast matching
* perform broad refactors

Do not ask command-by-command for routine allowed reads, searches, diffs, or edits.

---

## Default QUICK-task workflow

For a QUICK task:

1. Run `git status` once.
2. Inspect the named function or file.
3. Inspect only immediate call sites needed to understand the data flow.
4. Identify the existing canonical source of the requested value.
5. Make the smallest patch.
6. Review the complete diff yourself.
7. Check all call sites if a function signature changed.
8. Run:

   * `git diff --check`
   * `git diff --stat`
9. Stop and report.

Do not:

* write a long plan before editing
* narrate every search command
* repeatedly reread the same files
* scan the whole repository
* launch MATLAB
* regenerate output
* perform end-to-end validation
* commit or push

Only explain the plan before editing when the requested change affects numerical analysis, data inclusion, statistics, fitting, clustering, or saved outputs.

---

## Validation levels

Never escalate validation levels without explicit permission.

### Level 1 — Static verification

This is the default.

Includes:

* review the complete diff
* `git diff --check`
* inspect changed function signatures
* inspect all relevant call sites
* verify variable dimensions and indexing logically
* verify that canonical source fields are used
* check for silent fallback behavior
* check that unrelated run configuration was not changed

For most plotting and metadata tasks, Level 1 is sufficient because PK performs final MATLAB and visual validation.

### Level 2 — Targeted smoke test

Use only when:

* PK explicitly requests a test, or
* the task is STANDARD and a narrow test is clearly valuable

Requirements:

* test only the changed function or a small synthetic example
* no full dataset
* no PDF regeneration unless explicitly requested
* no raw-data writes
* expected runtime should be short
* one failed attempt is enough to stop if failure is environmental

Do not turn a failed smoke test into MATLAB installation or path troubleshooting.

### Level 3 — End-to-end validation

Use only when PK explicitly asks for language such as:

* run the full analysis
* regenerate the PDFs
* validate both chambers
* perform end-to-end validation
* run full `psycluster`
* audit the complete output

Level 3 is never implied merely because a code change affects a figure.

---

## MATLAB execution and environment rules

For QUICK tasks, do not launch MATLAB.

For STANDARD or DEEP tasks, launch MATLAB only when the requested validation requires it.

Never autonomously:

* run `restoredefaultpath`
* run `rehash toolboxcache`
* modify `startup.m`
* modify MATLAB preferences
* add the entire repository recursively with `addpath(genpath(...))`
* create temporary modified copies of `mainPipeline.m`
* alter run configuration to force another chamber or animal
* kill MATLAB processes
* troubleshoot MATLAB installation
* troubleshoot licenses
* repair toolbox paths
* repeatedly retry failed batch launches

If a MATLAB command fails because of path, toolbox, startup, mapped-drive, license, or environment problems:

1. Stop.
2. Report the exact error.
3. Leave the repository unchanged beyond the requested patch.
4. Let PK run the code in the normal configured MATLAB session.

A clear typo in the command may be corrected once. Environment failures should not be retried autonomously.

---

## Network-drive rules

This repository resides on a network drive.

Avoid expensive recursive operations.

For QUICK tasks:

* search specific files or narrow directories
* use exact symbols and function names
* avoid repository-wide `Get-ChildItem -Recurse`
* avoid repeated whole-tree `rg`
* avoid repeatedly reading large files in full
* do not scan raw-data directories

If a narrow search does not locate the implementation, report that rather than launching a broad multi-minute crawl.

---

## Scientific source-mapping rules

Scientific metadata and plotted values must come from the canonical analysis source already used by the pipeline.

### Prefer direct indexing

When the plotting loop already has a stable block or experiment index and the source vector is block-indexed, use direct indexing.

Example principle:

```matlab
value = canonicalVector(blockIdx);
```

Do not invent string-based date/run matching when a direct canonical index exists.

Do not create alternate mappings merely because field names differ across helper functions.

### No silent fallback

Never silently substitute:

* `N/A`
* zero
* empty text
* the previous experiment’s value
* the first matching row
* a guessed default
* a reconstructed label

when an expected scientific value cannot be mapped.

For eligible experiments with a missing canonical value:

* render or report `UNASSIGNED` only when that is explicitly useful
* emit a clear warning
* preserve the missing state
* do not misclassify it as intentionally excluded

Use `N/A` only when the canonical inclusion/exclusion mask explicitly says that the experiment is not applicable.

Never default an experiment to ineligible because optional metadata is absent.

### Canonical mapping checks

Before adding a scientific metadata label:

1. Identify the exact canonical source.
2. Identify its indexing convention.
3. Trace one representative experiment through the current loop.
4. Confirm that aggregate/distribution code uses the same assignment.
5. Check expected counts where possible.

For cluster labels, verify that the number of experiments assigned to each cluster matches the number used by the corresponding aggregate or distribution analysis.

### Pipeline-order dependency

If the requested value does not exist at the current plotting stage:

1. First look for an existing finalized saved field or context already passed downstream.
2. Prefer passing the canonical value forward with minimal changes.
3. Do not redesign pipeline order.
4. Do not add a second full rendering pass.
5. Do not regenerate or append PDFs differently without first telling PK why that is necessary.

If a local display request truly requires pipeline restructuring, stop and explain the dependency before implementing it.

---

## Current trusted scientific invariants

Do not change these unless PK explicitly requests a scientific change.

### Behavioral workflow

Current sequence:

`summary -> psyclusterPre -> psycluster`

Run separately for each animal and chamber.

### Power quantities

Spatial duty cycle:

```matlab
sDC = AreaON / AreaROI;
```

ROI-averaged power density:

```matlab
PDROI = PDDMD * sDC * tDC;
```

Units: `mW/mm^2`

Total delivered power:

```matlab
Ptotal = PDDMD * AreaON * tDC;
```

Units: `mW`

Current power clustering uses:

```matlab
Ptotal
```

It does not use `PDROI` as the clustering variable.

All calculations use full-precision numeric values. Rounding is display-only.

### Baseline mode

The current baseline-mode source is `baselineTS`.

* filled/nonempty `baselineTS` means separate baseline file
* empty/missing baseline file means combined baseline/opto block

### Delta definitions

Merged, horizontal, and vertical delta values use:

```matlab
deltaBias = con - incon;
```

```matlab
deltaMask = baseline - mean([con, incon]);
```

Zero contrast only matches exact zero.

Nonzero contrast matching uses a maximum discrepancy/spread of 5 contrast percentage points unless explicitly changed.

### Power-cluster assignment

For current power-cluster analyses, use the finalized assignment already produced by the clustering workflow.

When `powerEffectClusterByBlock` is the finalized block-indexed assignment in the active code path, map using the exact `blockIdx` for the plotted experiment.

Do not recompute cluster membership inside plotting functions.

Do not infer membership from rounded power boundaries.

---

## Working-tree and Git safety

Run `git status` once before editing.

Assume the working tree may already contain important uncommitted changes.

Never:

* reset
* restore unrelated files
* checkout over user changes
* stash without instruction
* rewrite run configuration
* overwrite unrelated edits
* clean untracked files
* commit or push without instruction

Directly relevant immediate call-site edits are allowed even when they were not explicitly named, but keep scope narrow.

Review the diff yourself after editing. PK generally does not manually inspect diffs, so the responsibility for checking the patch remains with Codex.

Do not dump a large diff into the final response unless requested.

If pre-existing changes are mixed into a touched file:

* preserve them
* distinguish them from the new patch in the summary
* do not attempt to clean them up

---

## Main entry point

The active entry point is:

`mainPipeline.m`

Do not assume `mainColStimPipeline.m` is current.

Older entry points and deprecated modes may be consulted for historical context, but do not copy old behavior into the current path without verifying it.

Do not change the top-level run configuration in `mainPipeline.m` unless the task explicitly asks for it.

---

## Trusted analysis modes

Current trusted modes:

* `expt`
* `summary`
* `psyclusterPre`
* `psycluster`
* `psyphidist`

Other modes may be deprecated or inconsistent.

Do not base a major change on an older mode without checking first.

---

## Documentation usage

Relevant documentation:

* `docs/psychometrics.md`
* `docs/DATA_DICTIONARY.md`
* `docs/STRUCT_FIELD_MAP.md`

For QUICK display-only changes, do not reread all documentation unless the needed field or invariant is unclear.

For changes affecting:

* calculations
* statistics
* fitting
* trial inclusion
* source fields
* saved outputs
* clustering
* scientific definitions

consult the relevant documentation before editing.

If documentation and active trusted code conflict, report the conflict rather than guessing.

---

## Code-editing principles

Prefer:

* minimal patches
* direct data flow
* existing helpers
* existing coding style
* clear MATLAB code
* explicit indexing
* fail-loud behavior for scientific metadata
* preserving existing function signatures

Avoid:

* speculative abstractions
* broad helper frameworks for one small task
* duplicated calculations
* parallel sources of truth
* unnecessary optional arguments
* fragile string matching
* silent exception swallowing
* broad refactors
* second-pass rendering
* unrelated cleanup

Do not change a function signature unless needed. If changed, inspect every call site before stopping.

Do not add defensive fallbacks that hide incorrect data flow.

---

## MATLAB style

Use clear, explicit MATLAB code.

Check:

* empty arrays
* missing fields
* finite values
* row/column orientation
* table versus struct access
* stable experiment/block indexing
* MATLAB version compatibility

Preserve existing coding style.

Avoid over-commenting obvious code.

Comments should explain non-obvious scientific or indexing logic.

Good:

```matlab
% Use the finalized block-indexed cluster assignment.
```

Avoid:

```matlab
% This code was added by Codex to improve robustness.
```

Never add comments referring to Codex, AI, prompts, or the user’s request.

---

## Debugging rules

When debugging:

1. Locate the exact failing output or line.
2. Trace the value back to its canonical source.
3. Reproduce with the smallest available example.
4. Patch the root cause minimally.
5. Perform static verification.
6. Suggest one focused runtime check for PK.

Do not begin with a full-dataset rerun.

Do not infer that an audit passed merely because a mathematical identity is internally consistent. Confirm that the correct source fields and display formatting are used when those are part of the reported problem.

---

## Progress-reporting rules

Do not narrate routine work command by command.

For QUICK tasks, remain quiet until:

* the patch is complete, or
* a genuine blocker requires PK’s decision

Do not send repeated messages such as:

* “I am now locating…”
* “Next I will inspect…”
* “I found another file…”
* “I am going to retry…”

For longer STANDARD or DEEP tasks, give a brief update only when there is meaningful progress or a blocker.

---

## Stop conditions

Stop immediately and report rather than continuing when:

* a QUICK task requires a pipeline redesign
* more than four files require meaningful edits for a supposedly local change
* the canonical scientific source cannot be identified
* direct indexing and current mappings disagree
* a MATLAB environment command fails
* validation requires changing run configuration
* a test would write unexpected generated output
* the task requires raw-data modification
* an unrelated bug is discovered
* a command is taking substantially longer than expected
* repeated searches are not narrowing the problem
* the requested change cannot be verified statically and runtime validation was not authorized

Do not autonomously turn the blocker into a new task.

---

## Expected final response

### QUICK task response

Keep it concise:

* files changed
* exact behavior changed
* self-review/static checks performed
* whether MATLAB was run
* one short manual verification step, if useful

Example structure:

```text
Changed:
- fileA.m: ...
- fileB.m: ...

Checks:
- reviewed full diff
- git diff --check passed
- MATLAB/full pipeline not run

Verify:
- run the normal analysis and inspect ...
```

### STANDARD or DEEP task response

Include:

* files changed
* scientific logic implemented
* tests actually run
* outputs actually generated
* unresolved risks or assumptions
* what was not verified

Never claim a run, test, visual check, or source validation occurred when it did not.

Do not instruct PK to inspect the diff as the primary verification. Codex must self-review the diff; PK performs the final scientific/output check.

---

## Project scientific context

The project supports macaque V1 experiments involving:

* visual orientation discrimination
* orientation columns
* optogenetic stimulation
* congruent and incongruent optostimulation
* psychometric curves
* behavioral biasing and masking
* widefield calcium imaging
* stimulation power
* power density
* stimulated-column counts
* projector/camera coordinate transforms

Common active data and metadata structures may include:

* `behavioralData`
* `behavioralData.optoTS`
* `baselineTS`
* `bitmapData`
* `mdlStruct`
* final experiment/block index arrays
* chamber-specific summary tables
* power-cluster assignments

Do not assume a field’s semantics solely from its name. Verify how the active trusted pipeline consumes it.
