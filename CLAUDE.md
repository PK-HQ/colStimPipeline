@AGENTS.md

# Agent Execution Rules

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


## Claude Code-specific instructions

- Begin substantial changes in Plan mode.
- Before editing, inspect git status, recent commits, relevant callers, and relevant saved-data structures.
- Keep changes compatible with MATLAB R2018b.
- Do not silently substitute legacy or folded analysis methods for signed-source models.
- Keep changes focused on the requested task.
- Report tests/checks performed and anything that could not be executed.
- Do not commit or push unless explicitly asked.