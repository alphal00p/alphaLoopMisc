# Gravity Numerator Evaluator

This repository is a path-finder study for evaluating quantum-gravity Feynman-graph numerators from DOT graphs. The current driver is [`evaluator.py`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/evaluator.py), which can process the split graph copies stored in [`dot_files/ScalarGravity_2L_processed.dot`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/dot_files/ScalarGravity_2L_processed.dot) and [`dot_files/ScalarGravity_3L_processed.dot`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/dot_files/ScalarGravity_3L_processed.dot).

Each DOT file contains several graph copies whose numerators sum back to the physical numerator. The split is technical and is used to keep the numerator at most linear in each EMR momentum component.

## Requirements

- Python environment with `symbolica`, `symbolica.community.spenso`, `symbolica.community.idenso`, `gammaloop`, `numpy`, `pydot`, `prettytable`, and `progressbar2`
- A valid Symbolica license in the environment
- The project virtualenv at [`.venv`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/.venv) if you want to reproduce the runs used in this repository

Typical shell setup:

```bash
cd /Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator
source .venv/bin/activate
export SYMBOLICA_LICENSE='dcec4a5e#6a95649c#7dca8216-8afe-57c8-975e-03eb5e68e4ee'
export SYMBOLICA_HIDE_BANNER=1
```

If `prettytable` or `progressbar2` is missing, run:

```bash
python3 -m pip install prettytable --upgrade
python3 -m pip install progressbar2 --upgrade
```

The helper script [update_dependencies.sh](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/update_dependencies.sh) also installs it after rebuilding the local `gammaloop` and `symbolica` wheels.

## Repository Layout

- [`evaluator.py`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/evaluator.py): main script
- [`dot_files/ScalarGravity_2L_processed.dot`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/dot_files/ScalarGravity_2L_processed.dot): 2-loop example input with 9 graph copies
- [`dot_files/ScalarGravity_3L_processed.dot`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/dot_files/ScalarGravity_3L_processed.dot): 3-loop example input with 27 graph copies
- [`ORIGINAL_OLD_BACKUP_VERSION/evaluator.py`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/ORIGINAL_OLD_BACKUP_VERSION/evaluator.py): old reference implementation

## What The Script Does

The current supported workflow is DOT-based:

1. Read one DOT file containing one or more graph copies.
2. For each graph, build a shared concrete kinematic point:
   - fixed external 4-momenta,
   - fixed loop-basis spatial momenta,
   - fixed OSE values.
3. Convert each edge `lmb_rep` into a concrete EMR 4-momentum.
4. Optionally evaluate the graph numerically with `spenso` by registering concrete `Q(edge)` tensors.
5. Optionally build a symbolic dot-product expression with `idenso`.
6. Sum the per-copy contributions.
7. If the dot-product form is enabled, build a compiled Symbolica evaluator.
8. Feed the evaluator with EMR dot products computed from the same concrete EMR momenta used by `spenso`.
9. Repeat that one aligned input point across a batch for timing.
10. Print timing and comparison summaries in colored `PrettyTable` output.

Two execution routes are available:

- `spenso` route: numerical tensor-network execution on concrete EMR `Q(edge)` tensors
- dot-product route: symbolic reduction to `spenso::dot(spenso::mink(4), Q(i), Q(j))`, followed by compiled Symbolica evaluation
- parametric `spenso` route: symbolic tensor-network execution on parametric EMR-component `Q(edge, spenso::cind(mu))` tensors, followed by compiled Symbolica evaluation

The comparison between the three routes is meaningful because they all now use the same physical EMR input point.

## General Structure Of `evaluator.py`

The code is organized around the `GraphProcessor` class.

- `GraphProcessor.__init__`: initializes the default kinematics, replacement rules, evaluator options, and parameter ordering
- `_build_edge_momenta`: evaluates each edge `lmb_rep` into a concrete EMR 4-vector
- `_build_emr_hep_library`: registers those concrete edge momenta as `Q(edge)` tensors for `spenso`
- `_build_graph_rules`: parses all node and edge numerator rules from the DOT graph
- `_get_graph_execution_plan`: provides the hardcoded contraction order and RAM cap for each supported graph family
- `execute_graph`: walks the graph in the chosen order and updates the `spenso` tensor network and/or symbolic dot-product expression
- `process_dot`: processes one DOT graph copy
- `process_dots`: processes all selected copies in a DOT file, sums them, runs the compiled evaluator, and prints the comparison report
- `_run_compiled_evaluator`: dedicated helper for evaluator build, compile, batched evaluation, and timing
- `_print_timing_report` and `_print_comparison_report`: terminal reporting helpers

## Current Kinematics

The script currently hardcodes one concrete test point in [`evaluator.py`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/evaluator.py):

- external momenta
  - `P0 = [1.0, 0.0, 1.0, 0.0]`
  - `P1 = [2.0, 1.0, 0.0, 1.0]`
  - `P2 = [1.5, 1.5, -1.0, 1.0]`
  - `P3 = [1.5, -0.5, 2.0, -0.0]`
- loop basis spatial momenta
  - loop `i` gets `[2/(i+1), 3/(i+1), 4/(i+1)]`
- OSE values
  - start from `(edge_id + 1) / n_edges`
  - overwrite the first entries with the energies of the external momenta

Those values define the loop-momentum basis point from which the EMR edge momenta and aligned dot products are derived.

## Usage

Basic syntax:

```bash
python3 evaluator.py <input.dot> [options]
```

Even though the positional argument is still called `expr_file`, the current maintained workflow is DOT input. If the input ends in `.txt`, the script dispatches to `process_expression`, which currently raises `NotImplementedError`.

## Command-Line Options

The table below describes every option currently parsed by [`evaluator.py`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/evaluator.py).

| Option | Meaning | Current effect in DOT workflow |
| --- | --- | --- |
| `expr_file` | Input path. Historically allowed text expressions. | Pass a `.dot` file here. `.txt` expression mode is currently not implemented. |
| `-v`, `--verbose` | Enable stepped execution mode with verbose logging. | Parsed, but currently not used by the DOT execution path. |
| `--n_loops <int>` | Number of loop basis vectors to seed. | Active. Must match the topology you want to evaluate. |
| `--n_edges <int>` | Number of edges used to define OSE values and dot-product parameter ordering. | Active. Should match the edge indexing used by the DOT input. |
| `--max_horner_scheme_variables <int>` | Limit the Horner-scheme variable count in the Symbolica evaluator builder. | Active when the compiled evaluator is built. |
| `--max_iterations_in_evaluator_building <int>` | Limit evaluator optimization rounds. | Active when the compiled evaluator is built. |
| `--n_cpe_iterations <int>` | Limit common-subexpression / CPE iterations in evaluator construction. | Active when the compiled evaluator is built. |
| `-ve`, `--verbose_evaluator` | Enable verbose output from evaluator construction. | Active for every compiled evaluator build in the DOT workflow, including both the dot-product evaluator and the parametric `spenso` evaluator. |
| `-ns`, `--no_spenso` | Disable the numerical `spenso` route. | Active. Useful for timing only the symbolic dot-product route. |
| `--n_graphs <int>` | Process only the first `n` graph copies in the DOT file. | Active. Good for quick tests. |
| `--n_core_for_evaluator_building <int>` | Number of cores given to the Symbolica evaluator builder. | Active when the compiled evaluator is built. |
| `--n_core_for_dot_processing <int>` | Number of worker processes used to process graph copies from the DOT file. | Active. `1` keeps the current serial behavior exactly. Values `>1` enable multicore DOT processing plus an aggregated live progress bar with total step progress, completed graphs, active components, and summed worker RAM. |
| `--sparse_momenta` | Use sparse tensor storage for momenta instead of dense 4-vectors. | Active for the current `spenso` EMR library and aligned dot-product input generation. |
| `-no_i`, `--no_complex_i` | Replace the explicit imaginary unit symbol `𝑖` by `1`. | Active during rule parsing. |
| `-c`, `--concrete_loop_momenta` | Execute with concrete loop-momentum values instead of symbolic loop-basis parameters. | Mostly legacy in the current DOT workflow. The maintained DOT route already derives a concrete aligned point. |
| `-fl`, `--float_constants` | Replace certain UFO constants by decimal floating-point numbers instead of rationals. | Active during rule parsing. |
| `-dots`, `--build_dot_products_form` | Build the symbolic dot-product expression. | Active. Required if you want the compiled evaluator path. |
| `--build_spenso_parametric_form` | Build the symbolic `spenso` scalar in terms of EMR momentum components instead of dot products. | Active. Builds a second compiled evaluator parameterized by `Q(edge, spenso::cind(mu))`. |
| `-flhep`, `--float_hep_lib` | Use the floating-point HEP tensor library instead of the atomic/symbolic one. | Active for `spenso` tensor registration and aligned dot-product input evaluation. |
| `-cose`, `--concrete_oses` | Execute with concrete OSE values instead of symbolic OSE parameters. | Mostly legacy in the current DOT workflow. The maintained DOT route already derives a concrete aligned point. |

## Output Files

During a run, the script may create:

- `dot_products_results__<dot-stem>__graphs_<count>__no_i_<0_or_1>__floatconst_<0_or_1>.txt`
  - cached summed dot-product expression for the current DOT input and option set
- `.symbolica_compiled/<label>.cpp`
  - generated C++ source for the compiled evaluator
- `.symbolica_compiled/<label>.so`
  - compiled shared object for the evaluator

## Example Runs

### 2-loop, dot-product route only

This skips `spenso`, builds the symbolic dot-product expression, compiles the evaluator, repeats one aligned input point across a batch, and prints timing information.

```bash
python3 evaluator.py dot_files/ScalarGravity_2L_processed.dot \
  -ve \
  --n_loops 2 \
  --n_edges 11 \
  --n_graphs 1 \
  --no_spenso \
  --build_dot_products_form \
  --n_core_for_evaluator_building 1
```

### 2-loop, aligned `spenso` vs compiled evaluator validation

This runs both routes on the same aligned EMR kinematics and prints a colored validation table with the absolute and relative mismatch.

```bash
python3 evaluator.py dot_files/ScalarGravity_2L_processed.dot \
  -ve \
  --n_loops 2 \
  --n_edges 11 \
  --n_graphs 1 \
  --build_dot_products_form \
  --n_core_for_evaluator_building 1
```

### 2-loop, direct `spenso` vs dot-product evaluator vs parametric `spenso` evaluator

This is the command for comparing all three paths in one run:

- `a)` dot-product symbolic construction plus compiled evaluator
- `b)` parametric `spenso` symbolic construction plus compiled evaluator
- `c)` direct numerical `spenso` evaluation

The parametric `spenso` scalar is much larger than the dot-product form, so this mode is useful for comparing build cost, evaluation speed, and numerical agreement.

```bash
python3 evaluator.py dot_files/ScalarGravity_2L_processed.dot \
  -ve \
  --n_loops 2 \
  --n_edges 11 \
  --n_graphs 1 \
  --build_dot_products_form \
  --build_spenso_parametric_form \
  --n_core_for_evaluator_building 1
```

### 2-loop, many graph copies

This is closer to the original study mode where the selected graph copies are summed before evaluator construction.

```bash
python3 evaluator.py dot_files/ScalarGravity_2L_processed.dot \
  -ve \
  --n_loops 2 \
  --n_edges 11 \
  --build_dot_products_form \
  --n_graphs 1000 \
  --n_core_for_evaluator_building 10
```

### 2-loop, multicore DOT-copy processing

This keeps the evaluator build single-core but parallelizes the per-copy graph processing phase. With `--n_core_for_dot_processing > 1`, the script switches from the serial per-step printout to a live aggregated progress bar showing total completed contraction steps, finished graphs, active worker components, and summed worker RAM usage.

```bash
python3 evaluator.py dot_files/ScalarGravity_2L_processed.dot \
  -ve \
  --n_loops 2 \
  --n_edges 11 \
  --n_graphs 4 \
  --build_dot_products_form \
  --n_core_for_dot_processing 2 \
  --n_core_for_evaluator_building 1
```

### 3-loop smoke test

Use a small subset first. The 3-loop graph family has a separate hardcoded contraction order and is heavier.

```bash
python3 evaluator.py dot_files/ScalarGravity_3L_processed.dot \
  -ve \
  --n_loops 3 \
  --n_edges 14 \
  --n_graphs 1 \
  --build_dot_products_form \
  --n_core_for_evaluator_building 1
```

## Reading The Output

The script typically prints:

- serial DOT processing output when `--n_core_for_dot_processing 1`
  - per-step progress lines
  - current step in the hardcoded contraction order
  - RAM usage
  - current tensor rank
- multicore DOT processing output when `--n_core_for_dot_processing > 1`
  - one aggregated progress bar
  - total completed steps versus the global target across all active graph components
  - completed graphs versus selected graphs
  - summed worker RAM usage
  - active worker/component labels with per-component step counters and current rank
- per-graph summary
  - `Numerical result` when `spenso` is enabled
  - symbolic expression size
- total summaries
  - `Total spenso result` when multiple graph copies are summed
  - green start/stop markers for each flavour
  - one final solid-border `PrettyTable` with columns `dot-products`, `spenso-numeric`, and `spenso-parametric`
  - the final scalar values
  - validation against the direct numerical `spenso` route when available
  - expression sizes for the evaluator-backed symbolic routes
  - tensor-network build time, symbolic construction time, evaluator build time, evaluator compile time, and runtime evaluation speed where applicable
  - explicit final printouts, immediately before the table, telling you where the dot-product and parametric-`spenso` evaluator expressions were saved

If `--no_spenso` is used, the `spenso-numeric` column is left empty in the final summary table.

## Important Caveats

- Only DOT input is currently supported in a maintained way.
- The contraction order is hardcoded in `_get_graph_execution_plan` and currently supports graph names starting with `h_diagram` and `hh_diagram`.
- The aligned numerical comparison only happens when both:
  - `spenso` is enabled
  - `--build_dot_products_form` is enabled
- Some CLI flags remain from older code paths and are currently legacy for the maintained DOT workflow, notably `--verbose`, `--concrete_loop_momenta`, and `--concrete_oses`.

## Reference

For the legacy implementation used as a performance sanity check, see [`ORIGINAL_OLD_BACKUP_VERSION/evaluator.py`](/Users/vjhirsch/Documents/Work/Projects/alphaLoop/git_alphaLoopMisc/gravity_evaluator/ORIGINAL_OLD_BACKUP_VERSION/evaluator.py).
