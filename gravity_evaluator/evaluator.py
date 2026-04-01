import argparse
import concurrent.futures
import math
import multiprocessing as mp
import os
from queue import Empty
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pydot
from numpy.typing import NDArray
try:
    from prettytable import PrettyTable, TableStyle
except ImportError as exc:
    raise ImportError(
        "prettytable is required for formatted comparison output. "
        "Run update_dependencies.sh or `python3 -m pip install prettytable --upgrade`."
    ) from exc
from symbolica import E, Expression, N, P, Polynomial, Replacement, S
from symbolica.community.idenso import (
    cook_indices,
    simplify_color,
    simplify_gamma,
    simplify_metrics,
    to_dots,
)
from symbolica.community.spenso import (
    ExecutionMode,
    LibraryTensor,
    Representation,
    Tensor,
    TensorLibrary,
    TensorNetwork,
    TensorStructure,
)


class Utils:
    @staticmethod
    def _current_ram_mb():
        try:
            rss_kb = int(
                subprocess.check_output(
                    ["ps", "-o", "rss=", "-p", str(os.getpid())],
                    text=True,
                ).strip()
            )
            return rss_kb / 1024.0
        except Exception:
            return None

    @staticmethod
    def stepped_execution(tn, hep_lib, max_steps=None, t_delta=0.1, verbose=True):
        i_step = 0
        t_start_fn = time.time()
        while True:
            i_step += 1
            if verbose:
                elapsed_s = time.time() - t_start_fn
                print("[{:.3f}s] Performing scalar step    # {}".format(
                    elapsed_s, i_step))
                ram_mb = Utils._current_ram_mb()
                if ram_mb is not None:
                    elapsed_s = time.time() - t_start_fn
                    print("[{:.3f}s] Current RAM usage: {:.2f} MiB".format(
                        elapsed_s, ram_mb))
            if t_delta is not None:
                time.sleep(t_delta)
            # print(tn)
            # print(i_step)
            tn.execute(n_steps=1, mode=ExecutionMode.Scalar, library=hep_lib)

            if verbose:
                elapsed_s = time.time() - t_start_fn
                print("[{:.3f}s] Performing reduction step # {}".format(
                    elapsed_s, i_step))

            # i_step += 1

            # print(tn)
            # print(i_step)
            tn.execute(n_steps=1, mode=ExecutionMode.Single, library=hep_lib)
            if verbose:
                elapsed_s = time.time() - t_start_fn
                print("[{:.3f}s] DONE. size at this point: {}".format(
                    elapsed_s, len(str(tn))))
                ram_mb = Utils._current_ram_mb()
                if ram_mb is not None:
                    elapsed_s = time.time() - t_start_fn
                    print("[{:.3f}s] Current RAM usage: {:.2f} MiB".format(
                        elapsed_s, ram_mb))

            if max_steps and i_step > max_steps:
                break
            try:
                _ = tn.result_tensor(hep_lib)
            except:
                continue
            break

    @staticmethod
    def get_dot_node_id_from_name(node_name: str) -> int | None:
        name = node_name.strip('"').split(":")[0]
        if name.startswith("ext"):
            return None
        return int(name)  # type: ignore


class Ansi:
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    CYAN = "\033[96m"
    RESET = "\033[0m"


@dataclass
class EvaluatorRunResult:
    label: str
    batch_size: int
    first_result: complex
    build_time_ms: float
    compile_time_ms: float
    eval_time_us: float


@dataclass
class DotProcessingResult:
    spenso_result: Expression | complex | None
    spenso_parametric_result: Expression | None
    dot_result: Expression | None
    edge_momenta: dict[int, NDArray[np.float64]]
    spenso_timings: "SpensoTimingBreakdown | None"
    spenso_parametric_timings: "SpensoParametricTimingBreakdown | None"
    dot_product_timings: "DotProductTimingBreakdown | None"


@dataclass
class ParallelDotWorkerResult:
    graph_index: int
    graph_name: str
    edge_momenta: dict[int, tuple[float, float, float, float]]
    spenso_result: complex | None
    dot_result_path: str | None
    spenso_parametric_result_path: str | None
    spenso_timings: "SpensoTimingBreakdown | None"
    spenso_parametric_timings: "SpensoParametricTimingBreakdown | None"
    dot_product_timings: "DotProductTimingBreakdown | None"
    symbolic_expression_size_mb: float


@dataclass
class SpensoTimingBreakdown:
    graph_count: int
    edge_momenta_ms: float
    emr_library_ms: float
    tensor_network_build_ms: float
    execution_ms: float

    @property
    def total_ms(self) -> float:
        return (
            self.edge_momenta_ms
            + self.emr_library_ms
            + self.tensor_network_build_ms
            + self.execution_ms
        )

    def __add__(self, other: "SpensoTimingBreakdown") -> "SpensoTimingBreakdown":
        return SpensoTimingBreakdown(
            graph_count=self.graph_count + other.graph_count,
            edge_momenta_ms=self.edge_momenta_ms + other.edge_momenta_ms,
            emr_library_ms=self.emr_library_ms + other.emr_library_ms,
            tensor_network_build_ms=self.tensor_network_build_ms +
            other.tensor_network_build_ms,
            execution_ms=self.execution_ms + other.execution_ms,
        )


@dataclass
class DotProductTimingBreakdown:
    graph_count: int
    construction_ms: float
    per_graph_construction_samples: tuple[tuple[int, float], ...]

    @property
    def min_construction_ms(self) -> float:
        entry = self.min_construction_entry
        return 0.0 if entry is None else entry[1]

    @property
    def max_construction_ms(self) -> float:
        entry = self.max_construction_entry
        return 0.0 if entry is None else entry[1]

    @property
    def avg_construction_ms(self) -> float:
        if not self.per_graph_construction_samples:
            return 0.0
        values = np.asarray(
            [value for _, value in self.per_graph_construction_samples],
            dtype=np.float64,
        )
        return float(np.mean(values))

    @property
    def median_construction_ms(self) -> float:
        if not self.per_graph_construction_samples:
            return 0.0
        values = np.asarray(
            [value for _, value in self.per_graph_construction_samples],
            dtype=np.float64,
        )
        return float(np.median(values))

    @property
    def min_construction_entry(self) -> tuple[int, float] | None:
        if not self.per_graph_construction_samples:
            return None
        return min(self.per_graph_construction_samples, key=lambda sample: sample[1])

    @property
    def max_construction_entry(self) -> tuple[int, float] | None:
        if not self.per_graph_construction_samples:
            return None
        return max(self.per_graph_construction_samples, key=lambda sample: sample[1])

    def __add__(self, other: "DotProductTimingBreakdown") -> "DotProductTimingBreakdown":
        return DotProductTimingBreakdown(
            graph_count=self.graph_count + other.graph_count,
            construction_ms=self.construction_ms + other.construction_ms,
            per_graph_construction_samples=(
                self.per_graph_construction_samples + other.per_graph_construction_samples
            ),
        )


@dataclass
class SpensoParametricTimingBreakdown:
    graph_count: int
    emr_library_ms: float
    tensor_network_build_ms: float
    construction_ms: float
    per_graph_construction_samples: tuple[tuple[int, float], ...]

    @property
    def total_ms(self) -> float:
        return self.emr_library_ms + self.tensor_network_build_ms + self.construction_ms

    @property
    def min_construction_ms(self) -> float:
        entry = self.min_construction_entry
        return 0.0 if entry is None else entry[1]

    @property
    def max_construction_ms(self) -> float:
        entry = self.max_construction_entry
        return 0.0 if entry is None else entry[1]

    @property
    def avg_construction_ms(self) -> float:
        if not self.per_graph_construction_samples:
            return 0.0
        values = np.asarray(
            [value for _, value in self.per_graph_construction_samples],
            dtype=np.float64,
        )
        return float(np.mean(values))

    @property
    def median_construction_ms(self) -> float:
        if not self.per_graph_construction_samples:
            return 0.0
        values = np.asarray(
            [value for _, value in self.per_graph_construction_samples],
            dtype=np.float64,
        )
        return float(np.median(values))

    @property
    def min_construction_entry(self) -> tuple[int, float] | None:
        if not self.per_graph_construction_samples:
            return None
        return min(self.per_graph_construction_samples, key=lambda sample: sample[1])

    @property
    def max_construction_entry(self) -> tuple[int, float] | None:
        if not self.per_graph_construction_samples:
            return None
        return max(self.per_graph_construction_samples, key=lambda sample: sample[1])

    def __add__(self, other: "SpensoParametricTimingBreakdown") -> "SpensoParametricTimingBreakdown":
        return SpensoParametricTimingBreakdown(
            graph_count=self.graph_count + other.graph_count,
            emr_library_ms=self.emr_library_ms + other.emr_library_ms,
            tensor_network_build_ms=self.tensor_network_build_ms +
            other.tensor_network_build_ms,
            construction_ms=self.construction_ms + other.construction_ms,
            per_graph_construction_samples=(
                self.per_graph_construction_samples + other.per_graph_construction_samples
            ),
        )


@dataclass
class GraphExecutionResult:
    spenso_result: Expression | None
    dot_result: Expression | None
    spenso_execution_ms: float
    dot_product_construction_ms: float


class GraphProcessor:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.use_float_hep_lib = args.float_hep_lib
        self.use_sparse_momenta = args.sparse_momenta
        self.batch_size_for_timing = 1000
        self.mink_rep = Representation.mink(4)
        self.peak_ram_usage_mb: float | None = Utils._current_ram_mb()

        self.loop_momenta: list[NDArray[np.float64]] = []
        for i_loop in range(args.n_loops):
            self.loop_momenta.append(np.array(
                [2.0 / float(i_loop + 1), 3.0 / float(i_loop + 1), 4.0 / float(i_loop + 1)], dtype=np.float64))

        self.external_momenta = [
            np.array([1.0, 0.0, 1.0, 0.0], dtype=np.float64),
            np.array([2.0, 1.0, 0.0, 1.0], dtype=np.float64),
            -np.array([-1.5, -1.5, 1.0, -1.0], dtype=np.float64),
            -np.array([-1.5, 0.5, -2.0, 0.0], dtype=np.float64),
        ]

        if args.no_complex_i:
            self.string_replacements = [("𝑖", "")]
        else:
            self.string_replacements = []

        self.oses = [(ie + 1) / float(args.n_edges)
                     for ie in range(args.n_edges)]
        for i, p in enumerate(self.external_momenta):
            if i < len(self.oses):
                self.oses[i] = p[0]

        if args.float_constants:
            self.global_replacement_rules = [
                Replacement(E("UFO::dim"), E("4.0")),
                Replacement(E("UFO::TTT"), E("0.5")),
                Replacement(E("UFO::SSTmpart2"), E("0.3333333333333333")),
                Replacement(E("UFO::SSTmpart1"), E("0.25")),
                Replacement(E("UFO::SST"), E("0.2")),
                Replacement(E("UFO::SSTTmpart2"), E("0.16666666666666666")),
            ]
        else:
            self.global_replacement_rules = [
                Replacement(E("UFO::dim"), E("4")),
                Replacement(E("UFO::TTT"), E("1/2")),
                Replacement(E("UFO::SSTmpart2"), E("1/3")),
                Replacement(E("UFO::SSTmpart1"), E("1/4")),
                Replacement(E("UFO::SST"), E("1/5")),
                Replacement(E("UFO::SSTTmpart2"), E("1/6")),
            ]

        if self.use_float_hep_lib:
            self.hep_lib_with_momenta = TensorLibrary.hep_lib()
        else:
            self.hep_lib_with_momenta = TensorLibrary.hep_lib_atom()

        self.params = []
        self.param_values = []
        self.basis_component_values: dict[Expression, float] = {}
        vector_name_p = S("P")
        for i in range(len(self.external_momenta)):
            p_structure = TensorStructure(self.mink_rep, N(
                i), name=vector_name_p)  # type: ignore
            if args.concrete_loop_momenta:
                if self.use_sparse_momenta:
                    p = LibraryTensor.sparse(p_structure, type_info=float)
                    p[1], p[2], p[3] = (
                        self.external_momenta[i][1],
                        self.external_momenta[i][2],
                        self.external_momenta[i][3],
                    )
                else:
                    p = LibraryTensor.dense(
                        p_structure,
                        [
                            0,
                        ]
                        + [vi for vi in self.external_momenta[i][1:]],
                    )
            else:
                if self.use_sparse_momenta:
                    p = LibraryTensor.sparse(p_structure, type_info=Expression)
                    for j in range(1, 4):
                        p[j] = E("P({},spenso::cind({}))".format(i, j))
                else:
                    p = LibraryTensor.dense(
                        p_structure,
                        [
                            E("0"),
                        ]
                        + [E("P({},spenso::cind({}))".format(i, j))
                           for j in range(1, 4)],
                    )
                self.params.extend(
                    E("P({},spenso::cind({}))".format(i, j)) for j in range(1, 4))
                self.param_values.extend(
                    self.external_momenta[i][j] for j in range(1, 4))
            for j in range(1, 4):
                self.basis_component_values[E(
                    f"P({i},spenso::cind({j}))")] = float(self.external_momenta[i][j])
            self.hep_lib_with_momenta.register(p)

        vector_name_k = S("K")
        if args.concrete_loop_momenta:
            for i, k in enumerate(self.loop_momenta):
                ki_structure = TensorStructure(self.mink_rep, N(
                    i), name=vector_name_k)  # pyright: ignore
                if self.use_sparse_momenta:
                    ki = LibraryTensor.sparse(ki_structure, type_info=float)
                    ki[1], ki[2], ki[3] = k[0], k[1], k[2]
                else:
                    ki = LibraryTensor.dense(
                        ki_structure,
                        [
                            0,
                        ]
                        + [vi for vi in k],
                    )
                self.hep_lib_with_momenta.register(ki)
        else:
            for i, k in enumerate(self.loop_momenta):
                ki_structure = TensorStructure(self.mink_rep, N(
                    i), name=vector_name_k)  # pyright: ignore
                if self.use_sparse_momenta:
                    ki = LibraryTensor.sparse(
                        ki_structure, type_info=Expression)
                    for j in range(1, 4):
                        ki[j] = E("K({},spenso::cind({}))".format(i, j))
                else:
                    ki = LibraryTensor.dense(
                        ki_structure,
                        [
                            E("0"),
                        ]
                        + [E("K({},spenso::cind({}))".format(i, j))
                           for j in range(1, 4)],
                    )
                self.hep_lib_with_momenta.register(ki)
            for i, k in enumerate(self.loop_momenta):
                for j in range(1, 4):
                    self.params.append(
                        E("K({},spenso::cind({}))".format(i, j)))
                    self.param_values.append(k[j - 1])
        for i, k in enumerate(self.loop_momenta):
            for j in range(1, 4):
                self.basis_component_values[E(
                    f"K({i},spenso::cind({j}))")] = float(k[j - 1])

        self.OSE_vars = []
        for ie, ose in enumerate(self.oses):
            ose_structure = TensorStructure(
                self.mink_rep, N(ie), name=S("OSE"))  # type: ignore
            ose_var = E(f"OSE({ie})")
            self.OSE_vars.append(ose_var)
            if args.concrete_oses:
                ose_tensor = LibraryTensor.sparse(
                    ose_structure, type_info=float)
                ose_tensor[0] = ose
            else:
                ose_tensor = LibraryTensor.sparse(
                    ose_structure, type_info=Expression)
                ose_tensor[0] = ose_var
                self.params.append(ose_var)
                self.param_values.append(ose)
            self.basis_component_values[ose_var] = float(ose)
            self.hep_lib_with_momenta.register(ose_tensor)

        self.dot_variables_params = []
        for i_e in range(args.n_edges):
            for j_e in range(args.n_edges):
                if j_e >= i_e:
                    self.dot_variables_params.append(
                        E(f"spenso::dot(spenso::mink(4),Q({i_e}),Q({j_e}))"))
        self.edge_basis_symbol_replacements = [
            Replacement(E("P(i_,a___)"), E("PBASIS(i_)")),
            Replacement(E("K(i_,a___)"), E("KBASIS(i_)")),
            Replacement(E("OSE(i_,a___)"), E("OSEBASIS(i_)")),
        ]
        self.edge_basis_vectors = self._build_edge_basis_vectors()

        self.evaluator_base_options = {
            "constants": {},
            "functions": {},
            "verbose": args.verbose_evaluator,
            "cpe_iterations": args.n_cpe_iterations,
            "n_cores": args.n_core_for_evaluator_building,
        }
        if args.max_horner_scheme_variables is not None:
            self.evaluator_base_options["max_horner_scheme_variables"] = args.max_horner_scheme_variables
        if args.max_iterations_in_evaluator_building is not None:
            self.evaluator_base_options["iterations"] = args.max_iterations_in_evaluator_building

    def _get_evaluator_options(self, params: list[Expression]) -> dict[str, Any]:
        evaluator_options = dict(self.evaluator_base_options)
        evaluator_options["params"] = params
        return evaluator_options

    @staticmethod
    def _colorize(text: str, color: str) -> str:
        return f"{color}{text}{Ansi.RESET}"

    @staticmethod
    def _format_complex(value: complex) -> str:
        return f"{value.real:+.16e}{value.imag:+.16e}j"

    @staticmethod
    def _format_duration_ms(duration_ms: float) -> str:
        if duration_ms >= 1000.0:
            return f"{duration_ms / 1000.0:.2f} s"
        if duration_ms >= 1.0:
            return f"{duration_ms:.2f} ms"
        return f"{duration_ms * 1000.0:.2f} us"

    def _coerce_scalar_to_complex(self, value: Any) -> complex:
        if isinstance(value, complex):
            return value
        if isinstance(value, np.ndarray):
            if value.size != 1:
                raise TypeError(
                    f"Expected a scalar-compatible array, got shape {value.shape}")
            return self._coerce_scalar_to_complex(value.reshape(-1)[0])
        if isinstance(value, np.generic):
            return complex(value.item())
        if isinstance(value, Expression):
            try:
                return complex(str(value))
            except Exception:
                return complex(value.evaluate_complex({}, {}))
        return complex(value)

    def _extract_first_evaluator_result(self, evaluation_output: Any) -> complex:
        if isinstance(evaluation_output, np.ndarray):
            return self._coerce_scalar_to_complex(evaluation_output.reshape(-1)[0])
        if isinstance(evaluation_output, (list, tuple)):
            return self._coerce_scalar_to_complex(evaluation_output[0])
        if hasattr(evaluation_output, "__getitem__"):
            return self._coerce_scalar_to_complex(evaluation_output[0])
        return self._coerce_scalar_to_complex(evaluation_output)

    def _compiled_artifact_stem(self, label: str) -> Path:
        sanitized = "".join(ch if ch.isalnum()
                            else "_" for ch in label).strip("_")
        artifact_dir = Path(".symbolica_compiled")
        artifact_dir.mkdir(exist_ok=True)
        return artifact_dir / (sanitized or "compiled")

    def _announce_flavour(self, flavour: str, action: str) -> None:
        print(self._colorize(f"[{action}] {flavour}", Ansi.GREEN))

    def _record_ram_usage(self, ram_usage_mb: float | None = None) -> None:
        if ram_usage_mb is None:
            ram_usage_mb = Utils._current_ram_mb()
        if ram_usage_mb is None:
            return
        if self.peak_ram_usage_mb is None or ram_usage_mb > self.peak_ram_usage_mb:
            self.peak_ram_usage_mb = float(ram_usage_mb)

    def _expression_output_path(
        self,
        dot_file_path: str,
        args: argparse.Namespace,
        label: str,
    ) -> Path:
        artifact_dir = Path(".symbolica_expressions")
        artifact_dir.mkdir(exist_ok=True)
        graph_count = args.n_graphs if args.n_graphs is not None else "all"
        return artifact_dir / (
            f"{label}__{Path(dot_file_path).stem}__"
            f"graphs_{graph_count}__"
            f"no_i_{int(args.no_complex_i)}__"
            f"floatconst_{int(args.float_constants)}.txt"
        )

    def _save_expression_artifact(
        self,
        expr: Expression,
        dot_file_path: str,
        args: argparse.Namespace,
        label: str,
    ) -> tuple[Path, float]:
        artifact_path = self._expression_output_path(
            dot_file_path, args, label)
        artifact_path.write_text(expr.to_canonical_string(), encoding="utf-8")
        formatted_artifact_path = artifact_path.with_name(
            f"{artifact_path.stem}_formatted{artifact_path.suffix}"
        )
        formatted_artifact_path.write_text(expr.format(), encoding="utf-8")

        return artifact_path, expr.get_byte_size() / 1_000_000.0

    @staticmethod
    def _format_optional_complex(value: complex | None) -> str:
        return "" if value is None else GraphProcessor._format_complex(value)

    def _format_optional_duration_ms(self, value: float | None) -> str:
        return "" if value is None else self._format_duration_ms(value)

    def _format_optional_duration_with_graph(
        self,
        entry: tuple[int, float] | None,
    ) -> str:
        if entry is None:
            return ""
        graph_index, duration_ms = entry
        return f"{self._format_duration_ms(duration_ms)} (g{graph_index})"

    @staticmethod
    def _format_optional_path(path: Path | None) -> str:
        return "" if path is None else str(path)

    @staticmethod
    def _format_optional_size_mb(size_mb: float | None) -> str:
        return "" if size_mb is None else f"{size_mb:.2f} MB"

    @staticmethod
    def _format_optional_eval_runtime(result: EvaluatorRunResult | None) -> str:
        if result is None:
            return ""
        return f"{result.eval_time_us:.2f} us/sample x {result.batch_size}"

    @staticmethod
    def _minkowski_dot(
        left: NDArray[np.float64],
        right: NDArray[np.float64],
    ) -> complex:
        return complex(
            float(left[0] * right[0] - np.dot(left[1:], right[1:])),
            0.0,
        )

    def _colorize_summary_value(
        self,
        value: str,
        *,
        color: str | None = None,
    ) -> str:
        if value == "":
            return ""
        if color is not None:
            return self._colorize(value, color)
        if value == "baseline":
            return self._colorize(value, Ansi.YELLOW)
        if value == "PASS":
            return self._colorize(value, Ansi.GREEN)
        if value == "FAIL":
            return self._colorize(value, Ansi.RED)
        return self._colorize(value, Ansi.CYAN)

    @staticmethod
    def _comparison_cells(
        candidate: complex | None,
        reference: complex | None,
    ) -> tuple[str, str, str]:
        if candidate is None or reference is None:
            return ("", "", "")
        abs_diff = abs(candidate - reference)
        rel_diff = abs_diff / max(abs(reference), 1e-30)
        is_close = np.isclose(candidate, reference, rtol=1e-10, atol=1e-12)
        status = "PASS" if bool(is_close) else "FAIL"
        return (status, f"{abs_diff:.16e}", f"{rel_diff:.16e}")

    def _print_final_summary_table(
        self,
        dot_result: complex | None,
        dot_evaluator_run: EvaluatorRunResult | None,
        dot_product_timings: DotProductTimingBreakdown | None,
        dot_expression_path: Path | None,
        dot_expression_size_mb: float | None,
        spenso_numeric_result: complex | None,
        spenso_timings: SpensoTimingBreakdown | None,
        spenso_parametric_result: complex | None,
        spenso_parametric_evaluator_run: EvaluatorRunResult | None,
        spenso_parametric_timings: SpensoParametricTimingBreakdown | None,
        spenso_parametric_expression_path: Path | None,
        spenso_parametric_expression_size_mb: float | None,
        peak_ram_usage_mb: float | None,
    ) -> None:
        dot_status, dot_abs_diff, dot_rel_diff = self._comparison_cells(
            dot_result,
            spenso_numeric_result,
        )
        param_status, param_abs_diff, param_rel_diff = self._comparison_cells(
            spenso_parametric_result,
            spenso_numeric_result,
        )

        table = PrettyTable()
        title_suffix_parts = [
            "cores: "
            f"dot={self.args.n_core_for_dot_processing}, "
            f"evaluator={self.args.n_core_for_evaluator_building}"
        ]
        if peak_ram_usage_mb is not None:
            title_suffix_parts.insert(
                0,
                f"peak RAM usage: {peak_ram_usage_mb:.2f} MiB",
            )
        table.title = (
            "Evaluation Strategy Summary"
            if not title_suffix_parts
            else f"Evaluation Strategy Summary ({' | '.join(title_suffix_parts)})"
        )
        table.set_style(TableStyle.SINGLE_BORDER)
        table.field_names = [
            "Metric",
            self._colorize("dot-products", Ansi.CYAN),
            self._colorize("spenso-numeric", Ansi.CYAN),
            self._colorize("spenso-parametric", Ansi.CYAN),
        ]
        for field_name in table.field_names:
            table.align[field_name] = "l"

        rows = [
            ("Final scalar",
             self._format_optional_complex(dot_result),
             self._format_optional_complex(spenso_numeric_result),
             self._format_optional_complex(spenso_parametric_result)),
            ("Validation vs spenso-numeric",
             dot_status,
             "baseline" if spenso_numeric_result is not None else "",
             param_status),
            ("Abs diff vs spenso-numeric",
             dot_abs_diff,
             "",
             param_abs_diff),
            ("Rel diff vs spenso-numeric",
             dot_rel_diff,
             "",
             param_rel_diff),
            ("Expression size",
             self._format_optional_size_mb(dot_expression_size_mb),
             "",
             self._format_optional_size_mb(spenso_parametric_expression_size_mb)),
            ("Q library build",
             "",
             self._format_optional_duration_ms(
                 None if spenso_timings is None else spenso_timings.emr_library_ms
             ),
             self._format_optional_duration_ms(
                 None if spenso_parametric_timings is None else spenso_parametric_timings.emr_library_ms
             )),
            ("Tensor network build",
             "",
             self._format_optional_duration_ms(
                 None if spenso_timings is None else spenso_timings.tensor_network_build_ms
             ),
             self._format_optional_duration_ms(
                 None if spenso_parametric_timings is None else spenso_parametric_timings.tensor_network_build_ms
             )),
            ("Symbolic construction total",
             self._format_optional_duration_ms(
                 None if dot_product_timings is None else dot_product_timings.construction_ms
             ),
             "",
             self._format_optional_duration_ms(
                 None if spenso_parametric_timings is None else spenso_parametric_timings.construction_ms
             )),
            ("Symbolic construction avg/graph",
             self._format_optional_duration_ms(
                 None if dot_product_timings is None else dot_product_timings.avg_construction_ms
             ),
             "",
             self._format_optional_duration_ms(
                 None if spenso_parametric_timings is None else spenso_parametric_timings.avg_construction_ms
             )),
            ("Symbolic construction min/graph",
             self._format_optional_duration_with_graph(
                 None if dot_product_timings is None else dot_product_timings.min_construction_entry
             ),
             "",
             self._format_optional_duration_with_graph(
                 None if spenso_parametric_timings is None else spenso_parametric_timings.min_construction_entry
             )),
            ("Symbolic construction median/graph",
             self._format_optional_duration_ms(
                 None if dot_product_timings is None else dot_product_timings.median_construction_ms
             ),
             "",
             self._format_optional_duration_ms(
                 None if spenso_parametric_timings is None else spenso_parametric_timings.median_construction_ms
             )),
            ("Symbolic construction max/graph",
             self._format_optional_duration_with_graph(
                 None if dot_product_timings is None else dot_product_timings.max_construction_entry
             ),
             "",
             self._format_optional_duration_with_graph(
                 None if spenso_parametric_timings is None else spenso_parametric_timings.max_construction_entry
             )),
            ("Evaluator build",
             self._format_optional_duration_ms(
                 None if dot_evaluator_run is None else dot_evaluator_run.build_time_ms
             ),
             "",
             self._format_optional_duration_ms(
                 None if spenso_parametric_evaluator_run is None else spenso_parametric_evaluator_run.build_time_ms
             )),
            ("Evaluator compile",
             self._format_optional_duration_ms(
                 None if dot_evaluator_run is None else dot_evaluator_run.compile_time_ms
             ),
             "",
             self._format_optional_duration_ms(
                 None if spenso_parametric_evaluator_run is None else spenso_parametric_evaluator_run.compile_time_ms
             )),
            ("Evaluation time",
             self._format_optional_eval_runtime(dot_evaluator_run),
             self._format_optional_duration_ms(
                 None if spenso_timings is None else spenso_timings.execution_ms
             ),
             self._format_optional_eval_runtime(spenso_parametric_evaluator_run)),
        ]

        for metric, dot_cell, numeric_cell, param_cell in rows:
            table.add_row([
                metric,
                self._colorize_summary_value(dot_cell),
                self._colorize_summary_value(numeric_cell),
                self._colorize_summary_value(param_cell),
            ])
        print(table)

    def _print_saved_expression_locations(
        self,
        dot_expression_path: Path | None,
        spenso_parametric_expression_path: Path | None,
    ) -> None:
        if dot_expression_path is not None:
            print(self._colorize(
                f"Saved dot-products evaluator expression to {
                    dot_expression_path.resolve()}",
                Ansi.GREEN,
            ))
        if spenso_parametric_expression_path is not None:
            print(self._colorize(
                "Saved spenso-parametric evaluator expression to "
                f"{spenso_parametric_expression_path.resolve()}",
                Ansi.GREEN,
            ))

    def _run_compiled_evaluator(
        self,
        expr: Expression,
        params: list[Expression],
        input_row: NDArray[np.complex128],
        batch_size: int,
        label: str,
        compile_mode: str = "complex",
        verbose: bool | None = None,
    ) -> EvaluatorRunResult:
        evaluator_options = self._get_evaluator_options(params)
        effective_verbose = self.args.verbose_evaluator if verbose is None else verbose
        evaluator_options["verbose"] = effective_verbose
        self._announce_flavour(label, "START")
        self._record_ram_usage()
        if effective_verbose:
            print(
                self._colorize(
                    f"[VERBOSE] Symbolica evaluator optimization enabled for {
                        label}",
                    Ansi.YELLOW,
                )
            )
        t_start = time.time()
        eager_eval = expr.evaluator(**evaluator_options)
        build_time_ms = (time.time() - t_start) * 1000.0
        self._record_ram_usage()

        artifact_stem = self._compiled_artifact_stem(label)
        t_start = time.time()
        compiled_eval = eager_eval.compile(
            artifact_stem.name,
            str(artifact_stem.with_suffix(".cpp")),
            str(artifact_stem.with_suffix(".so")),
            compile_mode,
            inline_asm="default",
        )
        compile_time_ms = (time.time() - t_start) * 1000.0
        self._record_ram_usage()

        batch = np.repeat(input_row[np.newaxis, :], batch_size, axis=0)
        t_start = time.time()
        evaluation_output = compiled_eval.evaluate(batch)
        eval_time_us = (time.time() - t_start) * \
            1_000_000.0 / float(batch_size)
        self._record_ram_usage()
        first_result = self._extract_first_evaluator_result(evaluation_output)

        result = EvaluatorRunResult(
            label=label,
            batch_size=batch_size,
            first_result=first_result,
            build_time_ms=build_time_ms,
            compile_time_ms=compile_time_ms,
            eval_time_us=eval_time_us,
        )
        self._announce_flavour(label, "STOP")
        return result

    def _build_edge_basis_vectors(self) -> dict[Expression, NDArray[np.float64]]:
        basis_vectors: dict[Expression, NDArray[np.float64]] = {}
        for i, momentum in enumerate(self.external_momenta):
            basis_vectors[E(f"PBASIS({i})")] = np.array(
                momentum, dtype=np.float64)
        for i, momentum in enumerate(self.loop_momenta):
            basis_vectors[E(f"KBASIS({i})")] = np.array(
                [0.0] + [float(value) for value in momentum],
                dtype=np.float64,
            )
        for ie, ose in enumerate(self.oses):
            basis_vectors[E(f"OSEBASIS({ie})")] = np.array(
                [float(ose), 0.0, 0.0, 0.0],
                dtype=np.float64,
            )
        return basis_vectors

    def _build_edge_momenta(
        self,
        dot_graph: pydot.Dot,
    ) -> dict[int, NDArray[np.float64]]:
        edge_momenta: dict[int, NDArray[np.float64]] = {}
        for edge in dot_graph.get_edges():
            edge_id = int(edge.get("id").strip('"'))
            lmb_rep = edge.get("lmb_rep")
            assert lmb_rep is not None, f"Missing lmb_rep for edge {edge_id}"
            emr_expr = E(f"OSE({edge_id},a___)") + E(lmb_rep.strip('"'))
            basis_expr = emr_expr.replace_multiple(
                self.edge_basis_symbol_replacements,
                repeat=True,
            )
            momentum = np.zeros(4, dtype=np.float64)
            for basis_symbol, basis_vector in self.edge_basis_vectors.items():
                coeff = self._coerce_scalar_to_complex(
                    basis_expr.coefficient(basis_symbol))
                if coeff != 0:
                    momentum += float(coeff.real) * basis_vector
            edge_momenta[edge_id] = momentum
        return edge_momenta

    def _build_emr_hep_library(
        self,
        edge_momenta: dict[int, NDArray[np.float64]],
    ) -> TensorLibrary:
        if self.use_float_hep_lib:
            hep_lib = TensorLibrary.hep_lib()
        else:
            hep_lib = TensorLibrary.hep_lib_atom()

        vector_name_q = S("Q")
        for edge_id, momentum in sorted(edge_momenta.items()):
            q_structure = TensorStructure(self.mink_rep, N(
                edge_id), name=vector_name_q)  # type: ignore
            if self.use_sparse_momenta:
                q_tensor = LibraryTensor.sparse(q_structure, type_info=float)
                for component, value in enumerate(momentum):
                    if value != 0.0:
                        q_tensor[component] = float(value)
            else:
                q_tensor = LibraryTensor.dense(
                    q_structure,
                    [float(value) for value in momentum],
                )
            hep_lib.register(q_tensor)
        return hep_lib

    @staticmethod
    def _get_sorted_edge_ids(edge_momenta: dict[int, NDArray[np.float64]]) -> list[int]:
        return sorted(edge_momenta)

    def _build_emr_component_params(
        self,
        edge_ids: list[int],
    ) -> list[Expression]:
        return [
            E(f"Q({edge_id},spenso::cind({component}))")
            for edge_id in edge_ids
            for component in range(4)
        ]

    def _build_parametric_emr_hep_library(
        self,
        edge_ids: list[int],
    ) -> TensorLibrary:
        hep_lib = TensorLibrary.hep_lib_atom()
        vector_name_q = S("Q")
        for edge_id in edge_ids:
            q_structure = TensorStructure(self.mink_rep, N(
                edge_id), name=vector_name_q)  # type: ignore
            component_params = [
                E(f"Q({edge_id},spenso::cind({component}))")
                for component in range(4)
            ]
            if self.use_sparse_momenta:
                q_tensor = LibraryTensor.sparse(
                    q_structure, type_info=Expression)
                for component, component_param in enumerate(component_params):
                    q_tensor[component] = component_param
            else:
                q_tensor = LibraryTensor.dense(q_structure, component_params)
            hep_lib.register(q_tensor)
        return hep_lib

    def _build_emr_component_input_row(
        self,
        edge_momenta: dict[int, NDArray[np.float64]],
        edge_ids: list[int],
    ) -> NDArray[np.complex128]:
        component_values: list[complex] = []
        for edge_id in edge_ids:
            momentum = edge_momenta[edge_id]
            component_values.extend(
                complex(float(value), 0.0) for value in momentum
            )
        return np.asarray(component_values, dtype=np.complex128)

    def _build_dot_product_input_row(
        self,
        edge_momenta: dict[int, NDArray[np.float64]],
    ) -> NDArray[np.complex128]:
        dot_values = []
        zero_momentum = np.zeros(4, dtype=np.float64)
        for left_edge_id in range(self.args.n_edges):
            left_momentum = edge_momenta.get(left_edge_id, zero_momentum)
            for right_edge_id in range(left_edge_id, self.args.n_edges):
                right_momentum = edge_momenta.get(right_edge_id, zero_momentum)
                dot_values.append(self._minkowski_dot(
                    left_momentum, right_momentum))
        assert len(dot_values) == len(self.dot_variables_params), (
            f"Built {len(dot_values)} aligned dot-product inputs but expected "
            f"{len(self.dot_variables_params)} parameters."
        )
        return np.asarray(dot_values, dtype=np.complex128)

    def _validate_dot_parameter_coverage(
        self,
        edge_momenta: dict[int, NDArray[np.float64]],
    ) -> None:
        if not edge_momenta:
            return
        max_edge_id = max(edge_momenta)
        if max_edge_id >= self.args.n_edges:
            raise ValueError(
                f"--n_edges {self.args.n_edges} is too small for this graph: "
                f"encountered edge id {
                    max_edge_id}, so at least {max_edge_id + 1} "
                "edge slots are required."
            )

    def _assert_matching_edge_momenta(
        self,
        reference: dict[int, NDArray[np.float64]],
        candidate: dict[int, NDArray[np.float64]],
        graph_name: str,
    ) -> None:
        if set(reference) != set(candidate):
            raise ValueError(
                f"EMR edge sets differ for graph {graph_name}: "
                f"{sorted(reference)} != {sorted(candidate)}"
            )
        for edge_id in sorted(reference):
            if not np.allclose(reference[edge_id], candidate[edge_id], atol=1e-12, rtol=0.0):
                raise ValueError(
                    f"Inconsistent EMR momentum for edge {
                        edge_id} in graph {graph_name}: "
                    f"{reference[edge_id]} != {candidate[edge_id]}"
                )

    def _dot_products_cache_path(
        self,
        dot_file_path: str,
        args: argparse.Namespace,
    ) -> Path:
        graph_count = args.n_graphs if args.n_graphs is not None else "all"
        return Path(
            "dot_products_results__"
            f"{Path(dot_file_path).stem}__"
            f"graphs_{graph_count}__"
            f"no_i_{int(args.no_complex_i)}__"
            f"floatconst_{int(args.float_constants)}.txt"
        )

    def _parse_rule_expression(self, raw_expression: str) -> Expression:
        str_expr = raw_expression.strip('"')
        for old, new in self.string_replacements:
            str_expr = str_expr.replace(old, new)
        rule = E(str_expr)
        return rule.replace_multiple(self.global_replacement_rules, repeat=True)

    def _build_graph_rules(
        self,
        dot_graph: pydot.Dot,
        hep_lib: TensorLibrary | None,
        build_tensor_networks: bool,
    ) -> tuple[
        dict[int, tuple[TensorNetwork | None, Expression]],
        dict[tuple[int | None, int | None],
             tuple[TensorNetwork | None, Expression]],
    ]:
        vertex_rules: dict[int, tuple[TensorNetwork | None, Expression]] = {}
        edge_rules: dict[tuple[int | None, int | None],
                         tuple[TensorNetwork | None, Expression]] = {}

        for node in dot_graph.get_nodes():
            node_id = Utils.get_dot_node_id_from_name(node.get_name())
            if node_id is None:
                continue
            rule = self._parse_rule_expression(
                node.get_attributes().get("num", "1"))
            tn: TensorNetwork | None = None
            if build_tensor_networks:
                assert hep_lib is not None
                tn = TensorNetwork(cook_indices(rule), hep_lib)
                tn.execute(hep_lib)
            vertex_rules[node_id] = (tn, rule)

        for edge in dot_graph.get_edges():
            src = Utils.get_dot_node_id_from_name(edge.get_source())
            dst = Utils.get_dot_node_id_from_name(edge.get_destination())
            rule = self._parse_rule_expression(
                edge.get_attributes().get("num", "1"))
            tn = None
            if build_tensor_networks:
                assert hep_lib is not None
                tn = TensorNetwork(cook_indices(rule), hep_lib)
                tn.execute(hep_lib)
            edge_rules[tuple(
                sorted([src, dst], key=lambda a: -1 if a is None else a))] = (tn, rule)

        return vertex_rules, edge_rules

    def _get_graph_execution_plan(
        self,
        graph_name: str,
    ) -> tuple[list[int | tuple[int | None, int | None]], float]:
        if graph_name.startswith("h_diagram"):
            return ([
                (None, 0),
                0,
                (0, 1),
                (0, 4),
                (None, 2),
                2,
                (2, 3),
                (2, 4),
                4,
                (4, 5),
                5,
                (5, 1),
                1,
                (None, 1),
                (5, 3),
                3,
                (None, 3),
            ], 15000.0)
        if graph_name.startswith("hh_diagram"):
            return ([
                (None, 0),
                0,
                (0, 1),
                (0, 5),
                (None, 3),
                3,
                (3, 4),
                (3, 5),
                5,
                (5, 6),
                1,
                (1, 2),
                (1, 6),
                6,
                (6, 7),
                4,
                (None, 4),
                (4, 7),
                7,
                (7, 2),
                2,
                (2, None),
            ], 40000.0)
        raise NotImplementedError(
            f"DOT processing is not implemented yet for diagram '{graph_name}'.")

    def process_expression(self, expr: Expression, args: argparse.Namespace):
        raise NotImplementedError(
            "Expression processing is no longer up to date, only DOT graph processing is supported for now.")

    def _collect_dot_results_serial(
        self,
        selected_dot_graphs: list[pydot.Dot],
        args: argparse.Namespace,
    ) -> tuple[list[tuple[str, DotProcessingResult]], dict[int, NDArray[np.float64]] | None]:
        all_dot_results: list[tuple[str, DotProcessingResult]] = []
        if args.build_dot_products_form:
            self._announce_flavour("dot-products construction", "START")
        if not args.no_spenso:
            self._announce_flavour("spenso-numeric direct evaluation", "START")
        if args.build_spenso_parametric_form:
            self._announce_flavour("spenso-parametric construction", "START")
        self._record_ram_usage()
        reference_edge_momenta: dict[int, NDArray[np.float64]] | None = None
        for dot_graph in selected_dot_graphs:
            print("Processing graph: {}".format(dot_graph.get_name()))
            res = self.process_dot(
                dot_graph, args, graph_display_index=len(all_dot_results) + 1)
            graph_name = dot_graph.get_name() or "<unnamed>"
            if reference_edge_momenta is None:
                reference_edge_momenta = res.edge_momenta
            else:
                self._assert_matching_edge_momenta(
                    reference_edge_momenta,
                    res.edge_momenta,
                    graph_name,
                )

            memory = 0
            if isinstance(res.spenso_result, Expression):
                memory += res.spenso_result.get_byte_size()
            if res.spenso_parametric_result is not None:
                memory += res.spenso_parametric_result.get_byte_size()
            if res.dot_result is not None:
                memory += res.dot_result.get_byte_size()
            print(f"Symbolic expression size: {memory / 1000000.0:.2f} MB")
            self._record_ram_usage()
            all_dot_results.append((graph_name, res))

        if args.build_dot_products_form:
            self._announce_flavour("dot-products construction", "STOP")
        if not args.no_spenso:
            self._announce_flavour("spenso-numeric direct evaluation", "STOP")
        if args.build_spenso_parametric_form:
            self._announce_flavour("spenso-parametric construction", "STOP")
        return all_dot_results, reference_edge_momenta

    def _format_parallel_active_components(
        self,
        component_states: dict[tuple[int, str], dict[str, Any]],
        max_width: int = 160,
    ) -> str:
        active_parts: list[str] = []
        for key in sorted(component_states):
            graph_index, _ = key
            status = component_states[key]
            if status["finished"]:
                continue
            active_parts.append(
                f"g{graph_index + 1}:{status['component']} "
                f"{status['current_step']}/{status['total_steps']} "
                f"r{status['rank']} {status['ram_usage_mb']:.0f}MiB"
            )
        if not active_parts:
            return "idle"
        active_text = " | ".join(active_parts)
        if len(active_text) <= max_width:
            return active_text
        return active_text[: max_width - 3] + "..."

    @staticmethod
    def _parallel_progress_term_width(max_width: int = 300) -> int:
        try:
            terminal_width = os.get_terminal_size().columns
        except OSError:
            terminal_width = max_width
        return max(120, min(max_width, terminal_width))

    def _collect_dot_results_parallel(
        self,
        dot_file_path: str,
        selected_dot_graphs: list[pydot.Dot],
        args: argparse.Namespace,
    ) -> tuple[list[tuple[str, DotProcessingResult]], dict[int, NDArray[np.float64]] | None]:
        try:
            import progressbar
        except ImportError as exc:
            raise ImportError(
                "progressbar2 is required for multicore DOT processing. "
                "Run `python3 -m pip install progressbar2 --upgrade`."
            ) from exc

        all_dot_results: list[tuple[str, DotProcessingResult]] = []
        if args.build_dot_products_form:
            self._announce_flavour("dot-products construction", "START")
        if not args.no_spenso:
            self._announce_flavour("spenso-numeric direct evaluation", "START")
        if args.build_spenso_parametric_form:
            self._announce_flavour("spenso-parametric construction", "START")

        total_target_steps = 0
        for dot_graph in selected_dot_graphs:
            graph_name = (dot_graph.get_name() or "<unnamed>").strip('"')
            order, _ = self._get_graph_execution_plan(graph_name)
            total_target_steps += len(order) * (1 +
                                                int(args.build_spenso_parametric_form))

        component_states: dict[tuple[int, str], dict[str, Any]] = {}
        completed_results: dict[int, ParallelDotWorkerResult] = {}
        progress_term_width = self._parallel_progress_term_width()
        active_max_width = max(40, progress_term_width - 140)
        widgets = [
            progressbar.Percentage(),
            " ",
            progressbar.Bar(marker="█", left="│", right="│"),
            " ",
            progressbar.Timer(),
            " ",
            progressbar.ETA(),
            " | ",
            progressbar.DynamicMessage("steps"),
            " | ",
            progressbar.DynamicMessage("graphs"),
            " | ",
            progressbar.DynamicMessage("ram"),
            " | ",
            progressbar.DynamicMessage("active"),
        ]

        with tempfile.TemporaryDirectory(prefix=".parallel_dot_results_", dir=".") as temp_dir:
            with mp.Manager() as manager:
                progress_queue = manager.Queue()
                ctx = mp.get_context("spawn")
                with concurrent.futures.ProcessPoolExecutor(
                    max_workers=args.n_core_for_dot_processing,
                    mp_context=ctx,
                ) as executor:
                    future_to_index = {
                        executor.submit(
                            _parallel_process_dot_worker,
                            args,
                            dot_graph.to_string(),
                            graph_index,
                            temp_dir,
                            progress_queue,
                        ): graph_index
                        for graph_index, dot_graph in enumerate(selected_dot_graphs)
                    }
                    pending = set(future_to_index)
                    with progressbar.ProgressBar(
                        max_value=max(total_target_steps, 1),
                        widgets=widgets,
                        redirect_stdout=True,
                        term_width=progress_term_width,
                    ) as bar:
                        last_bar_state: tuple[int, int,
                                              float, str] | None = None
                        while pending:
                            saw_event = False
                            while True:
                                try:
                                    event = progress_queue.get(
                                        timeout=0.1 if not saw_event else 0.0)
                                except Empty:
                                    break
                                saw_event = True
                                component_states[(event["graph_index"], event["component"])] = {
                                    "graph_name": event["graph_name"],
                                    "component": event["component"],
                                    "current_step": event["current_step"],
                                    "total_steps": event["total_steps"],
                                    "ram_usage_mb": event["ram_usage_mb"],
                                    "rank": event["rank"],
                                    "finished": event["event"] == "component_finished",
                                }

                            done_futures = {
                                future for future in pending if future.done()}
                            for future in done_futures:
                                payload = future.result()
                                completed_results[payload.graph_index] = payload
                            pending -= done_futures

                            current_steps = sum(
                                min(status["current_step"],
                                    status["total_steps"])
                                for status in component_states.values()
                            )
                            current_ram_mb = sum(
                                status["ram_usage_mb"]
                                for status in component_states.values()
                                if not status["finished"]
                            )
                            self._record_ram_usage(current_ram_mb)
                            active_text = self._format_parallel_active_components(
                                component_states,
                                max_width=active_max_width,
                            )
                            bar_state = (
                                min(current_steps, max(total_target_steps, 1)),
                                len(completed_results),
                                round(current_ram_mb, 1),
                                active_text,
                            )
                            if bar_state != last_bar_state:
                                bar.update(
                                    bar_state[0],
                                    steps=f"{
                                        current_steps}/{total_target_steps}",
                                    graphs=f"{len(completed_results)
                                              }/{len(selected_dot_graphs)}",
                                    ram=f"{current_ram_mb:.1f} MiB",
                                    active=active_text,
                                )
                                last_bar_state = bar_state

            reference_edge_momenta: dict[int,
                                         NDArray[np.float64]] | None = None
            for graph_index in sorted(completed_results):
                payload = completed_results[graph_index]
                edge_momenta = {
                    edge_id: np.asarray(values, dtype=np.float64)
                    for edge_id, values in payload.edge_momenta.items()
                }
                res = DotProcessingResult(
                    spenso_result=payload.spenso_result,
                    spenso_parametric_result=(
                        None
                        if payload.spenso_parametric_result_path is None
                        else E(Path(payload.spenso_parametric_result_path).read_text(encoding="utf-8").strip())
                    ),
                    dot_result=(
                        None
                        if payload.dot_result_path is None
                        else E(Path(payload.dot_result_path).read_text(encoding="utf-8").strip())
                    ),
                    edge_momenta=edge_momenta,
                    spenso_timings=payload.spenso_timings,
                    spenso_parametric_timings=payload.spenso_parametric_timings,
                    dot_product_timings=payload.dot_product_timings,
                )
                if reference_edge_momenta is None:
                    reference_edge_momenta = res.edge_momenta
                else:
                    self._assert_matching_edge_momenta(
                        reference_edge_momenta,
                        res.edge_momenta,
                        payload.graph_name,
                    )
                all_dot_results.append((payload.graph_name, res))

        if args.build_dot_products_form:
            self._announce_flavour("dot-products construction", "STOP")
        if not args.no_spenso:
            self._announce_flavour("spenso-numeric direct evaluation", "STOP")
        if args.build_spenso_parametric_form:
            self._announce_flavour("spenso-parametric construction", "STOP")
        return all_dot_results, reference_edge_momenta

    def process_dots(self, dot_file_path: str, args: argparse.Namespace):
        self._record_ram_usage()
        dot_graphs = pydot.graph_from_dot_file(dot_file_path)
        assert dot_graphs is not None, "Failed to parse DOT file: {}".format(
            dot_file_path)
        selected_dot_graphs = dot_graphs if args.n_graphs is None else dot_graphs[
            : args.n_graphs]
        if args.n_core_for_dot_processing == 1:
            all_dot_results, reference_edge_momenta = self._collect_dot_results_serial(
                selected_dot_graphs,
                args,
            )
        else:
            all_dot_results, reference_edge_momenta = self._collect_dot_results_parallel(
                dot_file_path,
                selected_dot_graphs,
                args,
            )

        if not all_dot_results:
            print("No graphs selected for processing.")
            return

        total_spenso_result: complex | None = None
        total_spenso_timings: SpensoTimingBreakdown | None = None
        total_spenso_parametric_result: Expression | None = None
        total_spenso_parametric_timings: SpensoParametricTimingBreakdown | None = None
        total_dot_product_timings: DotProductTimingBreakdown | None = None
        dot_expression_path: Path | None = None
        dot_expression_size_mb: float | None = None
        spenso_parametric_expression_path: Path | None = None
        spenso_parametric_expression_size_mb: float | None = None
        if all_dot_results[0][1].spenso_result is not None:
            total_spenso_value = sum(
                result.spenso_result for _, result in all_dot_results
            )  # type: ignore[arg-type]
            total_spenso_result = self._coerce_scalar_to_complex(
                total_spenso_value)
            total_spenso_timings = sum(
                (
                    result.spenso_timings
                    for _, result in all_dot_results
                    if result.spenso_timings is not None
                ),
                start=SpensoTimingBreakdown(
                    graph_count=0,
                    edge_momenta_ms=0.0,
                    emr_library_ms=0.0,
                    tensor_network_build_ms=0.0,
                    execution_ms=0.0,
                ),
            )
        if all_dot_results[0][1].spenso_parametric_result is not None:
            total_spenso_parametric_result = sum(
                result.spenso_parametric_result for _, result in all_dot_results
            )  # type: ignore[arg-type]
            assert isinstance(
                total_spenso_parametric_result,
                Expression,
            ), "Expected total parametric spenso result to be an Expression, got {}".format(type(total_spenso_parametric_result))
            total_spenso_parametric_timings = sum(
                (
                    result.spenso_parametric_timings
                    for _, result in all_dot_results
                    if result.spenso_parametric_timings is not None
                ),
                start=SpensoParametricTimingBreakdown(
                    graph_count=0,
                    emr_library_ms=0.0,
                    tensor_network_build_ms=0.0,
                    construction_ms=0.0,
                    per_graph_construction_samples=(),
                ),
            )
        if all_dot_results[0][1].dot_product_timings is not None:
            total_dot_product_timings = sum(
                (
                    result.dot_product_timings
                    for _, result in all_dot_results
                    if result.dot_product_timings is not None
                ),
                start=DotProductTimingBreakdown(
                    graph_count=0,
                    construction_ms=0.0,
                    per_graph_construction_samples=(),
                ),
            )

        cache_path = self._dot_products_cache_path(dot_file_path, args)
        total_dot_result: Expression | None = None
        if all_dot_results[0][1].dot_result is not None:
            total_dot_result = sum(
                result.dot_result for _, result in all_dot_results
            )  # type: ignore[arg-type]
            assert isinstance(
                total_dot_result,
                Expression,
            ), "Expected total dot product result to be an Expression, got {}".format(type(total_dot_result))
            cache_path.write_text(
                total_dot_result.to_canonical_string(),
                encoding="utf-8",
            )
        elif cache_path.is_file():
            total_dot_result = E(
                cache_path.read_text(encoding="utf-8").strip())

        if total_dot_result is not None:
            dot_expression_path, dot_expression_size_mb = self._save_expression_artifact(
                total_dot_result,
                dot_file_path,
                args,
                "dot_products",
            )
            self._record_ram_usage()
        if total_spenso_parametric_result is not None:
            spenso_parametric_expression_path, spenso_parametric_expression_size_mb = self._save_expression_artifact(
                total_spenso_parametric_result,
                dot_file_path,
                args,
                "spenso_parametric",
            )
            self._record_ram_usage()

        if total_dot_result is None and args.build_dot_products_form:
            print(
                "No dot-product expression available; skipping dot-product evaluator build.")

        assert reference_edge_momenta is not None
        self._validate_dot_parameter_coverage(reference_edge_momenta)
        edge_ids = self._get_sorted_edge_ids(reference_edge_momenta)
        dot_evaluator_run: EvaluatorRunResult | None = None
        if total_dot_result is not None:
            aligned_input_row = self._build_dot_product_input_row(
                reference_edge_momenta)
            dot_evaluator_run = self._run_compiled_evaluator(
                expr=total_dot_result,
                params=self.dot_variables_params,
                input_row=aligned_input_row,
                batch_size=self.batch_size_for_timing,
                label="dot-products",
                verbose=args.verbose_evaluator,
            )

        spenso_parametric_evaluator_run: EvaluatorRunResult | None = None
        if total_spenso_parametric_result is not None:
            aligned_emr_component_row = self._build_emr_component_input_row(
                reference_edge_momenta,
                edge_ids,
            )
            spenso_parametric_evaluator_run = self._run_compiled_evaluator(
                expr=total_spenso_parametric_result,
                params=self._build_emr_component_params(edge_ids),
                input_row=aligned_emr_component_row,
                batch_size=self.batch_size_for_timing,
                label="spenso-parametric",
                verbose=args.verbose_evaluator,
            )

        self._print_saved_expression_locations(
            dot_expression_path=dot_expression_path,
            spenso_parametric_expression_path=spenso_parametric_expression_path,
        )
        self._print_final_summary_table(
            dot_result=None if dot_evaluator_run is None else dot_evaluator_run.first_result,
            dot_evaluator_run=dot_evaluator_run,
            dot_product_timings=total_dot_product_timings,
            dot_expression_path=dot_expression_path,
            dot_expression_size_mb=dot_expression_size_mb,
            spenso_numeric_result=total_spenso_result,
            spenso_timings=total_spenso_timings,
            spenso_parametric_result=None if spenso_parametric_evaluator_run is None else spenso_parametric_evaluator_run.first_result,
            spenso_parametric_evaluator_run=spenso_parametric_evaluator_run,
            spenso_parametric_timings=total_spenso_parametric_timings,
            spenso_parametric_expression_path=spenso_parametric_expression_path,
            spenso_parametric_expression_size_mb=spenso_parametric_expression_size_mb,
            peak_ram_usage_mb=self.peak_ram_usage_mb,
        )

    def execute_graph(
        self,
        vertex_rules: dict[int, tuple[TensorNetwork | None, Expression]],
        edge_rules: dict[tuple[int | None, int | None], tuple[TensorNetwork | None, Expression]],
        order: list[int | tuple[int | None, int | None]],
        max_RAM: float | None,
        hep_lib: TensorLibrary | None,
        do_dot_products: bool = False,
        do_spenso: bool = True,
        progress_callback: Callable[[dict[str, Any]], None] | None = None,
        graph_name: str | None = None,
        component_label: str = "main",
        emit_step_logs: bool = True,
    ) -> GraphExecutionResult:
        muncher: TensorNetwork | None = None
        muncher_dot: Expression | None = None
        spenso_execution_s = 0.0
        dot_product_construction_s = 0.0

        def get_tn_for_step(step: int | tuple[int | None, int | None]) -> TensorNetwork:
            if isinstance(step, int):
                tn = vertex_rules[step][0]
            else:
                tn = edge_rules[tuple(
                    sorted([step[0], step[1]], key=lambda a: -1 if a is None else a))][0]
            if tn is None:
                raise ValueError(f"Missing TensorNetwork for step {step}")
            return tn

        def get_rule_for_step(step: int | tuple[int | None, int | None]) -> Expression:
            if isinstance(step, int):
                return vertex_rules[step][1]
            return edge_rules[tuple(sorted([step[0], step[1]], key=lambda a: -1 if a is None else a))][1]

        t_start = time.time()
        self._record_ram_usage()
        if progress_callback is not None:
            progress_callback({
                "event": "component_started",
                "graph_name": graph_name or "<unnamed>",
                "component": component_label,
                "current_step": 0,
                "total_steps": len(order),
                "ram_usage_mb": 0.0,
                "rank": 0,
            })
        for i_step, step in enumerate(order):
            ram_usage = Utils._current_ram_mb()
            self._record_ram_usage(ram_usage)
            if max_RAM is not None and ram_usage is not None and ram_usage > max_RAM:
                raise MemoryError(
                    f"RAM usage {ram_usage:.2f} MiB exceeded the specified "
                    f"limit of {max_RAM:.2f} MiB."
                )
            if muncher is None:
                current_rank = 0
            else:
                assert hep_lib is not None
                tensor_structure = muncher.result_tensor(
                    hep_lib).structure()
                current_rank = int(math.log(len(tensor_structure), 4) if len(
                    tensor_structure) > 0 else 0)
            ram_display = "n/a" if ram_usage is None else f"{ram_usage:.2f}"
            if emit_step_logs:
                print(
                    f"t={int(time.time() - t_start)
                         }s | Step {i_step + 1}: {step} "
                    f"[RAM: {ram_display} MiB, Rank: {current_rank}]"
                )

            spenso_step_start = time.time()
            if muncher is None:
                if do_spenso:
                    muncher = get_tn_for_step(step)
            else:
                if do_spenso:
                    muncher *= get_tn_for_step(step)
            if do_spenso:
                spenso_execution_s += time.time() - spenso_step_start

            if muncher_dot is None:
                if do_dot_products:
                    dot_step_start = time.time()
                    muncher_dot = get_rule_for_step(step)
                    dot_product_construction_s += time.time() - dot_step_start
            else:
                if do_dot_products:
                    dot_step_start = time.time()
                    muncher_dot *= get_rule_for_step(step)
                    dot_product_construction_s += time.time() - dot_step_start

            if do_spenso:
                spenso_step_start = time.time()
                muncher.execute()
                spenso_execution_s += time.time() - spenso_step_start
            if do_dot_products:
                dot_step_start = time.time()
                muncher_dot = to_dots(simplify_metrics(muncher_dot))
                dot_product_construction_s += time.time() - dot_step_start
                # print(muncher_dot)
            if progress_callback is not None:
                progress_callback({
                    "event": "step",
                    "graph_name": graph_name or "<unnamed>",
                    "component": component_label,
                    "current_step": i_step + 1,
                    "total_steps": len(order),
                    "ram_usage_mb": 0.0 if ram_usage is None else float(ram_usage),
                    "rank": current_rank,
                })

        if emit_step_logs:
            print(f"Execution completed in {time.time() - t_start:.3f}s")
        if do_spenso:
            assert muncher is not None
            assert hep_lib is not None
            n_entries = len(muncher.result_tensor(
                hep_lib).structure())
        else:
            n_entries = 1
        if n_entries > 1:
            raise ValueError(
                f"Warning: Final result has {n_entries} entries, expected 1")
        else:
            if do_spenso:
                assert muncher is not None
                spenso_step_start = time.time()
                muncher = muncher.result_scalar()
                spenso_execution_s += time.time() - spenso_step_start
                self._record_ram_usage()
        if progress_callback is not None:
            progress_callback({
                "event": "component_finished",
                "graph_name": graph_name or "<unnamed>",
                "component": component_label,
                "current_step": len(order),
                "total_steps": len(order),
                "ram_usage_mb": 0.0,
                "rank": 0,
            })

        return GraphExecutionResult(
            spenso_result=muncher,
            dot_result=muncher_dot,
            spenso_execution_ms=spenso_execution_s * 1000.0,
            dot_product_construction_ms=dot_product_construction_s * 1000.0,
        )

    def process_dot(
        self,
        dot_graph: pydot.Dot,
        args: argparse.Namespace,
        graph_display_index: int = 1,
        progress_callback: Callable[[dict[str, Any]], None] | None = None,
        emit_logs: bool = True,
    ) -> DotProcessingResult:
        t_start = time.time()
        self._record_ram_usage()
        edge_momenta = self._build_edge_momenta(dot_graph)
        edge_momenta_ms = (time.time() - t_start) * 1000.0
        edge_ids = self._get_sorted_edge_ids(edge_momenta)
        self._record_ram_usage()

        emr_library_ms = 0.0
        emr_hep_lib = None
        if not args.no_spenso:
            t_start = time.time()
            emr_hep_lib = self._build_emr_hep_library(edge_momenta)
            emr_library_ms = (time.time() - t_start) * 1000.0
            self._record_ram_usage()

        t_start = time.time()
        vertex_rules, edge_rules = self._build_graph_rules(
            dot_graph,
            emr_hep_lib,
            build_tensor_networks=not args.no_spenso,
        )
        tensor_network_build_ms = (
            time.time() - t_start) * 1000.0 if not args.no_spenso else 0.0
        self._record_ram_usage()
        graph_name = dot_graph.get_name()
        assert graph_name is not None
        graph_name = graph_name.strip('"')
        order, max_ram = self._get_graph_execution_plan(graph_name)
        res = self.execute_graph(
            vertex_rules=vertex_rules,
            edge_rules=edge_rules,
            order=order,
            do_dot_products=args.build_dot_products_form,
            do_spenso=not args.no_spenso,
            hep_lib=emr_hep_lib,
            max_RAM=max_ram,
            progress_callback=progress_callback,
            graph_name=graph_name,
            component_label="main",
            emit_step_logs=emit_logs,
        )
        self._record_ram_usage()

        spenso_parametric_result: Expression | None = None
        spenso_parametric_timings = None
        if args.build_spenso_parametric_form:
            t_start = time.time()
            parametric_emr_hep_lib = self._build_parametric_emr_hep_library(
                edge_ids)
            parametric_emr_library_ms = (time.time() - t_start) * 1000.0
            self._record_ram_usage()

            t_start = time.time()
            parametric_vertex_rules, parametric_edge_rules = self._build_graph_rules(
                dot_graph,
                parametric_emr_hep_lib,
                build_tensor_networks=True,
            )
            parametric_tn_build_ms = (time.time() - t_start) * 1000.0
            self._record_ram_usage()
            parametric_res = self.execute_graph(
                vertex_rules=parametric_vertex_rules,
                edge_rules=parametric_edge_rules,
                order=order,
                do_dot_products=False,
                do_spenso=True,
                hep_lib=parametric_emr_hep_lib,
                max_RAM=max_ram,
                progress_callback=progress_callback,
                graph_name=graph_name,
                component_label="parametric",
                emit_step_logs=emit_logs,
            )
            spenso_parametric_result = parametric_res.spenso_result
            self._record_ram_usage()
            spenso_parametric_timings = SpensoParametricTimingBreakdown(
                graph_count=1,
                emr_library_ms=parametric_emr_library_ms,
                tensor_network_build_ms=parametric_tn_build_ms,
                construction_ms=parametric_res.spenso_execution_ms,
                per_graph_construction_samples=(
                    (graph_display_index, parametric_res.spenso_execution_ms),),
            )

        if emit_logs:
            print(f"Returning result for graph {graph_name}")
        spenso_timings = None
        dot_product_timings = None
        if not args.no_spenso:
            spenso_timings = SpensoTimingBreakdown(
                graph_count=1,
                edge_momenta_ms=edge_momenta_ms,
                emr_library_ms=emr_library_ms,
                tensor_network_build_ms=tensor_network_build_ms,
                execution_ms=res.spenso_execution_ms,
            )
        if args.build_dot_products_form:
            dot_product_timings = DotProductTimingBreakdown(
                graph_count=1,
                construction_ms=res.dot_product_construction_ms,
                per_graph_construction_samples=(
                    (graph_display_index, res.dot_product_construction_ms),),
            )
        return DotProcessingResult(
            spenso_result=res.spenso_result,
            spenso_parametric_result=spenso_parametric_result,
            dot_result=res.dot_result,
            edge_momenta=edge_momenta,
            spenso_timings=spenso_timings,
            spenso_parametric_timings=spenso_parametric_timings,
            dot_product_timings=dot_product_timings,
        )


def _parallel_process_dot_worker(
    args: argparse.Namespace,
    dot_graph_data: str,
    graph_index: int,
    temp_dir: str,
    progress_queue: Any,
) -> ParallelDotWorkerResult:
    dot_graphs = pydot.graph_from_dot_data(dot_graph_data)
    assert dot_graphs is not None and len(
        dot_graphs) == 1, "Expected a single DOT graph payload"
    dot_graph = dot_graphs[0]
    processor = GraphProcessor(args)

    def progress_callback(event: dict[str, Any]) -> None:
        payload = dict(event)
        payload["graph_index"] = graph_index
        progress_queue.put(payload)

    res = processor.process_dot(
        dot_graph,
        args,
        graph_display_index=graph_index + 1,
        progress_callback=progress_callback,
        emit_logs=False,
    )

    dot_result_path: str | None = None
    if res.dot_result is not None:
        dot_path = Path(temp_dir) / f"dot_result__graph_{graph_index}.txt"
        dot_path.write_text(
            res.dot_result.to_canonical_string(), encoding="utf-8")
        dot_result_path = str(dot_path)

    spenso_parametric_result_path: str | None = None
    if res.spenso_parametric_result is not None:
        param_path = Path(temp_dir) / \
            f"spenso_parametric_result__graph_{graph_index}.txt"
        param_path.write_text(
            res.spenso_parametric_result.to_canonical_string(), encoding="utf-8")
        spenso_parametric_result_path = str(param_path)

    symbolic_expression_size_mb = 0.0
    if isinstance(res.spenso_result, Expression):
        symbolic_expression_size_mb += res.spenso_result.get_byte_size() / 1_000_000.0
    if res.spenso_parametric_result is not None:
        symbolic_expression_size_mb += res.spenso_parametric_result.get_byte_size() / \
            1_000_000.0
    if res.dot_result is not None:
        symbolic_expression_size_mb += res.dot_result.get_byte_size() / 1_000_000.0

    return ParallelDotWorkerResult(
        graph_index=graph_index,
        graph_name=(dot_graph.get_name() or "<unnamed>"),
        edge_momenta={
            edge_id: tuple(float(value) for value in momentum)
            for edge_id, momentum in res.edge_momenta.items()
        },
        spenso_result=(
            None
            if res.spenso_result is None
            else processor._coerce_scalar_to_complex(res.spenso_result)
        ),
        dot_result_path=dot_result_path,
        spenso_parametric_result_path=spenso_parametric_result_path,
        spenso_timings=res.spenso_timings,
        spenso_parametric_timings=res.spenso_parametric_timings,
        dot_product_timings=res.dot_product_timings,
        symbolic_expression_size_mb=symbolic_expression_size_mb,
    )


parser = argparse.ArgumentParser()
parser.add_argument("expr_file", nargs="?", default="expr_h.txt",
                    help="Path to the expression text file")
parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help="Enable stepped execution mode with verbose logging",
)
parser.add_argument(
    "--n_loops",
    type=int,
    default=2,
    help="number of loops",
)
parser.add_argument(
    "--n_edges",
    type=int,
    default=20,
    help="number of edges",
)
parser.add_argument(
    "--max_horner_scheme_variables",
    type=int,
    default=None,
    help="Max number of horner scheme variable in evaluator",
)
parser.add_argument(
    "--max_iterations_in_evaluator_building",
    type=int,
    default=None,
    help="Max number of evaluator optimizations rounds",
)
parser.add_argument(
    "--n_cpe_iterations",
    type=int,
    default=None,
    help="Max number of CPE iterations to consider",
)
parser.add_argument(
    "-ve",
    "--verbose_evaluator",
    action="store_true",
    default=False,
    help="Enable verbose evaluator construction",
)
parser.add_argument(
    "-ns",
    "--no_spenso",
    action="store_true",
    default=False,
    help="Disable spenso numerator construction",
)
parser.add_argument(
    "--n_graphs",
    type=int,
    default=None,
    help="number of graphs to process",
)
parser.add_argument(
    "--n_core_for_dot_processing",
    type=int,
    default=1,
    help="Number of cores for processing DOT graphs; set to 1 to keep the current serial behavior",
)
parser.add_argument(
    "--n_core_for_evaluator_building",
    type=int,
    default=1,
    help="Number of cores for building the evaluator",
)
parser.add_argument(
    "--sparse_momenta",
    action="store_true",
    help="Use sparse representation for loop momenta (only the non-zero components) instead of dense 4-vectors",
)
parser.add_argument(
    "-no_i",
    "--no_complex_i",
    action="store_true",
    default=False,
    help="Set complex i to 1",
)
parser.add_argument(
    "-c",
    "--concrete_loop_momenta",
    action="store_true",
    default=False,
    help="Execute with concrete loop momenta values (instead of symbolic parameters)",
)
parser.add_argument(
    "-fl",
    "--float_constants",
    action="store_true",
    default=False,
    help="Use floating point constants",
)
parser.add_argument(
    "-dots",
    "--build_dot_products_form",
    action="store_true",
    default=False,
    help="Build dot products expression",
)
parser.add_argument(
    "--build_spenso_parametric_form",
    action="store_true",
    default=False,
    help="Build a symbolic spenso scalar in terms of EMR momentum components",
)
parser.add_argument(
    "-flhep",
    "--float_hep_lib",
    action="store_true",
    default=False,
    help="Use floating point hep lib",
)
parser.add_argument(
    "-cose",
    "--concrete_oses",
    action="store_true",
    default=False,
    help="Execute with concrete ose values (instead of symbolic parameters)",
)
if __name__ == "__main__":
    args = parser.parse_args()

    processor = GraphProcessor(args)
    if args.expr_file.endswith(".txt"):
        processor.process_expression(
            E(Path(args.expr_file).read_text(encoding="utf-8").strip()), args)
    else:
        processor.process_dots(args.expr_file, args)
