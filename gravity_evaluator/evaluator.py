import argparse
import math
import os
import subprocess
import time
from asyncio.unix_events import SelectorEventLoop
from multiprocessing import Array
from pathlib import Path

import numpy as np
import pydot
from numpy._core.numerictypes import float64
from numpy.typing import NDArray
from python_utils.types import Number
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
                print("[{:.3f}s] Performing scalar step    # {}".format(elapsed_s, i_step))
                ram_mb = Utils._current_ram_mb()
                if ram_mb is not None:
                    elapsed_s = time.time() - t_start_fn
                    print("[{:.3f}s] Current RAM usage: {:.2f} MiB".format(elapsed_s, ram_mb))
            if t_delta is not None:
                time.sleep(t_delta)
            # print(tn)
            # print(i_step)
            tn.execute(n_steps=1, mode=ExecutionMode.Scalar, library=hep_lib)

            if verbose:
                elapsed_s = time.time() - t_start_fn
                print("[{:.3f}s] Performing reduction step # {}".format(elapsed_s, i_step))

            # i_step += 1

            # print(tn)
            # print(i_step)
            tn.execute(n_steps=1, mode=ExecutionMode.Single, library=hep_lib)
            if verbose:
                elapsed_s = time.time() - t_start_fn
                print("[{:.3f}s] DONE. size at this point: {}".format(elapsed_s, len(str(tn))))
                ram_mb = Utils._current_ram_mb()
                if ram_mb is not None:
                    elapsed_s = time.time() - t_start_fn
                    print("[{:.3f}s] Current RAM usage: {:.2f} MiB".format(elapsed_s, ram_mb))

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


class GraphProcessor:
    def __init__(self, args: argparse.Namespace):
        self.loop_momenta: list[NDArray[np.float64]] = []
        for i_loop in range(args.n_loops):
            self.loop_momenta.append(np.array([2.0 / float(i_loop + 1), 3.0 / float(i_loop + 1), 4.0 / float(i_loop + 1)], dtype=np.float64))

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

        self.oses = [(ie + 1) / float(args.n_edges) for ie in range(args.n_edges)]
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

        mink_rep = Representation.mink(4)

        # External momenta
        if args.float_hep_lib:
            self.hep_lib_with_momenta = TensorLibrary.hep_lib()
        else:
            self.hep_lib_with_momenta = TensorLibrary.hep_lib_atom()
        mink_rep = Representation.mink(4)

        # External momenta
        self.params = []
        self.param_values = []
        vector_name_p = S("P")
        for i in range(len(self.external_momenta)):
            p_structure = TensorStructure(mink_rep, N(i), name=vector_name_p)  # type: ignore
            if args.concrete_loop_momenta:
                if args.sparse_momenta:
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
                if args.sparse_momenta:
                    p = LibraryTensor.sparse(p_structure, type_info=Expression)
                    for j in range(1, 4):
                        p[j] = E("P({},spenso::cind({}))".format(i, j))
                else:
                    p = LibraryTensor.dense(
                        p_structure,
                        [
                            E("0"),
                        ]
                        + [E("P({},spenso::cind({}))".format(i, j)) for j in range(1, 4)],
                    )
                self.params.extend(E("P({},spenso::cind({}))".format(i, j)) for j in range(1, 4))
                self.param_values.extend(self.external_momenta[i][j] for j in range(1, 4))
            self.hep_lib_with_momenta.register(p)

        # Loop momenta
        vector_name_k = S("K")
        if args.concrete_loop_momenta:
            for i, k in enumerate(self.loop_momenta):
                ki_structure = TensorStructure(mink_rep, N(i), name=vector_name_k)  # pyright: ignore
                if args.sparse_momenta:
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
                ki_structure = TensorStructure(mink_rep, N(i), name=vector_name_k)  # pyright: ignore
                if args.sparse_momenta:
                    ki = LibraryTensor.sparse(ki_structure, type_info=Expression)
                    for j in range(1, 4):
                        ki[j] = E("K({},spenso::cind({}))".format(i, j))
                else:
                    ki = LibraryTensor.dense(
                        ki_structure,
                        [
                            E("0"),
                        ]
                        + [E("K({},spenso::cind({}))".format(i, j)) for j in range(1, 4)],
                    )
                self.hep_lib_with_momenta.register(ki)
            for i, k in enumerate(self.loop_momenta):
                for j in range(1, 4):
                    self.params.append(E("K({},spenso::cind({}))".format(i, j)))
                    self.param_values.append(k[j - 1])

        # OSEs
        self.OSE_vars = []
        for ie, ose in enumerate(self.oses):
            ose_structure = TensorStructure(mink_rep, N(ie), name=S("OSE"))  # type: ignore
            ose_var = E(f"OSE({ie})")
            self.OSE_vars.append(ose_var)
            if args.concrete_oses:
                ose_tensor = LibraryTensor.sparse(ose_structure, type_info=float)
                ose_tensor[0] = ose
            else:
                ose_tensor = LibraryTensor.sparse(ose_structure, type_info=Expression)
                ose_tensor[0] = ose_var
                self.params.append(ose_var)
                self.param_values.append(ose)
            self.hep_lib_with_momenta.register(ose_tensor)

        self.dot_variables_params = []
        for i_e in range(args.n_edges):
            for j_e in range(args.n_edges):
                if j_e >= i_e:
                    self.dot_variables_params.append(E(f"spenso::dot(spenso::mink(4),Q({i_e}),Q({j_e}))"))

        self.evaluator_options = {
            "constants": {},
            "params": self.params,
            "functions": {},
            "verbose": args.verbose_evaluator,
            "cpe_iterations": args.n_cpe_iterations,
            "n_cores": args.n_core_for_evaluator_building,
        }
        if args.max_horner_scheme_variables is not None:
            self.evaluator_options["max_horner_scheme_variables"] = args.max_horner_scheme_variables
        if args.max_iterations_in_evaluator_building is not None:
            self.evaluator_options["iterations"] = args.max_iterations_in_evaluator_building

    def process_expression(self, expr: Expression, args: argparse.Namespace):
        raise NotImplementedError("Expression processing is no longer up to date, only DOT graph processing is supported for now.")
        numerator = expr.replace_multiple(self.global_replacement_rules, repeat=True)

        n_random = len(self.params)
        batch_size = 1000000
        rng = np.random.default_rng(12345)

        # Random rank-2 array (batch_size, n_random) of float64 values.
        sample = rng.random((batch_size, n_random), dtype=np.float64)

        # Random rank-2 array (batch_size, n_random) of complex128 values,
        # with independent random real/imaginary parts for every batch sample.
        sample_complex = rng.random((batch_size, n_random), dtype=np.float64) + 1j * rng.random((batch_size, n_random), dtype=np.float64)

        if args.zero_last_sample:
            sample[-1, :] = 0.0
            sample_complex[-1, :] = 0.0 + 0.0j

        tn_numerical_to_execute = TensorNetwork(cook_indices(numerator), self.hep_lib_with_momenta)

        t_start = time.time()
        if args.verbose:
            Utils.stepped_execution(
                tn_numerical_to_execute,
                self.hep_lib_with_momenta,
                t_delta=0.0,
                verbose=True,
            )
        else:
            tn_numerical_to_execute.execute(self.hep_lib_with_momenta)

        numerical_evaluation = tn_numerical_to_execute.result_scalar()
        t_end = time.time()
        print("Execution time: {:.2f} seconds".format(t_end - t_start))
        t_start = time.time()
        eval = numerical_evaluation.evaluator(constants={}, params=self.params, functions={})
        print("Evaluation time evaluator: {:.2f} ms".format((time.time() - t_start) * 1000))

        t_start = time.time()
        eval.evaluate(sample * 100.0)
        print("Evaluation time real: {:.2f} mus".format((time.time() - t_start) * 1000_000 / float(batch_size)))

        t_start = time.time()
        eval.evaluate_complex(sample_complex)
        print("Evaluation time complex: {:.2f} mus".format((time.time() - t_start) * 1000_000 / float(batch_size)))

    def process_dots(self, dot_file_path: str, args: argparse.Namespace):
        dot_graphs = pydot.graph_from_dot_file(dot_file_path)
        assert dot_graphs is not None, "Failed to parse DOT file: {}".format(dot_file_path)
        all_dot_results: list[tuple[pydot.Dot, tuple[Expression | None, Expression | None]]] = []
        selected_dot_graphs = dot_graphs if args.n_graphs is None else dot_graphs[: args.n_graphs]
        # print(self.params)
        # stopme
        for dot_graph in selected_dot_graphs:
            print("Processing graph: {}".format(dot_graph.get_name()))
            res = self.process_dot(dot_graph, args)
            if res[0] is not None:
                try:
                    numerical_res = complex(str(res[0]))
                    print("Numerical result: {}".format(numerical_res))
                except Exception:
                    pass
            memory = 0
            if res[0] is not None:
                memory += res[0].get_byte_size()
            if res[1] is not None:
                memory += res[1].get_byte_size()
            print(f"Symbolic expression size: {memory / 1000000.0:.2f} MB")
            all_dot_results.append((dot_graph, res))

        if all_dot_results[0][1][0] is not None:
            total_res = sum(res[0] for _, res in all_dot_results)  # type: ignore
            assert isinstance(total_res, Expression), "Expected total result to be an Expression, got {}".format(type(total_res))

            DO_POLYNOMIALIZE = False
            if DO_POLYNOMIALIZE:
                print("Computing polynomial....")
                total_res_poly = total_res.rationalize(relative_error=1e-10).to_polynomial(vars=self.OSE_vars[:1], minimal_poly=P("y^2+1"))
                print("Done!....")
                print([shape for (shape, coeff) in total_res_poly.coefficient_list()])

            TEST_OSE_RAISINGS = False
            if TEST_OSE_RAISINGS:
                # Test OSE raisings
                print("Building OSE coefficients...")
                loop_ose_vars = self.OSE_vars[len(self.external_momenta) :]
                loop_ose_vars = self.OSE_vars[len(self.external_momenta) : len(self.external_momenta) + 1]
                total_res_coefficients = total_res.coefficient_list(*loop_ose_vars)
                shapes = [shape.to_polynomial(vars=loop_ose_vars).coefficient_list()[0][0] for (shape, coeff) in total_res_coefficients]
                print(shapes)
                max_ose_power = max(max(s) for s in shapes)
                if max_ose_power <= 1:
                    print("SUCCESS! No OSE raisings found!")

            DO_SPENSO_EXPRESSION_EVALUATOR = False
            if DO_SPENSO_EXPRESSION_EVALUATOR:
                print("Starting evaluator building...")
                t_start = time.time()
                eager_eval = total_res.evaluator(**self.evaluator_options)
                print("Evaluation time evaluator: {:.2f} ms".format((time.time() - t_start) * 1000))

                print("Starting evaluator compilation...")
                t_start = time.time()
                eval = eager_eval.compile("test", "test.cpp", "test.so", "complex", inline_asm="default")
                print("Evaluator compilation time: {:.2f} ms".format((time.time() - t_start) * 1000))

                n_random = len(self.params)
                batch_size = 1000
                rng = np.random.default_rng(12345)

                # Random rank-2 array (batch_size, n_random) of float64 values.
                sample = rng.random((batch_size, n_random), dtype=np.float64)

                # Random rank-2 array (batch_size, n_random) of complex128 values,
                # with independent random real/imaginary parts for every batch sample.
                sample_complex = rng.random((batch_size, n_random), dtype=np.float64) + 1j * rng.random((batch_size, n_random), dtype=np.float64)

                t_start = time.time()
                eval.evaluate(sample_complex)
                print("Evaluation time complex: {:.2f} mus".format((time.time() - t_start) * 1000_000 / float(batch_size)))

        IS_DOT_PRODUCT_RESULT_IS_PRESENT = False
        if all_dot_results[0][1][1] is not None:
            IS_DOT_PRODUCT_RESULT_IS_PRESENT = True
            total_res = sum(res[1] for _, res in all_dot_results)  # type: ignore
            assert isinstance(total_res, Expression), "Expected total dot product result to be an Expression, got {}".format(type(total_res))
            if not os.path.isfile("./dot_products_results.txt"):
                with open("./dot_products_results.txt", "w") as f:
                    f.write(total_res.to_canonical_string())

        elif os.path.isfile("./dot_products_results.txt"):
            IS_DOT_PRODUCT_RESULT_IS_PRESENT = True
            with open("./dot_products_results.txt", "r") as f:
                total_res = E(f.read())

        if IS_DOT_PRODUCT_RESULT_IS_PRESENT:
            self.evaluator_options["params"] = self.dot_variables_params
            print("Starting evaluator building...")
            t_start = time.time()
            eager_eval = total_res.evaluator(**self.evaluator_options)
            print("Evaluation time evaluator: {:.2f} ms".format((time.time() - t_start) * 1000))

            print("Starting evaluator compilation...")
            t_start = time.time()
            eval = eager_eval.compile("test", "test.cpp", "test.so", "complex", inline_asm="default")
            print("Evaluator compilation time: {:.2f} ms".format((time.time() - t_start) * 1000))

            n_random = len(self.dot_variables_params)
            batch_size = 1000
            rng = np.random.default_rng(12345)

            # Random rank-2 array (batch_size, n_random) of float64 values.
            sample = rng.random((batch_size, n_random), dtype=np.float64)

            # Random rank-2 array (batch_size, n_random) of complex128 values,
            # with independent random real/imaginary parts for every batch sample.
            sample_complex = rng.random((batch_size, n_random), dtype=np.float64) + 1j * rng.random((batch_size, n_random), dtype=np.float64)

            t_start = time.time()
            eval.evaluate(sample_complex)
            print("Evaluation time complex: {:.2f} mus".format((time.time() - t_start) * 1000_000 / float(batch_size)))

    def execute_graph(
        self, vertex_rules, edge_rules, order, max_RAM: float | None, do_dot_products=False, do_spenso=True
    ) -> tuple[Expression | None, Expression | None]:
        muncher = None
        muncher_dot = None

        def get_tn_for_step(step):
            if isinstance(step, int):
                return vertex_rules[step][0]
            else:
                return edge_rules[tuple(sorted([step[0], step[1]], key=lambda a: -1 if a is None else a))][0]

        def get_rule_for_step(step):
            if isinstance(step, int):
                return vertex_rules[step][1]
            else:
                return edge_rules[tuple(sorted([step[0], step[1]], key=lambda a: -1 if a is None else a))][1]

        t_start = time.time()
        for i_step, step in enumerate(order):
            ram_usage = Utils._current_ram_mb()
            if max_RAM is not None and ram_usage is not None and ram_usage > max_RAM:
                raise MemoryError(f"RAM usage {ram_usage:.2f} MiB exceeded the specified limit of {max_RAM:.2f} MiB.")
            if muncher is None:
                current_rank = 0
            else:
                tensor_structure = muncher.result_tensor(self.hep_lib_with_momenta).structure()
                current_rank = int(math.log(len(tensor_structure), 4) if len(tensor_structure) > 0 else 0)
            print(f"t={int(time.time() - t_start)}s | Step {i_step + 1}: {step} [RAM: {ram_usage:.2f} MiB, Rank: {current_rank}]")

            if muncher is None:
                if do_spenso:
                    muncher = get_tn_for_step(step)
            else:
                if do_spenso:
                    muncher *= get_tn_for_step(step)

            if muncher_dot is None:
                if do_dot_products:
                    muncher_dot = get_rule_for_step(step)
            else:
                if do_dot_products:
                    muncher_dot *= get_rule_for_step(step)

            if do_spenso:
                muncher.execute()
            if do_dot_products:
                muncher_dot = to_dots(simplify_metrics(muncher_dot))
                # print(muncher_dot)

        print(f"Execution completed in {time.time() - t_start:.3f}s")
        if do_spenso:
            n_entries = len(muncher.result_tensor(self.hep_lib_with_momenta).structure())
        else:
            n_entries = 1
        if n_entries > 1:
            raise ValueError(f"Warning: Final result has {n_entries} entries, expected 1")
        else:
            if do_spenso:
                muncher = muncher.result_scalar()

        return (muncher, muncher_dot)

    def process_dot(self, dot_graph: pydot.Dot, args: argparse.Namespace) -> tuple[Expression | None, Expression | None]:
        emr_to_lmb_replacements = []
        head = "Q"
        emr_to_ose_replacements = []
        for edge in dot_graph.get_edges():
            e_id = int(edge.get("id").strip('"'))
            emr_to_ose_replacements.append(Replacement(E(f"{head}({e_id},a___)"), E(f"OSE({e_id},a___)+{head}({e_id},a___)")))
            emr_to_lmb_replacements.append(Replacement(E(f"{head}({e_id},a___)"), E(edge.get("lmb_rep").strip('"'))))

        self.vertex_rules = {}
        self.edge_rules = {}
        for node in dot_graph.get_nodes():
            node_id = Utils.get_dot_node_id_from_name(node.get_name())
            if node_id is None:
                continue
            str_expr = node.get_attributes().get("num", "1").strip('"')
            for old, new in self.string_replacements:
                str_expr = str_expr.replace(old, new)
            rule = E(str_expr)
            rule = rule.replace_multiple(self.global_replacement_rules, repeat=True)
            rule_for_tn = rule.replace_multiple(emr_to_ose_replacements)
            rule_for_tn = rule_for_tn.replace_multiple(emr_to_lmb_replacements, repeat=True)
            tn = TensorNetwork(cook_indices(rule_for_tn), self.hep_lib_with_momenta)
            tn.execute(self.hep_lib_with_momenta)
            self.vertex_rules[node_id] = (tn, rule)
        for edge in dot_graph.get_edges():
            src = Utils.get_dot_node_id_from_name(edge.get_source())
            dst = Utils.get_dot_node_id_from_name(edge.get_destination())
            str_expr = edge.get_attributes().get("num", "1").strip('"')
            for old, new in self.string_replacements:
                str_expr = str_expr.replace(old, new)
            rule = E(str_expr)
            rule = rule.replace_multiple(self.global_replacement_rules, repeat=True)
            rule_for_tn = rule.replace_multiple(emr_to_ose_replacements)
            rule_for_tn = rule_for_tn.replace_multiple(emr_to_lmb_replacements, repeat=True)
            tn = TensorNetwork(cook_indices(rule_for_tn), self.hep_lib_with_momenta)
            tn.execute(self.hep_lib_with_momenta)
            self.edge_rules[tuple(sorted([src, dst], key=lambda a: -1 if a is None else a))] = (tn, rule)
        # print(self.vertex_rules)
        # print(self.edge_rules)
        # stop

        graph_name = dot_graph.get_name()
        assert graph_name is not None
        graph_name = graph_name.strip('"')
        if graph_name.startswith("h_diagram"):
            # h_diagram processing
            res = self.execute_graph(
                vertex_rules=self.vertex_rules,
                edge_rules=self.edge_rules,
                # order=[
                #     # First the bottom blackhole line
                #     (None, 0),
                #     0,
                #     (0, 1),
                #     1,
                #     (None, 1),
                #     # Then the second blackhole line
                #     (None, 2),
                #     2,
                #     (2, 3),
                #     3,
                #     (None, 3),
                #     # Then the graviton propagators leaving the first blackhole line
                #     (0, 4),
                #     (1, 5),
                #     # Then the graviton propagators leaving the second blackhole line
                #     (0, 4),
                #     (1, 5),
                #     # Then the first graviton feynman rule
                #     5,
                #     # Then the internal graviton propagator
                #     (4, 5),
                #     # Finally the second graviton feynman rule
                #     4,
                # ],
                order=[
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
                ],
                do_dot_products=args.build_dot_products_form,
                do_spenso=not args.no_spenso,
                max_RAM=15000.0,  # Set a RAM limit of 8000 MiB for the execution of the graph. Adjust this value based on your system's capabilities and requirements
            )
        elif graph_name.startswith("hh_diagram"):
            # hh_diagram processing
            res = self.execute_graph(
                vertex_rules=self.vertex_rules,
                edge_rules=self.edge_rules,
                order=[
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
                ],
                do_dot_products=args.build_dot_products_form,
                do_spenso=not args.no_spenso,
                max_RAM=40000.0,  # Set a RAM limit of 8000 MiB for the execution of the graph. Adjust this value based on your system's capabilities and requirements
            )
        else:
            raise NotImplementedError(f"DOT processing is not implemented yet for diagram '{graph_name}'.")

        print(f"Returning result for graph {graph_name}")
        return res


parser = argparse.ArgumentParser()
parser.add_argument("expr_file", nargs="?", default="expr_h.txt", help="Path to the expression text file")
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
        processor.process_expression(E(Path(args.expr_file).read_text(encoding="utf-8").strip()), args)
    else:
        processor.process_dots(args.expr_file, args)
