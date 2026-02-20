#!/usr/bin/env python

from __future__ import annotations

import os
from time import perf_counter

import gamma_chain
import numpy as np
from symbolica import E, Expression
from symbolica.community.idenso import cook_indices
from symbolica.community.spenso import (
    LibraryTensor,
    Representation,
    TensorLibrary,
    TensorNetwork,
    TensorStructure,
)

MAX_LENGTH = int(gamma_chain.gamma_chain_max_length.maxlength)
BENCHMARK_EVALUATIONS = 2000
HEAVY_BENCHMARK_EVALUATIONS = 3
RESULT_TOL = 1e-9
RUN_HEAVY_SPENSO = os.environ.get("RUN_HEAVY_SPENSO", "0") == "1"
HEAVY_UNIQUE_DUMMIES = 7  # 10
HEAVY_POSITIVE_INSERTIONS = 5  # 64

vbar = [3.0, 3.1, 3.2, 3.3]
u = [4.0, 4.1, 4.2, 4.3]
indices = [1, -1, -2, 2, -1, -2, -3, -4, -5, -6, -5, -4, -6, -3]
vectors = [
    [1.0, 1.1, 1.2, 1.3],
    [2.0, 2.1, 2.2, 2.3],
]


def _as_symbolica_coeff(value: complex) -> Expression:
    if abs(value.imag) < 1e-12:
        return E(str(int(round(value.real))))
    if abs(value.real) < 1e-12:
        if abs(value.imag - 1.0) < 1e-12:
            return E("1i")
        if abs(value.imag + 1.0) < 1e-12:
            return E("-1i")
    raise ValueError(f"Unsupported gamma coefficient: {value}")


def _build_fortran_gamma_library() -> TensorLibrary:
    bis = Representation.bis(4)
    mink = Representation.mink(4)
    gamma_structure = TensorStructure(bis, bis, mink, name="spenso::gamma_custom")
    gamma_tensor = LibraryTensor.sparse(gamma_structure, Expression)

    ci = 1j
    gamma_matrices = np.zeros((4, 4, 4), dtype=np.complex128)
    gamma_matrices[0] = np.array(
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]],
        dtype=np.complex128,
    )
    gamma_matrices[1] = np.array(
        [[0, 0, 0, -1], [0, 0, -1, 0], [0, 1, 0, 0], [1, 0, 0, 0]],
        dtype=np.complex128,
    )
    gamma_matrices[2] = np.array(
        [[0, 0, 0, -ci], [0, 0, ci, 0], [0, ci, 0, 0], [-ci, 0, 0, 0]],
        dtype=np.complex128,
    )
    gamma_matrices[3] = np.array(
        [[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0]],
        dtype=np.complex128,
    )

    for mu in range(4):
        for row in range(4):
            for col in range(4):
                value = gamma_matrices[mu, row, col]
                if abs(value) > 1e-12:
                    gamma_tensor[row, col, mu] = _as_symbolica_coeff(value)

    library = TensorLibrary.hep_lib_atom()
    library.register(gamma_tensor)
    return library


def _build_spenso_chain_expression(chain_indices: list[int]) -> Expression:
    bis_ids = [1000 + i for i in range(len(chain_indices) + 1)]
    repeated_mink_ids: dict[int, int] = {}
    next_repeated_id = 1
    next_external_id = 200

    factors = [f"vbar(spenso::bis(4,idx({bis_ids[0]})))"]
    for position, mu_index in enumerate(chain_indices):
        if mu_index < 0:
            if mu_index not in repeated_mink_ids:
                repeated_mink_ids[mu_index] = next_repeated_id
                next_repeated_id += 1
            mink_id = repeated_mink_ids[mu_index]
        else:
            mink_id = next_external_id + position

        factors.append(
            "spenso::gamma_custom("
            f"spenso::bis(4,idx({bis_ids[position]})),"
            f"spenso::bis(4,idx({bis_ids[position + 1]})),"
            f"spenso::mink(4,idx({mink_id}))"
            ")"
        )
        if mu_index > 0:
            factors.append(f"vector({mu_index},spenso::mink(4,idx({mink_id})))")

    factors.append(f"u(spenso::bis(4,idx({bis_ids[-1]})))")
    return E(" * ".join(factors))


def _build_spenso_params(n_vectors: int) -> list[Expression]:
    params: list[Expression] = []
    for i in range(4):
        params.append(E(f"vbar(spenso::cind({i}))"))
    for i in range(4):
        params.append(E(f"u(spenso::cind({i}))"))
    for vector_id in range(1, n_vectors + 1):
        for i in range(4):
            params.append(E(f"vector({vector_id},spenso::cind({i}))"))
    return params


def _to_spenso_vector_components(vector: list[float]) -> np.ndarray:
    # This component map is required to reproduce compute_chain's gamma convention.
    return np.array([vector[0], -vector[1], vector[2], -vector[3]], dtype=np.complex128)


def _build_spenso_input_values(
    vbar_values: list[float],
    u_values: list[float],
    vector_values: list[list[float]],
) -> np.ndarray:
    total_size = 4 + 4 + 4 * len(vector_values)
    values = np.zeros(total_size, dtype=np.complex128)
    cursor = 0
    values[cursor : cursor + 4] = np.array(vbar_values, dtype=np.complex128)
    cursor += 4
    values[cursor : cursor + 4] = np.array(u_values, dtype=np.complex128)
    cursor += 4
    for vector in vector_values:
        values[cursor : cursor + 4] = _to_spenso_vector_components(vector)
        cursor += 4
    return values


def _benchmark(label: str, evaluations: int, fn) -> tuple[complex, float]:
    start = perf_counter()
    out = 0j
    for _ in range(evaluations):
        out = complex(fn())
    elapsed = perf_counter() - start
    print(f"{label:<40} total={elapsed:.6f}s per_eval={elapsed / evaluations:.6e}s")
    return out, elapsed


def _setup_spenso_for_case(
    chain_indices: list[int],
    vector_values: list[list[float]],
    vbar_values: list[float],
    u_values: list[float],
    spenso_library: TensorLibrary,
):
    setup_start = perf_counter()
    spenso_expression = _build_spenso_chain_expression(chain_indices)
    # print("Spenso expression:\n", spenso_expression)
    spenso_network = TensorNetwork(cook_indices(spenso_expression), spenso_library)
    spenso_network.execute(spenso_library)
    spenso_scalar_expression = spenso_network.result_scalar()
    # print("spenso_scalar_expression:\n", spenso_scalar_expression)
    spenso_params = _build_spenso_params(len(vector_values))
    spenso_evaluator = spenso_scalar_expression.evaluator(
        constants={},
        functions={},
        params=spenso_params,
        iterations=100,
        n_cores=1,
        external_functions=None,
        conditionals=[],
        verbose=False,
        cpe_iterations=1000,
    )
    spenso_evaluator = spenso_evaluator.compile(
        "gamma_chain", "test_gamama_chain.cpp", "test_gamama_chain.so", inline_asm="default", optimization_level=3, number_type="complex"
    )
    spenso_inputs = _build_spenso_input_values(vbar_values, u_values, vector_values)
    spenso_setup_elapsed = perf_counter() - setup_start
    return spenso_evaluator, spenso_inputs, spenso_setup_elapsed


def _build_heavy_case() -> tuple[list[float], list[float], list[int], list[list[float]]]:
    negative_indices = list(range(-1, -HEAVY_UNIQUE_DUMMIES - 1, -1))
    negative_indices.extend(list(range(-HEAVY_UNIQUE_DUMMIES, 0)))
    positive_indices = [(i % 3) + 1 for i in range(HEAVY_POSITIVE_INSERTIONS)]

    heavy_indices: list[int] = []
    neg_pointer = 0
    pos_pointer = 0
    while neg_pointer < len(negative_indices) or pos_pointer < len(positive_indices):
        if pos_pointer < len(positive_indices):
            heavy_indices.append(positive_indices[pos_pointer])
            pos_pointer += 1
        if neg_pointer < len(negative_indices):
            heavy_indices.append(negative_indices[neg_pointer])
            neg_pointer += 1
        if neg_pointer < len(negative_indices):
            heavy_indices.append(negative_indices[neg_pointer])
            neg_pointer += 1

    heavy_indices = heavy_indices[:MAX_LENGTH]

    heavy_vectors = []
    for vector_id in range(1, 4):
        base = 0.25 * vector_id
        heavy_vectors.append([base + 0.13, base + 0.27, base + 0.41, base + 0.53])

    return [3.0, 3.1, 3.2, 3.3], [4.0, 4.1, 4.2, 4.3], heavy_indices, heavy_vectors


print("Testing benchmark chain of 14 gamma matrices:")
gamma_chain_result = gamma_chain.compute_chain(vbar, u, indices, vectors)
target_result = complex(-78354.8416 + 1312.256j)
print(gamma_chain_result)
print("Target benchmark result:")
print(target_result)
print("gamma_chain-target abs diff:", abs(gamma_chain_result - target_result))

spenso_lib = _build_fortran_gamma_library()
spenso_evaluator, spenso_inputs, spenso_setup_elapsed = _setup_spenso_for_case(indices, vectors, vbar, u, spenso_lib)
spenso_result = complex(spenso_evaluator.evaluate(spenso_inputs).item())
print("\nSpenso + Symbolica result:")
print(spenso_result)
print("spenso-gamma_chain abs diff:", abs(spenso_result - gamma_chain_result))
print("Results match:", abs(spenso_result - gamma_chain_result) < RESULT_TOL)

print("\nPerformance comparison:")
print(f"Spenso setup (library + TN + evaluator): {spenso_setup_elapsed:.6f}s")

_ = gamma_chain.compute_chain(vbar, u, indices, vectors)
_ = spenso_evaluator.evaluate(spenso_inputs)

_, gamma_eval_total = _benchmark(
    "gamma_chain.compute_chain",
    BENCHMARK_EVALUATIONS,
    lambda: gamma_chain.compute_chain(vbar, u, indices, vectors),
)
_, spenso_eval_total = _benchmark(
    "spenso evaluator.evaluate_complex",
    BENCHMARK_EVALUATIONS,
    lambda: spenso_evaluator.evaluate(spenso_inputs).item(),
)

print(f"speedup (gamma/spenso eval only): {gamma_eval_total / spenso_eval_total:.2f}x")
# print(
#     "spenso total including setup for this run:",
#     f"{spenso_setup_elapsed + spenso_eval_total:.6f}s",
# )

heavy_vbar, heavy_u, heavy_indices, heavy_vectors = _build_heavy_case()
heavy_unique_dummies = len({idx for idx in heavy_indices if idx < 0})
print("\nHeavy internal-contraction chain:")
print(f"length={len(heavy_indices)} unique_dummy_indices={heavy_unique_dummies} n_vectors={len(heavy_vectors)}")

single_start = perf_counter()
heavy_result = gamma_chain.compute_chain(heavy_vbar, heavy_u, heavy_indices, heavy_vectors)
single_elapsed = perf_counter() - single_start
print("gamma_chain heavy result:", heavy_result)
print(f"gamma_chain heavy single eval: {single_elapsed:.6f}s")

_, heavy_gamma_eval_total = _benchmark(
    "gamma_chain.compute_chain (heavy)",
    HEAVY_BENCHMARK_EVALUATIONS,
    lambda: gamma_chain.compute_chain(heavy_vbar, heavy_u, heavy_indices, heavy_vectors),
)

if RUN_HEAVY_SPENSO:
    heavy_spenso_evaluator, heavy_spenso_inputs, heavy_spenso_setup_elapsed = _setup_spenso_for_case(
        heavy_indices, heavy_vectors, heavy_vbar, heavy_u, spenso_lib
    )
    heavy_spenso_result = complex(heavy_spenso_evaluator.evaluate(heavy_spenso_inputs).item())
    print("heavy spenso result:", heavy_spenso_result)
    print("heavy spenso-gamma_chain abs diff:", abs(heavy_spenso_result - heavy_result))
    print(f"heavy spenso setup time: {heavy_spenso_setup_elapsed:.6f}s")
    _, heavy_spenso_eval_total = _benchmark(
        "spenso evaluator.evaluate_complex (heavy)",
        BENCHMARK_EVALUATIONS,
        lambda: heavy_spenso_evaluator.evaluate(heavy_spenso_inputs).item(),
    )
    print(
        "heavy speedup (gamma/spenso eval only):",
        f"{heavy_gamma_eval_total / HEAVY_BENCHMARK_EVALUATIONS / (heavy_spenso_eval_total / BENCHMARK_EVALUATIONS):.2f}x",
    )
else:
    print("Skipping heavy Spenso build/eval by default (set RUN_HEAVY_SPENSO=1 to enable).")
