from __future__ import annotations

import itertools as it
import multiprocessing as mp
from pathlib import Path
from typing import Union

import more_itertools as itx
import networkx as nx
import numpy as np
import pandas as pd
from networkx.utils import py_random_state
from snakeboost import snakemake_args

from notebooks.adjacency_matrix import AdjacencyMatrix
from notebooks.bn_metadata import read_metadata
from notebooks.utils import underscore as __


def filter_logile(bin: int, num_bins: int = 10):
    def inner(matrix):
        if bin >= num_bins:
            raise ValueError("bin must be less then num_bins")
        masked = np.ma.masked_equal(matrix, 0)
        log = np.ma.log10(masked)
        threshold = 10 ** (((log.max() - log.min()) * bin / num_bins) + log.min())
        return matrix < threshold

    return inner


@py_random_state(3)
def double_edge_swap(G, nswap=1, max_tries=100, seed=None):
    """Swap two edges in the graph while keeping the node degrees fixed.

    A double-edge swap removes two randomly chosen edges u-v and x-y
    and creates the new edges u-x and v-y::

     u--v            u  v
            becomes  |  |
     x--y            x  y

    If either the edge u-x or v-y already exist no swap is performed
    and another attempt is made to find a suitable edge pair.

    Parameters
    ----------
    G : graph
       An undirected graph

    nswap : integer (optional, default=1)
       Number of double-edge swaps to perform

    max_tries : integer (optional)
       Maximum number of attempts to swap edges

    seed : integer, random_state, or None (default)
        Indicator of random number generation state.
        See :ref:`Randomness<randomness>`.

    Returns
    -------
    G : graph
       The graph after double edge swaps.

    Notes
    -----
    Does not enforce any connectivity constraints.

    The graph G is modified in place.
    """
    if G.is_directed():
        raise nx.NetworkXError("double_edge_swap() not defined for directed graphs.")
    if nswap > max_tries:
        raise nx.NetworkXError("Number of swaps > number of tries allowed.")
    if len(G) < 4:
        raise nx.NetworkXError("Graph has less than four nodes.")
    # Instead of choosing uniformly at random from a generated edge list,
    # this algorithm chooses nonuniformly from the set of nodes with
    # probability weighted by degree.
    n = 0
    swapcount = 0
    keys, degrees = zip(*G.degree())  # keys, degree
    cdf = nx.utils.cumulative_distribution(degrees)  # cdf of degree
    discrete_sequence = nx.utils.discrete_sequence
    while swapcount < nswap:
        #        if random.random() < 0.5: continue # trick to avoid periodicities?
        # pick two random edges without creating edge list
        # choose source node indices from discrete distribution
        (ui, xi) = discrete_sequence(2, cdistribution=cdf, seed=seed)
        if ui == xi:
            continue  # same source, skip
        u = keys[ui]  # convert index to label
        x = keys[xi]
        # choose target uniformly from neighbors
        v = seed.choice(list(G[u]))
        y = seed.choice(list(G[x]))
        if v == y:
            continue  # same target, skip
        swaps = [
            ((x, y), (u, x)),
            ((u, v), (v, y)),
            ((u, x), (x, y)),
            ((v, y), (u, v)),
        ]
        attrs = [G.edges[src] if src in G.edges else None for src, _ in swaps]
        if (x in G[u]) ^ (y in G[v]):
            continue
        for (_, dest), attr in zip(swaps, attrs):
            if dest in G.edges:
                G.remove_edge(*dest)
            if attr is None:
                continue
            G.add_edge(*dest, **attr)
        swapcount += 1
        # if (x not in G[u]) and (y not in G[v]):  # don't create parallel edges
        #     G.add_edge(u, x, **G.edges[x, y])
        #     G.add_edge(v, y, **G.edges[u, v])
        #     G.remove_edge(u, v)
        #     G.remove_edge(x, y)
        #     swapcount += 1

        if n >= max_tries:
            e = (
                f"Maximum number of swap attempts ({n}) exceeded "
                f"before desired swaps achieved ({nswap})."
            )
            raise nx.NetworkXAlgorithmError(e)
        n += 1
    return G


get_k_values = lambda G: __.pipe(
    [degree for n, degree in G.degree],
    np.sort,
    np.unique,
)


def _weighted_rich_club(G):
    w_ranked = np.sort([edge["weight"] for edge in G.edges.values()])[::-1]
    k_values = get_k_values(G)
    results = np.empty(len(k_values))
    for rank, degree_thresh in enumerate(k_values):
        sG = G.subgraph([node for node, degree in G.degree if degree > degree_thresh])
        rich_weight = sum(edge["weight"] for edge in sG.edges.values())
        global_weight = w_ranked[: len(sG.edges)].sum()
        if global_weight:
            results[rank] = rich_weight / global_weight
        else:
            results[rank] = np.nan
    return results


def _randomize_graph(G, Q=100, seed=None):
    R = G.copy()
    E = R.number_of_edges()
    double_edge_swap(R, Q * E, max_tries=Q * E * 10, seed=seed)
    return _weighted_rich_club(R)


def weighted_rich_club(G, normalized=True, m=1000, Q=100, seed=None, threads=None):
    result = _weighted_rich_club(G)
    k_values = get_k_values(G)
    if normalized and m:
        randoms = np.empty((m, len(k_values)))
        with mp.Pool(threads) as pool:
            random_Gs = [
                pool.apply_async(_randomize_graph, [_G]) for _G in it.repeat(G, m)
            ]
            for i, R in enumerate(random_Gs):
                randoms[i] = R.get()
        result = result / randoms.mean(axis=0)

    return np.vstack([k_values, result]).T


def run(
    src: Path, filter_level: int, num_randoms: int, threads: Union[int, None] = None
):
    adj = (
        AdjacencyMatrix(
            raw=np.genfromtxt(src, delimiter=","),
            metadata=read_metadata(),
        )
        .mask_diagonal()
        .mask_equal(0)
        .mask_where(filter_logile(filter_level))
    )
    adj.props["distance"] = np.ma.filled(1 / adj.raw, np.NaN)
    if num_randoms:
        print(f"Normalizing against {num_randoms} random graphs")
    return pd.DataFrame(
        weighted_rich_club(adj.graph, normalized=True, m=num_randoms, threads=threads),
        columns=["k", "phi"],
    )

def main():
    args = snakemake_args(
        input=["in"],
        output=["out"],
        params={"filter_level": "--filter-level", "normalization": "--normalization"},
    )

    if isinstance(args.input, dict):
        raise TypeError("Inputs must be specified as a single item")
    if isinstance(args.output, dict):
        raise TypeError("Outputs must be specified as a single item")
    input = itx.one(args.input)
    output = itx.one(args.output)
    threads = args.threads or None
    params = args.params if isinstance(args.params, dict) else {}
    filter_levels = list(map(int, params.pop("filter_level", "").split(",")))
    num_randoms = int(params.pop("normalization", 1000))

    if len(params):
        raise TypeError(f"Unrecognized parameters: {params}")

    (
        pd.concat(
            list(map(lambda _: run(input, _, num_randoms, threads), filter_levels)),
            keys=filter_levels,
            names=["threshold"]
        )
        .reset_index()
        .drop(columns="level_1")
        .to_csv(output)
    )


if __name__ == "__main__":
    main()

