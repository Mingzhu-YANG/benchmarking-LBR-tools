"""Shared helpers: split-based edge extraction that is independent of node numbering
and of rooting.  All trees are treated as UNROOTED; edges are identified by the
bipartition (split) they induce on a given taxon set S.

Pruning is done implicitly: for a taxon subset S, every edge of the full reference tree
induces a bipartition of S; edges inducing the same bipartition (because the taxa between
them were pruned away, or because the root splits one unrooted edge in two) are SUMMED.
This is exactly what pruning + suppressing degree-2 nodes + unrooting does.
"""
import re
import dendropy


def popcount(x):
    return bin(x).count("1")


def read_newick(text, strip_re=None):
    t = dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True)
    if strip_re:
        for lf in t.leaf_node_iter():
            if lf.taxon is not None:
                lf.taxon.label = re.sub(strip_re, "", lf.taxon.label)
    return t


def raw_edges(tree, index):
    """Return ([(taxon_mask_below_edge, length), ...] for every non-root node, root_mask).
    `index` maps taxon label -> bit position. Tips absent from `index` get mask 0."""
    masks = {}
    out = []
    seed = tree.seed_node
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            lab = nd.taxon.label if nd.taxon is not None else None
            b = index.get(lab)
            m = (1 << b) if b is not None else 0
        else:
            m = 0
            for c in nd.child_node_iter():
                m |= masks[id(c)]
        masks[id(nd)] = m
        if nd is not seed:
            out.append((m, float(nd.edge_length) if nd.edge_length is not None else 0.0))
    return out, masks[id(seed)]


def restrict(edges, S, labels):
    """Project edges onto taxon set S (bitmask).
    Returns (terminal {label: length}, internal {split_mask: length}).
    Internal split key = side NOT containing the lowest-index taxon of S (unrooted-safe)."""
    n = popcount(S)
    anchor = S & -S
    term, inte = {}, {}
    for m, L in edges:
        sub = m & S
        k = popcount(sub)
        if k == 0 or k == n:
            continue
        comp = S ^ sub
        single = sub if k == 1 else (comp if n - k == 1 else 0)
        if single:
            lab = labels[single.bit_length() - 1]
            term[lab] = term.get(lab, 0.0) + L
        else:
            key = comp if (sub & anchor) else sub
            inte[key] = inte.get(key, 0.0) + L
    return term, inte


def mask_to_labels(mask, labels):
    return [labels[i] for i in range(mask.bit_length()) if (mask >> i) & 1]
