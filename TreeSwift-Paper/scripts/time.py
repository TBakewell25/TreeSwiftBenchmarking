#!/usr/bin/env python3
import collections; from collections.abc import MutableMapping; collections.MutableMapping = MutableMapping
from Bio import Phylo
from itertools import combinations
from treeswift import read_tree_newick
import dendropy
import ete3
import numpy
from CompactTree import compact_tree
NA = "NA" # what to print when function is not implemented

# get memory usage
def memory():
    from os import getpid; from psutil import Process
    return Process(getpid()).memory_info().rss

# distance matrix in ETE Toolkit
# obtained from https://github.com/linsalrob/EdwardsLab/blob/master/trees/tree_to_cophenetic_matrix.py
def distance_matrix_ete(tree):
    leaves = tree.get_leaves()
    paths = {x:set() for x in leaves}
    for n in leaves:
        if n.is_root():
            continue
        movingnode = n
        while not movingnode.is_root():
            paths[n].add(movingnode); movingnode = movingnode.up
    leaf_distances = {x.name:{} for x in leaves}
    for (leaf1, leaf2) in combinations(leaves, 2):
        uniquenodes = paths[leaf1] ^ paths[leaf2]
        distance = sum(x.dist for x in uniquenodes)
        leaf_distances[leaf1.name][leaf2.name] = distance
        leaf_distances[leaf2.name][leaf1.name] = distance
    return leaf_distances

# main code
from io import StringIO
from sys import argv
from time import time
if len(argv) != 4 or argv[2] not in {'treeswift','dendropy','biophylo','ete3', 'compacttree'}:
    print("USAGE: %s <tree_file> <treeswift_or_dendropy_or_biophylo_or_ete3> <task>"%argv[0]); exit(1)
if argv[1].lower().endswith('.gz'):
    from gzip import open as gopen
    treestr = gopen(argv[1]).read().decode().strip()
else:
    treestr = open(argv[1]).read().strip()
treeio = StringIO(treestr) # for Bio.Phylo

# postorder traversal
def postorder(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        for node in tree.postorder_node_iter():
            pass
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        for node in tree.find_clades(order='postorder'):
            pass
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_postorder():
            pass
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        for node in tree.traverse(strategy='postorder'):
            pass
        t_end = time()

    elif m == 'compacttree':
        tree = compact_tree.compact_tree(argv[1][:-3])
        total = 0.
        t_start = time()
        for node in compact_tree.traverse_postorder(tree):
            total += tree.get_edge_length(node)
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# preorder traversal
def preorder(m):
    if m == 'dendropy':
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_start = time()
        for node in tree.preorder_node_iter():
            pass
        t_end = time()
    elif m == 'biophylo':
        tree = Phylo.read(treeio, 'newick')
        t_start = time()
        for node in tree.find_clades(order='preorder'):
            pass
        t_end = time()
    elif m == 'treeswift':
        tree = read_tree_newick(treestr)
        t_start = time()
        for node in tree.traverse_preorder():
            pass
        t_end = time()
    elif m == 'ete3':
        tree = ete3.Tree(treestr,format=1)
        t_start = time()
        for node in tree.traverse(strategy='preorder'):
            pass
        t_end = time()
    elif m == 'compacttree':
        tree = compact_tree.compact_tree(argv[1][:-3])
        total = 0.
        t_start = time()
        for node in compact_tree.traverse_preorder(tree):
            total += tree.get_edge_length(node)
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# load tree
def load_tree(m):
    if m == 'dendropy':
        t_start = time()
        tree = dendropy.Tree.get(data=treestr, schema='newick')
        t_end = time()
    elif m == 'biophylo':
        t_start = time()
        tree = Phylo.read(treeio, 'newick')
        t_end = time()
    elif m == 'treeswift':
        t_start = time()
        tree = read_tree_newick(treestr)
        t_end = time()
    elif m == 'ete3':
        t_start = time()
        tree = ete3.Tree(treestr,format=1)
        t_end = time()
    elif m == 'compacttree':
        t_start = time()
        tree = compact_tree.compact_tree(argv[1][:-3])
        t_end = time()
    else:
        assert False, "Invalid tool: %s"%m
    return t_end-t_start

# memory
def measure_memory(m):
    if m == 'dendropy':
        m_start = memory()
        t = dendropy.Tree.get(data=treestr, schema='newick')
        t.encode_bipartitions()
        m_end = memory()
    elif m == 'biophylo':
        m_start = memory()
        t = Phylo.read(treeio, 'newick')
        m_end = memory()
    elif m == 'treeswift':
        m_start = memory()
        t = read_tree_newick(treestr)
        m_end = memory()
    elif m == 'ete3':
        m_start = memory()
        t = ete3.Tree(treestr,format=1)
        m_end = memory()
    elif m == 'compacttree':
        m_start = memory()
        t = compact_tree.compact_tree(argv[1][:-3])
        m_end = memory()
    else:
        assert False, "Invalid tool: %s"%m
    return m_end-m_start

TASKS = {
    'load_tree':load_tree,
    'memory':measure_memory,
    'postorder':postorder,
    'preorder':preorder,
}

# run
if argv[3] not in TASKS:
    print("Invalid task. Valid options: %s"%', '.join(sorted(TASKS.keys()))); exit(1)
print(TASKS[argv[3]](argv[2]))
