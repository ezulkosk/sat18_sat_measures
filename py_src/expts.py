import itertools
import numpy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import statsmodels.formula.api as sm
import sys
from tabulate import tabulate

from old_expts.data_analysis import type_map, all_subsets
from phdtk.db import create_sat_db_csv, init_cp2017_setup
import phdtk.latex_gen as latex_gen

comps = ["sc14-app", "sc13-app", "sc11-app", "sc09-app"]
solvers = ["maplecomsps", "glucose", "minisat"]

FSTR = "{0:.2f}"
FINE_FSTR = "{0:.5f}"
BIG_FSTR = "{0:.0f}"

def mmul(x, c):
    res = 1
    for i in c:
        res *= x[i]
    return res


def add_higher_order_features(df, subset, sz, heterogenous=False):
    """

    :param df:
    :param sz:
    :param heterogenous: if True, don't allow higher order features that use the same variable twice (e.g. V*V).

    :return:
    """
    fnames = []
    combs = []
    if sz == -1:
        upper = len(subset)
    else:
        upper = sz
    for i in range(2, upper + 1):
        for c in itertools.combinations_with_replacement(subset, i):
            flag = False
            if heterogenous:
                if len(set(c)) != len(c):
                    continue
            combs.append(c)

    for c in combs:
        fname = "*".join(list(c))
        fnames.append(fname)
        df[fname] = df.apply(lambda x : mmul(x, c), axis=1)
    # df['log_' + t] = df[t].apply(lambda x: 0.01 if x == 0 else x).apply(numpy.log)
    #x.append(reduce(operator.mul, c, 1))
    #print(fnames)
    return fnames


def data_summary(df, benchmarks, out_file, columns=None, names=None, caption_prefix=None):
    """
    Summarizes the average values of base features across each benchmark
    """
    if not caption_prefix:
        caption_prefix = ""
    df = df[df.benchmark.isin(benchmarks)]
    df['SAT'] = df['result']
    df.loc[df.SAT != "SAT", ['SAT']] = np.nan
    g = df.groupby("benchmark")
    res = g.aggregate('count')
    index = res.index
    res = res.append(res.sum(numeric_only=True), ignore_index=True)
    res.index = list(index) + ["Total"]
    # use simp_num_vars for # time instances, since maplecomsps contains some instances that get simplified away
    if not names:
        names = ["simp_num_vars", "simp_lsr_size", "simp_weak_size", "simp_q", "simp_backbones", "simp_tw_upper"]
    if not columns:
        columns = ["Instances", "LSR", "Weak", "Cmty", "Bones", "TW"]

    out = res[names]
    out.columns = columns

    with open(out_file, 'w') as o:
        latex_gen.insert_table(o, out.to_latex(), tabular=True, precomputed=True, tiny=False,
                               caption=caption_prefix + " The number of instances for which" +
                               " we were able to successfully compute each parameter. " +
                               "``Cmty'' refers to the community parameters; " +
                               "``Alpha+Dim'' refers to popularity and fractal dimension parameters;" +
                               "``TW'' denotes the treewidth upper bound; " +
                               "``Bones'' denotes backbone size; ``Merge+Res'' denote mergeability and resolvability.")
    # print(out.to_latex())


def average_metric_values(df, benchmarks, out_file, columns=None, names=None, caption_prefix=None, fstr=None):
    """
    Do metrics look better for app as opposed to random/crafted?
    """
    if not fstr:
        fstr = FSTR
    if not caption_prefix:
        caption_prefix = ""
    if "family" in benchmarks:
        df = df[df.benchmark.isin(["Application"])]
    else:
        df = df[df.benchmark.isin(benchmarks)]
        g = df.groupby("benchmark")

    if not columns:
        columns = ['benchmark', 'simp_lvr', 'simp_wvr', 'simp_q', 'simp_backbonesvr', 'simp_tw_uppervr']
    if not names:
        names = ['LSR/V', 'Weak/V', 'Q', 'Bones/V', 'TW/V']
    # print(columns)
    df = df[columns]

    if "family" in benchmarks:
        g = df.groupby("family")
    else:
        g = df.groupby("benchmark")


    res = g.aggregate('mean')
    res2 = g.aggregate('std')
    res3 = res.combine(res2, lambda x, y: [fstr.format(i) + " (" + fstr.format(j) + ")" for i, j in zip(x, y)])
    # print(res3)
    res3.columns = names

    with open(out_file, 'w') as o:
        latex_gen.insert_table(o, res3.to_latex(), tabular=True, precomputed=True, tiny=False,
                               caption=caption_prefix + " Mean (std. dev.) of several parameter values. ",
                               label="tab-meanstd")


def lsr_all_decs_comparison(df, benchmarks, out_file, caption_prefix=None):
    if not caption_prefix:
        caption_prefix = ""
    df = df[df.benchmark.isin(benchmarks)]
    df['simp_lsr_all_decs_overlap_ratio'] = df['simp_lsr_all_decs_intervr'] / df['simp_lsr_all_decs_unionvr']
    df = df[['benchmark', 'simp_lvr', 'simp_all_decsvr', 'simp_lsr_all_decs_overlap_ratio']]

    g = df.groupby("benchmark")
    res = g.aggregate('mean')
    res2 = g.aggregate('std')

    '''
    for col in res:
        res[col] = [np.nan if (not isinstance(val, str) and np.isnan(val)) else
                   (val if isinstance(val, str) else str(float(val)))
                   for val in res[col].tolist()]

    for col in res2:
        res2[col] = [np.nan if (not isinstance(val, str) and np.isnan(val)) else
                   (val if isinstance(val, str) else str(float(val)))
                   for val in res2[col].tolist()]
    '''
    # print(res)
    # print(res2)
    # print(res.combine(res2, lambda x, y: str(x) + " (" + str(y) + ")"))
    res3 = res.combine(res2, lambda x, y: [FSTR.format(i) + " (" + FSTR.format(j) + ")" for i, j in zip(x, y)])
    res3.columns = ["Laser", "All Decisions", "Overlap Ratio"]

    # print(res3)
    # print(res3.to_latex())

    with open(out_file, 'w') as o:
        latex_gen.insert_table(o, res3.to_latex(), precomputed=True, tiny=False, tabular=True,
                               caption=caption_prefix + " Mean (std. dev.) of Laser produced backdoor sizes " +
                               "versus all decision variables. " +
                               "Overlap Ratio is the size of the set " +
                               "$(Laser \\cap All Decisions) / (Laser \\cup All Decisions)$.",
                               label="lsr_vs_all_decs_table")


def structure_logging_summary(df, benchmarks, out_file, full=False):
    """
    do metrics look better for app as opposed to random/crafted?
    """
    print("Structure logging")
    out_str = ""

    df = df[df.benchmark.isin(benchmarks)]

    # 'struct_gini_normalized_picks', 'struct_ar_gini_normalized_picks', 'struct_nr_gini_normalized_picks',
    # 'struct_gini_normalized_clauses', 'struct_ar_gini_normalized_clauses', 'struct_nr_gini_normalized_clauses',
    df = df[['benchmark', 'name', 'struct_lsr', 'struct_ar_lsr', 'struct_nr_lsr',
             'simp_maplesat_time', 'simp_maplesat_ar_time', 'simp_maplesat_nr_time',
             'simp_maplesat_conflicts', 'simp_maplesat_ar_conflicts', 'simp_maplesat_nr_conflicts',
             'struct_avg_clause_lsr', 'struct_ar_avg_clause_lsr', 'struct_nr_avg_clause_lsr']]
    df = df.dropna()
    if full:
        # ['benchmark',
        # 'struct_gini_normalized_picks',
        # 'struct_ar_gini_normalized_picks',
        # 'struct_nr_gini_normalized_picks'],
        # ['benchmark',
        # 'struct_gini_normalized_clauses',
        # 'struct_ar_gini_normalized_clauses',
        # 'struct_nr_gini_normalized_clauses'],
        feature_lists = [
            ['benchmark', 'struct_lsr', 'struct_ar_lsr', 'struct_nr_lsr'],
            ['benchmark', 'struct_avg_clause_lsr', 'struct_ar_avg_clause_lsr', 'struct_nr_avg_clause_lsr'],
            ['benchmark', 'simp_maplesat_conflicts', 'simp_maplesat_ar_conflicts', 'simp_maplesat_nr_conflicts'],
            ['benchmark', 'simp_maplesat_time', 'simp_maplesat_ar_time', 'simp_maplesat_nr_time']
        ]

        # 'P1: Community-based Spatial Locality of Decisions',
        # 'P2: Community-based Spatial Locality of Learnt Clauses',
        expt_name_list = [
            'LSR Size',
            'Avg. Clause LSR',
            'Num Conflicts',
            'Solving Time (s)']
        best = ["min", "min", "min", "min"]
    else:
        feature_lists = [
            ['benchmark', 'struct_lsr', 'struct_ar_lsr', 'struct_nr_lsr'],
            ['benchmark', 'struct_avg_clause_lsr', 'struct_ar_avg_clause_lsr', 'struct_nr_avg_clause_lsr'],
            ['benchmark', 'simp_maplesat_conflicts', 'simp_maplesat_ar_conflicts', 'simp_maplesat_nr_conflicts'],
            ['benchmark', 'simp_maplesat_time', 'simp_maplesat_ar_time', 'simp_maplesat_nr_time']
        ]

        expt_name_list = [
            'LSR Size',
            'Avg. Clause LSR',
            'Num Conflicts',
            'Solving Time (s)']
        best = ["min", "min", "min", "min"]

    end_row = " \\\\ \\hline"

    out_str += "\\begin{center}\n"
    out_str += "\\begin{tabular}{ |l|c|c|c| }\n"

    # header
    out_str += "\\hline\n"
    out_str += " & ".join(
        ["\\textbf{" + i + "}" for i in ["Property", "Luby", "Always Restart", "Never Restart"]]) + end_row + "\n"

    for l, e, b in zip(feature_lists, expt_name_list, best):
        df2 = df[l]
        g = df2.groupby("benchmark")
        res = g.aggregate('mean')
        res2 = g.aggregate('std')
        res3 = res.combine(res2, lambda x, y: [FSTR.format(i) + " (" + FSTR.format(j) + ")"
                                               if i <= 1000
                                               else BIG_FSTR.format(i) + " (" + BIG_FSTR.format(j) + ")"
                                               for i, j in zip(x, y)])
        out_str += e + "& "
        for index, row in res3.iterrows():
            pre = ""
            post = "\\\\"
            nums = [float(row[fname].split()[0]) for fname in l[1:]]
            nums_and_std = [row[fname] for fname in l[1:]]
            high = -1
            low = 9999999
            high_index = -1
            low_index = -1

            for i in range(len(nums)):
                if str(nums[i]) == "nan":
                    continue
                else:
                    if nums[i] > high:
                        high = nums[i]
                        high_index = i
                    if nums[i] < low:
                        low = nums[i]
                        low_index = i

            for i in range(len(nums)):
                if str(nums[i]) == "nan":
                    continue

                if b == "min" and low_index == i:
                    nums_and_std[low_index] = "\\textbf{" + nums_and_std[low_index] + "}"
                elif b == "max" and high_index == i:
                    nums_and_std[high_index] = "\\textbf{" + nums_and_std[high_index] + "}"

            out_str += pre + " & ".join(str(i) for i in nums_and_std) + post + "\n"
    out_str += "\\hline\n"
    out_str += "\\end{tabular}\n"
    out_str += "\\end{center}\n"
    with open(out_file, 'w') as o:
        latex_gen.insert_table(o, out_str, tabular=True, precomputed=True, tiny=False, label="tab_lens",
                               caption="Comparison of LSR measures and solving time for various restart policies" +
                                       " on the Agile benchmark. LSR sizes are normalized by the number of variables.")

def pocr_structure_logging_summary(df, benchmark_names, identifiers, headers, caption, label, branching_heuristics, out_file, sat=None, fstr=None):

    # TODO how to filter instances (maybe take all that vsids + lrb finish, and then note in the paper the difference
    # for random

    # print("Structure logging")
    out_str = ""

    if not fstr:
        fstr = FSTR

    HOW = "any"

    b_map = {"LRB": "lrb", "VSIDS": "vsids", "Random": "random"}
    r_map = {"Luby": "luby", "Always": "ar", "Never": "nr"}

    main_cols = ['benchmark', 'name', 'bucket', 'simp_num_vars', 'simp_backbones', 'simp_backbonesvr', 'result',
             'pocr_lrb_ar_all_decs',
             'pocr_lrb_ar_bb_flips',
             'pocr_lrb_ar_bb_subsumed',
             'pocr_lrb_ar_gini_clauses',
             'pocr_lrb_ar_gini_picks',
             'pocr_lrb_ar_lsr',
             'pocr_lrb_luby_all_decs',
             'pocr_lrb_luby_bb_flips',
             'pocr_lrb_luby_bb_subsumed',
             'pocr_lrb_luby_gini_clauses',
             'pocr_lrb_luby_gini_picks',
             'pocr_lrb_luby_lsr',
             'pocr_lrb_nr_all_decs',
             'pocr_lrb_nr_bb_flips',
             'pocr_lrb_nr_bb_subsumed',
             'pocr_lrb_nr_gini_clauses',
             'pocr_lrb_nr_gini_picks',
             'pocr_lrb_nr_lsr',
             'pocr_random_ar_all_decs',
             'pocr_random_ar_bb_flips',
             'pocr_random_ar_bb_subsumed',
             'pocr_random_ar_gini_clauses',
             'pocr_random_ar_gini_picks',
             'pocr_random_ar_lsr',
             'pocr_random_luby_all_decs',
             'pocr_random_luby_bb_flips',
             'pocr_random_luby_bb_subsumed',
             'pocr_random_luby_gini_clauses',
             'pocr_random_luby_gini_picks',
             'pocr_random_luby_lsr',
             'pocr_random_nr_all_decs',
             'pocr_random_nr_bb_flips',
             'pocr_random_nr_bb_subsumed',
             'pocr_random_nr_gini_clauses',
             'pocr_random_nr_gini_picks',
             'pocr_random_nr_lsr',
             'pocr_vsids_ar_all_decs',
             'pocr_vsids_ar_bb_flips',
             'pocr_vsids_ar_bb_subsumed',
             'pocr_vsids_ar_gini_clauses',
             'pocr_vsids_ar_gini_picks',
             'pocr_vsids_ar_lsr',
             'pocr_vsids_luby_all_decs',
             'pocr_vsids_luby_bb_flips',
             'pocr_vsids_luby_bb_subsumed',
             'pocr_vsids_luby_gini_clauses',
             'pocr_vsids_luby_gini_picks',
             'pocr_vsids_luby_lsr',
             'pocr_vsids_nr_all_decs',
             'pocr_vsids_nr_bb_flips',
             'pocr_vsids_nr_bb_subsumed',
             'pocr_vsids_nr_gini_clauses',
             'pocr_vsids_nr_gini_picks',
             'pocr_vsids_nr_lsr',
             'pocr_lrb_ar_time',
             'pocr_lrb_luby_time',
             'pocr_lrb_nr_time',
             'pocr_random_ar_time',
             'pocr_random_luby_time',
             'pocr_random_nr_time',
             'pocr_vsids_ar_time',
             'pocr_vsids_luby_time',
             'pocr_vsids_nr_time',
             'pocr_lrb_ar_bb_conflicts',
             'pocr_lrb_luby_bb_conflicts',
             'pocr_lrb_nr_bb_conflicts',
             'pocr_random_ar_bb_conflicts',
             'pocr_random_luby_bb_conflicts',
             'pocr_random_nr_bb_conflicts',
             'pocr_vsids_ar_bb_conflicts',
             'pocr_vsids_luby_bb_conflicts',
             'pocr_vsids_nr_bb_conflicts',
             'pocr_lrb_ar_bb_propagations',
             'pocr_lrb_luby_bb_propagations',
             'pocr_lrb_nr_bb_propagations',
             'pocr_random_ar_bb_propagations',
             'pocr_random_luby_bb_propagations',
             'pocr_random_nr_bb_propagations',
             'pocr_vsids_ar_bb_propagations',
             'pocr_vsids_luby_bb_propagations',
             'pocr_vsids_nr_bb_propagations',
             'pocr_lrb_ar_bb_subsumed_raw',
             'pocr_lrb_luby_bb_subsumed_raw',
             'pocr_lrb_nr_bb_subsumed_raw',
             'pocr_random_ar_bb_subsumed_raw',
             'pocr_random_luby_bb_subsumed_raw',
             'pocr_random_nr_bb_subsumed_raw',
             'pocr_vsids_ar_bb_subsumed_raw',
             'pocr_vsids_luby_bb_subsumed_raw',
             'pocr_vsids_nr_bb_subsumed_raw',
             'pocr_lrb_ar_decisions',
             'pocr_lrb_luby_decisions',
             'pocr_lrb_nr_decisions',
             'pocr_random_ar_decisions',
             'pocr_random_luby_decisions',
             'pocr_random_nr_decisions',
             'pocr_vsids_ar_decisions',
             'pocr_vsids_luby_decisions',
             'pocr_vsids_nr_decisions',
             'pocr_lrb_ar_conflicts',
             'pocr_lrb_luby_conflicts',
             'pocr_lrb_nr_conflicts',
             'pocr_random_ar_conflicts',
             'pocr_random_luby_conflicts',
             'pocr_random_nr_conflicts',
             'pocr_vsids_ar_conflicts',
             'pocr_vsids_luby_conflicts',
             'pocr_vsids_nr_conflicts',
             'pocr_lrb_ar_lsr_cmty_spread',
             'pocr_lrb_luby_lsr_cmty_spread',
             'pocr_lrb_nr_lsr_cmty_spread',
             'pocr_random_ar_lsr_cmty_spread',
             'pocr_random_luby_lsr_cmty_spread',
             'pocr_random_nr_lsr_cmty_spread',
             'pocr_vsids_ar_lsr_cmty_spread',
             'pocr_vsids_luby_lsr_cmty_spread',
             'pocr_vsids_nr_lsr_cmty_spread',
             'pocr_lrb_ar_lsr_cmty_largest_ratio',
             'pocr_lrb_luby_lsr_cmty_largest_ratio',
             'pocr_lrb_nr_lsr_cmty_largest_ratio',
             'pocr_random_ar_lsr_cmty_largest_ratio',
             'pocr_random_luby_lsr_cmty_largest_ratio',
             'pocr_random_nr_lsr_cmty_largest_ratio',
             'pocr_vsids_ar_lsr_cmty_largest_ratio',
             'pocr_vsids_luby_lsr_cmty_largest_ratio',
             'pocr_vsids_nr_lsr_cmty_largest_ratio',
             'proof_lrb_ar_avg_deps',
             'proof_lrb_ar_avg_merges',
             'proof_lrb_ar_avg_merges_normalized_by_deps',
             'proof_lrb_ar_gini_cmty_clauses',
             'proof_lrb_ar_gini_cmty_clauses_normalized',
             'proof_lrb_ar_gini_cmty_occs',
             'proof_lrb_ar_gini_cmty_occs_normalized',
             'proof_lrb_ar_num_merges',
             'proof_lrb_ar_pf_lit_mse',
             'proof_lrb_ar_pf_var_mse',
             'proof_lrb_ar_q_original_partition',
             'proof_lrb_luby_avg_deps',
             'proof_lrb_luby_avg_merges',
             'proof_lrb_luby_avg_merges_normalized_by_deps',
             'proof_lrb_luby_gini_cmty_clauses',
             'proof_lrb_luby_gini_cmty_clauses_normalized',
             'proof_lrb_luby_gini_cmty_occs',
             'proof_lrb_luby_gini_cmty_occs_normalized',
             'proof_lrb_luby_num_merges',
             'proof_lrb_luby_pf_lit_mse',
             'proof_lrb_luby_pf_var_mse',
             'proof_lrb_luby_pf_lit_rmse',
             'proof_lrb_luby_pf_var_rmse',
             'proof_lrb_ar_pf_lit_rmse',
             'proof_lrb_ar_pf_var_rmse',
             'proof_lrb_nr_pf_lit_rmse',
             'proof_lrb_nr_pf_var_rmse',
             'proof_lrb_luby_q_original_partition',
             'proof_lrb_nr_avg_deps',
             'proof_lrb_nr_avg_merges',
             'proof_lrb_nr_avg_merges_normalized_by_deps',
             'proof_lrb_nr_gini_cmty_clauses',
             'proof_lrb_nr_gini_cmty_clauses_normalized',
             'proof_lrb_nr_gini_cmty_occs',
             'proof_lrb_nr_gini_cmty_occs_normalized',
             'proof_lrb_nr_num_merges',
             'proof_lrb_nr_pf_lit_mse',
             'proof_lrb_nr_pf_var_mse',
             'proof_lrb_nr_q_original_partition'
                 ]

    if not branching_heuristics:
        branching_heuristics = ["LRB", "VSIDS", "Random"]

    df = df[main_cols]
    for c in main_cols:
        if ('lsr' in c or 'all_decs' in c) and 'lsr_cmty' not in c:
            df[c + "_vr"] = df[c] / df['simp_num_vars']

    bases = ['pocr_lrb_ar',
             'pocr_lrb_luby',
             'pocr_lrb_nr',
             'pocr_random_ar',
             'pocr_random_luby',
             'pocr_random_nr',
             'pocr_vsids_ar',
             'pocr_vsids_luby',
             'pocr_vsids_nr']

    #for i in bases:
    #    df[i + "_glr"] = df[i + "_conflicts"] / df[i + "_decisions"]

    # Generic Tables
    all_dfs = []
    for b in benchmark_names:
        d = df[df.benchmark.isin([b])]
        if sat:
            d = d.drop(d[d.result != sat].index)
        all_dfs.append(d)
        print(len(d))

    avgs = []
    stds = []
    #print(identifiers)
    all_cols = []
    for i in identifiers:
        l = [col for col in df.columns if i in col and "_bb_" not in col] # TODO NOTE removing backbone stuff for now
        all_cols += l
    temp_cols = []
    for col in all_cols:
        found = False
        for b in branching_heuristics:
            k = b_map[b]
            if k in col:
                found = True
                break
        if found:
            temp_cols.append(col)
    all_cols = temp_cols
    #all_cols = ["name", "bucket"] + all_cols
    all_cols = ["simp_num_vars"] + all_cols
    # print(all_cols)
    print("ALLDFS:", benchmark_names)
    for d in all_dfs:
        d = d[all_cols]
        # sys.exit()
        #d = d.drop(d[d.bucket != "hardware"].index) # TODO warning
        d = d.dropna(how="all") # this is only for the print, uses any below
        #print(d)
        d = d.dropna(how=HOW)
        #print(d)
        # print(d.mean())
        d = d.drop_duplicates()
        d = d.reset_index(drop=True)
        # print(len(d))
        a = np.mean(d)
        std = np.std(d)
        avgs.append(a)
        stds.append(std)
        #print(a)
        print("Num Instances ", sat, len(d))
        if len(d) < 2:
            print("Not enough instances, returning")
            return
    # produce table
    end_row = " \\\\ \\hline"

    out_str = ""
    out_str += "\\begin{table}[t]\n"
    out_str += "\\begin{center}\n"
    out_str += "\\begin{tabular}{ |l|l||c|c||c|c| }\n"



    out_str += "\\hline\n"
    out_str += "\\multicolumn{2}{ |c|| }{\\textbf{Heuristic}} " + \
        "".join(["& \\multicolumn{" + str(len(identifiers)) + "}{ c"
                 + (" || " if j != len(benchmark_names) - 1 else " | ") + "}{\\textbf{" + i + "}} "
                 for i,j in zip(benchmark_names, range(len(benchmark_names)))]) + end_row + "\n"
    out_str += " & ".join(
        ["\\textbf{" + i + "}" for i in
         ["Branching", "Restart"] + (headers * len(benchmark_names))]) + end_row + "\n"
    for b in branching_heuristics:
        for r in ["Luby", "Always" , "Never"]:
            row = ""
            if r == "Luby":
                row += "\multirow{3}{4em}{" + b + "} &"
            else:
                row += "& "
            row += r + " & "
            key = "pocr_" + b_map[b] + "_" + r_map[r] + "_"
            # if key not in main_cols and not key.endswith("_vr"):
            #    print(key)
            #    key = "proof_" + b_map[b] + "_" + r_map[r] + "_"

            #for a in avgs:
            #    print(a)
            try:
                row += " & ".join(
                    fstr.format(i) + " (" + fstr.format(j)  + ")" for i,j in [(a[key + ident], std[key + ident]) for a, std, ident in zip(avgs * len(identifiers) * len(benchmark_names), stds * len(identifiers) * len(benchmark_names), identifiers * len(benchmark_names))])
            except:
                key = "proof_" + b_map[b] + "_" + r_map[r] + "_"
                row += " & ".join(
                    fstr.format(i) + " (" + fstr.format(j) + ")" for i, j in
                    [(a[key + ident], std[key + ident]) for a, std, ident in
                     zip(avgs * len(identifiers) * len(benchmark_names), stds * len(identifiers) * len(benchmark_names),
                         identifiers * len(benchmark_names))])

            out_str += row
            if r == "Never":
                out_str += end_row + "\n"
            else:
                out_str += "\\\\\n"
    out_str += "\\end{tabular}\n"
    out_str += "\\end{center}\n"
    out_str += "\\caption{" + (sat + " ONLY. " if sat else "") + caption + "}\n"
    out_str += "\\label{" + label + "}\n"
    out_str += "\\end{table}\n\n"

    # print(out_str)
    out_file.write(out_str)
    return


def create_extra_ratios(df):
    # df['qcor'] = df['q'] / df['num_cmtys']
    df['simp_qcor'] = df['simp_q'] / df['simp_num_cmtys']
    # df['cvr'] = df['num_clauses'] / df['num_vars']
    df['simp_cvr'] = df['simp_num_clauses'] / df['simp_num_vars']
    # df['tw_upper_vr'] = df['tw_upper'] / df['num_vars']
    df['simp_tw_uppervr'] = df['simp_tw_upper'] / df['simp_num_vars']
    df['simp_backbonesvr'] = df['simp_backbones'] / df['simp_num_vars']
    df['simp_wvr'] = df['simp_weak_size'] / df['simp_num_vars']
    df['simp_lvr'] = df['simp_lsr_size'] / df['simp_num_vars']
    df['simp_unionwvr'] = df['simp_num_vars_in_any_weak'] / df['simp_num_vars']
    df['simp_cmtysvr'] = df['simp_num_cmtys'] / df['simp_num_vars']

    df['simp_all_decsvr'] = df['simp_all_decs'] / df['simp_num_vars']
    df['simp_lsr_all_decs_unionvr'] = df['simp_lsr_all_decs_union'] / df['simp_num_vars']
    df['simp_lsr_all_decs_intervr'] = df['simp_lsr_all_decs_inter'] / df['simp_num_vars']


def regression_helper(df, benchmarks=None, times=None, subsets=None, subset_size_filter=None, rotate=False,
                      filter_under_second=False,
                      scale_features=True, log_time=True, heterogeneous=True, highest_order_features=-1,
                      grab_all=False, ridge=True):
    """
    Performs Ridge regression on subsets of features vs time.

    :param df: the DataFrame
    :param benchmarks: list of string IDs of the considered benchmarks
    :param times: list of solvers to be considered
    :param subsets: if given, only perform regression on these subsets of features
    :param subset_size_filter: if given, tries all subsets of features of this size
    :param rotate: used with subset_size_filter and multiple benchmarks (only dumps the best feature sets of each)
    :param filter_under_second: if True, remove any instances that solved under 1 second
    :param scale_features: currently a mean zero, std. dev. of 1 (time is not scaled by this)
    :param log_time: take log of time before regression
    :param heterogeneous: (Should probably set to True). Whenever we create a higher order feature, don't allow e.g. Q*Q
    :param highest_order_features: only allow subsets of features of given maximal size
    :param grab_all: output results of all regressions, rather than just the best for each benchmark
    :param ridge: if true, use ridge regression, else use linear
    :return:
    """

    if not benchmarks:
        benchmarks = ["app", "random", "crafted", "agile"]
    if not times:
        times = ['time']
    if not subsets:
        data_types = ["simp_num_vars", "simp_num_clauses", "simp_cvr",  # basic
                      "simp_num_cmtys", "simp_q", "simp_qcor",  # cmty
                      # "simp_lvr",  # lsr
                      "simp_tw_uppervr",  # tw
                      "dimacs_merges", "dimacs_resolutions",  # merge
                      "simp_fractal_dim_vig", "simp_fractal_dim_cvig", "simp_alpha_var",  # fractal / pop
                      "inter_cmty_merge_over_res", "avg_intra_cmty_merge_over_res"
                      ]

        df = df[data_types + ['benchmark', 'time']]
        df = df.dropna()

        subsets = [list(i) for i in all_subsets(data_types, 1)]

    # filter out subsets of wrong size
    if subset_size_filter:
        subsets = [i for i in subsets if len(i) == subset_size_filter]
        # print("WARNING FILTERING FOR RES REMOVE LATER!!")
        # subsets = [i for i in subsets if 'dimacs_resolutions' in i or 'dimacs_merges' in i]
        print("SUBSETS: ", len(subsets))

    if log_time:
        for t in times:
            df['log_' + t] = df[t].apply(lambda x: 0.01 if x <= 0 else x)
            df['log_' + t] = df['log_' + t].apply(numpy.log)
    else:
        # TODO clean
        for t in times:
            df['log_' + t] = df[t]

    features = []
    r2_values = []
    num_instances = []
    subset_count = 0
    for s in subsets:
        subset_count += 1
        print("Subset count: ", subset_count)
        r2_inst = []
        num_insts_inst = []
        for c in benchmarks:
            for t in times:
                # print(r2_values)
                if benchmarks == ["all"]:
                    curr_df = df
                else:
                    curr_df = df.loc[df['benchmark'] == c]
                if filter_under_second:
                    curr_df = curr_df.loc[df[t] > 1]
                curr_df = curr_df[s + ['log_' + t]].dropna()
                if len(curr_df) < 4:
                    print("df too small", c, s)
                    r2_inst.append("N/A")
                    num_insts_inst.append("N/A")
                    continue

                # scale features, but not time
                if scale_features:
                    # use this for 0-1 scaling
                    # curr_df[curr_df.columns.difference(['log_' + t])] -=
                    #        curr_df[curr_df.columns.difference(['log_' + t])].min()
                    # curr_df[curr_df.columns.difference(['log_' + t])]
                    #        /= curr_df[curr_df.columns.difference(['log_' + t])].max()

                    # use this for mean 0 std 1 scaling
                    curr_df[curr_df.columns.difference(['log_' + t])] \
                        -= curr_df[curr_df.columns.difference(['log_' + t])].mean()
                    curr_df[curr_df.columns.difference(['log_' + t])] \
                        /= curr_df[curr_df.columns.difference(['log_' + t])].std()

                fnames = add_higher_order_features(curr_df, s, highest_order_features, heterogenous=heterogeneous)

                model = sm.ols(data=curr_df, formula="log_" + t + " ~ " + "+".join(s + fnames))

                if ridge:
                    # if L1_wt = 0, then it's ridge; if it's 1, then it's lasso
                    res = model.fit_regularized(L1_wt=0)
                else:
                    res = model.fit()
                # TODO this is the summary line
                print(res.summary())
                print(res.pvalues)

                r2 = res.rsquared_adj

                # print(r2)
                #print(curr_df)
                #sys.exit()
                instances = len(curr_df.index)
                if instances < 1:
                    r2 = "N/A"
                    instances = "N/A"
                r2_inst.append(r2)
                num_insts_inst.append(instances)

        na_test = set(r2_inst)
        if list(na_test) == ["N/A"]:
            print("in here")
            continue
        features.append("$" + "\oplus{}".join([type_map[t] for t in s]) + "$")
        r2_values.append(r2_inst)
        num_instances.append(num_insts_inst)

    print("counts", len(features), len(r2_values), len(num_instances))
    items = list(zip(features, r2_values, num_instances))
    print("ITEMS:", items)
    big_list = []
    if rotate:
        for i in range(len(benchmarks)):
            curr_items = sorted(items, key=lambda x: x[1][i] if x[1][i] != "N/A" else -1, reverse=True)
            big_list.append(curr_items)
    else:
        big_list.append(items)

    out_table_rows = []
    index = 0
    print("Big list:", big_list)
    for l in big_list:
        rows = []
        if not l:
            continue
        if not grab_all:
            f, r2s, insts = l[0]
            out_row = [f]
            count = 0
            for i, j in zip(r2s, insts):
                value = str(FSTR.format(i)) + " (" + str(j) + ")" if i != "N/A" else "N/A"
                if count == index:
                    value = "\\textbf{" + value + "}"
                out_row.append(value)
                count += 1
            out_table_rows.append(out_row)
            index += 1
        for f, r2s, insts in l:
            row = [f] + [str(FSTR.format(i)) + " (" + str(j) + ")" if i != "N/A" else "N/A" for i, j in zip(r2s, insts)]
            rows.append(row)
            if grab_all:
                out_table_rows.append(row)
    return out_table_rows


def regression(df, benchmarks, out_file, caption_prefix=None):
    """
    Tests if subsets of features correlate with solving time.
    """

    #print(df['dimacs_merges'])
    #sys.exit()


    heterogeneous_r2 = regression_helper(df, benchmarks=benchmarks, subsets=[
        ["simp_num_vars", "simp_num_clauses", "intra_ratio", "inter_cmty_merge_over_res", "avg_intra_cmty_merge_over_res"],
        ["simp_num_vars", "simp_num_clauses", "dimacs_merges", "dimacs_resolutions", "intra_ratio"],
        ["simp_num_vars", "simp_num_clauses", "dimacs_merges", "dimacs_resolutions", "avg_intra_cmty_merge_over_res"],
        ["simp_num_vars", "simp_num_clauses", "dimacs_merges", "dimacs_resolutions", "inter_cmty_merge_over_res"],
        ["simp_num_vars", "simp_num_clauses", "dimacs_merges", "dimacs_resolutions"],
        ["simp_num_vars", "simp_num_clauses", "dimacs_merges", "dimacs_resolutions", "mvr"],
        ["simp_num_vars", "simp_num_clauses", "dimacs_merges", "mvr"],
        ["simp_num_vars", "simp_num_clauses", "mvr"],
        ["simp_num_vars", "simp_num_clauses", "simp_cvr"],
        ["simp_num_vars", "simp_num_clauses", "simp_num_cmtys", "simp_q"],
        ["simp_num_vars", "simp_num_clauses", "simp_num_cmtys", "simp_q", "simp_qcor"],
        ["simp_num_vars", "simp_num_clauses", "simp_lsr_size", "simp_lvr"],
        ["simp_num_vars", "simp_num_clauses", "simp_num_min_weak", "simp_weak_size"],
        ["simp_num_vars", "simp_num_clauses", "simp_backbones", "simp_backbonesvr"],
        ["simp_num_vars", "simp_num_clauses", "simp_tw_upper", "simp_tw_uppervr"],
        ["simp_num_vars", "simp_num_clauses", "dimacs_merges", "dimacs_resolutions", "dimacs_merge_resolutions_cmty_avg"],
        ["simp_num_vars", "simp_num_clauses", "dimacs_merges", "dimacs_merge_resolutions_cmty_avg"],
        ["simp_num_vars", "simp_num_clauses", "dimacs_resolutions", "dimacs_merges"],
        ["simp_num_vars", "simp_num_clauses", "dimacs_resolutions", "dimacs_merge_resolutions_cmty_avg"],
        ["simp_num_vars", "simp_num_clauses", "simp_fractal_dim_vig", "simp_fractal_dim_cvig", "simp_alpha_var"]
    ],
                                         rotate=False, grab_all=True, ridge=False)

    """
    data_types = ["simp_num_vars", "simp_num_clauses", "simp_cvr",  # basic
                  "simp_num_cmtys", "simp_q", "simp_qcor",  # cmty
                  "simp_lsr_size", "simp_lvr",  # lsr
                  "simp_tw_upper", "simp_tw_uppervr",  # tw
                  "dimacs_merges", "dimacs_resolutions", "dimacs_merge_resolutions_cmty_avg",  # merge
                  "simp_fractal_dim_vig", "simp_fractal_dim_cvig", "simp_alpha_var" # fractal / pop
                  ]
    """
    data_types = ["simp_num_vars", "simp_num_clauses", "simp_cvr",  # basic
                  "simp_num_cmtys", "simp_q", "simp_qcor",  # cmty
                  # "simp_lvr",  # lsr
                  "simp_tw_uppervr",  # tw
                  "dimacs_merges", "dimacs_resolutions", # merge
                  "simp_fractal_dim_vig", "simp_fractal_dim_cvig", "simp_alpha_var",  # fractal / pop
                  "inter_cmty_merge_over_res", "avg_intra_cmty_merge_over_res"
                  ]

    df = df[data_types + ['benchmark', 'time']]
    df = df.dropna()


    # NOTE: ordered based on significance values
    #best_combined_r2 = regression_helper(df, benchmarks=benchmarks, subsets=[
    #    ["simp_q", "simp_cvr", "simp_lvr", "simp_qcor", "simp_num_clauses"],
    #    ["simp_tw_uppervr", "simp_q", "simp_num_cmtys", "simp_tw_upper", "simp_lvr"],
    #    ["simp_qcor", "simp_lvr", "simp_num_clauses", "simp_lsr_size", "simp_q"],
    #    ["simp_num_cmtys", "simp_tw_uppervr", "simp_cvr", "simp_tw_upper", "simp_q"]
    #],
    #rotate=False, grab_all=True, ridge=False)

    #best_combined_r2 = regression_helper(df, benchmarks=benchmarks, subsets=[
    #    ["simp_num_vars", "simp_cvr", "simp_lvr", "simp_qcor", "simp_num_clauses", "simp_q"],
    #    ["simp_tw_uppervr", "simp_q", "simp_num_clauses", "simp_qcor", "simp_lvr", "dimacs_merge_resolutions_cmty_avg"],
    #    ["dimacs_resolutions", "simp_alpha_var", "simp_lvr", "simp_num_cmtys", "simp_cvr", "simp_num_clauses"],
    #    ["simp_lsr_size", "dimacs_resolutions", "simp_num_clauses", "simp_tw_upper", "simp_q", "simp_lvr"]
    #],
    #rotate=False, grab_all=True, ridge=False)

    best_combined_r2 = regression_helper(df, benchmarks=benchmarks, subsets=[
       ["simp_cvr", "simp_q", "dimacs_merges", "dimacs_resolutions", "simp_fractal_dim_vig"],
       ["simp_num_clauses", "simp_cvr", "simp_q", "simp_tw_uppervr", "simp_fractal_dim_cvig"],
       ["simp_num_vars", "simp_cvr", "simp_num_cmtys", "simp_q", "simp_fractal_dim_vig"],
       ["simp_num_vars", "simp_cvr", "simp_q", "simp_tw_uppervr", "dimacs_merges"]
    ],
    rotate=False, grab_all=True, ridge=False)

    #best_combined_r2 = regression_helper(df, benchmarks=benchmarks, subset_size_filter=5, rotate=True, ridge=False)

    #best_combined_r2 = regression_helper(df, benchmarks=benchmarks, subset_size_filter=6, highest_order_features=3, rotate=True, ridge=False)
    # best_combined_r2 = regression_helper(df, benchmarks=benchmarks, subset_size_filter=6, highest_order_features=-1, rotate=True, ridge=False)
    #best_combined_r2 = regression_helper(df, benchmarks=benchmarks, subset_size_filter=5, highest_order_features=-1, rotate=True, ridge=False)

    rows = heterogeneous_r2 + [["\\hline"]] + best_combined_r2

    with open(out_file, 'w') as o:
        latex_gen.insert_table(o, rows, tiny=False, headers=["Feature Set"] + benchmarks,
                               caption=caption_prefix + " Adjusted R$^2$ values for the given features, "
                               + "compared to log of MapleCOMSPS' solving time. "
                               + "The number in parentheses indicates the number of instances "
                               + "that were considered in each case. The lower section considers "
                               + "heterogeneous sets of features across different parameter types.",
                               label="tab-regressions", tabular=True)



def collect_data(df, benchmarks, tables_dir, prefix_label, tex):
    if not prefix_label:
        caption_prefix = ""
    else:
        caption_prefix = prefix_label.replace("_", " ").upper() + ":"

    data_summary(df, benchmarks, tables_dir + prefix_label + "datasummary.tex", caption_prefix=caption_prefix)
    average_metric_values(df, benchmarks, tables_dir + prefix_label + "averagemetrics.tex",
                          caption_prefix=caption_prefix)
    # structure_logging_summary(d, ["Agile"], tables_dir + l + "structure_logging.tex",
    #                           full=True, caption_prefix=caption_prefix)
    lsr_all_decs_comparison(df,
                            benchmarks,
                            tables_dir + prefix_label + "lsr_all_decs_comparison.tex",
                            caption_prefix=caption_prefix)

    # print(df)
    print("REGRESSION")
    regression(df, benchmarks, tables_dir + prefix_label + "regression.tex", caption_prefix=caption_prefix)
    tex.write("\\begin{table}[t]\n")
    for i in ["datasummary.tex", "averagemetrics.tex", "lsr_all_decs_comparison.tex",
              "regression.tex"]:  # "regression.tex", "same_insts_regression.tex", "bridgebd.tex"]:
        filename = tables_dir + prefix_label + i
        tex.write("\\input{" + filename + "}")
        tex.write("\n")
    tex.write("\\end{table}\n")


def collect_refined_app_data(df, benchmarks, tables_dir, tex, split_sat=True, split_buckets=True):
    # to add a refinement, create a tuple (label, constraint)
    refinements = []
    if split_sat:
        sat_refinement = ("sat_", lambda x: x.result != 'SAT')
        unsat_refinement = ("unsat_", lambda x: x.result != 'UNSAT')
        refinements.append((sat_refinement, unsat_refinement))
    if split_buckets:
        hardware_refinement = ("hardware_", lambda x: x.bucket != 'hardware')
        software_refinement = ("software_", lambda x: x.bucket != 'software')
        crypto_refinement = ("crypto_", lambda x: x.bucket != 'crypto')
        refinements.append((hardware_refinement, software_refinement, crypto_refinement))

    all_refinements = itertools.product(*refinements)
    for refine_list in all_refinements:

        curr_df = df.copy(deep=True)
        full_label = ""
        for (label, refinement) in refine_list:
            full_label += label
            print(curr_df['bucket'])
            curr_df = curr_df.drop(curr_df[refinement].index)
        collect_data(curr_df, benchmarks, tables_dir, full_label, tex)


def main():
    base_dir, data_dir, case_studies, fields = init_cp2017_setup()
    benchmarks = ["Application", "Crafted", "Random", "Agile"]

    families = [
        '2d-strip-packing',
        'argumentation',
        'bio',
        'crypto-aes',
        'crypto-des',
        'crypto-gos',
        'crypto-md5',
        'crypto-sha',
        'crypto-vmpc',
        'diagnosis',
        'fpga-routing',
        'hardware-bmc',
        'hardware-bmc-ibm',
        'hardware-cec',
        'hardware-manolios',
        'hardware-manolius',
        'hardware-velev',
        'planning',
        'scheduling',
        'scheduling-pesp',
        'software-bit-verif',
        'software-bmc',
        'termination'
    ]


    pd.set_option('display.max_rows', 1500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

    tex_dir = "/home/ezulkosk/cp2017_benchmarks/analysis/"
    tables_dir = tex_dir + "tables/"
    tex_file = tex_dir + "/db.tex"
    tex = latex_gen.gen_tex_article(tex_file)
    # tex.write('\\tiny\n')

    if True:
        create_sat_db_csv(data_dir, fields)
    df = pd.read_csv(data_dir + 'db.csv')

    # TODO replacing sat lsr with sat min
    # print(df['simp_sat_min_lsr'])

    df['simp_lsr_size'] = df['simp_lsr_opt_size']
    #df['simp_lsr_size2'] = df['simp_sat_min_lsr']
    #df['simp_lsr_size2'][df['simp_sat_min_lsr'].isnull()] = df['simp_lsr_size']
    #df['simp_lsr_size'] = df['simp_lsr_size2']
    #todo make sure this is desired
    print("TODO check here!!!")

    df['mvr'] = df['dimacs_merges'] / df['dimacs_resolutions']

    df['dimacs_merges'] /= df['simp_num_clauses'] * df['simp_num_clauses']
    df['dimacs_resolutions'] /= df['simp_num_clauses'] * df['simp_num_clauses']


    df['simp_num_vars_recip'] = 1 / df['simp_num_vars']
    #df["pocr_lrb_luby_var_mse"] = df["pocr_lrb_luby_var_mse"] * df["pocr_lrb_luby_var_mse"] / df['simp_num_vars']
    #df["pocr_lrb_luby_var_rmse"] = np.sqrt(df["pocr_lrb_luby_var_mse"])

    #df["pocr_lrb_luby_lit_mse"] = df["pocr_lrb_luby_lit_mse"] * df["pocr_lrb_luby_lit_mse"] / df['simp_num_vars']
    #df["pocr_lrb_luby_lit_rmse"] = np.sqrt(df["pocr_lrb_luby_lit_mse"])


    # df['relative_var_picks_change'] = df["pocr_lrb_luby_var_mse"] * df['simp_num_vars_recip']

    for i in ["luby", "ar", "nr"]:
        df["proof_lrb_" + i + "_pf_var_mse"] = df["proof_lrb_" + i + "_pf_var_mse"] * df["proof_lrb_" + i + "_pf_var_mse"] / df['simp_num_vars']
        df["proof_lrb_" + i + "_pf_var_rmse"] = np.sqrt(df["proof_lrb_" + i + "_pf_var_mse"])

        df["proof_lrb_" + i + "_pf_lit_mse"] = df["proof_lrb_" + i + "_pf_lit_mse"] * df["proof_lrb_" + i + "_pf_lit_mse"] / (2*df['simp_num_vars'])
        df["proof_lrb_" + i + "_pf_lit_rmse"] = np.sqrt(df["proof_lrb_" + i + "_pf_lit_mse"])
    print(df['proof_lrb_luby_pf_var_rmse'])
    # df['relative_var_proof_occs_change'] = df["proof_lrb_luby_pf_var_mse"] / df['simp_num_vars_recip']

    # add extra ratios
    create_extra_ratios(df)

    df['time'] = df['simp_maplecomsps'] # - df['simp_maplecomsps_pp']
    # df['time'] = df['simp_lingeling'] # - df['simp_maplecomsps_pp']


    df['benchmark'] = df['benchmark'].map({'agile': 'Agile',
                                           'crafted': 'Crafted',
                                           'random': 'Random',
                                           'app': 'Application'
                                           })

    df.drop(df[df.simp_num_vars < 1].index, inplace=True)
    # print(len(df))
    df.drop(df[df.sub_benchmark== 'app'].index, inplace=True)
    # print(len(df))
    # sys.exit()

    # drop any duplicate entries
    df.drop_duplicates(inplace=True)
    #print(df)
    #sys.exit()

    # collect_refined_app_data(df, benchmarks, tables_dir, tex)

    collect_data(df, benchmarks, tables_dir, "", tex)
    data_summary(df, benchmarks, tables_dir + "" + "datasummary.tex",
                 names=["simp_lsr_size",
                        "simp_q",
                        "simp_alpha_var",
                        "simp_weak_size",
                        "simp_tw_upper",
                        "simp_backbones",
                        "dimacs_resolutions"],
                 columns=["LSR",
                          "Cmty",
                          "Alpha+Dim",
                          "Weak",
                          "TW",
                          "Bones",
                          "Merge+Res"], caption_prefix="")

    data_summary(df, benchmarks, tables_dir + "proof_" + "datasummary.tex",
                 names=["proof_lrb_luby_pf_var_mse",
                        "proof_lrb_luby_pf_lit_mse",
                        "proof_lrb_luby_avg_merges",
                        "proof_lrb_luby_gini_cmty_clauses",
                        "proof_lrb_luby_q_original_partition"],
                 columns=["Proof Var MSE",
                          "Proof Lit MSE",
                          "Avg Merges",
                          "Gini Cmty Clauses",
                          "Proof Q Orig"], caption_prefix="Proof ")

    average_metric_values(df, benchmarks, tables_dir + "proof_" + "averagemetrics.tex",
                          columns=["benchmark",
                                   "proof_lrb_luby_pf_var_mse",
                                 "proof_lrb_luby_pf_lit_mse",
                                 "proof_lrb_luby_avg_merges",
                                 "proof_lrb_luby_gini_cmty_clauses",
                                 "proof_lrb_luby_q_original_partition"],
                          names=[ "Proof Var MSE",
                                   "Proof Lit MSE",
                                   "Avg Merges",
                                   "Gini Cmty Clauses",
                                   "Proof Q Orig"],
                          caption_prefix="Proof ")

    average_metric_values(df, benchmarks, tables_dir + "" + "averagemetrics.tex",
                          columns=["benchmark",
                                   "simp_lvr",
                                   "simp_q",
                                   "simp_fractal_dim_vig",
                                   "simp_fractal_dim_cvig",
                                   "simp_wvr",
                                   "simp_backbonesvr",
                                   "simp_tw_uppervr"
                                   ],
                          names=["LSR/V",
                                 "Q",
                                 "DimVIG",
                                 "DimCVIG",
                                 "Weak/V",
                                 "Bones/V",
                                 "TW/V"
                                 ],
                          caption_prefix="")

    average_metric_values(df, benchmarks, tables_dir + "" + "averagemetrics2.tex",
                          columns=["benchmark",
                                   "dimacs_resolutions",
                                   "dimacs_merges",
                                   "dimacs_merge_resolutions_cmty_avg"],
                          names=[
                                 "Res",
                                 "Merge",
                                 "CM"],
                          caption_prefix="",
                          fstr=FINE_FSTR)

    average_metric_values(df, benchmarks, tables_dir + "proof_" + "popularity_vars.tex",
                          columns=["benchmark",
                                   "simp_num_vars_recip",
                                   #"relative_var_picks_change",
                                   #"relative_var_proof_occs_change",
                                   "pocr_lrb_luby_var_rmse",
                                   "proof_lrb_luby_pf_var_rmse"],
                          names=["1/#Vars",
                                 #"Picks Var MSE / AvgPop",
                                 #"Proof Var MSE / AvgPop",
                                 "Picks Var RMSE",
                                 "Proof Var RMSE"],
                          caption_prefix="POPULARITY VARS. ",
                          fstr=FINE_FSTR)

    average_metric_values(df, benchmarks, tables_dir + "proof_" + "popularity_lits.tex",
                          columns=["benchmark",
                                   "simp_num_vars_recip",
                                   # "relative_var_picks_change",
                                   # "relative_var_proof_occs_change",
                                   "pocr_lrb_luby_lit_rmse",
                                   "proof_lrb_luby_pf_lit_rmse"],
                          names=["1/#Vars",
                                 # "Picks Var MSE / AvgPop",
                                 # "Proof Var MSE / AvgPop",
                                 "Picks Lit RMSE",
                                 "Proof Lit RMSE"],
                          caption_prefix="POPULARITY LITS. ",
                          fstr=FINE_FSTR)

    average_metric_values(df, ["family"], tables_dir + "proof_" + "family_popularity_vars.tex",
                          columns=["family",
                                   "simp_alpha_var",
                                   "pocr_lrb_luby_var_rmse",
                                   "proof_lrb_luby_pf_var_mse"],
                          names=["Alpha Var",
                                 "Picks Var MSE",
                                 "Proof Var MSE"],
                          caption_prefix="POPULARITY VARS BY FAMILY. ",
                          fstr=FINE_FSTR)

    """
    average_metric_values(df, benchmarks, tables_dir + "proof_" + "popularity_lits.tex",
                          columns=["benchmark",
                                   "simp_alpha_var",
                                   "pocr_lrb_luby_lit_mse",
                                   "proof_lrb_luby_pf_lit_mse"],
                          names=["Alpha Var",
                                 "Picks Lit MSE",
                                 "Proof Lit MSE"],
                          caption_prefix="POPULARITY LITS (GRAPH BY ALPHA VAR). ",
                          fstr=FINE_FSTR)
    """

    average_metric_values(df, benchmarks, tables_dir + "proof_" + "cmty.tex",
                          columns=["benchmark",
                                   "simp_q",
                                   "pocr_lrb_luby_gini_picks",
                                   "pocr_lrb_luby_gini_clauses",
                                   "proof_lrb_luby_gini_cmty_clauses_normalized",
                                   "proof_lrb_luby_gini_cmty_occs_normalized",
                                   "proof_lrb_luby_q_original_partition"],
                          names=["Q",
                                 "Gini Picks",
                                 "Gini Clauses",
                                 "Gini Proof Occs",
                                 "Gini Proof Clauses",
                                 "Proof Q"],
                          caption_prefix="COMMUNITY. ",
                          fstr=None)

    tex.write("\\begin{table}[t]\n")
    for i in ["datasummary.tex", "averagemetrics.tex", "popularity_vars.tex", "popularity_lits.tex"]:  # "regression.tex", "same_insts_regression.tex", "bridgebd.tex"]:
        filename = tables_dir + "proof_" + i
        tex.write("\\input{" + filename + "}")
        tex.write("\n")
    tex.write("\\end{table}\n")

    tex.write("\\begin{table}[t]\n")
    for i in ["family_popularity_vars", "cmty.tex"]:  # "regression.tex", "same_insts_regression.tex", "bridgebd.tex"]:
        filename = tables_dir + "proof_" + i
        tex.write("\\input{" + filename + "}")
        tex.write("\n")
    tex.write("\\end{table}\n")

    tex.write("\\begin{table}[t]\n")

    for i in ["datasummary.tex", "averagemetrics.tex", "averagemetrics2.tex"]:  # "regression.tex", "same_insts_regression.tex", "bridgebd.tex"]:
        filename = tables_dir + "" + i
        tex.write("\\input{" + filename + "}")
        tex.write("\n")
    tex.write("\\end{table}\n")

    # sys.exit()

    b = ["Agile"]
    pocr_identifiers = ["gini_picks", "gini_clauses"]
    pocr_headers = ["Gini Picks"]
    pocr_captions = "Gini Picks"
    pocr_label = "tab:picks"
    sat = "UNSAT"

    # print(df.head())
    pocr_expt_params = [#(b, ["gini_picks", "gini_clauses"], ["Gini Picks", "Gini Clauses"],
                        # "Cmty Gini Expt.", "tab:gini"),

                        #(b, ["lsr_cmty_largest_ratio", "lsr_cmty_spread"], ["Largest Ratio", "Spread"], "LSR Cmty Expt.",
                        #"tab:lsr_cmty"),

                        #(b, ["lsr_vr", "all_decs", "time", "conflicts"], ["LSR", "All Decisions", "Time", "Conflicts"], "LSR Expt.",
                        #"tab:lsr_all_decs")
                        #print("need pocr lrb apps instances asap")
                        (["Agile"], ["lsr_vr", "all_decs_vr", "time", "conflicts"], ["LSR", "All Decisions", "Time (s)", "Conflicts"], "LSR Expt.",
                             "tab:lsr_agile", None),

                        (["Application"], ["lsr_vr", "all_decs_vr", "time", "conflicts"], ["LSR", "All Decisions", "Time (s)", "Conflicts"], "LSR Expt.",
                         "tab:lsr_app", ["LRB"]),

                        (["Agile", "Application", "Crafted", "Random"], ["gini_picks"], ["Gini Picks"],
                         "Cmty Expt.",
                         "tab:cmty_picks_all", ["LRB"]),

                        (["Agile", "Application", "Crafted", "Random"], ["gini_clauses"], ["Gini Clauses"],
                         "Cmty Expt.",
                         "tab:cmty_clauses_all", ["LRB"]),

                        (["Agile", "Application", "Crafted", "Random"], ["pf_var_rmse"], ["Pf Var RMSE"],
                         "Pop Expt.",
                         "tab:pf_pop_var_all", ["LRB"]),

                        (["Agile", "Application", "Crafted", "Random"], ["pf_lit_rmse"], ["Pf Lit RMSE"],
                         "Pop Expt.",
                         "tab:pf_pop_lit_all", ["LRB"]),

                        (["Agile", "Application", "Crafted", "Random"], ["gini_cmty_occs_normalized"],
                         ["Pf Gini Picks"],
                         "Cmty Proof Expt.",
                         "tab:pf_cmty_occs_all", ["LRB"]),

                        (["Agile", "Application", "Crafted", "Random"], ["gini_cmty_clauses_normalized"],
                         ["Pf Gini Clauses"],
                         "Cmty Proof Expt.",
                         "tab:pf_cmty_clauses_all", ["LRB"]),

                        #q_original_partition
                        (["Agile", "Application", "Crafted", "Random"], ["q_original_partition"],
                         ["Q"],
                         "Cmty Proof Q Expt.",
                         "tab:pf_cmty_q_all", ["LRB"]),




                        #(["Agile"], ["lsr_vr", "time", "conflicts"], ["LSR", "Time", "Conflicts"], "LSR Expt.",
                        # "tab:lsr_agile", ["LRB"]),

                        #(b, ["lsr_vr", "time", "conflicts"], ["LSR", "Time", "Conflicts"], "LSR Expt.",
                        #     "tab:lsr_agile")

    ]


    for p in pocr_expt_params:
        for s in [None]:# , "SAT", "UNSAT"]:
            print(*p)
            if "rmse" in str(p):
                pocr_structure_logging_summary(df, *p, tex, sat=s, fstr=FINE_FSTR)
            else:
                pocr_structure_logging_summary(df, *p, tex, sat=s)

    latex_gen.end_tex_article(tex)
    latex_gen.open_tex_pdf(tex_file)


if __name__ == '__main__':
    main()
