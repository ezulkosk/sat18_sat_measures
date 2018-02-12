import os
import random
import sys
import matplotlib.pyplot as plt
import numpy
from matplotlib.backends.backend_pdf import PdfPages
import math

from numpy.ma import arctanh, mean
from scipy.stats import spearmanr, pearsonr
from tabulate import tabulate

FSTR = "{0:.2f}"



def produce_scatter_plot_graph(graph, pdf):
    # print("Drawing ", graph[0])
    name, points = graph
    X, Y, XSAT, YSAT, XUNSAT, YUNSAT = [], [], [], [], [], []
    for i in points:
        if i[2] == "SAT":
            XSAT.append(i[0])
            YSAT.append(i[1])
        else:
            XUNSAT.append(i[0])
            YUNSAT.append(i[1])
        X.append(i[0])
        Y.append(i[1])

    plt.plot(X, Y, linestyle="-", marker="", linewidth=2.0)
    plt.plot(XSAT, YSAT, linestyle="", marker="+", color='r', markersize=7, label="SAT", linewidth=2.0)
    plt.plot(XUNSAT, YUNSAT, linestyle="", marker="o", color='b', markersize=7, label="UNSAT", linewidth=2.0)
    plt.xlabel("Num Increased Merges")
    plt.ylabel("Time (s)")
    plt.title(name)
    plt.legend()
    pdf.savefig(figsize=(18,18))

    plt.close()

def main(base_dir):
    cnf_dir = base_dir + "cnf/"
    maplesat_dir = base_dir + "maplesat/"
    merges_dir = base_dir + "merges/"
    proof_dir = base_dir + "proof_out/"

    time_map = {}
    res_map = {}
    clause_length_map = {}
    proof_map = {}
    names = []
    actual_merges = []

    with open(base_dir + "actual_merges") as stream:
        for line in stream:
            actual_merges.append(line.strip())

    merge_files = os.listdir(merges_dir)
    merges_map = {}
    # print("SKIPPING")
    for f in merge_files:
        arr = f.split(".")
        if 'merge' in arr:
            name = ".".join(arr[:-3])
            num_merge = int(arr[-2])
            full_name = name + ".merge." + str(num_merge)
            main = False
        else:
            #if "T100" in f or "T3" in f:
            #    name = arr[0]
            #else:
            if "unif" in f:
                name = ".".join(arr[:-1])
            else:
                name = arr[0] + "." + arr[1]
            num_merge = 0
            full_name = name
            main = True


        merge_file = merges_dir + full_name + ".merges"
        # print(name, merge_file, full_name)
        # if main:
        # print(merge_file)
        with open(merge_file) as stream:
            for line in stream:
                if line.startswith("NumMerge"):
                    arr = line.strip().split(",")
                    merges_map[(name, num_merge)] = int(arr[1])

    proof_files = os.listdir(proof_dir)
    proof_map = {}
    # print("SKIPPING")
    for f in proof_files:
        arr = f.split(".")
        if 'merge' in arr:
            name = ".".join(arr[:-3])
            num_merge = int(arr[-2])
            full_name = name + ".merge." + str(num_merge)
            main = False
        else:
            # if "T100" in f or "T3" in f:
            #    name = arr[0]
            # else:
            if "unif" in f:
                name = ".".join(arr[:-1])
            else:
                name = arr[0] + "." + arr[1]
            num_merge = 0
            full_name = name
            main = True

        proof_file = proof_dir + full_name + ".proof_analyses"
        # print(name, merge_file, full_name)
        # if main:
        # print(merge_file)
        useful_clauses, total_clauses, useless_input, merges_proof, merges_locality_avg, merges_normalized = None, None, None, None, None, None
        with open(proof_file) as stream:
            for line in stream:
                if line.startswith("Useful"):
                    arr = line.strip().split()
                    useful_clauses = arr[1]
                    total_clauses = arr[3]
                elif line.startswith("Merges,"):
                    arr = line.strip().split(",")
                    merges_proof = arr[1]
                elif line.startswith("MergeLocalityAverage,"):
                    arr = line.strip().split(",")
                    merges_locality_avg = arr[1]
                elif line.startswith("MergeLocalityNormalizedByNumDeps,"):
                    arr = line.strip().split(",")
                    merges_normalized = arr[1]
                elif line.startswith("UselessInput"):
                    arr = line.strip().split()
                    useless_input = arr[1]
        if useful_clauses and total_clauses and useless_input and merges_proof and merges_locality_avg and merges_normalized:
            proof_map[(name, num_merge)] = (int(useful_clauses), int(total_clauses), int(useless_input),
                                            int(merges_proof), float(merges_locality_avg),
                                            float(merges_normalized))


    maplesat_files = os.listdir(maplesat_dir)
    main_file_times = []
    for f in maplesat_files:
        if "merge" in f and f[:-len(".out")] not in actual_merges:
            continue
        # print(f)

        arr = f.split(".")
        #for i in arr:
        #    if i.startswith('s'):
        #        name += i
        #        break
        #    else:
        #        name += i + "."
        if 'merge' in arr:
            name = ".".join(arr[:-3])
            num_merge = int(arr[-2])
            full_name = name + ".merge." + str(num_merge)
            if name not in names:
                names.append(name)
            main = False
        else:
            name = arr[0] + "." + arr[1]
            num_merge = 0
            full_name = name
            if name not in names:
                names.append(name)
            main = True

        maplesat_file = maplesat_dir + full_name + ".out"
        res = "UNKNOWN"
        time = ""
        if os.path.exists(maplesat_file):
            with open(maplesat_file) as stream:
                clause_length_list = []
                for line in stream:
                    if line.startswith("SAT"):
                        res = "SAT"
                    elif line.startswith("UNSAT"):
                        res = "UNSAT"
                    elif line.startswith("CPU "):
                        time = line.strip().split()[3]
                    elif "% |" in line:
                        clause_length_list.append(int(line.split()[9]))


        # print(name, num_merge, res, time)
        if time != "":
            time_map[(name, num_merge)] = time
            res_map[(name, num_merge)] = res
            clause_length_map[(name, num_merge)] = clause_length_list


    names = sorted(names)
    graphs = []
    # print(len(names))
    max_merges = {}
    for name in names:
        # print("------")
        # print(name)
        points = []
        for i in range(0, 501, 1):
            if (name, i) in time_map.keys():
                # print(i, time_map[(name, i)], res_map[(name,i)])
                points.append((i, float(time_map[(name, i)]), res_map[(name,i)], clause_length_map[(name, i)],
                               proof_map.get((name, i), None)))
                max_merges[name] = i
        # print([float(i[1]) for i in points])
        if len(points) == 0 or max([float(i[1]) for i in points]) < 1:
            continue
        graphs.append((name, points))

    # print(" ".join(["merge_cnf/" + k + ".merge." + str(v) for k, v in max_merges.items()]))

    spearman_corrs = {}
    pearson_corrs = {}
    times = {}
    start_times = {}
    end_times = {}
    merges_initial = {}
    merges_final = {}
    if "random" in base_dir:
        ts = ["SAT", "UNSAT"]
    else:
        ts = [
            "1.8",
            "1.9",
            "2.0",
            "2.1",
            "2.2",
            "2.3",
            "2.4",
            "2.5",
            # "3.0",
            "5.0",
            "10.0",
            "100.0"
        ]

    for t in ts:
        times[t] = []
        spearman_corrs[t] = []
        pearson_corrs[t] = []
        start_times[t] = []
        end_times[t] = []
        merges_initial[t] = []
        merges_final[t] = []

    sat_count = 0
    unsat_count = 0
    with PdfPages(base_dir + "scaling.pdf") as pdf:
        # each figure has specified graph and transform, one page for each graph
        for graph in graphs:
            print("---------------")
            if len(graph[1]) > 1:
                X = [i[0] for i in graph[1]]
                Y = [i[1] for i in graph[1]]
                LENS = [i[3] for i in graph[1]]
                merge_data = [i[4] for i in graph[1]]
                #for l in LENS:
                #    print(l)
                for l in merge_data:
                    if l and l[3] and l[1]:
                        print(l[3] / l[1])
                        # print(l[2]/l[1])
                # print(LENS)
                RES = [i[2] for i in graph[1]]
                if len(set(RES)) > 1:
                    #print("HIT")
                    #print(graph[0])
                    continue
                t = graph[0]
                #print(t)
                #print(graph)
                if "unif" in t:
                    if graph[1][0][2] == "UNSAT":
                        t = "UNSAT"
                        unsat_count += 1
                    else:
                        t = "SAT"
                        sat_count += 1
                else:
                    t = t.split("T")[1]

                    t = str(t.split("_")[0])

                # print(t)
                start_times[t] = start_times[t] + [graph[1][0][1]]
                end_times[t] = end_times[t] + [graph[1][-1][1]]
                # print(graph[0])
                merges_initial[t] = merges_initial[t] + [merges_map[(graph[0],0)]]
                for k in merges_map.keys():
                    if k[0].startswith(graph[0]) and k[1] != 0:
                        merges_final[t] = merges_final[t] + [merges_map[k]]

                # print(t)
                spearman_corrs[t] = spearman_corrs[t] + [spearmanr(X, Y)[0]]
                pearson_corrs[t] = pearson_corrs[t] + [pearsonr(X, Y)[0]]

                # produce_scatter_plot_graph(graph, pdf)
                for i in graph[1]:
                    # print(i)
                    if i[2] == "SAT":
                        # produce_scatter_plot_graph(graph, pdf)
                        break

    print("SAT", sat_count, "UNSAT", unsat_count)
    rows = []
    for k in ts:
        fisher = []
        for v in spearman_corrs[k]:
            fisher.append(arctanh(v))
        spearman_avg = mean(fisher)
        fisher = []
        for v in pearson_corrs[k]:
            fisher.append(arctanh(v))
        pearson_avg = mean(fisher)
        # print(k, math.tanh(spearman_avg), math.tanh(pearson_avg), start_times[k], end_times[k])
        start = FSTR.format(mean(start_times[k])) +" (" + FSTR.format(numpy.std(start_times[k])) +")"
        end = FSTR.format(mean(end_times[k])) +" (" + FSTR.format(numpy.std(end_times[k])) +")"
        rows.append([k, FSTR.format(math.tanh(spearman_avg)), FSTR.format(math.tanh(pearson_avg)),
                     FSTR.format(mean(merges_initial[k])), start, FSTR.format(mean(merges_final[k])), end ])
    print(tabulate(rows, headers=["T", "Spearman", "Pearson", "Initial Merges", "Initial Time (s)",
                                  "Final Merges", "Final Time (s)"], tablefmt="latex"))
    # print("Base File Avg Time:", mean(main_file_times))

def create_popsim_cnf():
    ts = [
        "1.8",
        "1.9",
        "2.0",
        "2.1",
        "2.2",
        "2.3",
        "2.4",
        "2.5",
        # "3.0",
        "5.0",
        "10.0",
        "100.0"
    ]
    ts.reverse()
    B=0
    K=3
    k=0
    for i in ts:
        for j in range(0, 10):
            if float(i) < 2:
                n = "5000"
                m = "21250"
            elif float(i) < 2.2:
                n = "2000"
                m = "8500"
            #elif float(i) < 3:
            #    #n = "2500"
            #    #m = "10625"
            #    n = "1000"
            #    m = "4250"
            #    continue
            #elif float(i) < 5:
            #    n = "600"
            #    m = "2550"
            else:
                n = "300"
                m = "1275"
            s = str(random.randint(1, 5000))
            cmd = "timeout 120 ./psc3.1 -n " + n + " -m " + m + " -K 3 -B 0 -b 0.5 -k 0 -T " + i + " -s " + s + " > " + \
                "cnf/ps_n" + n + "_m" + m + "_T" + i + "_s" + s + ".cnf;"
            print(cmd)




if __name__ == '__main__':
    # create base cnf files
    #create_popsim_cnf()
    # gen merges of base cnf files: use genexpt_create_merge_files.sh
    # get times of all files: genexpt_create_maplesat_jobs.sh, genexpt_create_maplesat_merge_jobs.sh
    # get #merges of first and last: TODO
    # run this script's main
    base_dir = "/home/ezulkosk/cp2017_benchmarks/popsim/"
    # base_dir = "/home/ezulkosk/cp2017_benchmarks/popsim_random/"

    main(base_dir)