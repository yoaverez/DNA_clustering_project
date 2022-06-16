import argparse
import os
import shutil
import threading
import random
from datetime import datetime
from collections import Counter

import numpy as np

# from simulator import *
from metrics import *
from evaluation import *
from gen_input import *
import time
import math
# random.seed(datetime.now())

Number_of_steps = 220
Windows_size = 4
Similarity_threshold = 15


def file_to_clustering(file_path):
    """
    Return a tuple of all the attributes for ClusteringInfo object from a givan file in the evyat format.
    The result tuple contains:
    |   clustering: A dict in the form: {cluster_id (int): list of the cluster reads (list of strings)}
    |   original_strand_dict: {strand_id (int): the actual read (string)}
    |   reads_err_original_strand_dict: {read_id (int): the id of the origin strand of the read (int)}
    |   reads_err: A list of all the reads. the i'th element is the read with read_id = i.
    |   reads_err_dict: {read_id (int): the read itself (string)}
    |   C_dict: { cluster rep (full str) : List of all the reads that belong to that cluster }
    |   C_reps: [(Read, Cluster rep of the cluster to which the read belongs to)] must be sorted lexicographic.
    """
    reads_err = []  # will have all the reads from the sequencing stage
    reads_err_dict = {}
    strand_id = 0
    clustering = {}
    original_strand_dict = {}  # map from orig strand id to the actual strand
    reads_err_original_strand_dict = {}  # map from read_err to it's orig strand id
    C_dict = {}  # { cluster rep (full str) : List of all the reads that belong to that cluster }
    C_reps = []  # [(Read, Cluster rep of the cluster to which the read belongs to)] must be sorted lexicographic by read.
    with open(file_path, 'r') as evyat_f:
        line = "j"
        while line:
            line = evyat_f.readline()
            if not line:
                break
            original_strand_dict.update({strand_id: line.strip()})
            clustering[strand_id] = []
            line = evyat_f.readline()  # line == '*****************************\n':
            line = evyat_f.readline()
            cluster = []
            while line != '\n':
                striped_line = line.strip()
                reads_err.append(striped_line)
                clustering[strand_id].append(striped_line)
                reads_err_dict.update({len(reads_err)-1: striped_line})
                reads_err_original_strand_dict.update({len(reads_err)-1: strand_id})
                C_reps.append((striped_line, original_strand_dict[strand_id]))
                cluster.append(striped_line)
                line = evyat_f.readline()
            C_dict.update({original_strand_dict[strand_id]: cluster})
            strand_id = strand_id + 1
            line = evyat_f.readline()
    C_reps = sorted(C_reps, key=lambda x: x[0])
    return clustering, original_strand_dict, reads_err_original_strand_dict, reads_err, reads_err_dict, C_dict, C_reps


def hash_fun(x, a, w, l):
    """
    An implementation of the hash function from the article section 3.2.
    It finds the first occurrence of substring a in x and return the substring in x that starts
    with the starting index of a and its length is w + l.
    """
    # assume that x has substring a
    # but what happens if it return -1???
    ind = x.find(a)
    return x[ind:min(len(x), ind + w + l)]


def dna_str_to_decimal(st):
    """
    return the value of the DNA string in base 10.
    We interpreted DNA string as a number in base 4
    where A = 0 | C = 1 | G = 2 | T = 3
    """
    # All 3-grams: AAA,AAG,AAC,AAT,AGA,AGC,...,TTA,TTG,TTC,TTT
    # A-0, G-1, C-2, T-3
    # index of CAT = Decimal representaion of (203)
    # it's actually 2*4^2 + 0*4^1 + 3*4^0 its quad not decimal
    N_q = {"A": 0, "C": 1, "G": 2, "T": 3}
    dec = 0
    for i in range(0, len(st)):
        dec += N_q[st[i]] * (4 ** (len(st) - i - 1))
    return dec


def bin_sig(x, q):
    """ An implementation of the binary signature from the article section 3.3"""
    bs = [0] * (4 ** q)  # if bs[i] = 1 it means that there is at least 1 substring
                         # of x that in quad base means i
                         # for example: if x is CAT... then CAT is substring in x
                         # and it's value is 203(in base 4) = i(in base 10)
    for i in range(0, len(x) - q + 1):
        st = x[i:i + q]
        bs[dna_str_to_decimal(st)] = 1
    bs_str = ''.join(str(e) for e in bs)
    return bs_str


def rep_find(inp, reads_leaders):
    """
    | Args:
    |    inp is a read id
    |    parent[i] = the minimum id of all the read in the cluster that contains read i
    |    e.g. if read i is in the cluster {k, i, j} where k < i < j
    |    than parent[k], parent[i], parent[j] = k
    | return:
    |    the id of the cluster that contains the read with id=ind
    """
    temp = inp
    cnt = 0
    while reads_leaders[temp] != temp and cnt < 10:
        cnt += 1
        temp = reads_leaders[temp]
    return temp


def create_clustering(reads_leaders):
    """
    | Args:
    |     reads_leaders[i] = the minimum id of all the read in the cluster that contains read i
    |     e.g. if read i is in the cluster {k, i, j} where k < i < j
    |     than reads_leaders[k], reads_leaders[i], reads_leaders[j] = k

    | return:
    |    Clustering that is construct from the parent list.
    |    The output is a dict where the keys are the id of a cluster
    |    and the values are the cluster - list of reads id.
    """
    clstr = {}
    for i in range(0, len(reads_leaders)):
        clstr[i] = []
    for i in range(0, len(reads_leaders)):
        clstr[rep_find(i, reads_leaders)].append(i)
    return clstr


def min_max(val1, val2):
    """ return the tuple: (min{val1, val2}, max{val1, val2})"""
    min_val = min(val1, val2)
    max_val = max(val1, val2)
    return min_val, max_val


def condition0(bin_sig_arr, hash_C_til, id1, read1, read2, th_low, th_high, r, index_size=0):
    """
    The original condition in the algorithm for merging two clusters
    """
    return (((ham_dis(bin_sig_arr[hash_C_til[id1][0]], bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_low) or
             ((ham_dis(bin_sig_arr[hash_C_til[id1][0]],
                       bin_sig_arr[hash_C_til[id1 + 1][0]]) <= th_high) and
              edit_dis(read1, read2) <= r)))


def hash_based_cluster(reads_err, number_of_steps=Number_of_steps, windows_size=Windows_size,
                       similarity_threshold=Similarity_threshold, index_size=0, cond_func=condition0):
    """
    Implementation of the microsoft hash based clustering algorithm.
    | Args:
    |   reads_err: reads_err: A list in which the i'th element is the string DNA of the read with id = i
    |   number_of_steps: The number of iterations of the algorithm. Default is 220.
    |   windows_size: A parameter for calculate the binary signature. Default is 4.
    |   similarity_threshold: A bound for the algorithm. Default is 15.
    |   index_size: the number of chars in each index. e.g. for AAA,AAC,AAG,AAT,.... the size is 3.
    Returns a clustering. A list of clusters, each cluster in a sorted list if all the cluster reads ids.
    """
    # reads_err = clustering_info.reads_err  # will have all the reads from the sequencing stage

    reads_err_ind = [0] * (len(reads_err))
    read_leaders = [0] * (len(reads_err))
    bin_sig_arr = []

    local_comm = number_of_steps

    q = windows_size
    # computing the binary signature for all the reads
    for i in range(0, len(reads_err)):
        reads_err_ind[i] = (i, reads_err[i][index_size:])
        read_leaders[i] = i
        bin_sig_arr.append(bin_sig(reads_err[i][index_size:], q))

    C_til = create_clustering(read_leaders)  # short for C^(~) i.e. C tilda
    dist_arr = []

    # th_low, th_high:
    # finding Hamming distance from the first read
    for i in range(1, len(reads_err)):
        dist_arr.append(ham_dis(bin_sig_arr[i], bin_sig_arr[0]))

    # sort the dist array
    dist_arr.sort()
    # the weird definition of theta_low and theta_high
    # I should probably ask why this specific theta's
    for i in range(0, 1000):
        if dist_arr[i + 1] - dist_arr[i] > 10:
            k = i
            break
    th_low = min(dist_arr[i] + 5, dist_arr[i + 1] - 8)
    th_high = min(dist_arr[i + 1] - 5, th_low + 10)

    r = similarity_threshold
    # set parameters from the article section 5.1:
    w = math.ceil(math.log(len(reads_err[0][index_size:]), 4))
    l = math.ceil(math.log(len(reads_err), 4))
    #######################################
    # Add shuffle on reads_err_ind
    reads_err_ind_shuffle = random.sample(reads_err_ind, k=len(reads_err_ind))
    reads_err_ind_shuffle_ids = [0] * len(reads_err_ind_shuffle)
    for index, (read_id, read) in enumerate(reads_err_ind_shuffle):
        reads_err_ind_shuffle_ids[read_id] = index
    # print(reads_err_ind[0:10])
    # print(reads_err_ind_shuffle[0:10])
    #######################################
    for lcl_step in range(0, local_comm):
        # report_func(total_bar_size, (2 * tmp_bar) + lcl_step)
        hash_select = [0] * (len(reads_err))
        # picking random representatives samples from each cluster:
        for i in range(0, len(C_til)):
            if len(C_til[i]) != 0:
                hash_select[random.choice(C_til[i])] = 1

        # computing the hash value for each representative:
        a = rand_perm(w)
        hash_C_til = [0] * (len(reads_err))
        for i in range(0, len(reads_err_ind_shuffle)):
            j = reads_err_ind_shuffle_ids[i]
            if hash_select[i] == 0:
                hash_C_til[i] = (reads_err_ind_shuffle[j][0], "")
            else:
                hash_C_til[i] = (reads_err_ind_shuffle[j][0], hash_fun(reads_err_ind_shuffle[j][1], a, w, l))
        # sort the hash values by the strings values
        hash_C_til.sort(key=lambda x: x[1])

        cnt = 0
        # implementation of lines 8-10 in the algorithm in the article:
        for i in range(0, len(hash_C_til) - 1):
            if hash_C_til[i][1] == "":
                continue
            else:
                if hash_C_til[i][1] == hash_C_til[i + 1][1]:
                    x = reads_err[hash_C_til[i][0]][index_size:]
                    y = reads_err[hash_C_til[i + 1][0]][index_size:]
                    # (edit_dis(x[:5], y[:5]) <= 3) and
                    # if ((index_size == 0 or (edit_dis(x[:index_size], y[:index_size]) <= 3)) and
                    #     ((ham_dis(bin_sig_arr[hash_C_til[i][0]], bin_sig_arr[hash_C_til[i + 1][0]]) <= th_low) or
                    #         ((ham_dis(bin_sig_arr[hash_C_til[i][0]],
                    #                        bin_sig_arr[hash_C_til[i + 1][0]]) <= th_high) and
                    #          edit_dis(x, y) <= r))):
                    if cond_func(bin_sig_arr, hash_C_til, i, x, y, th_low, th_high, r, index_size):
                        cnt += 1
                        min_temp, max_temp = min_max(rep_find(hash_C_til[i][0], read_leaders),
                                                          rep_find(hash_C_til[i + 1][0], read_leaders))
                        # merging x and y clusters
                        C_til[min_temp].extend(C_til[max_temp])
                        C_til[max_temp] = []
                        read_leaders[max_temp] = min_temp

    return [sorted(x) for x in list(C_til.values()) if x != []], bin_sig_arr


def str_clustering_to_ids(cluster_info):
    """
    Make a the evyat file like the output of the algorithm.
    | Args:
    |   cluster_info: An instance of class ClusteringInfo
    Return:
        A list of clusters. Each cluster is a list of reads_ids
    """
    clustering = []
    for i in range(len(cluster_info.original_strand_dict())):
        clustering.append([])
    for key, value in cluster_info.reads_err_original_strand_dict.items():
        clustering[value].append(key)
    return [sorted(cluster) for cluster in clustering]


def find_best_cluster(singleton_id, singleton_cluster_id, candidates_ids, reads_err,
                      clusters_reps, bin_sign_arr, index_size=6):
    """
    An auxiliary for the handle_singletons function.

    :return: The cluster id of the most likely cluster to be the origin of the given singleton.
    :param singleton_id: The singleton whom we want to get his origin cluster.
    :param singleton_cluster_id: The singleton cluster id in case we don't find a good origin cluster.
    :param candidates_ids: List of all the clusters ids which are candidates to be the singleton origin cluster.
    :param reads_err: mapping from read_id (int) to the actual read (string).
    :param clusters_reps: A dict of the form: {cluster_id: list of representatives ids for this cluster}
    :param bin_sign_arr: The binary signatures of all the reads. bin_sign_arr[i] = the signature of the read with id=i.
    :param index_size: The number of symbols dedicated for the read's index.
    """
    best_cluster_id = singleton_cluster_id
    # First filter
    candidates_ids = [c for c in candidates_ids if c != singleton_cluster_id]
    compare = []
    true_candidates_sizes = []
    for candidate_id in candidates_ids:
        reps = clusters_reps[candidate_id][0]
        ham = 0
        for i, rep in enumerate(reps):
            ham += ham_dis(bin_sign_arr[rep], bin_sign_arr[singleton_id])
        compare.append(ham / len(reps))
    if len(compare) != 0:
        # true candidates:
        min_ham = min(compare)
        min_ham_cluster_id = candidates_ids[np.argmin(compare)]
        # condition for not good enough candidates:
        if not (min_ham <= 95
                or (edit_dis(clusters_reps[min_ham_cluster_id][1], reads_err[singleton_id][: index_size]) <= 1
                    and min_ham <= 90)):
            return singleton_cluster_id

        # otherwise there are good candidates:
        true_candidates = [cluster_id for index, cluster_id in enumerate(candidates_ids) if compare[index] <= min_ham + 6]
        true_candidates_compare = [compare[index] for index, cluster_id in enumerate(candidates_ids) if compare[index] <= min_ham + 6]
        true_candidates_sizes.append(len(true_candidates))
        best_cluster_id = min_ham_cluster_id
        inside_ham_dist = {}
        if len(true_candidates) > 2:
            abs_diffs = {}
            for index, true_can in enumerate(true_candidates):
                ham = 0
                for i, rep1 in enumerate(clusters_reps[true_can][0]):
                    for j, rep2 in enumerate(clusters_reps[true_can][0]):
                        if i < j:
                            ham += ham_dis(bin_sign_arr[rep1], bin_sign_arr[rep2])
                inside_ham_dist[true_can] = 0 if ham == 0 else ham / ((len(clusters_reps[true_can][0]) ** 2 - len(clusters_reps[true_can][0])) / 2)
                abs_diffs[true_can] = abs(inside_ham_dist[true_can] - true_candidates_compare[index])
            best_cluster_id = min(abs_diffs.items(), key=lambda pair: pair[1])[0]

    return best_cluster_id


def handle_singletons_with_index_ver5_5(algo_clustering, orig_cluster_info, bin_sign_arr, index_size,
                                        threshold=10, num_of_hashes=20, num_epochs=2,
                                        converge=False, convergence_param=10,
                                        convergence_ratio=0.8, return_stats=False):
    """
    :return: clustering after handling singletons for a given clustering.
    :param algo_clustering: The clustering we want to improve. Should be list of lists.
    :param orig_cluster_info: Used both for stats and for maps between read ids to the actual reads.
    :param bin_sign_arr: The binary signature array that relates to the reads and is constant for those reads.
    :param index_size: The number of symbols dedicate for indexing.
    :param threshold: The maximum number of representatives we sample from each cluster. Default to 10.
    :param num_of_hashes: The number of different hash functions to use. Default to 20.
    :param num_epochs: The number of iterations to preform. Default to 2.
    :param converge: True if you want to run iterations until convergence (we don't handle a lot singletons).
        Default to False
    :param convergence_param: Used only if converge=True. Stop the iterations if we handled less than convergence_param
        singletons in the last iteration. Default to 10.
    :param convergence_ratio: Used only if converge=True. Stop the iterations if after the current iteration we remains
        with more than convergence_ratio singletons from the previous iteration. Default to 0.8.
    :param return_stats: If true, we also return the number of unwanted singletons at each epoch and the time it took
        for each epoch in minutes. Default to False.
    """
    algo_clustering_copy = copy.deepcopy(algo_clustering)
    times_for_each_epoch = [0]
    ########################################################################
    # can not use this vars because they have knowledge about the real clustering:
    num_of_remaining_singletons = []
    stat1, stat2 = find_clusters_stats(algo_clustering_copy, orig_cluster_info)
    num_of_remaining_singletons.append(len(find_unwanted_singletons(stat1)))
    ########################################################################
    number_of_singletons = None  # will represent both wanted and unwanted singletons.
                                 # Meaning sees only the algo clustering and not the real clustering.
    print(f"########################### start running handle singletons ###########################")
    epoch_id = 0
    while epoch_id < num_epochs or converge:
        start = time.perf_counter_ns()
        reads_err = orig_cluster_info.reads_err
        stat1, stat2 = find_clusters_stats(algo_clustering_copy, orig_cluster_info)
        w = math.ceil(math.log(len(reads_err[0][index_size:]), 4))
        l = math.ceil(math.log(len(reads_err), 4))-3
        hashing_a = []
        for i in range(num_of_hashes):
            a = rand_perm(w)
            hashing_a.append(a)

        candidates = {}
        hash_dict = {}
        # find singletons
        singletons = []
        for cluster_id, cluster in enumerate(algo_clustering_copy):
            if len(cluster) == 1:
                singleton_id = cluster[0]
                candidates[singleton_id] = {}
                singletons.append((singleton_id, cluster_id))
                for a in hashing_a:
                    hash_value = hash_fun(reads_err[singleton_id][index_size:], a, w, l)
                    if hash_value == '':
                        continue
                    if hash_value in hash_dict:
                        hash_dict[hash_value]['singletons_ids'].append(singleton_id)
                    else:
                        hash_dict[hash_value] = {'cluster_ids': [],
                                                 'singletons_ids': [singleton_id]}

        # meaning it's the first iteration.
        if number_of_singletons is None:
            number_of_singletons = len(singletons)
            print(f'number of singletons before: {number_of_singletons}')
        # convergence:
        elif converge and (number_of_singletons - len(singletons) <= convergence_param
                           or len(singletons)/number_of_singletons > convergence_ratio):
            number_of_singletons = len(singletons)
            break
        # regular epoch:
        else:
            number_of_singletons = len(singletons)
            print(f'number of singletons after {epoch_id} epochs: {number_of_singletons}')

        # find clusters reps and their index:
        clusters_reps = {}
        for cluster_id, cluster in enumerate(algo_clustering_copy):
            true_index = []
            representatives = random.sample(cluster, min([threshold, len(cluster)]))
            index_mat = np.array([list(reads_err[read][:index_size]) for read in representatives]).T
            for row in index_mat:
                row_list = list(row)
                true_index.append(max(row_list, key=row_list.count))
            str_index = ''.join(true_index)
            clusters_reps[cluster_id] = (representatives, str_index)

            for a in hashing_a:
                for rep in representatives:
                    hash_value = hash_fun(reads_err[rep][index_size:], a, w, l)
                    if hash_value in hash_dict:
                        hash_dict[hash_value]['cluster_ids'].append(cluster_id)
                        for singleton_id in hash_dict[hash_value]['singletons_ids']:
                            if cluster_id in candidates[singleton_id]:
                                candidates[singleton_id][cluster_id] += 1
                            else:
                                candidates[singleton_id][cluster_id] = 1

        # iterate over all the candidates for each singleton:
        for singleton_id, singleton_cluster_id in singletons:
            cluster_id = find_best_cluster(singleton_id, singleton_cluster_id,
                                           list(candidates[singleton_id].keys()), reads_err, clusters_reps,
                                           bin_sign_arr, index_size)
            if cluster_id != singleton_cluster_id:
                algo_clustering_copy[cluster_id].extend(algo_clustering_copy[singleton_cluster_id])
                algo_clustering_copy[singleton_cluster_id] = []

        algo_clustering_copy = [x for x in algo_clustering_copy if x != []]
        end = time.perf_counter_ns()
        elapsed_in_ns = end - start
        times_for_each_epoch.append((elapsed_in_ns * math.pow(10, -9))/60)
        stat1, stat2 = find_clusters_stats(algo_clustering_copy, orig_cluster_info)
        num_of_remaining_singletons.append(len(find_unwanted_singletons(stat1)))

        epoch_id += 1

    print(f'number of singletons after {epoch_id} epochs: {number_of_singletons}')
    if converge:
        print(f'converge after {epoch_id} epochs')
    print(f"########################### finish running handle singletons ###########################")
    if return_stats:
        return [sorted(x) for x in algo_clustering_copy if x != []], '',\
               num_of_remaining_singletons, times_for_each_epoch
    return [sorted(x) for x in algo_clustering_copy if x != []], ''


def separate_cluster(cluster, bin_sign_arr, threshold=10):
    """
    An auxiliary function for the handle_unions function.

    :return: If we think the given cluster is union of clusters then tuple of two disjoint clusters
        s.t. C1 union C2 = the given cluster. Otherwise, tuple of two None i.e. (None, None)
    :param cluster: The cluster id of the suspect to be a union of clusters.
    :param bin_sign_arr: The binary signature array that relates to the reads and is constant for those reads.
    :param threshold: The maximum number of representatives we sample from each cluster. Default to 10.
    """
    # finding cluster representatives:
    n = min([threshold, len(cluster)])
    reps = random.sample(cluster, n)
    # calculate the reps hamming distance from one another:
    ham_dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ham_dist_matrix[i][j] = ham_dis(bin_sign_arr[reps[i]], bin_sign_arr[reps[j]])

    # find the two reads that are farthest from each other
    hami, hamj = np.unravel_index(ham_dist_matrix.argmax(), ham_dist_matrix.shape)

    # finding the minimum hamming distance from the two read we found:
    ham_within_cluster_1 = min([ham_j for j, ham_j in enumerate(ham_dist_matrix[hami]) if hami != j])
    ham_within_cluster_2 = min([ham_i for i, ham_i in enumerate(ham_dist_matrix[hamj]) if hamj != i])

    diff1 = ham_dist_matrix[hami][hamj] - ham_within_cluster_1
    diff2 = ham_dist_matrix[hami][hamj] - ham_within_cluster_2

    # if the condition is true it means that the cluster is probably not a union of clusters.
    if ham_dist_matrix[hami][hamj] < 115 or diff1 < 30 or diff2 < 30:
        return None, None

    cluster1 = [reps[hami]]
    cluster2 = [reps[hamj]]
    # separate the cluster to two clusters based on the hamming distances we computed earlier:
    for read_id in cluster:
        if read_id != reps[hami] and read_id != reps[hamj]:
            ham1 = ham_dis(bin_sign_arr[cluster1[0]], bin_sign_arr[read_id])
            ham2 = ham_dis(bin_sign_arr[cluster2[0]], bin_sign_arr[read_id])
            if ham1 < ham2:
                cluster1.append(read_id)
            else:
                cluster2.append(read_id)

    # check that the two cluster we created are really not suppose to be in the same cluster:
    dists_inside_cluster1 = []
    for i in range(len(cluster1)):
        for j in range(i + 1, len(cluster1)):
            dists_inside_cluster1.append(ham_dis(bin_sign_arr[cluster1[i]], bin_sign_arr[cluster1[j]]))
    dists_inside_cluster2 = []
    for i in range(len(cluster2)):
        for j in range(i + 1, len(cluster2)):
            dists_inside_cluster2.append(ham_dis(bin_sign_arr[cluster2[i]], bin_sign_arr[cluster2[j]]))

    mean_dists_inside_cluster1 = 0 if len(dists_inside_cluster1) == 0 else np.mean(dists_inside_cluster1)
    mean_dists_inside_cluster2 = 0 if len(dists_inside_cluster2) == 0 else np.mean(dists_inside_cluster2)
    if mean_dists_inside_cluster1 > 100 and mean_dists_inside_cluster2 > 100:
        return None, None
    dists_between_clusters = []
    for i in range(len(cluster1)):
        for j in range(len(cluster2)):
            dists_between_clusters.append(ham_dis(bin_sign_arr[cluster1[i]], bin_sign_arr[cluster2[j]]))

    mean_dists_between_clusters = 0 if len(dists_between_clusters) == 0 else np.mean(dists_between_clusters)
    if mean_dists_between_clusters-min([mean_dists_inside_cluster1, mean_dists_inside_cluster2]) <= 25:
        return None, None
    return cluster1, cluster2


def handle_unions(algo_clustering, orig_cluster_info, bin_sign_arr, index_size, threshold=10, log=True):
    """
    :return: clustering after handling unions for a given clustering.
    :param algo_clustering: The clustering we want to improve. Should be list of lists.
    :param orig_cluster_info: Used both for stats and for maps between read ids to the actual reads.
    :param bin_sign_arr: The binary signature array that relates to the reads and is constant for those reads.
    :param index_size: The number of symbols dedicate for indexing.
    :param threshold: The maximum number of representatives we sample from each cluster. Default to 10.
    :param log: If True return also a string containing the log of this function. Default to True.
    """
    avg_cluster_size = sum([len(cluster) for cluster in algo_clustering])/len(algo_clustering)
    count_wrong = 0
    count_right = 0
    sep_dict = {}
    log_str = f'#################################################\n' \
              f'               handle_unions_log                 \n' \
              f'#################################################\n'
    # for checks only
    stats1, stats2 = find_clusters_stats(algo_clustering, orig_cluster_info)
    unwanted_unions = find_unwanted_unions(stats1)
    for cluster_id, cluster in enumerate(algo_clustering):
        if len(cluster) >= avg_cluster_size:
            cluster1, cluster2 = separate_cluster(cluster, bin_sign_arr, threshold=threshold)
            if cluster1 is not None:
                sep_dict[cluster_id] = (cluster1, cluster2)
                if cluster_id in unwanted_unions:
                    count_right += 1
                else:
                    log_str += f"{stats1[cluster_id]}\n"
                    log_str += f"{len(cluster1)=}\n{len(cluster2)=}\n"
                    count_wrong += 1
    log_str += f"{count_wrong=}\n"
    log_str += f"{count_right=}\n"
    for cluster_id, (cluster1, cluster2) in sep_dict.items():
        algo_clustering[cluster_id] = cluster1
        algo_clustering.append(cluster2)
    if log:
        return algo_clustering, log_str
    return algo_clustering


def arrange_clustering(algo_clustering, orig_cluster_info):
    """
    Arrange the output of the clustering algorithm so that the order of the clusters
    and the original strands will be similar to the original (true) clustering.

    :param algo_clustering: The output from the clustering algorithm.
    :param orig_cluster_info: all the known information on the true clustering.
        It should be ClusteringInfo object.
    :return: the evyat format string of the arranged clustering.
    """
    res = ''
    for cluster in algo_clustering:
        orig_strand_candidates = []
        for cluster_element in cluster:
            orig_strand_candidates.append(orig_cluster_info.reads_err_original_strand_dict.get(cluster_element))
        orig_strand_id = max(orig_strand_candidates, key=orig_strand_candidates.count)
        res += str(orig_cluster_info.original_strand_dict.get(orig_strand_id)) + '\n'
        res += '*****************************\n'
        for cluster_element in cluster:
            res += orig_cluster_info.read_err_dict.get(cluster_element) + '\n'
        res += '\n\n'
    return res


def file_to_cluster(file_path):
    """
    :return: A clustering from a file as a dict
        where the keys are cluster id and the value in a list of all the reads (strings) that are in the cluster.
    """
    strand_id = 0
    cluster = {}
    with open(file_path, 'r') as evyat_f:
        line = "j"
        while line:
            line = evyat_f.readline()
            if not line:
                break

            cluster[strand_id] = []
            line = evyat_f.readline()  # line == '*****************************\n':
            line = evyat_f.readline()
            while line != '\n':
                cluster[strand_id].append(line.strip())
                line = evyat_f.readline()
            strand_id = strand_id + 1
            line = evyat_f.readline()
    return cluster


def algo_clustering_to_file_aux(input_path, index_size):
    """
    An auxiliary function for the algo_clustering_to_file function.

    :param input_path: The full or relative path to an inputs reads file.
    :param index_size: The number of symbols dedicate for indexing.
    """
    input_file_name = input_path[input_path.rfind("/"):]
    output_path = "files/minion_idt/6000 strands in size 150 with x2 errors and cluster avg of 40/algo_results" + input_file_name.replace(".txt", "_algo_result.txt")
    clustering_info = ClusteringInfo(file_path=input_path)
    C_til, bin_sig_arr = hash_based_cluster(clustering_info.reads_err, index_size=index_size)
    print(C_til[0])
    with open(output_path, 'w', newline='\n') as f:
        for cluster in C_til:
            for read_id in cluster:
                f.write(f"{read_id}\n")
            f.write("***\n")
        f.write("bin_sign:\n")
        for i in range(len(bin_sig_arr)):
            f.write(f"{bin_sig_arr[i]}\n")


def algo_clustering_to_file(index_size):
    """
    Saving clustering into file.

    :param index_size: The number of symbols dedicate for indexing.
    """
    input_path = "files/minion_idt/6000 strands in size 150 with x2 errors and cluster avg of 40/evyat files/evyat0_index.txt"
    for i in range(5):
        # if i != 4:
        if True:
            curr_input_path = input_path.replace("_index.txt", f"{i}_index.txt")
            algo_clustering_to_file_aux(curr_input_path, index_size)


def file_to_algo_clustering(path):
    """
    :return: Tuple containing clustering (list of lists) and binary signature array.
    :param path: The full or relative path to a clustering file.
    """
    clustering = []
    bin_sig_arr = []
    cluster_id = 0
    clustering.append([])
    with open(path, 'r') as f:
        # find the clustering
        line = f.readline().strip()
        while line != "bin_sign:":
            if line[0] != '*':
                clustering[cluster_id].append(int(line))
                line = f.readline()
            else:
                line = f.readline().strip()
                if line != "bin_sign:":  # means new cluster
                    clustering.append([])
                    cluster_id += 1

        # find the bin_signs
        line = f.readline().strip()
        while line:
            bin_sig_arr.append(line)
            line = f.readline().strip()
    return clustering, bin_sig_arr


def from_no_index_to_index_via_indices_file(indices_file_path, input_strands_file_path):
    """
    Merge input file with indices file to compose indexed input.

    :param indices_file_path: The full or relative path to an indices file.
    :param input_strands_file_path: The full or relative path to an inputs reads file.
    """
    with open(indices_file_path, 'r', newline='\n') as ind_file:
        index_line = ind_file.readline().rstrip()
        index_size = len(index_line)
        with open(input_strands_file_path, 'r', newline='\n') as input_file:
            input_line = input_file.readline().rstrip()
            output_path = input_strands_file_path.replace('.txt', f'_index_{index_size}.txt')
            with open(output_path, 'w', newline='\n') as out:
                while input_line:
                    output_line = index_line + input_line + '\n'
                    out.write(output_line)
                    index_line = ind_file.readline().rstrip()
                    input_line = input_file.readline().rstrip()


class ClusteringInfo:
    """
    | Attributes:
    |   clustering: A dict in the form: {cluster_id (int): list of the cluster reads (list of strings)}
    |   original_strand_dict: {strand_id (int): the actual read (string)}
    |   reads_err_original_strand_dict: {read_id (int): the id of the origin strand of the read (int)}
    |   reads_err: A list of all the reads. the i'th element is the read with read_id = i.
    |   reads_err_dict: {read_id (int): the read itself (string)}
    |   C_dict: { cluster rep (full str) : List of all the reads that belong to that cluster }
    |   C_reps: [(Read, Cluster rep of the cluster to which the read belongs to)] must be sorted lexicographic.
    | methods:
    |   __init__: Initialize a class object via a given object of the same class
    |             or via a file in the format of evyat file.
    |   __str__: Return an evyat format string contains the output of the clustering algorithm.
    |   str_to_file: Writes the output of __str__ to the given file.
    """
    def __init__(self, clustering_info=None, file_path=None):
        clusters_info = None
        if clustering_info is not None:
            clusters_info = clustering_info
        elif file_path is not None:
            clusters_info = file_to_clustering(file_path)
        if clusters_info is not None:
            self.clustering = clusters_info[0]
            self.original_strand_dict = clusters_info[1]
            self.reads_err_original_strand_dict = clusters_info[2]
            self.reads_err = clusters_info[3]
            self.reads_err_dict = clusters_info[4]
            self.C_dict = clusters_info[5]
            self.C_reps = clusters_info[6]
        else:
            self.clustering = None
            self.original_strand_dict = None
            self.reads_err_original_strand_dict = None
            self.reads_err = None
            self.reads_err_dict = None
            self.C_dict = None
            self.C_reps = None

    def __str__(self, number_of_steps=Number_of_steps, windows_size=Windows_size,
                similarity_threshold=Similarity_threshold):
        """
        return a string that represent the output of the clustering algorithm.
        | Args for the clustering algorithm:
        |   number_of_steps: The number of iterations of the algorithm. Default is 220.
        |   windows_size: A parameter for calculate the binary signature. Default is 4.
        |   similarity_threshold: A bound for the algorithm. Default is 15.

        """
        algo_clustering = hash_based_cluster(self.reads_err, number_of_steps, windows_size, similarity_threshold)
        return arrange_clustering(algo_clustering, self)

    def str_to_file(self, file_path, number_of_steps=Number_of_steps, windows_size=Windows_size,
                    similarity_threshold=Similarity_threshold):
        """
        Write the output of the __str__ method to a given file.
        Meaning, write the result of the clustering algorithm to the given file.
        | Args for the clustering algorithm:
        |   number_of_steps: The number of iterations of the algorithm. Default is 220.
        |   windows_size: A parameter for calculate the binary signature. Default is 4.
        |   similarity_threshold: A bound for the algorithm. Default is 15.
        """
        with open(file_path, 'w', newline='\n') as f:
            f.write(self.__str__(number_of_steps, windows_size, similarity_threshold))
