#!/usr/bin/env python3

###############################################################################
# Program: TransmissionCluster.py
# Type: Python Script
# Version: 1.0
# Author: Steven J. Clipman
# Description: Efficient distance-free method for defining clusters from phylogenetic trees
# License: MIT
###############################################################################

from queue import Queue
from treeswift import read_tree_newick
import PySimpleGUI as sg
import os
import matplotlib.pyplot as plt
import math
import statistics

NUM_THRESH = 1000  # number of thresholds to calculate genetic distance over


# cut out the current node's subtree (by setting all nodes' DELETED to True) and return list of leaves
def cut(node):
    cluster = list()
    descendants = Queue()
    descendants.put(node)
    while not descendants.empty():
        descendant = descendants.get()
        if descendant.DELETED:
            continue
        descendant.DELETED = True
        descendant.left_dist = 0
        descendant.right_dist = 0
        descendant.edge_length = 0
        if descendant.is_leaf():
            cluster.append(str(descendant))
        else:
            for c in descendant.children:
                descendants.put(c)
    return cluster


# initialize properties of input tree and return set containing taxa of leaves
def prep(tree, support):
    tree.resolve_polytomies()
    tree.suppress_unifurcations()
    leaves = set()
    for node in tree.traverse_postorder():
        if node.edge_length is None:
            node.edge_length = 0
        node.DELETED = False
        if node.is_leaf():
            leaves.add(str(node))
        else:
            try:
                node.confidence = float(str(node))
            except:
                node.confidence = 100.  # give edges without support values support 100
            if node.confidence < support:  # don't allow low-support edges
                node.edge_length = float('inf')
    return leaves


# min_clusters_threshold_max, but all clusters must define a clade
def min_clusters_threshold_max_clade(tree, threshold, support):
    leaves = prep(tree, support)
    clusters = list()
    for node in tree.traverse_postorder():
        # if I've already been handled, ignore me
        if node.DELETED:
            continue

        # find my undeleted max distances to leaf
        if node.is_leaf():
            node.left_dist = 0
            node.right_dist = 0
        else:
            children = list(node.children)
            if children[0].DELETED and children[1].DELETED:
                cut(node)
                continue
            if children[0].DELETED:
                node.left_dist = 0
            else:
                node.left_dist = max(children[0].left_dist, children[0].right_dist) + children[0].edge_length
            if children[1].DELETED:
                node.right_dist = 0
            else:
                node.right_dist = max(children[1].left_dist, children[1].right_dist) + children[1].edge_length

            # if my kids are screwing things up, cut both
            if node.left_dist + node.right_dist > threshold:
                cluster_l = cut(children[0])
                node.left_dist = 0
                cluster_r = cut(children[1])
                node.right_dist = 0

                # add cluster
                for cluster in (cluster_l, cluster_r):
                    if len(cluster) != 0:
                        clusters.append(cluster)
                        for leaf in cluster:
                            leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters


# pick the threshold between 0 and "distance threshold" that maximizes number of (non-singleton) clusters
def argmax_clusters(method, tree, threshold, support, display_fig):
    if display_fig is True:
        distfile = open("TransmissionCluster_PlotData_NumClusters_by_DistanceThreshold.txt", 'w')
        distfile.write("Distance\tNumClusters\n")
    from copy import deepcopy
    thresholds = [i*threshold/NUM_THRESH for i in range(NUM_THRESH+1)]
    best = None
    best_num = -1
    best_t = -1
    distv = []
    for i, t in enumerate(thresholds):
        sg.OneLineProgressMeter('TransmissionCluster', i+1, len(thresholds)-1, 'key', 'Computing best genetic distance threshold...', orientation='h')
        clusters = method(deepcopy(tree), t, support)
        num_non_singleton = len([c for c in clusters if len(c) > 1])
        if display_fig is True:
            distfile.write("%s\t%s\n" % (t, num_non_singleton))
        distv.extend([t]*num_non_singleton)
        if num_non_singleton > best_num:
            best = clusters
            best_num = num_non_singleton
            best_t = round(t, 3)
    outfile.write("Genetic Distance Uperbound: %f\n" % threshold)
    outfile.write("Best Distance Threshold: %f\n" % best_t)

    if display_fig is True:
        distfile.close()
        bin_size = round(math.sqrt(len(distv)))
        plt.figure(2)
        plt.hist(distv, bins=bin_size)
        plt.ylabel('Number of Clusters')
        plt.xlabel('Genetic Distance Threshold')

    return best


# plot distance histogram
def gen_hist(tree, display_fig):
    if display_fig is True:
        histfile = open("TransmissionCluster_PlotData_Pairwise_Distance_Histogram.txt", 'w')
    pw_dists = []
    distance_matrix = tree.distance_matrix(leaf_labels=True)
    for u in distance_matrix.keys():
        for v in distance_matrix[u].keys():
            pw_dists.append(distance_matrix[u][v])
            if display_fig is True:
                histfile.write("%s\t%s\t%s\n" % (u, v, distance_matrix[u][v]))

    bin_size = int(math.ceil(math.sqrt(len(pw_dists)) / 10.0)) * 10
    plt.figure(1)
    plt.hist(pw_dists, bins=bin_size)
    plt.ylabel('Count')
    plt.xlabel('Sample Pairwise Genetic Distance')
    histarray = plt.hist(pw_dists, bins=bin_size)[0]
    binsarray = plt.hist(pw_dists, bins=bin_size)[1]
    if display_fig is True:
        histfile.close()
    return histarray, binsarray


# get upper limit for computing genetic distance thresholds
def get_dist_limit(hist_plot):
    histarray = hist_plot[0]
    binsarray = hist_plot[1]
    ff = histarray[:5]
    meanff = statistics.mean(ff)
    maxarray = []
    for i in range(5, len(histarray)):
        curSet = histarray[i-5:i]
        if statistics.mean(curSet) < meanff:
            maxarray.append(binsarray[i])

    d = round(float(maxarray[1]), 3)

    return d


# generate edge list to visualize clusters in gephi
def generateEdgeList(distfilename, logfilename):
    outname = "network_diagram_edge_list.txt"
    logfile = open(logfilename, 'r')
    outfile = open(outname, 'w')
    distfile = open(distfilename, 'r')

    dline = distfile.readline()
    while dline[:7] != 'Species':
        dline = distfile.readline()
    dline = distfile.readline()

    distdic = {}
    while dline[:5] != 'Table':
        data = dline.split()
        id = data[0] + '-' + data[1]
        distdic[id] = data[2]
        dline = distfile.readline()

    distfile.close()

    line = logfile.readline()
    while line != 'ClusterNum\tNumberOfSamples\tSampleNames\n':
        line = logfile.readline()

    line = logfile.readline()

    while line != '':
        clusteredSamples = line[line.index("[")+1:line.index("]")]
        clusteredSamples = clusteredSamples.split(',')
        if len(clusteredSamples) == 2:
            writeline = clusteredSamples[0] + '\t' + clusteredSamples[1] + '\n'
        else:
            for i in range(len(clusteredSamples)):
                for j in range(i+1, len(clusteredSamples)):
                    id1 = clusteredSamples[i] + '-' + clusteredSamples[j]
                    id2 = clusteredSamples[j] + '-' + clusteredSamples[i]
                    if id1 in distdic:
                        distance = distdic[id1]
                    elif id2 in distdic:
                        distance = distdic[id2]

                    writeline = clusteredSamples[i] + '\t' + clusteredSamples[j] + '\t' + distance + '\n'
                    outfile.write(writeline)

        line = logfile.readline()

    logfile.close()
    outfile.close()

    rawedgelist = open(outname, 'r')
    newoutname = outname + "_filtered"
    newout = open(newoutname, 'w')
    newout.write("Source\tTarget\n")
    line = rawedgelist.readline()
    cdic = {}
    while line != '':
        data = line.split()
        if len(data) == 2:
            newout.write(line)
        else:
            id1 = data[0]
            id2 = data[1]
            dist = float(data[2])
            infocols = [id2, dist]
            if id1 not in cdic:
                cdic[id1] = infocols
            else:
                curdata = cdic[id1]
                curdist = float(curdata[1])
                if curdist > dist:
                    cdic[id1] = infocols

        line = rawedgelist.readline()

    rawedgelist.close()

    for key in cdic.keys():
        ol = key + '\t' + cdic[key][0] + '\n'
        newout.write(ol)

    newout.close()
    os.system("rm network_diagram_edge_list.txt")
    os.system("mv network_diagram_edge_list.txt_filtered network_diagram_edgeList.txt")


if __name__ == "__main__":
    # Render GUI window
    passing = False
    window = ''
    while passing is not True:
        if window != '':
            window.Close()
        layout = [[sg.Image('resources/logo.png')],
                    [sg.Text('Newick Tree File:', font=('Helvetica', 13, 'bold')), sg.InputText(font=('Helvetica 13'), key='infilename'), sg.FileBrowse(font=('Helvetica 13'))],
                    [sg.Text('Output Filename:', font=('Helvetica', 13, 'bold')), sg.InputText(font=('Helvetica 13'), default_text='TransmissionCluster_Results.txt', text_color='gray', key='outfilename')],
                    [sg.Text('Genetic Distance Threshold (optional):', font=('Helvetica 13')), sg.InputText(font=('Helvetica 13'), key='dist'), sg.Checkbox('Compute Best Distance Threshold', font=('Helvetica 13'), default=False, key='df')],
                    [sg.Text('Support Threshold (optional):', font=('Helvetica 13')), sg.InputText(font=('Helvetica 13'), key='support')],
                    [sg.Checkbox('Plot Histograms', font=('Helvetica 13'), default=False, key='plothist')],
                    [sg.OK('Analyze', font=('Helvetica 13'))]]

        window = sg.Window('TransmissionCluster', layout)
        event, values = window.Read()

        # parse user arguments
        try:
            float(values['dist'])
            if float(values['dist']) > 1 or float(values['dist']) < 0:
                sg.Popup("Error: Genetic distance threshold must be between 0 and 1.", font=('Helvetica', 13, 'bold'))
                passing = False
            else:
                passing = True
        except ValueError:
            if values['df'] is not True:
                sg.Popup("Error: Genetic distance threshold must be between 0 and 1.", font=('Helvetica', 13, 'bold'))
                passing = False
            else:
                passing = True

        if values['support'] != '':
            try:
                float(values['support'])
                if float(values['support']) > 1 or float(values['support']) < 0:
                    sg.Popup("Error: Support threshold must be between 0 and 1.", font=('Helvetica', 13, 'bold'))
                    passing = False
                else:
                    passing = True
            except ValueError:
                sg.Popup("Error: Support threshold must be between 0 and 1.", font=('Helvetica', 13, 'bold'))
                passing = False
        else:
            passing = True

        if os.path.exists(values['infilename']) is not True:
            sg.Popup("Error: Input tree not found.", font=('Helvetica', 13, 'bold'))
            passing = False

    infile = open(values['infilename'], 'r')
    outfile = open(values['outfilename'], 'w')
    if values['support'] == '':
        values['support'] = '-inf'
    trees = list()
    for line in infile:
        if isinstance(line, bytes):
            l = line.decode().strip()
        else:
            l = line.strip()
        trees.append(read_tree_newick(l))

    # run algorithm
    outfile.write("** TransmissionCluster Results **\n")
    outfile.write("Input File: %s\n" % values['infilename'])
    outfile.write("Support Threshold: %s\n" % values['support'])
    for t, tree in enumerate(trees):
        # plot pairwise distances
        visable = False
        if values['plothist'] is True:
            visable = True
        if values['df'] is False:
            outfile.write("Genetic Distance Threshold: %s\n" % values['dist'])
            if visable is True:
                gen_hist(tree, visable)
            clusters = min_clusters_threshold_max_clade(tree, float(values['dist']), float(values['support']))
        else:
            histarray = gen_hist(tree, visable)
            d = get_dist_limit(histarray)
            clusters = argmax_clusters(min_clusters_threshold_max_clade, tree, float(d), float(values['support']), visable)
        cluster_num = 1
        clust_members = {}
        for cluster in clusters:
            if len(cluster) > 1:
                for l in cluster:
                    if cluster_num in clust_members:
                        samplenames = clust_members[cluster_num]
                        samplenames.append(l)
                        clust_members[cluster_num] = samplenames
                    else:
                        samplenames = [l]
                        clust_members[cluster_num] = samplenames
                cluster_num += 1
        totalclusters = clust_members
        cluster_num -= 1
        outfile.write('Found %s clusters\n\n' % cluster_num)
        header = "ClusterNum\tNumberOfSamples\tSampleNames\n"
        outfile.write(header)
        for k in clust_members.keys():
            outfile.write("%s\t%s\t[%s]\n" % (k, len(clust_members[k]), (','.join(clust_members[k]))))
    outfile.close()
    sg.PopupOK('Process Complete!',
        'Results have been written to the output file:\n%s' % values['outfilename'],
        'Plots will now be displayed (if option checked)...', font=('Helvetica', 13))

    if visable is True:
        plt.show()
