#!/usr/bin/env python3
from queue import Queue
from treeswift import read_tree_newick
import PySimpleGUI as sg
import os
NUM_THRESH = 1000  # number of thresholds for the threshold-free methods to use


# cut out the current node's subtree (by setting all nodes' DELETED to True) and return list of leaves
def cut(node):
    cluster = list()
    descendants = Queue(); descendants.put(node)
    while not descendants.empty():
        descendant = descendants.get()
        if descendant.DELETED:
            continue
        descendant.DELETED = True
        descendant.left_dist = 0; descendant.right_dist = 0; descendant.edge_length = 0
        if descendant.is_leaf():
            cluster.append(str(descendant))
        else:
            for c in descendant.children:
                descendants.put(c)
    return cluster


# initialize properties of input tree and return set containing taxa of leaves
def prep(tree, support):
    tree.resolve_polytomies(); tree.suppress_unifurcations()
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
            node.left_dist = 0; node.right_dist = 0
        else:
            children = list(node.children)
            if children[0].DELETED and children[1].DELETED:
                cut(node); continue
            if children[0].DELETED:
                node.left_dist = 0
            else:
                node.left_dist = max(children[0].left_dist,children[0].right_dist) + children[0].edge_length
            if children[1].DELETED:
                node.right_dist = 0
            else:
                node.right_dist = max(children[1].left_dist,children[1].right_dist) + children[1].edge_length

            # if my kids are screwing things up, cut both
            if node.left_dist + node.right_dist > threshold:
                cluster_l = cut(children[0])
                node.left_dist = 0
                cluster_r = cut(children[1])
                node.right_dist = 0

                # add cluster
                for cluster in (cluster_l,cluster_r):
                    if len(cluster) != 0:
                        clusters.append(cluster)
                        for leaf in cluster:
                            leaves.remove(leaf)

    # add all remaining leaves to a single cluster
    if len(leaves) != 0:
        clusters.append(list(leaves))
    return clusters


# pick the threshold between 0 and "distance threshold" that maximizes number of (non-singleton) clusters
def argmax_clusters(method, tree, threshold, support):
    distfile = open("distance_permutations.txt", 'w')
    distfile.write("Distance\tNumClusters\n")
    from copy import deepcopy
    assert threshold > 0, "Threshold must be positive"
    thresholds = [i*threshold/NUM_THRESH for i in range(NUM_THRESH+1)]
    best = None; best_num = -1; best_t = -1
    for i, t in enumerate(thresholds):
        sg.OneLineProgressMeter('TransmissionCluster', i+1, len(thresholds)-1, 'key', 'Computing best genetic distance threshold...', orientation='h')
        clusters = method(deepcopy(tree), t, support)
        num_non_singleton = len([c for c in clusters if len(c) > 1])
        distfile.write("%s\t%s\n" % (t, num_non_singleton))
        if num_non_singleton > best_num:
            best = clusters; best_num = num_non_singleton; best_t = t
    outfile.write("Best Distance Threshold: %f\n" % best_t)
    distfile.close()
    return best


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
        distdic[id]=data[2]
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

    rawedgelist = open(outname,'r')
    newoutname = outname + "_filtered"
    newout = open(newoutname,'w')
    newout.write("Source\tTarget\n")
    line = rawedgelist.readline()
    cdic = {}
    while line != '':
        data = line.split()
        if len(data)==2:
            newout.write(line)
        else:
            id1 = data[0]
            id2 = data[1]
            dist = float(data[2])
            infocols = [id2,dist]
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
    #GUI
    layout = [ [sg.Image('resources/logo.png')],
                [sg.Text('Newick Tree File:', font=('Helvetica 12')), sg.InputText(key='infilename'), sg.FileBrowse(font=('Helvetica 12'))],
                [sg.Text('Output Filename:', font=('Helvetica 12')), sg.InputText(key='outfilename')],
                [sg.Text('Genetic Distance Threshold:', font=('Helvetica 12')), sg.InputText(key='dist'), sg.Checkbox('Compute Best Distance Threshold', font=('Helvetica 12'), default=False,key='df')],
                [sg.Text('Support Threshold (optional):', font=('Helvetica 12')), sg.InputText(key='support')],
                [sg.Text('Pairwise distance file (optional, see documentation):', font=('Helvetica 12')), sg.InputText(key = 'distfile'), sg.FileBrowse(font=('Helvetica 12'))],
                [sg.OK(font=('Helvetica 12')), sg.Cancel('Quit',font=('Helvetica 12'))] ]

    window = sg.Window('TransmissionCluster', layout)
    event, values = window.Read()

    # parse user arguments
    # assert float(values['dist']) >= 0, "ERROR: Genertic distance threshold must be at least 0"
    # assert float(values['dist']) <= 1, "ERROR: Genertic distance threshold must be less than 1"
    # assert float(values['support']) >= 0 or float(values['support']) == float('-inf'), "ERROR: Branch support threshold must be at least 0"
    # assert float(values['support']) <= 1, "ERROR: Branch support threshold must be less than 1"

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
        if values['df'] is False:
            clusters = min_clusters_threshold_max_clade(tree, float(values['dist']), float(values['support']))
        else:
            clusters = argmax_clusters(min_clusters_threshold_max_clade, tree, float(values['dist']), float(values['support']))
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

    if values['distfile'] != '':
        generateEdgeList(values['distfile'], values['outfilename'])

    sg.PopupOK('Process Complete!\nResults have been written to the output file.')
