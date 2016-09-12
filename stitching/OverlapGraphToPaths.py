#!/usr/bin/env python
import networkx as nx
import sys
import Queue
import argparse
import Overlap

ap = argparse.ArgumentParser(description="Find a path through overlap graphs")
ap.add_argument("ovp", help="Overlap file")
ap.add_argument("graph",help="Graph file")
ap.add_argument("paths", help="Output paths")
ap.add_argument("--tgraph", help="Graph with transitive edge removal", default=None)
args = ap.parse_args()


def OrderDestByOverlap(n,g,o):
    adj = []
    for dest in g[n]:
        dist = o[g.edge[n][dest]['index']].DistanceToStart()
        if dist != -1:
            adj.append((dist, dest))
    return sorted(adj)

def RemoveTips(g,l):
    nTips = 0
    for n in g.nodes():
        if g.in_degree(n) == 0 or g.out_degree(n) == 0:
            # Found a source or sink
            src = n
            srcLen = 0
            srcOut = n
            while g.in_degree(src) == 1 and srcLen <= l:
                srcOut = src
                src = g.in_edges(src)[0]
                srcLen += 1

            dest = n
            destLen = 0
            destOut = n
            while g.out_degree(dest) == 1 and destLen + srcLen <= l:
                destOut = dest
                dest = g[dest].keys()[0]
                destLen += 1
            toRemove = []
            if srcLen + destLen <= l:
                # Found a tip, remove the path from the src node to the
                # end of the dest
                if srcOut != src:
                    toRemove.append((src, srcOut))
                    while srcOut != n:
                        src = srcOut
                        srcOut = g.in_edges(src)[0]
                        toRemove.append((src,srcOut))
                dest = n
                while g.out_degree(dest) == 1:
                    destOut = g[dest].keys()[0]
                    toRemove.append((dest, destOut))
                    dest = destOut
                nTips += len(toRemove)
                g.remove_edges_from(toRemove)
    print "removed " + str(nTips) + " tips"
    return nTips

def SoleOut(g,n):
    return g[n].keys()[0]

def PathToBranch(g, e):
    path = [e[0], e[1]]
    d = e[1]
    while g.in_degree(d) == 1 and g.out_degree(d) == 1:
        s = d
        d = SoleOut(g,d)
        path.append(d)
    return path


def GetSourceNodes(g):
    #
    # Return a list of source nodes in topological order.
    #
    nodes = []
    for n in nx.topological_sort(g):
        if g.in_degree(n) == 0:
            nodes.append(n)
    return nodes

BRANCH_OUT=0
BRANCH_IN=1
BRANCH_INOUT=2
def GetBranchingNodes(g, branch=BRANCH_OUT):
    forkNodes = []
    for n in g.nodes():
        if branch == BRANCH_OUT or branch == BRANCH_INOUT:
            if g.out_degree(n) == 2:
                forkNodes.append(n)
        if branch == BRANCH_IN or branch == BRANCH_INOUT:
            if g.in_degree(n) == 2:
                forkNodes.append(n)
    return forkNodes



def RemoveSimpleBulges(g, maxBulgeLength):
    forkNodes = GetBranchingNodes(g, BRANCH_OUT)
    for source in forkNodes:
        dests = g[source].keys()
        if len(dests) != 2:
            print "ERROR, should be on a branch"
        paths = [PathToBranch(g, (source, dest)) for dest in dests]
        if paths[0][-1] != paths[1][-1]:
            print "Path from " + source + " is not a bulge, " + str(len(paths[0])) + "\t" + str(len(paths[1]))
            continue
        print "Paths: "
        print str(paths[0])
        print str(paths[1])

def LinearPaths(g):
    srcNodes = GetSourceNodes(g)
    nx.set_node_attributes(g, 'visited', {n: False for n in g.nodes()})
    optPaths = []
    for src in srcNodes:
        branchStack = [(src, 1)]
        branchPrev = [src]
        terminalNodes = []
        terminalNodeLengths = []
        path = []
        optPath = [src]
        tmp = src

        while len(branchStack) > 0:
            (node, index) = branchStack.pop()
            print "Branching - node: {}, index: {}, cur path: {}".format(node, index, len(path))
            path = path[0:index]
            path.append(node)
            g.node[node]['visited'] = True
            # if this is a branching path, add dest to the Queue
            # process nonbranching paths.

            while g.out_degree(node) == 1:
                g.node[node]['visited'] = True
                node = SoleOut(g, node)
                path.append(node)
                if g.node[node]['visited'] == True:
                    break


            #
            # First case, the search ended at a branching node
            #
            if g.node[node]['visited'] == False and g.out_degree(node) > 1:
                for dest in g[node].keys():
                    branchStack.append((dest, len(path)))
                continue

            # Reached either the end of a path, or a branch.  Handle both separately
            if g.out_degree(node) == 0:
                if len(path) > len(optPath):
                    optPath = path
        print "source: " + src + " has path of length " + str(len(optPath))
        optPaths.append(optPath)
    return optPaths




def RemoveTransitiveEdges(g, o):
    #
    # perform simple a->b, b->c, a-/->c transitive edge removal.
    # g - digraph
    # o - overlaps, indexed by g[node]['index']
    #
    transitive = []
    for n in g.nodes():
        inPlay = { dest for dest in g[n].keys() }
        adj = OrderDestByOverlap(n, g, o)
        for (dist, neighbor) in adj:
            for second in g[neighbor].keys():
                if second in inPlay:
                    transitive.append((n,second))
    g.remove_edges_from(transitive)
    print "removed " + str(len(transitive)) + " edges"

g = nx.read_gml(args.graph)
overlapFile = open(args.ovp)
pathFile = open(args.paths,'w')
overlaps = [Overlap.Overlap(line) for line in overlapFile]

RemoveTransitiveEdges(g, overlaps)
RemoveTips(g,1)
paths = LinearPaths(g)
for p in paths:
    pathFile.write("\t".join(p) + "\n")

if args.tgraph is not None:
    nx.write_gml(g, args.tgraph)
