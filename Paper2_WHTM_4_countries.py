# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:28:14 2019

@author: Dr. Utkarsh Verma
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 18:17:32 2019

@author: Dr. Utkarsh Verma
"""

import random
from itertools import combinations
import math
import bisect
import sys
import matplotlib.pyplot as plt
import numpy as np



def power_set(List):
    PS = [list(j) for i in range(len(List)) for j in combinations(List, i+1)]
    return PS

rng1 = np.random.RandomState(10)
rng2 = np.random.RandomState(11)
rng3 = np.random.RandomState(12)
# If assigned, DEBUG(str) is called with lots of debug messages.
DEBUG = None
"""def DEBUG(s):
    from sys import stderr
    print('DEBUG:', s, file=stderr)
"""

# Check delta2/delta3 computation after every substage;
# only works on integer weights, slows down the algorithm to O(n^4).
CHECK_DELTA = False

# Check optimality of solution before returning; only works on integer weights.
CHECK_OPTIMUM = True

import rpy2.robjects as robjects

from rpy2.robjects.vectors import StrVector
from pandas import *

#base=importr('base')
#game=importr('gametheory')

robjects.r('x=c()')
robjects.r('x[1]=22')
print (robjects.r('x'))

from rpy2.robjects.packages import importr
utils = importr('utils')
utils.install_packages('GameTheory')
game=importr("GameTheory")

print (game)

definegame=robjects.r['DefineGame']
shapley_value=robjects.r['ShapleyValue']
summary=robjects.r['summary']
Nucleolus=robjects.r['Nucleolus']
dataframe=robjects.r['data.frame']
unlist=robjects.r['unlist']





def maxWeightMatching(edges, maxcardinality=False):
    """Compute a maximum-weighted matching in the general undirected
    weighted graph given by "edges".  If "maxcardinality" is true,
    only maximum-cardinality matchings are considered as solutions.

    Edges is a sequence of tuples (i, j, wt) describing an undirected
    edge between vertex i and vertex j with weight wt.  There is at most
    one edge between any two vertices; no vertex has an edge to itself.
    Vertices are identified by consecutive, non-negative integers.

    Return a list "mate", such that mate[i] == j if vertex i is
    matched to vertex j, and mate[i] == -1 if vertex i is not matched.

    This function takes time O(n ** 3)."""

    #
    # Vertices are numbered 0 .. (nvertex-1).
    # Non-trivial blossoms are numbered nvertex .. (2*nvertex-1)
    #
    # Edges are numbered 0 .. (nedge-1).
    # Edge endpoints are numbered 0 .. (2*nedge-1), such that endpoints
    # (2*k) and (2*k+1) both belong to edge k.
    #
    # Many terms used in the comments (sub-blossom, T-vertex) come from
    # the paper by Galil; read the paper before reading this code.
    #

    # Python 2/3 compatibility.
    from sys import version as sys_version
    if sys_version < '3':
        integer_types = (int, long)
    else:
        integer_types = (int,)

    # Deal swiftly with empty graphs.
    if not edges:
        return [ ]

    # Count vertices.
    nedge = len(edges)
    nvertex = 0
    for (i, j, w) in edges:
        assert i >= 0 and j >= 0 #and i != j
        if i >= nvertex:
            nvertex = i + 1
        if j >= nvertex:
            nvertex = j + 1

    # Find the maximum edge weight.
    maxweight = max(0, max([ wt for (i, j, wt) in edges ]))

    # If p is an edge endpoint,
    # endpoint[p] is the vertex to which endpoint p is attached.
    # Not modified by the algorithm.
    endpoint = [ edges[p//2][p%2] for p in range(2*nedge) ]

    # If v is a vertex,
    # neighbend[v] is the list of remote endpoints of the edges attached to v.
    # Not modified by the algorithm.
    neighbend = [ [ ] for i in range(nvertex) ]
    for k in range(len(edges)):
        (i, j, w) = edges[k]
        neighbend[i].append(2*k+1)
        neighbend[j].append(2*k)

    # If v is a vertex,
    # mate[v] is the remote endpoint of its matched edge, or -1 if it is single
    # (i.e. endpoint[mate[v]] is v's partner vertex).
    # Initially all vertices are single; updated during augmentation.
    mate = nvertex * [ -1 ]

    # If b is a top-level blossom,
    # label[b] is 0 if b is unlabeled (free);
    #             1 if b is an S-vertex/blossom;
    #             2 if b is a T-vertex/blossom.
    # The label of a vertex is found by looking at the label of its
    # top-level containing blossom.
    # If v is a vertex inside a T-blossom,
    # label[v] is 2 iff v is reachable from an S-vertex outside the blossom.
    # Labels are assigned during a stage and reset after each augmentation.
    label = (2 * nvertex) * [ 0 ]

    # If b is a labeled top-level blossom,
    # labelend[b] is the remote endpoint of the edge through which b obtained
    # its label, or -1 if b's base vertex is single.
    # If v is a vertex inside a T-blossom and label[v] == 2,
    # labelend[v] is the remote endpoint of the edge through which v is
    # reachable from outside the blossom.
    labelend = (2 * nvertex) * [ -1 ]

    # If v is a vertex,
    # inblossom[v] is the top-level blossom to which v belongs.
    # If v is a top-level vertex, v is itself a blossom (a trivial blossom)
    # and inblossom[v] == v.
    # Initially all vertices are top-level trivial blossoms.
    inblossom = list(range(nvertex))

    # If b is a sub-blossom,
    # blossomparent[b] is its immediate parent (sub-)blossom.
    # If b is a top-level blossom, blossomparent[b] is -1.
    blossomparent = (2 * nvertex) * [ -1 ]

    # If b is a non-trivial (sub-)blossom,
    # blossomchilds[b] is an ordered list of its sub-blossoms, starting with
    # the base and going round the blossom.
    blossomchilds = (2 * nvertex) * [ None ]

    # If b is a (sub-)blossom,
    # blossombase[b] is its base VERTEX (i.e. recursive sub-blossom).
    blossombase = list(range(nvertex)) + nvertex * [ -1 ]

    # If b is a non-trivial (sub-)blossom,
    # blossomendps[b] is a list of endpoints on its connecting edges,
    # such that blossomendps[b][i] is the local endpoint of blossomchilds[b][i]
    # on the edge that connects it to blossomchilds[b][wrap(i+1)].
    blossomendps = (2 * nvertex) * [ None ]

    # If v is a free vertex (or an unreached vertex inside a T-blossom),
    # bestedge[v] is the edge to an S-vertex with least slack,
    # or -1 if there is no such edge.
    # If b is a (possibly trivial) top-level S-blossom,
    # bestedge[b] is the least-slack edge to a different S-blossom,
    # or -1 if there is no such edge.
    # This is used for efficient computation of delta2 and delta3.
    bestedge = (2 * nvertex) * [ -1 ]

    # If b is a non-trivial top-level S-blossom,
    # blossombestedges[b] is a list of least-slack edges to neighbouring
    # S-blossoms, or None if no such list has been computed yet.
    # This is used for efficient computation of delta3.
    blossombestedges = (2 * nvertex) * [ None ]

    # List of currently unused blossom numbers.
    unusedblossoms = list(range(nvertex, 2*nvertex))

    # If v is a vertex,
    # dualvar[v] = 2 * u(v) where u(v) is the v's variable in the dual
    # optimization problem (multiplication by two ensures integer values
    # throughout the algorithm if all edge weights are integers).
    # If b is a non-trivial blossom,
    # dualvar[b] = z(b) where z(b) is b's variable in the dual optimization
    # problem.
    dualvar = nvertex * [ maxweight ] + nvertex * [ 0 ]

    # If allowedge[k] is true, edge k has zero slack in the optimization
    # problem; if allowedge[k] is false, the edge's slack may or may not
    # be zero.
    allowedge = nedge * [ False ]

    # Queue of newly discovered S-vertices.
    queue = [ ]

    # Return 2 * slack of edge k (does not work inside blossoms).
    def slack(k):
        (i, j, wt) = edges[k]
        return dualvar[i] + dualvar[j] - 2 * wt

    # Generate the leaf vertices of a blossom.
    def blossomLeaves(b):
        if b < nvertex:
            yield b
        else:
            for t in blossomchilds[b]:
                if t < nvertex:
                    yield t
                else:
                    for v in blossomLeaves(t):
                        yield v

    # Assign label t to the top-level blossom containing vertex w
    # and record the fact that w was reached through the edge with
    # remote endpoint p.
    def assignLabel(w, t, p):
        if DEBUG: DEBUG('assignLabel(%d,%d,%d)' % (w, t, p))
        b = inblossom[w]
        assert label[w] == 0 and label[b] == 0
        label[w] = label[b] = t
        labelend[w] = labelend[b] = p
        bestedge[w] = bestedge[b] = -1
        if t == 1:
            # b became an S-vertex/blossom; add it(s vertices) to the queue.
            queue.extend(blossomLeaves(b))
            if DEBUG: DEBUG('PUSH ' + str(list(blossomLeaves(b))))
        elif t == 2:
            # b became a T-vertex/blossom; assign label S to its mate.
            # (If b is a non-trivial blossom, its base is the only vertex
            # with an external mate.)
            base = blossombase[b]
            assert mate[base] >= 0
            assignLabel(endpoint[mate[base]], 1, mate[base] ^ 1)

    # Trace back from vertices v and w to discover either a new blossom
    # or an augmenting path. Return the base vertex of the new blossom or -1.
    def scanBlossom(v, w):
        if DEBUG: DEBUG('scanBlossom(%d,%d)' % (v, w))
        # Trace back from v and w, placing breadcrumbs as we go.
        path = [ ]
        base = -1
        while v != -1 or w != -1:
            # Look for a breadcrumb in v's blossom or put a new breadcrumb.
            b = inblossom[v]
            if label[b] & 4:
                base = blossombase[b]
                break
            assert label[b] == 1
            path.append(b)
            label[b] = 5
            # Trace one step back.
            assert labelend[b] == mate[blossombase[b]]
            if labelend[b] == -1:
                # The base of blossom b is single; stop tracing this path.
                v = -1
            else:
                v = endpoint[labelend[b]]
                b = inblossom[v]
                assert label[b] == 2
                # b is a T-blossom; trace one more step back.
                assert labelend[b] >= 0
                v = endpoint[labelend[b]]
            # Swap v and w so that we alternate between both paths.
            if w != -1:
                v, w = w, v
       # Remove breadcrumbs.
        for b in path:
            label[b] = 1
        # Return base vertex, if we found one.
        return base

    # Construct a new blossom with given base, containing edge k which
    # connects a pair of S vertices. Label the new blossom as S; set its dual
    # variable to zero; relabel its T-vertices to S and add them to the queue.
    def addBlossom(base, k):
        (v, w, wt) = edges[k]
        bb = inblossom[base]
        bv = inblossom[v]
        bw = inblossom[w]
        # Create blossom.
        b = unusedblossoms.pop()
        if DEBUG: DEBUG('addBlossom(%d,%d) (v=%d w=%d) -> %d' % (base, k, v, w, b))
        blossombase[b] = base
        blossomparent[b] = -1
        blossomparent[bb] = b
        # Make list of sub-blossoms and their interconnecting edge endpoints.
        blossomchilds[b] = path = [ ]
        blossomendps[b] = endps = [ ]
        # Trace back from v to base.
        while bv != bb:
            # Add bv to the new blossom.
            blossomparent[bv] = b
            path.append(bv)
            endps.append(labelend[bv])
            assert (label[bv] == 2 or
                    (label[bv] == 1 and labelend[bv] == mate[blossombase[bv]]))
            # Trace one step back.
            assert labelend[bv] >= 0
            v = endpoint[labelend[bv]]
            bv = inblossom[v]
        # Reverse lists, add endpoint that connects the pair of S vertices.
        path.append(bb)
        path.reverse()
        endps.reverse()
        endps.append(2*k)
        # Trace back from w to base.
        while bw != bb:
            # Add bw to the new blossom.
            blossomparent[bw] = b
            path.append(bw)
            endps.append(labelend[bw] ^ 1)
            assert (label[bw] == 2 or
                    (label[bw] == 1 and labelend[bw] == mate[blossombase[bw]]))
            # Trace one step back.
            assert labelend[bw] >= 0
            w = endpoint[labelend[bw]]
            bw = inblossom[w]
        # Set label to S.
        assert label[bb] == 1
        label[b] = 1
        labelend[b] = labelend[bb]
        # Set dual variable to zero.
        dualvar[b] = 0
        # Relabel vertices.
        for v in blossomLeaves(b):
            if label[inblossom[v]] == 2:
                # This T-vertex now turns into an S-vertex because it becomes
                # part of an S-blossom; add it to the queue.
                queue.append(v)
            inblossom[v] = b
        # Compute blossombestedges[b].
        bestedgeto = (2 * nvertex) * [ -1 ]
        for bv in path:
            if blossombestedges[bv] is None:
                # This subblossom does not have a list of least-slack edges;
                # get the information from the vertices.
                nblists = [ [ p // 2 for p in neighbend[v] ]
                            for v in blossomLeaves(bv) ]
            else:
                # Walk this subblossom's least-slack edges.
                nblists = [ blossombestedges[bv] ]
            for nblist in nblists:
                for k in nblist:
                    (i, j, wt) = edges[k]
                    if inblossom[j] == b:
                        i, j = j, i
                    bj = inblossom[j]
                    if (bj != b and label[bj] == 1 and
                        (bestedgeto[bj] == -1 or
                         slack(k) < slack(bestedgeto[bj]))):
                        bestedgeto[bj] = k
            # Forget about least-slack edges of the subblossom.
            blossombestedges[bv] = None
            bestedge[bv] = -1
        blossombestedges[b] = [ k for k in bestedgeto if k != -1 ]
        # Select bestedge[b].
        bestedge[b] = -1
        for k in blossombestedges[b]:
            if bestedge[b] == -1 or slack(k) < slack(bestedge[b]):
                bestedge[b] = k
        if DEBUG: DEBUG('blossomchilds[%d]=' % b + repr(blossomchilds[b]))

    # Expand the given top-level blossom.
    def expandBlossom(b, endstage):
        if DEBUG: DEBUG('expandBlossom(%d,%d) %s' % (b, endstage, repr(blossomchilds[b])))
        # Convert sub-blossoms into top-level blossoms.
        for s in blossomchilds[b]:
            blossomparent[s] = -1
            if s < nvertex:
                inblossom[s] = s
            elif endstage and dualvar[s] == 0:
                # Recursively expand this sub-blossom.
                expandBlossom(s, endstage)
            else:
                for v in blossomLeaves(s):
                    inblossom[v] = s
        # If we expand a T-blossom during a stage, its sub-blossoms must be
        # relabeled.
        if (not endstage) and label[b] == 2:
            # Start at the sub-blossom through which the expanding
            # blossom obtained its label, and relabel sub-blossoms untili
            # we reach the base.
            # Figure out through which sub-blossom the expanding blossom
            # obtained its label initially.
            assert labelend[b] >= 0
            entrychild = inblossom[endpoint[labelend[b] ^ 1]]
            # Decide in which direction we will go round the blossom.
            j = blossomchilds[b].index(entrychild)
            if j & 1:
                # Start index is odd; go forward and wrap.
                j -= len(blossomchilds[b])
                jstep = 1
                endptrick = 0
            else:
                # Start index is even; go backward.
                jstep = -1
                endptrick = 1
            # Move along the blossom until we get to the base.
            p = labelend[b]
            while j != 0:
                # Relabel the T-sub-blossom.
                label[endpoint[p ^ 1]] = 0
                label[endpoint[blossomendps[b][j-endptrick]^endptrick^1]] = 0
                assignLabel(endpoint[p ^ 1], 2, p)
                # Step to the next S-sub-blossom and note its forward endpoint.
                allowedge[blossomendps[b][j-endptrick]//2] = True
                j += jstep
                p = blossomendps[b][j-endptrick] ^ endptrick
                # Step to the next T-sub-blossom.
                allowedge[p//2] = True
                j += jstep
            # Relabel the base T-sub-blossom WITHOUT stepping through to
            # its mate (so don't call assignLabel).
            bv = blossomchilds[b][j]
            label[endpoint[p ^ 1]] = label[bv] = 2
            labelend[endpoint[p ^ 1]] = labelend[bv] = p
            bestedge[bv] = -1
            # Continue along the blossom until we get back to entrychild.
            j += jstep
            while blossomchilds[b][j] != entrychild:
                # Examine the vertices of the sub-blossom to see whether
                # it is reachable from a neighbouring S-vertex outside the
                # expanding blossom.
                bv = blossomchilds[b][j]
                if label[bv] == 1:
                    # This sub-blossom just got label S through one of its
                    # neighbours; leave it.
                    j += jstep
                    continue
                for v in blossomLeaves(bv):
                    if label[v] != 0:
                        break
                # If the sub-blossom contains a reachable vertex, assign
                # label T to the sub-blossom.
                if label[v] != 0:
                    assert label[v] == 2
                    assert inblossom[v] == bv
                    label[v] = 0
                    label[endpoint[mate[blossombase[bv]]]] = 0
                    assignLabel(v, 2, labelend[v])
                j += jstep
        # Recycle the blossom number.
        label[b] = labelend[b] = -1
        blossomchilds[b] = blossomendps[b] = None
        blossombase[b] = -1
        blossombestedges[b] = None
        bestedge[b] = -1
        unusedblossoms.append(b)

    # Swap matched/unmatched edges over an alternating path through blossom b
    # between vertex v and the base vertex. Keep blossom bookkeeping consistent.
    def augmentBlossom(b, v):
        if DEBUG: DEBUG('augmentBlossom(%d,%d)' % (b, v))
        # Bubble up through the blossom tree from vertex v to an immediate
        # sub-blossom of b.
        t = v
        while blossomparent[t] != b:
            t = blossomparent[t]
        # Recursively deal with the first sub-blossom.
        if t >= nvertex:
            augmentBlossom(t, v)
        # Decide in which direction we will go round the blossom.
        i = j = blossomchilds[b].index(t)
        if i & 1:
            # Start index is odd; go forward and wrap.
            j -= len(blossomchilds[b])
            jstep = 1
            endptrick = 0
        else:
            # Start index is even; go backward.
            jstep = -1
            endptrick = 1
        # Move along the blossom until we get to the base.
        while j != 0:
            # Step to the next sub-blossom and augment it recursively.
            j += jstep
            t = blossomchilds[b][j]
            p = blossomendps[b][j-endptrick] ^ endptrick
            if t >= nvertex:
                augmentBlossom(t, endpoint[p])
            # Step to the next sub-blossom and augment it recursively.
            j += jstep
            t = blossomchilds[b][j]
            if t >= nvertex:
                augmentBlossom(t, endpoint[p ^ 1])
            # Match the edge connecting those sub-blossoms.
            mate[endpoint[p]] = p ^ 1
            mate[endpoint[p ^ 1]] = p
            if DEBUG: DEBUG('PAIR %d %d (k=%d)' % (endpoint[p], endpoint[p^1], p//2))
        # Rotate the list of sub-blossoms to put the new base at the front.
        blossomchilds[b] = blossomchilds[b][i:] + blossomchilds[b][:i]
        blossomendps[b]  = blossomendps[b][i:]  + blossomendps[b][:i]
        blossombase[b] = blossombase[blossomchilds[b][0]]
        assert blossombase[b] == v

    # Swap matched/unmatched edges over an alternating path between two
    # single vertices. The augmenting path runs through edge k, which
    # connects a pair of S vertices.
    def augmentMatching(k):
        (v, w, wt) = edges[k]
        if DEBUG: DEBUG('augmentMatching(%d) (v=%d w=%d)' % (k, v, w))
        if DEBUG: DEBUG('PAIR %d %d (k=%d)' % (v, w, k))
        for (s, p) in ((v, 2*k+1), (w, 2*k)):
            # Match vertex s to remote endpoint p. Then trace back from s
            # until we find a single vertex, swapping matched and unmatched
            # edges as we go.
            while 1:
                bs = inblossom[s]
                assert label[bs] == 1
                assert labelend[bs] == mate[blossombase[bs]]
                # Augment through the S-blossom from s to base.
                if bs >= nvertex:
                    augmentBlossom(bs, s)
                # Update mate[s]
                mate[s] = p
                # Trace one step back.
                if labelend[bs] == -1:
                    # Reached single vertex; stop.
                    break
                t = endpoint[labelend[bs]]
                bt = inblossom[t]
                assert label[bt] == 2
                # Trace one step back.
                assert labelend[bt] >= 0
                s = endpoint[labelend[bt]]
                j = endpoint[labelend[bt] ^ 1]
                # Augment through the T-blossom from j to base.
                assert blossombase[bt] == t
                if bt >= nvertex:
                    augmentBlossom(bt, j)
                # Update mate[j]
                mate[j] = labelend[bt]
                # Keep the opposite endpoint;
                # it will be assigned to mate[s] in the next step.
                p = labelend[bt] ^ 1
                if DEBUG: DEBUG('PAIR %d %d (k=%d)' % (s, t, p//2))

    # Verify that the optimum solution has been reached.
    def verifyOptimum():
        if maxcardinality:
            # Vertices may have negative dual;
            # find a constant non-negative number to add to all vertex duals.
            vdualoffset = max(0, -min(dualvar[:nvertex]))
        else:
            vdualoffset = 0
        # 0. all dual variables are non-negative
        assert min(dualvar[:nvertex]) + vdualoffset >= 0
        assert min(dualvar[nvertex:]) >= 0
        # 0. all edges have non-negative slack and
        # 1. all matched edges have zero slack;
        for k in range(nedge):
            (i, j, wt) = edges[k]
            s = dualvar[i] + dualvar[j] - 2 * wt
            iblossoms = [ i ]
            jblossoms = [ j ]
            while blossomparent[iblossoms[-1]] != -1:
                iblossoms.append(blossomparent[iblossoms[-1]])
            while blossomparent[jblossoms[-1]] != -1:
                jblossoms.append(blossomparent[jblossoms[-1]])
            iblossoms.reverse()
            jblossoms.reverse()
            for (bi, bj) in zip(iblossoms, jblossoms):
                if bi != bj:
                    break
                s += 2 * dualvar[bi]
            #assert s >= 0
            if mate[i] // 2 == k or mate[j] // 2 == k:
                assert mate[i] // 2 == k and mate[j] // 2 == k
                assert s == 0
        # 2. all single vertices have zero dual value;
        for v in range(nvertex):
            assert mate[v] >= 0 or dualvar[v] + vdualoffset == 0
        # 3. all blossoms with positive dual value are full.
        for b in range(nvertex, 2*nvertex):
            if blossombase[b] >= 0 and dualvar[b] > 0:
                assert len(blossomendps[b]) % 2 == 1
                for p in blossomendps[b][1::2]:
                    assert mate[endpoint[p]] == p ^ 1
                    assert mate[endpoint[p ^ 1]] == p
        # Ok.

    # Check optimized delta2 against a trivial computation.
    def checkDelta2():
        for v in range(nvertex):
            if label[inblossom[v]] == 0:
                bd = None
                bk = -1
                for p in neighbend[v]:
                    k = p // 2
                    w = endpoint[p]
                    if label[inblossom[w]] == 1:
                        d = slack(k)
                        if bk == -1 or d < bd:
                            bk = k
                            bd = d
                if DEBUG and (bestedge[v] != -1 or bk != -1) and (bestedge[v] == -1 or bd != slack(bestedge[v])):
                    DEBUG('v=' + str(v) + ' bk=' + str(bk) + ' bd=' + str(bd) + ' bestedge=' + str(bestedge[v]) + ' slack=' + str(slack(bestedge[v])))
                assert (bk == -1 and bestedge[v] == -1) or (bestedge[v] != -1 and bd == slack(bestedge[v]))

    # Check optimized delta3 against a trivial computation.
    def checkDelta3():
        bk = -1
        bd = None
        tbk = -1
        tbd = None
        for b in range(2 * nvertex):
            if blossomparent[b] == -1 and label[b] == 1:
                for v in blossomLeaves(b):
                    for p in neighbend[v]:
                        k = p // 2
                        w = endpoint[p]
                        if inblossom[w] != b and label[inblossom[w]] == 1:
                            d = slack(k)
                            if bk == -1 or d < bd:
                                bk = k
                                bd = d
                if bestedge[b] != -1:
                    (i, j, wt) = edges[bestedge[b]]
                    assert inblossom[i] == b or inblossom[j] == b
                    assert inblossom[i] != b or inblossom[j] != b
                    assert label[inblossom[i]] == 1 and label[inblossom[j]] == 1
                    if tbk == -1 or slack(bestedge[b]) < tbd:
                        tbk = bestedge[b]
                        tbd = slack(bestedge[b])
        if DEBUG and bd != tbd:
            DEBUG('bk=%d tbk=%d bd=%s tbd=%s' % (bk, tbk, repr(bd), repr(tbd)))
        assert bd == tbd

    # Main loop: continue until no further improvement is possible.
    for t in range(nvertex):

        # Each iteration of this loop is a "stage".
        # A stage finds an augmenting path and uses that to improve
        # the matching.
        if DEBUG: DEBUG('STAGE %d' % t)

        # Remove labels from top-level blossoms/vertices.
        label[:] = (2 * nvertex) * [ 0 ]

        # Forget all about least-slack edges.
        bestedge[:] = (2 * nvertex) * [ -1 ]
        blossombestedges[nvertex:] = nvertex * [ None ]

        # Loss of labeling means that we can not be sure that currently
        # allowable edges remain allowable througout this stage.
        allowedge[:] = nedge * [ False ]

        # Make queue empty.
        queue[:] = [ ]
 
        # Label single blossoms/vertices with S and put them in the queue.
        for v in range(nvertex):
            if mate[v] == -1 and label[inblossom[v]] == 0:
                assignLabel(v, 1, -1)

        # Loop until we succeed in augmenting the matching.
        augmented = 0
        while 1:

            # Each iteration of this loop is a "substage".
            # A substage tries to find an augmenting path;
            # if found, the path is used to improve the matching and
            # the stage ends. If there is no augmenting path, the
            # primal-dual method is used to pump some slack out of
            # the dual variables.
            if DEBUG: DEBUG('SUBSTAGE')

            # Continue labeling until all vertices which are reachable
            # through an alternating path have got a label.
            while queue and not augmented:

                # Take an S vertex from the queue.
                v = queue.pop()
                if DEBUG: DEBUG('POP v=%d' % v)
                assert label[inblossom[v]] == 1

                # Scan its neighbours:
                for p in neighbend[v]:
                    k = p // 2
                    w = endpoint[p]
                    # w is a neighbour to v
                    if inblossom[v] == inblossom[w]:
                        # this edge is internal to a blossom; ignore it
                        continue
                    if not allowedge[k]:
                        kslack = slack(k)
                        if kslack <= 0:
                            # edge k has zero slack => it is allowable
                            allowedge[k] = True
                    if allowedge[k]:
                        if label[inblossom[w]] == 0:
                            # (C1) w is a free vertex;
                            # label w with T and label its mate with S (R12).
                            assignLabel(w, 2, p ^ 1)
                        elif label[inblossom[w]] == 1:
                            # (C2) w is an S-vertex (not in the same blossom);
                            # follow back-links to discover either an
                            # augmenting path or a new blossom.
                            base = scanBlossom(v, w)
                            if base >= 0:
                                # Found a new blossom; add it to the blossom
                                # bookkeeping and turn it into an S-blossom.
                                addBlossom(base, k)
                            else:
                                # Found an augmenting path; augment the
                                # matching and end this stage.
                                augmentMatching(k)
                                augmented = 1
                                break
                        elif label[w] == 0:
                            # w is inside a T-blossom, but w itself has not
                            # yet been reached from outside the blossom;
                            # mark it as reached (we need this to relabel
                            # during T-blossom expansion).
                            assert label[inblossom[w]] == 2
                            label[w] = 2
                            labelend[w] = p ^ 1
                    elif label[inblossom[w]] == 1:
                        # keep track of the least-slack non-allowable edge to
                        # a different S-blossom.
                        b = inblossom[v]
                        if bestedge[b] == -1 or kslack < slack(bestedge[b]):
                            bestedge[b] = k
                    elif label[w] == 0:
                        # w is a free vertex (or an unreached vertex inside
                        # a T-blossom) but we can not reach it yet;
                        # keep track of the least-slack edge that reaches w.
                        if bestedge[w] == -1 or kslack < slack(bestedge[w]):
                            bestedge[w] = k

            if augmented:
                break

            # There is no augmenting path under these constraints;
            # compute delta and reduce slack in the optimization problem.
            # (Note that our vertex dual variables, edge slacks and delta's
            # are pre-multiplied by two.)
            deltatype = -1
            delta = deltaedge = deltablossom = None

            # Verify data structures for delta2/delta3 computation.
            if CHECK_DELTA:
                checkDelta2()
                checkDelta3()

            # Compute delta1: the minumum value of any vertex dual.
            if not maxcardinality:
                deltatype = 1
                delta = min(dualvar[:nvertex])

            # Compute delta2: the minimum slack on any edge between
            # an S-vertex and a free vertex.
            for v in range(nvertex):
                if label[inblossom[v]] == 0 and bestedge[v] != -1:
                    d = slack(bestedge[v])
                    if deltatype == -1 or d < delta:
                        delta = d
                        deltatype = 2
                        deltaedge = bestedge[v]

            # Compute delta3: half the minimum slack on any edge between
            # a pair of S-blossoms.
            for b in range(2 * nvertex):
                if ( blossomparent[b] == -1 and label[b] == 1 and
                     bestedge[b] != -1 ):
                    kslack = slack(bestedge[b])
                    if isinstance(kslack, integer_types):
                        assert (kslack % 2) == 0
                        d = kslack // 2
                    else:
                        d = kslack / 2
                    if deltatype == -1 or d < delta:
                        delta = d
                        deltatype = 3
                        deltaedge = bestedge[b]

            # Compute delta4: minimum z variable of any T-blossom.
            for b in range(nvertex, 2*nvertex):
                if ( blossombase[b] >= 0 and blossomparent[b] == -1 and
                     label[b] == 2 and
                     (deltatype == -1 or dualvar[b] < delta) ):
                    delta = dualvar[b]
                    deltatype = 4
                    deltablossom = b

            if deltatype == -1:
                # No further improvement possible; max-cardinality optimum
                # reached. Do a final delta update to make the optimum
                # verifyable.
                assert maxcardinality
                deltatype = 1
                delta = max(0, min(dualvar[:nvertex]))

            # Update dual variables according to delta.
            for v in range(nvertex):
                if label[inblossom[v]] == 1:
                    # S-vertex: 2*u = 2*u - 2*delta
                    dualvar[v] -= delta
                elif label[inblossom[v]] == 2:
                    # T-vertex: 2*u = 2*u + 2*delta
                    dualvar[v] += delta
            for b in range(nvertex, 2*nvertex):
                if blossombase[b] >= 0 and blossomparent[b] == -1:
                    if label[b] == 1:
                        # top-level S-blossom: z = z + 2*delta
                        dualvar[b] += delta
                    elif label[b] == 2:
                        # top-level T-blossom: z = z - 2*delta
                        dualvar[b] -= delta

            # Take action at the point where minimum delta occurred.
            if DEBUG: DEBUG('delta%d=%f' % (deltatype, delta))
            if deltatype == 1: 
                # No further improvement possible; optimum reached.
                break
            elif deltatype == 2:
                # Use the least-slack edge to continue the search.
                allowedge[deltaedge] = True
                (i, j, wt) = edges[deltaedge]
                if label[inblossom[i]] == 0:
                    i, j = j, i
                assert label[inblossom[i]] == 1
                queue.append(i)
            elif deltatype == 3:
                # Use the least-slack edge to continue the search.
                allowedge[deltaedge] = True
                (i, j, wt) = edges[deltaedge]
                assert label[inblossom[i]] == 1
                queue.append(i)
            elif deltatype == 4:
                # Expand the least-z blossom.
                expandBlossom(deltablossom, False)

            # End of a this substage.

        # Stop when no more augmenting path can be found.
        if not augmented:
            break

        # End of a stage; expand all S-blossoms which have dualvar = 0.
        for b in range(nvertex, 2*nvertex):
            if ( blossomparent[b] == -1 and blossombase[b] >= 0 and
                 label[b] == 1 and dualvar[b] == 0 ):
                expandBlossom(b, True)

    # Verify that we reached the optimum solution.
    if CHECK_OPTIMUM:
        verifyOptimum()

    # Transform mate[] such that mate[v] is the vertex to which v is paired.
    for v in range(nvertex):
        if mate[v] >= 0:
            mate[v] = endpoint[mate[v]]
    for v in range(nvertex):
        assert mate[v] == -1 or mate[mate[v]] == v

    return mate


def generating_bg(n):
    x=[4.0/211,44.0/211,59.0/211,60.0/211,102.0/211,104.0/211,115.0/211,119.0/211,120.0/211,121.0/211,122.0/211,123.0/211,
       156.0/211,200.0/211,210.0/211,1]
    generated_group=[]
    generated_rec=[]
    generated_don=[]
    for i in range(n):
        random_num=rng1.uniform(0,1)
        if random_num>=0 and random_num <=x[0]:
            #print ('generated group is A,A')
            generated_group.append('A,A')
            generated_rec.append('A')
            generated_don.append('A')
        elif random_num>=x[0] and random_num <=x[1]:
            
            #print ('generated group is A,B')
            generated_group.append('A,B')
            generated_rec.append('A')
            generated_don.append('B')
            
        elif random_num>=x[1] and random_num <=x[2]:
            
            #print 'generated group is A,AB'
            generated_group.append('A,AB')
            generated_rec.append('A')
            generated_don.append('AB')
            
        elif random_num>=x[2] and random_num <=x[3]:
            
            #print 'generated group is A,O'
            generated_group.append('A,O')
            generated_rec.append('A')
            generated_don.append('O')
            
        elif random_num>=x[3] and random_num <=x[4]:
            
            #print 'generated group is B,A'
            generated_group.append('B,A')
            generated_rec.append('B')
            generated_don.append('A')
            
        elif random_num>=x[4] and random_num <=x[5]:
            
            #print 'generated group is B,B'
            generated_group.append('B,B')
            generated_rec.append('B')
            generated_don.append('B')
            
        elif random_num>=x[5] and random_num <=x[6]:
            #print 'generated group is B,AB'
            generated_group.append('B,AB')
            generated_rec.append('B')
            generated_don.append('A')
            
        elif random_num>=x[6] and random_num <=x[7]:
            
            #print 'generated group is B,O'
            generated_group.append('B,O')
            generated_rec.append('B')
            generated_don.append('O')
            
        elif random_num>=x[7] and random_num <=x[8]:
            
            #print 'generated group is AB,A'
            generated_group.append('AB,A')
            generated_rec.append('AB')
            generated_don.append('A')
            
        elif random_num>=x[8] and random_num <=x[9]:
            
            #print 'generated group is AB,B'
            generated_group.append('AB,B')
            generated_rec.append('AB')
            generated_don.append('B')
            
        elif random_num>=x[9] and random_num <=x[10]:
            
            #print 'generated group is AB,AB'
            generated_group.append('AB,AB')
            generated_rec.append('AB')
            generated_don.append('AB')
            
        elif random_num>=x[10] and random_num <=x[11]:
            
            #print 'generated group is A,A'
            generated_group.append('AB,O')
            generated_rec.append('AB')
            generated_don.append('O')
            
        elif random_num>=x[11] and random_num <=x[12]:
            
            #print 'generated group is O,A'
            generated_group.append('O,A')
            generated_rec.append('O')
            generated_don.append('A')
            
        elif random_num>=x[12] and random_num <=x[13]:
            
            #print 'generated group is O,B'
            generated_group.append('O,B')
            generated_rec.append('O')
            generated_don.append('B')
            
        elif random_num>=x[13] and random_num <=x[14]:
            
            #print 'generated group is O,AB'
            generated_group.append('O,AB')
            generated_rec.append('O')
            generated_don.append('AB')
        elif random_num>=x[14] and random_num <=1:
            
            #print 'generated group is O,O'
            generated_group.append('O,O')
            generated_rec.append('O')
            generated_don.append('O')
        else:
            print ('No pair is generated')
    return generated_rec,generated_don

def compatible_edges(recipient_bg, donor_bg):
    Edges=[]
    for i in range(0,len(donor_bg)):
        for j in range(0,len(recipient_bg)):
            if ((donor_bg[i]==recipient_bg[j] and donor_bg[j]==recipient_bg[i] and donor_bg[i]!=recipient_bg[i] and donor_bg[j]!=recipient_bg[j]) and (donor_bg[i]!="Matched" and recipient_bg[i]!="Matched" and donor_bg[j]!="Matched" and recipient_bg[j]!="Matched")):
                Edges.append((i,j,1))
                continue

    return Edges

def Shapley_value(n,recipient_1,donor_1,recipient_2,donor_2,recipient_3,donor_3,recipient_4,donor_4):

    n=4
    
    Country1_rec=recipient_1
    Country1_don=donor_1

    Country2_rec=recipient_2
    Country2_don=donor_2

    Country3_rec=recipient_3
    Country3_don=donor_3

    Country4_rec=recipient_4
    Country4_don=donor_4

    Country12_rec=recipient_1+recipient_2
    Country12_don=donor_1+donor_2

    Country13_rec=recipient_1+recipient_3
    Country13_don=donor_1+donor_3

    Country14_rec=recipient_1+recipient_4
    Country14_don=donor_1+donor_4

    Country23_rec=recipient_2+recipient_3
    Country23_don=donor_2+donor_3

    Country24_rec=recipient_2+recipient_4
    Country24_don=donor_2+donor_4

    Country34_rec=recipient_3+recipient_4
    Country34_don=donor_3+donor_4

    Country123_rec=recipient_1+recipient_2+recipient_3
    Country123_don=donor_1+donor_2+donor_3

    Country124_rec=recipient_1+recipient_2+recipient_4
    Country124_don=donor_1+donor_2+donor_4

    Country134_rec=recipient_1+recipient_3+recipient_4
    Country134_don=donor_1+donor_3+donor_4

    Country234_rec=recipient_2+recipient_3+recipient_4
    Country234_don=donor_2+donor_3+donor_4

    Country1234_rec=recipient_1+recipient_2+recipient_3+recipient_4
    Country1234_don=donor_1+donor_2+donor_3+donor_4

    #Value for country 1
    Country1=compatible_edges(Country1_rec,Country1_don)
    Country1_pairs=maxWeightMatching(Country1, True)

    count_1=0
    for i in Country1_pairs:
        if i!= -1:
            count_1=count_1+1

    #Value for country 2
    Country2=compatible_edges(Country2_rec,Country2_don)
    Country2_pairs=maxWeightMatching(Country2, True)

    count_2=0
    for i in Country2_pairs:
        if i!= -1:
            count_2=count_2+1

    #Value for country 3
    Country3=compatible_edges(Country3_rec,Country3_don)   
    Country3_pairs=maxWeightMatching(Country3, True)

    count_3=0
    for i in Country3_pairs:
        if i!= -1:
            count_3=count_3+1

    #Value for country 4
    Country4=compatible_edges(Country4_rec,Country4_don)        
    Country4_pairs=maxWeightMatching(Country4, True)

    count_4=0
    for i in Country4_pairs:
        if i!= -1:
            count_4=count_4+1

    #Value for country 1 and 2
    Country12=compatible_edges(Country12_rec,Country12_don)                
    Country12_pairs=maxWeightMatching(Country12, True)

    count_12=0
    for i in Country12_pairs:
        if i!= -1:
            count_12=count_12+1
    
    #Value for country 1 and 3
    Country13=compatible_edges(Country13_rec,Country13_don)
    Country13_pairs=maxWeightMatching(Country13, True)

    count_13=0
    for i in Country13_pairs:
        if i!= -1:
            count_13=count_13+1

    #Value for country 1 and 4
    Country14=compatible_edges(Country14_rec,Country14_don)
    Country14_pairs=maxWeightMatching(Country14, True)

    count_14=0
    for i in Country14_pairs:
        if i!= -1:
            count_14=count_14+1
    
    #Value for country 2 and 3
    Country23=compatible_edges(Country23_rec,Country23_don)
    Country23_pairs=maxWeightMatching(Country23, True)

    count_23=0
    for i in Country23_pairs:
        if i!= -1:
            count_23=count_23+1

    #Value for country 2 and 4
    Country24=compatible_edges(Country24_rec,Country24_don)
    Country24_pairs=maxWeightMatching(Country24, True)

    count_24=0
    for i in Country24_pairs:
        if i!= -1:
            count_24=count_24+1

    #Value for country 3 and 4
    Country34=compatible_edges(Country34_rec,Country34_don)    
    Country34_pairs=maxWeightMatching(Country34, True)

    count_34=0
    for i in Country34_pairs:
        if i!= -1:
            count_34=count_34+1

    #Value for country 1, 2 and 3
    Country123=compatible_edges(Country123_rec,Country123_don) 
    Country123_pairs=maxWeightMatching(Country123, True)

    count_123=0
    for i in Country123_pairs:
        if i!= -1:
            count_123=count_123+1

    #Value for country 1, 2 and 4
    Country124=compatible_edges(Country124_rec,Country124_don)
    Country124_pairs=maxWeightMatching(Country124, True)

    count_124=0
    for i in Country124_pairs:
        if i!= -1:
            count_124=count_124+1

    #Value for country 1, 3 and 4
    Country134=compatible_edges(Country134_rec,Country134_don)
    Country134_pairs=maxWeightMatching(Country134, True)

    count_134=0
    for i in Country134_pairs:
        if i!= -1:
            count_134=count_134+1

    #Value for country 2, 3 and 4
    Country234=compatible_edges(Country234_rec,Country234_don)
    Country234_pairs=maxWeightMatching(Country234, True)

    count_234=0
    for i in Country234_pairs:
        if i!= -1:
            count_234=count_234+1

    #Value for country 1, 2, 3 and 4
    Country1234=compatible_edges(Country1234_rec,Country1234_don)
    Country1234_pairs=maxWeightMatching(Country1234, True)

    count_1234=0
    for i in Country1234_pairs:
        if i!= -1:
            count_1234=count_1234+1
    
    characteristic_function = [count_1,count_2,count_3,count_4,count_12,count_13,count_14,count_23,count_24,count_34,count_123,count_124,count_134,count_234,count_1234]
    tempList = list([i for i in range(n)])
    N = power_set(tempList)
    shapley_values = []
    for i in range(n):
        shapley = 0
        for j in N:
            if i not in j:
                cmod = len(j)
                Cui = j[:]
                bisect.insort_left(Cui,i)
                l = N.index(j)
                k = N.index(Cui)
                temp = float(float(characteristic_function[k]) - float(characteristic_function[l])) *\
                           float(math.factorial(cmod) * math.factorial(n - cmod - 1)) / float(math.factorial(n))
                shapley += temp
                # if i is 0:
                #     print j, Cui, cmod, n-cmod-1, characteristic_function[k], characteristic_function[l], math.factorial(cmod), math.factorial(n - cmod - 1), math.factorial(n)

        cmod = 0
        Cui = [i]
        k = N.index(Cui)
        temp = float(characteristic_function[k]) * float(math.factorial(cmod) * math.factorial(n - cmod - 1)) / float(math.factorial(n))
        shapley += temp

        shapley_values.append(shapley)
        
    return shapley_values

def Nucleous(n,recipient_1,donor_1,recipient_2,donor_2,recipient_3,donor_3,recipient_4,donor_4):
    
    Country1_rec=recipient_1
    Country1_don=donor_1

    Country2_rec=recipient_2
    Country2_don=donor_2

    Country3_rec=recipient_3
    Country3_don=donor_3

    Country4_rec=recipient_4
    Country4_don=donor_4

    Country12_rec=recipient_1+recipient_2
    Country12_don=donor_1+donor_2

    Country13_rec=recipient_1+recipient_3
    Country13_don=donor_1+donor_3

    Country14_rec=recipient_1+recipient_4
    Country14_don=donor_1+donor_4

    Country23_rec=recipient_2+recipient_3
    Country23_don=donor_2+donor_3

    Country24_rec=recipient_2+recipient_4
    Country24_don=donor_2+donor_4

    Country34_rec=recipient_3+recipient_4
    Country34_don=donor_3+donor_4

    Country123_rec=recipient_1+recipient_2+recipient_3
    Country123_don=donor_1+donor_2+donor_3

    Country124_rec=recipient_1+recipient_2+recipient_4
    Country124_don=donor_1+donor_2+donor_4

    Country134_rec=recipient_1+recipient_3+recipient_4
    Country134_don=donor_1+donor_3+donor_4

    Country234_rec=recipient_2+recipient_3+recipient_4
    Country234_don=donor_2+donor_3+donor_4

    Country1234_rec=recipient_1+recipient_2+recipient_3+recipient_4
    Country1234_don=donor_1+donor_2+donor_3+donor_4

    #Value for country 1
    Country1=compatible_edges(Country1_rec,Country1_don)
    Country1_pairs=maxWeightMatching(Country1, True)

    count_1=0
    for i in Country1_pairs:
        if i!= -1:
            count_1=count_1+1

    #Value for country 2
    Country2=compatible_edges(Country2_rec,Country2_don)
    Country2_pairs=maxWeightMatching(Country2, True)

    count_2=0
    for i in Country2_pairs:
        if i!= -1:
            count_2=count_2+1

    #Value for country 3
    Country3=compatible_edges(Country3_rec,Country3_don)   
    Country3_pairs=maxWeightMatching(Country3, True)

    count_3=0
    for i in Country3_pairs:
        if i!= -1:
            count_3=count_3+1

    #Value for country 4
    Country4=compatible_edges(Country4_rec,Country4_don)        
    Country4_pairs=maxWeightMatching(Country4, True)

    count_4=0
    for i in Country4_pairs:
        if i!= -1:
            count_4=count_4+1

    #Value for country 1 and 2
    Country12=compatible_edges(Country12_rec,Country12_don)                
    Country12_pairs=maxWeightMatching(Country12, True)

    count_12=0
    for i in Country12_pairs:
        if i!= -1:
            count_12=count_12+1
    
    #Value for country 1 and 3
    Country13=compatible_edges(Country13_rec,Country13_don)
    Country13_pairs=maxWeightMatching(Country13, True)

    count_13=0
    for i in Country13_pairs:
        if i!= -1:
            count_13=count_13+1

    #Value for country 1 and 4
    Country14=compatible_edges(Country14_rec,Country14_don)
    Country14_pairs=maxWeightMatching(Country14, True)

    count_14=0
    for i in Country14_pairs:
        if i!= -1:
            count_14=count_14+1
    
    #Value for country 2 and 3
    Country23=compatible_edges(Country23_rec,Country23_don)
    Country23_pairs=maxWeightMatching(Country23, True)

    count_23=0
    for i in Country23_pairs:
        if i!= -1:
            count_23=count_23+1

    #Value for country 2 and 4
    Country24=compatible_edges(Country24_rec,Country24_don)
    Country24_pairs=maxWeightMatching(Country24, True)

    count_24=0
    for i in Country24_pairs:
        if i!= -1:
            count_24=count_24+1

    #Value for country 3 and 4
    Country34=compatible_edges(Country34_rec,Country34_don)    
    Country34_pairs=maxWeightMatching(Country34, True)

    count_34=0
    for i in Country34_pairs:
        if i!= -1:
            count_34=count_34+1

    #Value for country 1, 2 and 3
    Country123=compatible_edges(Country123_rec,Country123_don) 
    Country123_pairs=maxWeightMatching(Country123, True)

    count_123=0
    for i in Country123_pairs:
        if i!= -1:
            count_123=count_123+1

    #Value for country 1, 2 and 4
    Country124=compatible_edges(Country124_rec,Country124_don)
    Country124_pairs=maxWeightMatching(Country124, True)

    count_124=0
    for i in Country124_pairs:
        if i!= -1:
            count_124=count_124+1

    #Value for country 1, 3 and 4
    Country134=compatible_edges(Country134_rec,Country134_don)
    Country134_pairs=maxWeightMatching(Country134, True)

    count_134=0
    for i in Country134_pairs:
        if i!= -1:
            count_134=count_134+1

    #Value for country 2, 3 and 4
    Country234=compatible_edges(Country234_rec,Country234_don)
    Country234_pairs=maxWeightMatching(Country234, True)

    count_234=0
    for i in Country234_pairs:
        if i!= -1:
            count_234=count_234+1

    #Value for country 1, 2, 3 and 4
    Country1234=compatible_edges(Country1234_rec,Country1234_don)
    Country1234_pairs=maxWeightMatching(Country1234, True)

    count_1234=0
    for i in Country1234_pairs:
        if i!= -1:
            count_1234=count_1234+1
            
    COALITIONS = [count_1,count_2,count_3,count_4,count_12,count_13,count_14,count_23,count_24,count_34,count_123,count_124,count_134,count_234,count_1234]

    LEMAIRE=definegame(4,COALITIONS)
    #print ('DefineGame',LEMAIRE)
    LEMAIRENucleolus = Nucleolus(LEMAIRE)
    print ('Nucleolus',LEMAIRENucleolus)
    k=unlist(LEMAIRENucleolus)
    Nucleous =[float(k[4]),float(k[5]),float(k[6]),float(k[7])]
    return Nucleous


def modified_graph(BG_rec, BG_don, n1, a1, b1,n2, a2, b2, n3, a3, b3, n4, a4, b4):

    B1=abs(n1-b1)
    B2=abs(n2-b2)
    B3=abs(n3-b3)
    B4=abs(n4-b4)
    A1=abs(b1-a1)
    A2=abs(b2-a2)
    A3=abs(b3-a3)
    A4=abs(b4-a4)

    Pairs=[]
    for i in range(len(BG_don)):
        for j in range(len(BG_rec)):
            if BG_don[i]==BG_rec[j] and BG_don[j]==BG_rec[i] and BG_don[i]!=BG_rec[i] and BG_don[j]!=BG_rec[j]:
                Pairs.append((i,j,1))
                continue

    for i in range(int(B1)):
        BG_rec.append("U")
                
    for i in range(n1):
        for j in range(len(BG_rec)):
            if BG_rec[j]=="U":
                Pairs.append((i,j,0))
                continue

    for i in range(int(A1)):
        BG_don.append('A11')
        BG_rec.append('A1')

    for i in range(n1):
        for j in range(len(BG_rec)):
            if BG_rec[j]=="A1":
                Pairs.append((i,j,0))
                continue
    
    for i in range(int(B2)):
        BG_rec.append("V")

    for i in range(n1,n1+n2):
        for j in range(len(BG_rec)):
            if BG_rec[j]=="V":
                Pairs.append((i,j,0))
                continue

    for i in range(int(A2)):
        BG_don.append('A22')
        BG_rec.append('A2')

    for i in range(n1,n1+n2):
        for j in range(len(BG_rec)):
            if BG_rec[j]=="A2":
                Pairs.append((i,j,0))
                continue

    for i in range(int(B3)):
        BG_rec.append("W")

    for i in range(n1+n2,n1+n2+n3):
        for j in range(len(BG_rec)):
            if BG_rec[j]=="W":
                Pairs.append((i,j,0))
                continue

    for i in range(int(A3)):
        BG_don.append('A33')
        BG_rec.append('A3')

    for i in range(n1+n2,n1+n2+n3):
        for j in range(len(BG_rec)):
            if BG_rec[j]=="A3":
                Pairs.append((i,j,0))
                continue

    for i in range(int(B4)):
        BG_rec.append("X")
                
    for i in range(n1+n2+n3,n1+n2+n3+n4):
        for j in range(len(BG_rec)):
            if BG_rec[j]=="X":
                Pairs.append((i,j,0))
                continue

    for i in range(int(A4)):
        BG_don.append('A44')
        BG_rec.append('A4')

    for i in range(n1+n2+n3,n1+n2+n3+n4):
        for j in range(len(BG_rec)):
            if BG_rec[j]=="A4":
                Pairs.append((i,j,0))
                continue
		 
    if len(Pairs)%2 !=0:
        BG_don.append('Add')
        BG_rec.append('Add')

    for i in range(len(BG_don)):
        for j in range(len(BG_rec)):
            if (BG_don[i]=="A11" or BG_don[i]=="A22" or BG_don[i]=="A33" or BG_don[i]=="A44" or BG_don[i]=="Add" ) and (BG_rec[j]=="A1" or BG_rec[j]=="A2" or BG_rec[j]=="A3" or BG_rec[j]=="A4" or BG_rec[j]=="Add"):
                Pairs.append((i,j,0))
                continue

    #print ('Modified country set',Country1_mod)
    #print ("Pairs", Pairs)
    
    return Pairs


IN_sol_mul_rep_c1=[]
IN_sol_mul_rep_c2=[]
IN_sol_mul_rep_c3=[]
IN_sol_mul_rep_c4=[]

GR_sol_mul_rep_c1=[]
GR_sol_mul_rep_c2=[]
GR_sol_mul_rep_c3=[]
GR_sol_mul_rep_c4=[]

SV_sol_mul_rep_c1=[]
SV_sol_mul_rep_c2=[]
SV_sol_mul_rep_c3=[]
SV_sol_mul_rep_c4=[]

NU_sol_mul_rep_c1=[]
NU_sol_mul_rep_c2=[]
NU_sol_mul_rep_c3=[]
NU_sol_mul_rep_c4=[]

for x in range(1):
	
	print ('Replication',x)
	
	RC_C1=0 #Remaining Credit of country 1 in Shapley Value case
	RC_C2=0
	RC_C3=0
	RC_C4=0
	
	NU_RC_C1=0 #Remaining Credit of country 1 in Nucleous case
	NU_RC_C2=0
	NU_RC_C3=0
	NU_RC_C4=0
	
	Number_recipient_country1=[] #Number of recipients generated in each round
	Number_recipient_country2=[]
	Number_recipient_country3=[]
	Number_recipient_country4=[]
	
	Recipient_BG_country1_total=[] #To record recipient BG over all 20 replications for countries
	Recipient_BG_country2_total=[]
	Recipient_BG_country3_total=[]
	Recipient_BG_country4_total=[]
	
	Donor_BG_country1_total=[] #To record Donor BG over all 20 replications for countries
	Donor_BG_country2_total=[]
	Donor_BG_country3_total=[]
	Donor_BG_country4_total=[]
	
	IN_BG_Rec_country1=[] #Global Random Recipient country 1
	IN_BG_Rec_country2=[]
	IN_BG_Rec_country3=[]
	IN_BG_Rec_country4=[]
	
	IN_BG_Don_country1=[]
	IN_BG_Don_country2=[]
	IN_BG_Don_country3=[]
	IN_BG_Don_country4=[]
	
	IN_sol_country1=[] #To record the optimal solution for each country if countries are match individually
	IN_sol_country2=[]
	IN_sol_country3=[]
	IN_sol_country4=[]
	    
	GR_BG_Rec_country1=[] #Global Random Recipient country 1
	GR_BG_Rec_country2=[]
	GR_BG_Rec_country3=[]
	GR_BG_Rec_country4=[]
	
	GR_BG_Don_country1=[]
	GR_BG_Don_country2=[]
	GR_BG_Don_country3=[]
	GR_BG_Don_country4=[]
	    
	GR_BG_Rec_Total=[]
	GR_BG_Don_Total=[]
	
	GR_sol_country1=[] # To record the solution for each country if assigned randomly in Grand Coalition
	GR_sol_country2=[]
	GR_sol_country3=[]
	GR_sol_country4=[]
	
	SV_BG_Rec_country1=[] #Bg recipient of country 1 considered for global coalition using shapley value
	SV_BG_Rec_country2=[]
	SV_BG_Rec_country3=[]
	SV_BG_Rec_country4=[]
	
	SV_BG_Don_country1=[]
	SV_BG_Don_country2=[]
	SV_BG_Don_country3=[]
	SV_BG_Don_country4=[]
	    
	SV_BG_Rec_Total=[]
	SV_BG_Don_Total=[]
	
	SV_sol_country1=[] #To record the solution for each country in Grand Coalition if assigned through theorem 3
	SV_sol_country2=[]
	SV_sol_country3=[]
	SV_sol_country4=[]
	
	NU_BG_Rec_country1=[] #Bg recipient of country 1 considered for global coalition using Nucleous
	NU_BG_Rec_country2=[]
	NU_BG_Rec_country3=[]
	NU_BG_Rec_country4=[]
	
	NU_BG_Don_country1=[]
	NU_BG_Don_country2=[]
	NU_BG_Don_country3=[]
	NU_BG_Don_country4=[]
	    
	NU_BG_Rec_Total=[]
	NU_BG_Don_Total=[]
	
	NU_sol_country1=[] #To record the solution for each country in Grand Coalition using Nucleous
	NU_sol_country2=[]
	NU_sol_country3=[]
	NU_sol_country4=[]
	
	IN_dropout_country1=0
	IN_dropout_country2=0
	IN_dropout_country3=0
	IN_dropout_country4=0
	
	GR_dropout_country1=0
	GR_dropout_country2=0
	GR_dropout_country3=0
	GR_dropout_country4=0
	
	SV_dropout_country1=0
	SV_dropout_country2=0
	SV_dropout_country3=0
	SV_dropout_country4=0
	
	NU_dropout_country1=0
	NU_dropout_country2=0
	NU_dropout_country3=0
	NU_dropout_country4=0
	
	IN_Waiting_time_country1=0
	IN_Waiting_time_country2=0
	IN_Waiting_time_country3=0
	IN_Waiting_time_country4=0
	
	GR_Waiting_time_country1=0
	GR_Waiting_time_country2=0
	GR_Waiting_time_country3=0
	GR_Waiting_time_country4=0
	
	SV_Waiting_time_country1=0
	SV_Waiting_time_country2=0
	SV_Waiting_time_country3=0
	SV_Waiting_time_country4=0
	
	NU_Waiting_time_country1=0
	NU_Waiting_time_country2=0
	NU_Waiting_time_country3=0
	NU_Waiting_time_country4=0
	
	p=0.3
	prob_fail_country1=0
	prob_fail_country2=0
	prob_fail_country3=0
	prob_fail_country4=0
	
	for z in range(20):
	    print ('Round',z)
	    random1=int(rng1.uniform(20,25)) #Number of pairs needed to be generated for country 1
	    random2=int(rng1.uniform(10,15))
	    random3=int(rng1.uniform(20,25))
	    random4=int(rng1.uniform(10,15))
	    
	    print ('random 1',random1)
	    print ('random 2',random2)
	    print ('random 3',random3)
	    print ('random 4',random4)
	
	    Number_recipient_country1.append(random1) #Store the number of recipients generated at stage z
	    Number_recipient_country2.append(random2)
	    Number_recipient_country3.append(random3)
	    Number_recipient_country4.append(random4)
	    
	    recipient1, donor1 = generating_bg(random1) #generating BG of pairs
	    recipient2, donor2 = generating_bg(random2)
	    recipient3, donor3 = generating_bg(random3)
	    recipient4, donor4 = generating_bg(random4)
	    
	    #print ('BG recipient1 and donor 1',recipient1, donor1)
	    #print ('BG recipient2 and donor 2',recipient2, donor2)
	    #print ('BG recipient3 and donor 3',recipient3, donor3)
	    #print ('BG recipient4 and donor 4',recipient4, donor4)
	
	    Recipient_BG_country1_total=Recipient_BG_country1_total+recipient1 # BG of Total number of recipients available for z^th run
	    Recipient_BG_country2_total=Recipient_BG_country2_total+recipient2
	    Recipient_BG_country3_total=Recipient_BG_country3_total+recipient3
	    Recipient_BG_country4_total=Recipient_BG_country4_total+recipient4
	    
	    #print ('Recipient_BG_country1_total before',Recipient_BG_country1_total)
	
	    Donor_BG_country1_total=Donor_BG_country1_total+donor1 #BG of total number of donors available for z^th run
	    Donor_BG_country2_total=Donor_BG_country2_total+donor2
	    Donor_BG_country3_total=Donor_BG_country3_total+donor3
	    Donor_BG_country4_total=Donor_BG_country4_total+donor4
	
	    compatible_edges_country1= compatible_edges(Recipient_BG_country1_total, Donor_BG_country1_total) #Compatible edges for each country
	    compatible_edges_country2= compatible_edges(Recipient_BG_country2_total, Donor_BG_country2_total)
	    compatible_edges_country3= compatible_edges(Recipient_BG_country3_total, Donor_BG_country3_total)
	    compatible_edges_country4= compatible_edges(Recipient_BG_country4_total, Donor_BG_country4_total)
	
	    count1=0
	    count2=0
	    count3=0
	    count4=0
	    for i in range(len(Recipient_BG_country1_total)):
	        if Recipient_BG_country1_total[i]!="Matched":
	            count1=count1+1
	    for i in range(len(Recipient_BG_country2_total)):
	        if Recipient_BG_country2_total[i]!="Matched":
	            count2=count2+1
	    for i in range(len(Recipient_BG_country3_total)):
	        if Recipient_BG_country3_total[i]!="Matched":
	            count3=count3+1
	    for i in range(len(Recipient_BG_country4_total)):
	        if Recipient_BG_country4_total[i]!="Matched":
	            count4=count4+1
	
	    #Number of recipient matched in country 1 without any coalition
	    Country1_matched_index=maxWeightMatching(compatible_edges_country1, True)
	    #print ('Country1_matched_index',Country1_matched_index)
	    
	    
	    Country1_count=0
	    indx=0
	    for i in Country1_matched_index:
	        if i!= -1:
	            Country1_count=Country1_count+1
	            Recipient_BG_country1_total[i]="Matched"
	            Donor_BG_country1_total[i]="Matched"
	            continue
	    
	 
	             
	    IN_sol_country1.append(Country1_count)
	    #print ('Recipient_BG_country1_total',Recipient_BG_country1_total)
	    #Number of recipient matched in country 2 without any coalition        
	    Country2_matched_index=maxWeightMatching(compatible_edges_country2, True)
	
	    Country2_count=0
	    for i in Country2_matched_index:
	        if i!= -1:
	            Country2_count=Country2_count+1
	            Recipient_BG_country2_total[i]="Matched"
	            Donor_BG_country2_total[i]="Matched"
	
	    IN_sol_country2.append(Country2_count)
	
	    #Number of recipient matched in country 3 without any coalition
	    Country3_matched_index=maxWeightMatching(compatible_edges_country3, True)
	
	    Country3_count=0
	    for i in Country3_matched_index:
	        if i!= -1:
	            Country3_count=Country3_count+1
	            Recipient_BG_country3_total[i]="Matched"
	            Donor_BG_country3_total[i]="Matched"
	
	    IN_sol_country3.append(Country3_count)
	
	    #Number of recipient matched in country 4 without any coalition
	    Country4_matched_index=maxWeightMatching(compatible_edges_country4, True)
	
	    Country4_count=0
	    for i in Country4_matched_index:
	        if i!= -1:
	            Country4_count=Country4_count+1
	            Recipient_BG_country4_total[i]="Matched"
	            Donor_BG_country4_total[i]="Matched"
	
	    IN_sol_country4.append(Country4_count)
	
	    #Deleting the pairs who are matched in zth run
	    #country 1
	    index1_R=[]
	    index1_D=[]
	    for i in range(len(Recipient_BG_country1_total)):
		    k1=rng3.uniform(0,1)
		    #print ('k1',k1)
		    if Recipient_BG_country1_total[i]!='Matched' and k1>=p:
			    index1_R.append(Recipient_BG_country1_total[i])
			    IN_Waiting_time_country1=IN_Waiting_time_country1 + 1
		    elif Recipient_BG_country1_total[i]!='Matched' and k1 <=p:
	              IN_dropout_country1=IN_dropout_country1+1
		    if Donor_BG_country1_total[i]!='Matched' and k1>=p:
	              index1_D.append(Donor_BG_country1_total[i])
	            
	    Recipient_BG_country1_total=index1_R
	    Donor_BG_country1_total=index1_D
	    #print ('Recipient_BG_country1_total',Recipient_BG_country1_total)
	    #country 2
	    index2_R=[]
	    index2_D=[]
	    for i in range(len(Recipient_BG_country2_total)):
		    k2=rng3.uniform(0,1)
		    #print ('k2',k2)
		    if Recipient_BG_country2_total[i]!='Matched' and k2>=p:
	              index2_R.append(Recipient_BG_country2_total[i])
	              IN_Waiting_time_country2=IN_Waiting_time_country2 + 1
		    elif Recipient_BG_country2_total[i]!='Matched' and k2 <=p:
	              IN_dropout_country2=IN_dropout_country2+1
		    if Donor_BG_country2_total[i]!='Matched' and k2>=p:
	              index2_D.append(Donor_BG_country2_total[i])
	    Recipient_BG_country2_total=index2_R
	    Donor_BG_country2_total=index2_D
	    #country 3
	    index3_R=[]
	    index3_D=[]
	    
	    for i in range(len(Recipient_BG_country3_total)):
		    k3=rng3.uniform(0,1)
		    #print ('k3',k3)
		    if Recipient_BG_country3_total[i]!='Matched' and k3>=p:
			    index3_R.append(Recipient_BG_country3_total[i])
			    IN_Waiting_time_country3=IN_Waiting_time_country3 + 1
		    elif Recipient_BG_country3_total[i]!='Matched' and k3 <=p:
	              IN_dropout_country3=IN_dropout_country3+1
		    if Donor_BG_country3_total[i]!='Matched' and k3>=p:
	              index3_D.append(Donor_BG_country3_total[i])
			    
	    Recipient_BG_country3_total=index3_R
	    Donor_BG_country3_total=index3_D
	    #country 4
	    index4_R=[]
	    index4_D=[]
	    
	    for i in range(len(Recipient_BG_country4_total)):
		    k4=rng3.uniform(0,1)
		    #print ('k4',k4)
		    if Recipient_BG_country4_total[i]!='Matched' and k4>=p:
			    index4_R.append(Recipient_BG_country4_total[i])
			    IN_Waiting_time_country4=IN_Waiting_time_country4 + 1
		    elif Recipient_BG_country4_total[i]!='Matched' and k4 <=p:
			    IN_dropout_country4=IN_dropout_country4+1
		    if Donor_BG_country4_total[i]!='Matched' and k4>=p:
			    index4_D.append(Donor_BG_country4_total[i])
	    Recipient_BG_country4_total=index4_R
	    Donor_BG_country4_total=index4_D
	    
	    
	    #Program for Random matching in a Global Registry

	    GR_BG_Rec_country1=GR_BG_Rec_country1+recipient1 #Global Random Recipient country 1
	    GR_BG_Rec_country2=GR_BG_Rec_country2+recipient2
	    GR_BG_Rec_country3=GR_BG_Rec_country3+recipient3
	    GR_BG_Rec_country4=GR_BG_Rec_country4+recipient4
	
	    GR_BG_Don_country1=GR_BG_Don_country1+donor1
	    GR_BG_Don_country2=GR_BG_Don_country2+donor2
	    GR_BG_Don_country3=GR_BG_Don_country3+donor3
	    GR_BG_Don_country4=GR_BG_Don_country4+donor4
	    
	    GR_BG_Rec_Total=GR_BG_Rec_country1+GR_BG_Rec_country2+GR_BG_Rec_country3+GR_BG_Rec_country4
	    GR_BG_Don_Total=GR_BG_Don_country1+GR_BG_Don_country2+GR_BG_Don_country3+GR_BG_Don_country4
	
	    #print ('GR_BG_Rec_Total',GR_BG_Rec_Total,len(GR_BG_Rec_Total),'GR_BG_Don_Total',GR_BG_Don_Total,len(GR_BG_Don_Total))
	
	    Compatible_edges_total= compatible_edges(GR_BG_Rec_Total, GR_BG_Don_Total)
	
	    #Number of recipient matched in all countries with coalition and random allocation
	    Country_matched_index=maxWeightMatching(Compatible_edges_total, True)
	
	    n1=len(GR_BG_Rec_country1)
	    n2=len(GR_BG_Rec_country2)
	    n3=len(GR_BG_Rec_country3)
	    n4=len(GR_BG_Rec_country4)
	    
	    #print ('n1,n1+n2,n1+n2+n3,n1+n2+n3+n4',n1,n1+n2,n1+n2+n3,n1+n2+n3+n4)
	
	    Country1_count=0
	    Country2_count=0
	    Country3_count=0
	    Country4_count=0
	    
	    #print ('Country matched index before',Country_matched_index)
	    #print ('GR_BG_Rec_country1 before',GR_BG_Rec_country1,len(GR_BG_Rec_country1))
	    
	    for i in Country_matched_index:
	        if i!= -1 and i <=n1-1:
	            Country1_count=Country1_count+1
	            GR_BG_Rec_country1[i]="Matched"
	            GR_BG_Don_country1[i]="Matched"
	        if i!=-1 and i >=n1 and i <= n1+n2-1:
	            Country2_count=Country2_count+1
	            GR_BG_Rec_country2[i-n1]="Matched"
	            GR_BG_Don_country2[i-n1]="Matched"
	        if i!=-1 and i >=n1+n2 and i <= n1+n2+n3-1:
	            Country3_count=Country3_count+1
	            GR_BG_Rec_country3[i-(n1+n2)]="Matched"
	            GR_BG_Don_country3[i-(n1+n2)]="Matched"
	        if i!=-1 and i >=n1+n2+n3 and i <= n1+n2+n3+n4:
	            Country4_count=Country4_count+1
	            GR_BG_Rec_country4[i-(n1+n2+n3)]="Matched"
	            GR_BG_Don_country4[i-(n1+n2+n3)]="Matched"
		
	    GR_sol_country1.append(Country1_count)
	    GR_sol_country2.append(Country2_count)
	    GR_sol_country3.append(Country3_count)
	    GR_sol_country4.append(Country4_count)
	
	    #Deleting the pairs who are matched in zth run
	    #country 1
	    index1_R=[]
	    index1_D=[]
	    for i in range(len(GR_BG_Rec_country1)):
		    g1=rng3.uniform(0,1)
		    #print ('g1',g1)
		    if GR_BG_Rec_country1[i]!='Matched' and g1>=p:
			    index1_R.append(GR_BG_Rec_country1[i])
			    GR_Waiting_time_country1=GR_Waiting_time_country1 + 1
		    elif GR_BG_Rec_country1[i]!='Matched' and g1<=p:
			    GR_dropout_country1=GR_dropout_country1+1
		    if GR_BG_Don_country1[i]!='Matched' and g1>=p:
			    index1_D.append(GR_BG_Don_country1[i])
	            
	    GR_BG_Rec_country1=index1_R
	    GR_BG_Don_country1=index1_D
	    
	    #country 2
	    index2_R=[]
	    index2_D=[]
	    for i in range(len(GR_BG_Rec_country2)):
		    g2=rng3.uniform(0,1)
		    #print ('g2',g2)
		    if GR_BG_Rec_country2[i]!='Matched' and g2>=p:
			    index2_R.append(GR_BG_Rec_country2[i])
			    GR_Waiting_time_country2=GR_Waiting_time_country2 + 1
		    elif GR_BG_Rec_country2[i]!='Matched' and g2<=p:
			    GR_dropout_country2=GR_dropout_country2+1
		    if GR_BG_Don_country2[i]!='Matched' and g2>=p:
			    index2_D.append(GR_BG_Don_country2[i])
	            
	    GR_BG_Rec_country2=index2_R
	    GR_BG_Don_country2=index2_D
	    #country 3
	    index3_R=[]
	    index3_D=[]
	    for i in range(len(GR_BG_Rec_country3)):
		    g3=rng3.uniform(0,1)
		    #print ('g3',g3)
		    if GR_BG_Rec_country3[i]!='Matched' and g3>=p:
			    index3_R.append(GR_BG_Rec_country3[i])
			    GR_Waiting_time_country3=GR_Waiting_time_country3 + 1
		    elif GR_BG_Rec_country3[i]!='Matched' and g3<=p:
			    GR_dropout_country3=GR_dropout_country3+1
		    if GR_BG_Don_country3[i]!='Matched' and g3>=p:
			    index3_D.append(GR_BG_Don_country3[i])
	    GR_BG_Rec_country3=index3_R
	    GR_BG_Don_country3=index3_D
	    #country 4
	    index4_R=[]
	    index4_D=[]
	    for i in range(len(GR_BG_Rec_country4)):
		    g4=rng3.uniform(0,1)
		    #print ('g4',g4)
		    if GR_BG_Rec_country4[i]!='Matched' and g4>=p:
			    index4_R.append(GR_BG_Rec_country4[i])
			    GR_Waiting_time_country4=GR_Waiting_time_country4 + 1
		    elif GR_BG_Rec_country4[i]!='Matched' and g4<=p:
			    GR_dropout_country4=GR_dropout_country4+1
		    if GR_BG_Don_country4[i]!='Matched' and g4>=p:
			    index4_D.append(GR_BG_Don_country4[i])
	            
	    GR_BG_Rec_country4=index4_R
	    GR_BG_Don_country4=index4_D
	    
	    SV_BG_Rec_country1=SV_BG_Rec_country1+recipient1 #Shapley value Recipient country 1
	    SV_BG_Rec_country2=SV_BG_Rec_country2+recipient2
	    SV_BG_Rec_country3=SV_BG_Rec_country3+recipient3
	    SV_BG_Rec_country4=SV_BG_Rec_country4+recipient4
	
	    SV_BG_Don_country1=SV_BG_Don_country1+donor1
	    SV_BG_Don_country2=SV_BG_Don_country2+donor2
	    SV_BG_Don_country3=SV_BG_Don_country3+donor3
	    SV_BG_Don_country4=SV_BG_Don_country4+donor4
	    
	    SV_BG_Rec_Total=SV_BG_Rec_country1+SV_BG_Rec_country2+SV_BG_Rec_country3+SV_BG_Rec_country4
	    SV_BG_Don_Total=SV_BG_Don_country1+SV_BG_Don_country2+SV_BG_Don_country3+SV_BG_Don_country4
	
	    shapley_values=Shapley_value(4,SV_BG_Rec_country1,SV_BG_Don_country1,SV_BG_Rec_country2,SV_BG_Don_country2,SV_BG_Rec_country3,SV_BG_Don_country3,SV_BG_Rec_country4,SV_BG_Don_country4)

	    C1=shapley_values[0]+RC_C1
	    C2=shapley_values[1]+RC_C2
	    C3=shapley_values[2]+RC_C3
	    C4=shapley_values[3]+RC_C4
	
	    a1=int(C1)-1
	    a2=int(C2)-1
	    a3=int(C3)-1
	    a4=int(C4)-1
	    b1=a1+2
	    b2=a2+2
	    b3=a3+2
	    b4=a4+2
	    n1=len(SV_BG_Rec_country1)
	    n2=len(SV_BG_Rec_country2)
	    n3=len(SV_BG_Rec_country3)
	    n4=len(SV_BG_Rec_country4)
	
	    New_graph=modified_graph(SV_BG_Rec_Total,SV_BG_Don_Total,n1,a1,b1,n2,a2,b2,n3,a3,b3,n4,a4,b4)
        
	    SV_Country_matched_index=maxWeightMatching(New_graph, True)
	
	    indx=0
	    Country1_count=0
	    Country2_count=0
	    Country3_count=0
	    Country4_count=0
	    for i in SV_Country_matched_index:
	        if indx <=n1+n2+n3+n4-1:
	            if i!= -1 and i <=n1-1 :
	                Country1_count=Country1_count+1
	                SV_BG_Rec_country1[i]="Matched"
	                SV_BG_Don_country1[i]="Matched"
	            if i!=-1 and i >=n1 and i <= n1+n2-1:
	                Country2_count=Country2_count+1
	                SV_BG_Rec_country2[i-n1]="Matched"
	                SV_BG_Don_country2[i-n1]="Matched"
	            if i!=-1 and i >=n1+n2 and i <= n1+n2+n3-1:
	                Country3_count=Country3_count+1
	                SV_BG_Rec_country3[i-(n1+n2)]="Matched"
	                SV_BG_Don_country3[i-(n1+n2)]="Matched"
	            if i!=-1 and i >=n1+n2+n3 and i <= n1+n2+n3+n4-1:
	                Country4_count=Country4_count+1
	                SV_BG_Rec_country4[i-(n1+n2+n3)]="Matched"
	                SV_BG_Don_country4[i-(n1+n2+n3)]="Matched"
	        indx=indx+1    
		
	    SV_sol_country1.append(Country1_count)
	    SV_sol_country2.append(Country2_count)
	    SV_sol_country3.append(Country3_count)
	    SV_sol_country4.append(Country4_count)
	
	    RC_C1=C1-Country1_count
	    RC_C2=C2-Country2_count
	    RC_C3=C3-Country3_count
	    RC_C4=C4-Country4_count
	    
	    #Deleting the pairs who are matched in zth run
	    #country 1
	    index1_R=[]
	    index1_D=[]
	    for i in range(len(SV_BG_Rec_country1)):
	        g1=rng3.uniform(0,1)
	        if SV_BG_Rec_country1[i]!='Matched' and g1>=p:
	            index1_R.append(SV_BG_Rec_country1[i])
	            SV_Waiting_time_country1=SV_Waiting_time_country1 + 1
	        elif SV_BG_Rec_country1[i]!='Matched' and g1<=p:
	            SV_dropout_country1=SV_dropout_country1+1
	            
	        if SV_BG_Don_country1[i]!='Matched' and g1>=p:
	            index1_D.append(SV_BG_Don_country1[i])
	    SV_BG_Rec_country1=index1_R
	    SV_BG_Don_country1=index1_D
	    #country 2
	    index2_R=[]
	    index2_D=[]
	    for i in range(len(SV_BG_Rec_country2)):
	        g2=rng3.uniform(0,1)
	        if SV_BG_Rec_country2[i]!='Matched' and g2>=p:
	            index2_R.append(SV_BG_Rec_country2[i])
	            SV_Waiting_time_country2=SV_Waiting_time_country2 + 1
	        elif SV_BG_Rec_country2[i]!='Matched' and g2<=p:
	            SV_dropout_country2=SV_dropout_country2+1
	            
	        if SV_BG_Don_country2[i]!='Matched' and g2>=p:
	            index2_D.append(SV_BG_Don_country2[i])
	    SV_BG_Rec_country2=index2_R
	    SV_BG_Don_country2=index2_D
	    #country 3
	    index3_R=[]
	    index3_D=[]
	    for i in range(len(SV_BG_Rec_country3)):
	        g3=rng3.uniform(0,1)
	        if SV_BG_Rec_country3[i]!='Matched' and g3>=p:
	            index3_R.append(SV_BG_Rec_country3[i])
	            SV_Waiting_time_country3=SV_Waiting_time_country3 + 1
	        elif SV_BG_Rec_country3[i]!='Matched' and g3<=p:
	            SV_dropout_country3=SV_dropout_country3+1
	            
	        if SV_BG_Don_country3[i]!='Matched' and g3>=p:
	            index3_D.append(SV_BG_Don_country3[i])
	    SV_BG_Rec_country3=index3_R
	    SV_BG_Don_country3=index3_D
	    #country 4
	    index4_R=[]
	    index4_D=[]
	    for i in range(len(SV_BG_Rec_country4)):
	        g4=rng3.uniform(0,1)
	        if SV_BG_Rec_country4[i]!='Matched' and g4>=p:
	            index4_R.append(SV_BG_Rec_country4[i])
	            SV_Waiting_time_country4=SV_Waiting_time_country4 + 1
	        elif SV_BG_Rec_country4[i]!='Matched' and g4<=p:
	            SV_dropout_country4=SV_dropout_country4+1
	            
	        if SV_BG_Don_country4[i]!='Matched' and g4>=p:
	            index4_D.append(SV_BG_Don_country4[i])
	    SV_BG_Rec_country4=index4_R
	    SV_BG_Don_country4=index4_D  
	    #end of program for Global allocation using credit system
	
	    #Program for Global allocation using Nucleous as target solution
	    
	    NU_BG_Rec_country1=NU_BG_Rec_country1+recipient1 #Shapley value Recipient country 1
	    NU_BG_Rec_country2=NU_BG_Rec_country2+recipient2
	    NU_BG_Rec_country3=NU_BG_Rec_country3+recipient3
	    NU_BG_Rec_country4=NU_BG_Rec_country4+recipient4
	
	    NU_BG_Don_country1=NU_BG_Don_country1+donor1
	    NU_BG_Don_country2=NU_BG_Don_country2+donor2
	    NU_BG_Don_country3=NU_BG_Don_country3+donor3
	    NU_BG_Don_country4=NU_BG_Don_country4+donor4
	    
	    NU_BG_Rec_Total=NU_BG_Rec_country1+NU_BG_Rec_country2+NU_BG_Rec_country3+NU_BG_Rec_country4
	    NU_BG_Don_Total=NU_BG_Don_country1+NU_BG_Don_country2+NU_BG_Don_country3+NU_BG_Don_country4
	    
	    NU=Nucleous(4,NU_BG_Rec_country1,NU_BG_Don_country1,NU_BG_Rec_country2,NU_BG_Don_country2,NU_BG_Rec_country3,NU_BG_Don_country3,NU_BG_Rec_country4,NU_BG_Don_country4)
	
	    C1=NU[0]+NU_RC_C1
	    C2=NU[1]+NU_RC_C2
	    C3=NU[2]+NU_RC_C3
	    C4=NU[3]+NU_RC_C4
	
	    a1=int(C1)-1
	    a2=int(C2)-1
	    a3=int(C3)-1
	    a4=int(C4)-1
	    b1=a1+2
	    b2=a2+2
	    b3=a3+2
	    b4=a4+2
	    n1=len(NU_BG_Rec_country1)
	    n2=len(NU_BG_Rec_country2)
	    n3=len(NU_BG_Rec_country3)
	    n4=len(NU_BG_Rec_country4)
	
	    New_graph=modified_graph(NU_BG_Rec_Total,NU_BG_Don_Total,n1,a1,b1,n2,a2,b2,n3,a3,b3,n4,a4,b4)
	
	    NU_Country_matched_index=maxWeightMatching(New_graph, True)
	
	    indx=0
	    Country1_count=0
	    Country2_count=0
	    Country3_count=0
	    Country4_count=0
	    
	    for i in NU_Country_matched_index:
	        if indx <=n1+n2+n3+n4-1:
	            if i!= -1 and i <=n1-1:
	                Country1_count=Country1_count+1
	                NU_BG_Rec_country1[i]="Matched"
	                NU_BG_Don_country1[i]="Matched"
	            if i!=-1 and i >=n1 and i <= n1+n2-1:
	                Country2_count=Country2_count+1
	                NU_BG_Rec_country2[i-n1]="Matched"
	                NU_BG_Don_country2[i-n1]="Matched"
	            if i!=-1 and i >=n1+n2 and i <= n1+n2+n3-1:
	                Country3_count=Country3_count+1
	                NU_BG_Rec_country3[i-(n1+n2)]="Matched"
	                NU_BG_Don_country3[i-(n1+n2)]="Matched"
	            if i!=-1 and i >=n1+n2+n3 and i <= n1+n2+n3+n4-1:
	                Country4_count=Country4_count+1
	                NU_BG_Rec_country4[i-(n1+n2+n3)]="Matched"
	                NU_BG_Don_country4[i-(n1+n2+n3)]="Matched"
	        indx=indx+1
	    
	    
	# =============================================================================
	#     for i in NU_Country_matched_index:
	#         if indx <=n1+n2+n3+n4-1:
	#             if i!= -1 and i <=n1-1:
	#                 Country1_count=Country1_count+1
	#                 NU_BG_Rec_country1[i]="Matched"
	#                 NU_BG_Don_country1[i]="Matched"
	#             if i!=-1 and i >=n1 and i <= n1+n2-1:
	#                 Country2_count=Country2_count+1
	#                 NU_BG_Rec_country2[i-n1]="Matched"
	#                 NU_BG_Don_country2[i-n1]="Matched"
	#             if i!=-1 and i >=n1+n2 and i <= n1+n2+n3-1:
	#                 Country3_count=Country3_count+1
	#                 NU_BG_Rec_country3[i-(n1+n2)]="Matched"
	#                 NU_BG_Don_country3[i-(n1+n2)]="Matched"
	#             if i!=-1 and i >=n1+n2+n3 and i <= n1+n2+n3+n4-1:
	#                 Country4_count=Country4_count+1
	#                 NU_BG_Rec_country4[i-(n1+n2+n3)]="Matched"
	#                 NU_BG_Don_country4[i-(n1+n2+n3)]="Matched"
	#         indx=indx+1
	# 
	# =============================================================================
	    NU_sol_country1.append(Country1_count)
	    NU_sol_country2.append(Country2_count)
	    NU_sol_country3.append(Country3_count)
	    NU_sol_country4.append(Country4_count)
	
	    NU_RC_C1=C1-Country1_count
	    NU_RC_C2=C2-Country2_count
	    NU_RC_C3=C3-Country3_count
	    NU_RC_C4=C4-Country4_count
	    #Deleting the pairs who are matched in zth run
	    #country 1
	    index1_R=[]
	    index1_D=[]
	    for i in range(len(NU_BG_Rec_country1)):
	        g1=rng3.uniform(0,1)
	        if NU_BG_Rec_country1[i]!='Matched' and g1>=p:
	            index1_R.append(NU_BG_Rec_country1[i])
	            NU_Waiting_time_country1=NU_Waiting_time_country1 + 1
	        elif NU_BG_Rec_country1[i]!='Matched' and g1<=p:
	            NU_dropout_country1=NU_dropout_country1+1
	        if NU_BG_Don_country1[i]!='Matched' and g1>=p:
	            index1_D.append(NU_BG_Don_country1[i])
	    NU_BG_Rec_country1=index1_R
	    NU_BG_Don_country1=index1_D
	    #country 2
	    index2_R=[]
	    index2_D=[]
	    for i in range(len(NU_BG_Rec_country2)):
	        g2=rng3.uniform(0,1)
	        if NU_BG_Rec_country2[i]!='Matched' and g2>=p:
	            index2_R.append(NU_BG_Rec_country2[i])
	            NU_Waiting_time_country2=NU_Waiting_time_country2 + 1
	        elif NU_BG_Rec_country2[i]!='Matched' and g2<=p:
	            NU_dropout_country2=NU_dropout_country2+1
	            
	        if NU_BG_Don_country2[i]!='Matched' and g2>=p:
	            index2_D.append(NU_BG_Don_country2[i])
	    NU_BG_Rec_country2=index2_R
	    NU_BG_Don_country2=index2_D
	    #country 3
	    index3_R=[]
	    index3_D=[]
	    for i in range(len(NU_BG_Rec_country3)):
	        g3=rng3.uniform(0,1)
	        if NU_BG_Rec_country3[i]!='Matched' and g3>=p:
	            index3_R.append(NU_BG_Rec_country3[i])
	            NU_Waiting_time_country3=NU_Waiting_time_country3 + 1
	        elif NU_BG_Rec_country3[i]!='Matched' and g3<=p:
	            NU_dropout_country3=NU_dropout_country3+1
	            
	        if NU_BG_Don_country3[i]!='Matched' and g3>=p:
	            index3_D.append(NU_BG_Don_country3[i])
	    NU_BG_Rec_country3=index3_R
	    NU_BG_Don_country3=index3_D
	    #country 4
	    index4_R=[]
	    index4_D=[]
	    for i in range(len(NU_BG_Rec_country4)):
	        g4=rng3.uniform(0,1)
	        if NU_BG_Rec_country4[i]!='Matched' and g4>=p:
	            index4_R.append(NU_BG_Rec_country4[i])
	            NU_Waiting_time_country4=NU_Waiting_time_country4 + 1
	        elif NU_BG_Rec_country4[i]!='Matched' and g4<=p:
	            NU_dropout_country4=NU_dropout_country4+1
	        if NU_BG_Don_country4[i]!='Matched' and g4>=p:
	            index4_D.append(NU_BG_Don_country4[i])
	    NU_BG_Rec_country4=index4_R
	    NU_BG_Don_country4=index4_D  
	    #end of program for Global allocation using credit system

	IN_sol_mul_rep_c1.append(np.mean(IN_sol_country1))
	IN_sol_mul_rep_c2.append(np.mean(IN_sol_country2))
	IN_sol_mul_rep_c3.append(np.mean(IN_sol_country3))
	IN_sol_mul_rep_c4.append(np.mean(IN_sol_country4))
	
	GR_sol_mul_rep_c1.append(np.mean(GR_sol_country1))
	GR_sol_mul_rep_c2.append(np.mean(GR_sol_country2))
	GR_sol_mul_rep_c3.append(np.mean(GR_sol_country3))
	GR_sol_mul_rep_c4.append(np.mean(GR_sol_country4))

	SV_sol_mul_rep_c1.append(np.mean(SV_sol_country1))
	SV_sol_mul_rep_c2.append(np.mean(SV_sol_country2))
	SV_sol_mul_rep_c3.append(np.mean(SV_sol_country3))
	SV_sol_mul_rep_c4.append(np.mean(SV_sol_country4))
	
	NU_sol_mul_rep_c1.append(np.mean(NU_sol_country1))
	NU_sol_mul_rep_c2.append(np.mean(NU_sol_country2))
	NU_sol_mul_rep_c3.append(np.mean(NU_sol_country3))
	NU_sol_mul_rep_c4.append(np.mean(NU_sol_country4))
	
	
print ('Avg num of recipient matched in country 1 in IN case',np.mean(IN_sol_mul_rep_c1))
print ('Avg num of recipient matched in country 2 in IN case',np.mean(IN_sol_mul_rep_c2))
print ('Avg num of recipient matched in country 3 in IN case',np.mean(IN_sol_mul_rep_c3))
print ('Avg num of recipient matched in country 4 in IN case',np.mean(IN_sol_mul_rep_c4))

print ('Avg num of recipients matched in all countries in IN case',np.mean(IN_sol_mul_rep_c1)+np.mean(IN_sol_mul_rep_c2)+np.mean(IN_sol_mul_rep_c3)+np.mean(IN_sol_mul_rep_c4))

print ('Avg num of recipient matched in country 1 in GR case',np.mean(GR_sol_mul_rep_c1))
print ('Avg num of recipient matched in country 2 in GR case',np.mean(GR_sol_mul_rep_c2))
print ('Avg num of recipient matched in country 3 in GR case',np.mean(GR_sol_mul_rep_c3))
print ('Avg num of recipient matched in country 4 in GR case',np.mean(GR_sol_mul_rep_c4))

print ('Avg num of recipients matched in all countries in IN case',np.mean(GR_sol_mul_rep_c1)+np.mean(GR_sol_mul_rep_c2)+np.mean(GR_sol_mul_rep_c3)+np.mean(GR_sol_mul_rep_c4))

print ('Avg num of recipient matched in country 1 in SV case',np.mean(SV_sol_mul_rep_c1))
print ('Avg num of recipient matched in country 2 in SV case',np.mean(SV_sol_mul_rep_c2))
print ('Avg num of recipient matched in country 3 in SV case',np.mean(SV_sol_mul_rep_c3))
print ('Avg num of recipient matched in country 4 in SV case',np.mean(SV_sol_mul_rep_c4))

print ('Avg num of recipients matched in all countries in IN case',np.mean(SV_sol_mul_rep_c1)+np.mean(SV_sol_mul_rep_c2)+np.mean(SV_sol_mul_rep_c3)+np.mean(SV_sol_mul_rep_c4))

print ('Avg num of recipient matched in country 1 in NU case',np.mean(NU_sol_mul_rep_c1))
print ('Avg num of recipient matched in country 2 in NU case',np.mean(NU_sol_mul_rep_c2))
print ('Avg num of recipient matched in country 3 in NU case',np.mean(NU_sol_mul_rep_c3))
print ('Avg num of recipient matched in country 4 in NU case',np.mean(NU_sol_mul_rep_c4))

print ('Avg num of recipients matched in all countries in IN case',np.mean(NU_sol_mul_rep_c1)+np.mean(NU_sol_mul_rep_c2)+np.mean(NU_sol_mul_rep_c3)+np.mean(NU_sol_mul_rep_c4))

I=np.mean(IN_sol_mul_rep_c1)+np.mean(IN_sol_mul_rep_c2)+np.mean(IN_sol_mul_rep_c3)+np.mean(IN_sol_mul_rep_c4)
R=np.mean(GR_sol_mul_rep_c1)+np.mean(GR_sol_mul_rep_c2)+np.mean(GR_sol_mul_rep_c3)+np.mean(GR_sol_mul_rep_c4)
S=np.mean(SV_sol_mul_rep_c1)+np.mean(SV_sol_mul_rep_c2)+np.mean(SV_sol_mul_rep_c3)+np.mean(SV_sol_mul_rep_c4)
N=np.mean(NU_sol_mul_rep_c1)+np.mean(NU_sol_mul_rep_c2)+np.mean(NU_sol_mul_rep_c3)+np.mean(NU_sol_mul_rep_c4)


y1=[np.mean(IN_sol_mul_rep_c1),np.mean(IN_sol_mul_rep_c2),np.mean(IN_sol_mul_rep_c3),np.mean(IN_sol_mul_rep_c4),I]
y2=[np.mean(GR_sol_mul_rep_c1),np.mean(GR_sol_mul_rep_c2),np.mean(GR_sol_mul_rep_c3),np.mean(GR_sol_mul_rep_c4),R]
y3=[np.mean(SV_sol_mul_rep_c1),np.mean(SV_sol_mul_rep_c2),np.mean(SV_sol_mul_rep_c3),np.mean(SV_sol_mul_rep_c4),S]
y4=[np.mean(NU_sol_mul_rep_c1),np.mean(NU_sol_mul_rep_c2),np.mean(NU_sol_mul_rep_c3),np.mean(NU_sol_mul_rep_c4),N]


xpos=np.arange(5)
plt.xticks(xpos,('Country1','Country2','Country3','Country4','Total'))
plt.bar(xpos-0.30, y1, color = 'grey', width = 0.15, label='Individual allocation')
plt.bar(xpos-0.10, y2, color = 'red', width = 0.15, label='Global random allocation')
plt.bar(xpos+0.10, y3, color= 'blue',width=0.15, label='Credit based allocation using Shapley Value')
plt.bar(xpos+0.30, y4, color= 'green',width=0.15, label='Credit based allocation using Nucleous')
plt.xlabel('Countries and total, Dropout probability = 0.3 \n Arrival rate C1 - Uniform(20-25), C2 - Uniform(10-15), C3 - Uniform(20-25), C4 - Uniform(10-15)')
plt.ylabel('Total number of recipients matched over 20 rounds')
plt.title('Comparison between individual vs Global random vs credit based allocation with hard to match patients')
plt.legend()
plt.show()

print ('Total matches in IN case',I)
print ('Total matches in GR case',R)
print ('Total matches in SV case',S)
print ('Total matches in NU case',N)











