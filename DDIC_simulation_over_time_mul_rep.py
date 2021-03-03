# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 11:16:25 2019

@author: Dr. Utkarsh Verma
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 19:27:34 2019

@author: Dr. Utkarsh Verma
"""

import random
from itertools import combinations
import math
import bisect
import matplotlib.pyplot as plt
import numpy as np

DEBUG = None
"""def DEBUG(s):
    from sys import stderr
    print('DEBUG:', s, file=stderr)
"""
random.seed(11)

# Check delta2/delta3 computation after every substage;
# only works on integer weights, slows down the algorithm to O(n^4).
CHECK_DELTA = False

# Check optimality of solution before returning; only works on integer weights.
CHECK_OPTIMUM = True

#Program for Deceased donor chains implementation over time

#Assumption for simulations:
#1. BG distribution of pairs follow ASTRA registry data
#2. BG distribution for Deceased donor will follow general BG distribution in India
#3. Arrival rate for Pair will follow uniform/Poisson between 10 to 15.
#4. Arrival rate for DD will follow uniform/poisson between 5 to 10.
#5. We will assume that there is a large wait list for DD kidney so all unutilize DD kidney will be
#   consumed in the wait list and there won't be any DD kidney which will go unutilize
#6. We will simulate this for a period of 5 years.  
#7. Probability of failure for a unmatch pair to go to next round will be 0.1 or 0.3

#Algorithm steps:

#1. Generate Pairs according to the ASTRA BG distribution 
#2. Generate DD and for each DD add 4 different types of BG wait list patients for that each DD 
#   is consumed with the program
#3. Create compatible edges considering both DD and pairs
#4. Solve it for optimality using Edmond's algorithm 
#5. Repeat this experiment for 5 years (60 rounds)
#.6 Record the relative gain of DD chains over independent functioning registries


#Edmond's Algorithm
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
            assert s >= 0
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

#BG generations of pairs based on ASTRA data
def generating_bg(n):
    x=[4.0/211,44.0/211,59.0/211,60.0/211,102.0/211,104.0/211,115.0/211,119.0/211,120.0/211,121.0/211,122.0/211,123.0/211,
       156.0/211,200.0/211,210.0/211,1]
    generated_group=[]
    generated_rec=[]
    generated_don=[]
    for i in range(n):
        random_num=random.uniform(0,1)
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

#Generating BG for DD
def generating_bg_DD(n):
	x=[0.37,0.69,0.92,1] #BG probability - O,B,A,AB
	generated_DD=[]
	for i in range(n):
		random_num=random.uniform(0,1)
		if random_num <=x[0]:
			generated_DD.append('O')
		elif random_num >= x[0] and random_num <= x[1]:
			generated_DD.append('B')
		elif random_num >= x[1] and random_num <= x[2]:
			generated_DD.append('A')
		elif random_num >= x[2] and random_num <= 1:
			generated_DD.append('AB')
			
	return generated_DD   

#Creating compatible edges
def compatible_edges(recipient_bg, donor_bg):
    Edges=[]
    for i in range(0,len(donor_bg)):	    
        for j in range(0,len(recipient_bg)):
            if ((donor_bg[i]==recipient_bg[j] and donor_bg[j]==recipient_bg[i] and donor_bg[i]!=recipient_bg[i] and donor_bg[j]!=recipient_bg[j] or donor_bg[i]=='O' or recipient_bg[j]=='AB') and (donor_bg[i]!="Matched" and recipient_bg[i]!="Matched" and donor_bg[j]!="Matched" and recipient_bg[j]!="Matched")):
                Edges.append((i,j,2))
                continue

    return Edges

def compatible_edges_DD(recipient_bg,DD_bg):
	Edges=[]
	WL_bg=DD_bg
	k=len(DD_bg)
	l=len(recipient_bg)
	for i in range(len(DD_bg)):
		for j in range(len(recipient_bg)):
			#print ('i,j',i+len(recipient_bg),j)
			if (DD_bg[i]==recipient_bg[j] and DD_bg[i]!="Matched" and recipient_bg[j]!="Matched"):
				Edges.append((i+len(recipient_bg),j,2))
				continue
	
	for i in range(len(DD_bg)):
		for j in range(len(WL_bg)):
			if (DD_bg[i]==WL_bg[j] and DD_bg[i]!='Matched' and WL_bg[j]!='Matched'):
				Edges.append((l+i,k+l+j,1))
				continue
	return Edges

Rep_Additional_gain=[]
Rep_Swap_solution=[]
Rep_DDIC_solution=[]

Rep_Pair_waiting_time=0 # Waiting time in swap allocation
Rep_Dropout_pair=[] #Number of dropout pairs in swap allocation
Rep_DD_Pair_waiting_time=0
Rep_DD_Dropout_pair=[]

Rep_WT_BG_O_Swap=[] #Waiting time for O recipient in swap registry
Rep_WT_BG_A_Swap=[]
Rep_WT_BG_B_Swap=[]
Rep_WT_BG_AB_Swap=[]

Rep_WT_BG_O_DDIC=[] #Waiting time for O recipient in DDIC
Rep_WT_BG_A_DDIC=[]
Rep_WT_BG_B_DDIC=[]
Rep_WT_BG_AB_DDIC=[]

Rep_DO_BG_O_Swap=[] # Total number of dropouts for O recipients in swap registry
Rep_DO_BG_A_Swap=[]
Rep_DO_BG_B_Swap=[]
Rep_DO_BG_AB_Swap=[]

Rep_DO_BG_O_DDIC=[] #Total number of dropouts for O recipient in DDIC
Rep_DO_BG_A_DDIC=[]
Rep_DO_BG_B_DDIC=[]
Rep_DO_BG_AB_DDIC=[]

Rep_count_O=[] #Count total number of O recipient generated in swap registry
Rep_count_A=[]
Rep_count_B=[]
Rep_count_AB=[]

n=60 #total number of Simulation rounds
p=0.3 #probability of dropout for a pair
r=30 #total number of repitation of the process
	
#Define arrays for storing data 
for a in range(r):	    
	Num_pair_generated=[] #Number of pairs generated in each round
	Num_DD_generated=[] #Number of DD generated in each round
	
	Donor_BG_pair_total=[] #List of donor's BG for pairs
	Recipient_BG_pair_total=[] #List of recipient's BG for pairs
	
	DDIC_Donor_BG_pair_total=[] #List of donor's BG for pairs in DDIC allocation
	DDIC_Recipient_BG_pair_total=[] #List of recipient's BG for pairs in DDIC allocation
	
	DD_BG_total=[] #List of DD's BG
	
	Num_pair_matched=[] #Number of pairs matched in each round
	Num_WL_matched=[] #Number of wait list patients matched in each round
	
	Additional_gain=[]
	Swap_solution=[]
	DDIC_solution=[]
	
	Pair_waiting_time=0 # Waiting time in swap allocation
	Dropout_pair=[] #Number of dropout pairs in swap allocation
	DD_Pair_waiting_time=0
	DD_Dropout_pair=[]
	
	WT_BG_O_Swap=0 #Waiting time for O recipient in swap registry
	WT_BG_A_Swap=0
	WT_BG_B_Swap=0
	WT_BG_AB_Swap=0
	
	WT_BG_O_DDIC=0 #Waiting time for O recipient in DDIC
	WT_BG_A_DDIC=0
	WT_BG_B_DDIC=0
	WT_BG_AB_DDIC=0
	
	DO_BG_O_Swap=0 # Total number of dropouts for O recipients in swap registry
	DO_BG_A_Swap=0
	DO_BG_B_Swap=0
	DO_BG_AB_Swap=0
	
	DO_BG_O_DDIC=0 #Total number of dropouts for O recipient in DDIC
	DO_BG_A_DDIC=0
	DO_BG_B_DDIC=0
	DO_BG_AB_DDIC=0
	
	count_O=0 #Count total number of O recipient generated in swap registry
	count_A=0
	count_B=0
	count_AB=0
	
	
	for i in range(n):
		Random_pair=random.randint(20,25) #Number of pairs generated in each round
		Random_DD=random.randint(10,15) #Number of DD generated in each round
		print ('Random pair',Random_pair)
		print ('Random DD',Random_DD)
		
		Num_pair_generated.append(Random_pair) #Store number of pairs generated
		Num_DD_generated.append(Random_DD) #Store number of DD generated
		
		generated_rec, generated_don=generating_bg(Random_pair) #Generating BG of pairs
		generated_DD=generating_bg_DD(Random_DD) #Generating BG of DD
		
		#print ('BG_rec_pair, BG_donor_pair',generated_rec, generated_don)
		#print ('BG_donor_DD', generated_DD)
		
		for i in range(len(generated_rec)):
			if generated_rec[i]=='O':
				count_O=count_O+1
			if generated_rec[i]=='A':
				count_A=count_A+1
			if generated_rec[i]=='B':
				count_B=count_B+1
			if generated_rec[i]=='AB':
				count_AB=count_AB+1
		
		#print ('generated_rec,generated_don',generated_rec,generated_don)
		#print ('Generated DD',generated_DD)
		
		Donor_BG_pair_total=Donor_BG_pair_total+generated_don #Total number of recipients considered in the current round for swap allocation
		Recipient_BG_pair_total=Recipient_BG_pair_total+generated_rec #Total number of donors considered in the current round for swap allocation
		Compatible_edges_pair=compatible_edges(Recipient_BG_pair_total, Donor_BG_pair_total)
		#print ('Compatible_edges_pair',Compatible_edges_pair)
		
		Matched_index_pair=maxWeightMatching(Compatible_edges_pair, True)
		#print ('Matched_index_pair',Matched_index_pair)
		
		count_matched_pair=0
		for i in Matched_index_pair:
			if i!=-1:
				count_matched_pair=count_matched_pair+1
				continue
		swap_DD=count_matched_pair+2*Random_DD
		print ('Swap_DD',swap_DD)
		Swap_solution.append(swap_DD)
		
		for i in range(len(Matched_index_pair)):
			if Matched_index_pair[i]!=-1:
				Recipient_BG_pair_total[i]='Matched'
				Donor_BG_pair_total[i]='Matched'
				continue
			
		#print ('Recipient_BG_pair_total after',Recipient_BG_pair_total)
		
		indx_R=[]
		indx_D=[]
		Num_Dropout_pair=0
		for i in range(len(Recipient_BG_pair_total)):
			g1=random.uniform(0,1)
			if Recipient_BG_pair_total[i]!='Matched' and g1 >=p:
				indx_R.append(Recipient_BG_pair_total[i])
				indx_D.append(Donor_BG_pair_total[i])
				Pair_waiting_time=Pair_waiting_time+1
				if Recipient_BG_pair_total[i]=='O':
					WT_BG_O_Swap=WT_BG_O_Swap+1
				if Recipient_BG_pair_total[i]=='A':
					WT_BG_A_Swap=WT_BG_A_Swap+1
				if Recipient_BG_pair_total[i]=='B':
					WT_BG_B_Swap=WT_BG_B_Swap+1
				if Recipient_BG_pair_total[i]=='AB':
					WT_BG_AB_Swap=WT_BG_AB_Swap+1
			elif Recipient_BG_pair_total[i]!='Matched' and g1 <=p:
				Num_Dropout_pair=Num_Dropout_pair+1
				if Recipient_BG_pair_total[i]=='O':
					DO_BG_O_Swap=DO_BG_O_Swap+1
				if Recipient_BG_pair_total[i]=='A':
					DO_BG_A_Swap=DO_BG_A_Swap+1
				if Recipient_BG_pair_total[i]=='B':
					DO_BG_B_Swap=DO_BG_B_Swap+1
				if Recipient_BG_pair_total[i]=='AB':
					DO_BG_AB_Swap=DO_BG_AB_Swap+1
		
		Dropout_pair.append(Num_Dropout_pair)
		#print ('Indx_R',indx_R)
		#print ('Indx_D',indx_D)
		Recipient_BG_pair_total=indx_R
		Donor_BG_pair_total=indx_D
		#print ('Recipient BG total, Donor BG total',Recipient_BG_pair_total,Donor_BG_pair_total)
		#DD_BG_total=DD_BG_total+generated_DD
		
		DDIC_Donor_BG_pair_total=DDIC_Donor_BG_pair_total+generated_don #Total number of recipients considered in the current round for swap allocation
		DDIC_Recipient_BG_pair_total=DDIC_Recipient_BG_pair_total+generated_rec #Total number of donors considered in the current round for swap allocation
		
		Compatible_edges_pair_DD=compatible_edges(DDIC_Recipient_BG_pair_total, DDIC_Donor_BG_pair_total)
		#print ('Compatible_edges_pair',Compatible_edges_pair_DD)
		Compatible_edges_DD=compatible_edges_DD(DDIC_Recipient_BG_pair_total, generated_DD)
		#print ('Compatible_edges_DD',Compatible_edges_DD)
		Compatible_edges=Compatible_edges_pair_DD+Compatible_edges_DD
		#print ('Compatible_edges',Compatible_edges)
		
		Matched_index_DD=maxWeightMatching(Compatible_edges, True)
		#print ('Matched_index_DD',Matched_index_DD)
				
		count_matched_DD=0
		for i in Matched_index_DD:
			if i!=-1 and i<=len(DDIC_Recipient_BG_pair_total)-1:
				count_matched_DD=count_matched_DD+1
				continue
		DDIC=count_matched_DD+2*Random_DD
		print ('DDIC sol',DDIC)
		DDIC_solution.append(DDIC)
		
		#print ('Recipient_BG_pair_total',DDIC_Recipient_BG_pair_total)
		
		for i in range(len(DDIC_Recipient_BG_pair_total)):
			if Matched_index_DD[i]!=-1:
				DDIC_Recipient_BG_pair_total[i]='Matched'
				DDIC_Donor_BG_pair_total[i]='Matched'
				continue
			
		#print ('Recipient_BG_pair_total after',DDIC_Recipient_BG_pair_total)
		
		indx_R=[]
		indx_D=[]
		Num_DD_Dropout_pair=0
		for i in range(len(DDIC_Recipient_BG_pair_total)):
			k=random.uniform(0,1)
			if DDIC_Recipient_BG_pair_total[i]!='Matched' and k>=p:
				indx_R.append(DDIC_Recipient_BG_pair_total[i])
				indx_D.append(DDIC_Donor_BG_pair_total[i])
				DD_Pair_waiting_time=DD_Pair_waiting_time+1
				if DDIC_Recipient_BG_pair_total[i]=='O':
					WT_BG_O_DDIC=WT_BG_O_DDIC+1
				if DDIC_Recipient_BG_pair_total[i]=='A':
					WT_BG_A_DDIC=WT_BG_A_DDIC+1
				if DDIC_Recipient_BG_pair_total[i]=='B':
					WT_BG_B_DDIC=WT_BG_B_DDIC+1
				if DDIC_Recipient_BG_pair_total[i]=='AB':
					WT_BG_AB_DDIC=WT_BG_AB_DDIC+1
			elif DDIC_Recipient_BG_pair_total[i]!='Matched' and k<=p:
				Num_DD_Dropout_pair=Num_DD_Dropout_pair+1
				if DDIC_Recipient_BG_pair_total[i]=='O':
					DO_BG_O_DDIC=DO_BG_O_DDIC+1
				if DDIC_Recipient_BG_pair_total[i]=='A':
					DO_BG_A_DDIC=DO_BG_A_DDIC+1
				if DDIC_Recipient_BG_pair_total[i]=='B':
					DO_BG_B_DDIC=DO_BG_B_DDIC+1
				if DDIC_Recipient_BG_pair_total[i]=='AB':
					DO_BG_AB_DDIC=DO_BG_AB_DDIC+1
		
		DD_Dropout_pair.append(Num_DD_Dropout_pair) #Number of dropouts added in each round
		#print ('Indx_R',indx_R)
		#print ('Indx_D',indx_D)
		DDIC_Recipient_BG_pair_total=indx_R
		DDIC_Donor_BG_pair_total=indx_D
			
		print ('Pair matched',count_matched_pair)
		print ('Pair and DD matched',count_matched_DD)
		print ('Additional pair matched',count_matched_DD-count_matched_pair)
		Additional_gain.append(count_matched_DD-count_matched_pair)
		
	
# =============================================================================
# 	print ('Total gain of DDIC over long run',sum(Additional_gain))
# 	print ('Total waiting time in swap registry',Pair_waiting_time)
# 	print ('Total waiting time in merged registry',DD_Pair_waiting_time)
# 	print ('Total number of dropout in swap registry',sum(Dropout_pair))
# 	print ('Total number of dropout in merged registry',sum(DD_Dropout_pair))
# 	
# 	print ('Waiting time for O type recipient in swap registry',WT_BG_O_Swap/(count_O))
# 	print ('Waiting time for A type recipient in swap registry',WT_BG_A_Swap/(count_A))
# 	print ('Waiting time for B type recipient in swap registry',WT_BG_B_Swap/(count_B))
# 	print ('Waiting time for AB type recipient in swap registry',WT_BG_AB_Swap/(count_AB))
# 	
# 	print ('Waiting time for O type recipient in DDIC registry',WT_BG_O_DDIC/(count_O))
# 	print ('Waiting time for A type recipient in DDIC registry',WT_BG_A_DDIC/(count_A))
# 	print ('Waiting time for B type recipient in DDIC registry',WT_BG_B_DDIC/(count_B))
# 	print ('Waiting time for AB type recipient in DDIC registry',WT_BG_AB_DDIC/(count_AB))
# 	
# 	print ('Total number of dropouts for O type recipient in swap registry',DO_BG_O_Swap)
# 	print ('Total number of dropouts for A type recipient in swap registry',DO_BG_A_Swap)
# 	print ('Total number of dropouts for B type recipient in swap registry',DO_BG_B_Swap)
# 	print ('Total number of dropouts for AB type recipient in swap registry',DO_BG_AB_Swap)
# 	
# 	print ('Total number of dropouts for O type recipient in DDIC registry',DO_BG_O_DDIC)
# 	print ('Total number of dropouts for A type recipient in DDIC registry',DO_BG_A_DDIC)
# 	print ('Total number of dropouts for B type recipient in DDIC registry',DO_BG_B_DDIC)
# 	print ('Total number of dropouts for AB type recipient in DDIC registry',DO_BG_AB_DDIC)
# 	
# 	plt.plot(Swap_solution, label='Current Process solution', linestyle='--')
# 	plt.plot(DDIC_solution, label='DDIC solution', linestyle='-')
# 	plt.plot(Dropout_pair, label='Dropouts in Current Process solution', linestyle='--')
# 	plt.plot(DD_Dropout_pair, label='Dropouts in DDIC solution')
# 	plt.xlabel('Number of rounds \n Arrival rate - PKE~U(20,25), DD~U(10,15), Dropout Probability=0.3')
# 	plt.ylabel('Number of recipients')
# 	plt.style.use('ggplot')
# 	#plt.title("Comparison of Current Process with DDIC")
# 	#plt.legend(loc='upper right')
# 	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
# 	           ncol=2, mode="expand", borderaxespad=0.)
# 	plt.show()
# =============================================================================
	 
	Rep_Additional_gain.append(Additional_gain)
	Rep_Swap_solution.append(Swap_solution)
	Rep_DDIC_solution.append(DDIC_solution)
	
	Rep_Dropout_pair.append(Dropout_pair)
	Rep_DD_Dropout_pair.append(DD_Dropout_pair)
	
	Rep_WT_BG_O_Swap.append(WT_BG_O_Swap/(count_O))
	Rep_WT_BG_A_Swap.append(WT_BG_A_Swap/(count_A))
	Rep_WT_BG_B_Swap.append(WT_BG_B_Swap/(count_B))
	Rep_WT_BG_AB_Swap.append(WT_BG_AB_Swap/(count_AB))
	
	Rep_WT_BG_O_DDIC.append(WT_BG_O_DDIC/(count_O))
	Rep_WT_BG_A_DDIC.append(WT_BG_A_DDIC/(count_A))
	Rep_WT_BG_B_DDIC.append(WT_BG_B_DDIC/(count_B))
	Rep_WT_BG_AB_DDIC.append(WT_BG_AB_DDIC/(count_AB))
	
	Rep_DO_BG_O_Swap.append(DO_BG_O_Swap)
	Rep_DO_BG_A_Swap.append(DO_BG_A_Swap)
	Rep_DO_BG_B_Swap.append(DO_BG_B_Swap)
	Rep_DO_BG_AB_Swap.append(DO_BG_AB_Swap)
	
	Rep_DO_BG_O_DDIC.append(DO_BG_O_DDIC)
	Rep_DO_BG_A_DDIC.append(DO_BG_A_DDIC)
	Rep_DO_BG_B_DDIC.append(DO_BG_B_DDIC)
	Rep_DO_BG_AB_DDIC.append(DO_BG_AB_DDIC)


print ('Avg Waiting time for O type recipient in swap registry',np.mean(Rep_WT_BG_O_Swap))
print ('Avg Waiting time for A type recipient in swap registry',np.mean(Rep_WT_BG_A_Swap))
print ('Avg Waiting time for B type recipient in swap registry',np.mean(Rep_WT_BG_B_Swap))
print ('Avg Waiting time for AB type recipient in swap registry',np.mean(Rep_WT_BG_AB_Swap))

print ('Avg Waiting time for O type recipient in DDIC registry',np.mean(Rep_WT_BG_O_DDIC))
print ('Avg Waiting time for A type recipient in DDIC registry',np.mean(Rep_WT_BG_A_DDIC))
print ('Avg Waiting time for B type recipient in DDIC registry',np.mean(Rep_WT_BG_B_DDIC))
print ('Avg Waiting time for AB type recipient in DDIC registry',np.mean(Rep_WT_BG_AB_DDIC))

print ('Avg Dropouts for O type recipient in swap registry',np.mean(Rep_DO_BG_O_Swap))
print ('Avg Dropouts for A type recipient in swap registry',np.mean(Rep_DO_BG_A_Swap))
print ('Avg Dropouts for B type recipient in swap registry',np.mean(Rep_DO_BG_B_Swap))
print ('Avg Dropouts for AB type recipient in swap registry',np.mean(Rep_DO_BG_AB_Swap))

print ('Avg Dropouts for O type recipient in DDIC registry',np.mean(Rep_DO_BG_O_DDIC))
print ('Avg Dropouts for A type recipient in DDIC registry',np.mean(Rep_DO_BG_A_DDIC))
print ('Avg Dropouts for B type recipient in DDIC registry',np.mean(Rep_DO_BG_B_DDIC))
print ('Avg Dropouts for AB type recipient in DDIC registry',np.mean(Rep_DO_BG_AB_DDIC))

#print ('Rep Swap Solution', Rep_Swap_solution)
#print ('Rep DDIC Solution', Rep_DDIC_solution)

#To find the averages of each round in number of transplants
Avg_swap_solution=[]
for l in range(n):
	a=[]
	for j in range(len(Rep_Swap_solution)):
		a.append(Rep_Swap_solution[j][l])
	Avg_swap_solution.append(np.mean(a))

#print ('Avg Swap Solution',Avg_swap_solution)

Avg_DDIC_solution=[]
for l in range(n):
	a=[]
	for j in range(len(Rep_DDIC_solution)):
		a.append(Rep_DDIC_solution[j][l])
	Avg_DDIC_solution.append(np.mean(a))

#print ('Avg DDIC Solution',Avg_DDIC_solution)

#To find the averages of each round in number of dropouts
Avg_Dropout_pair=[]
for l in range(n):
	a=[]
	for j in range(len(Rep_Dropout_pair)):
		a.append(Rep_Dropout_pair[j][l])
	Avg_Dropout_pair.append(np.mean(a))

#print ('Avg Dropout_pair',Avg_Dropout_pair)

Avg_DD_Dropout_pair=[]
for l in range(n):
	a=[]
	for j in range(len(Rep_DD_Dropout_pair)):
		a.append(Rep_DD_Dropout_pair[j][l])
	Avg_DD_Dropout_pair.append(np.mean(a))

#print ('Avg DD_Dropout_pair',Avg_DD_Dropout_pair)

plt.plot(Avg_swap_solution, label='Current Process solution', linestyle='--')
plt.plot(Avg_DDIC_solution, label='DDIC solution', linestyle='-')
plt.plot(Avg_Dropout_pair, label='Avg Dropouts in Current Process solution', linestyle='--')
plt.plot(Avg_DD_Dropout_pair, label='Avg Dropouts in DDIC solution')
plt.xlabel('Number of rounds \n Arrival rate - PKE~U(20,25), DD~U(10,15), Dropout Probability=0.3')
plt.ylabel('Number of recipients')
plt.style.use('ggplot')
#plt.title("Comparison of Current Process with DDIC")
#plt.legend(loc='upper right')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
	           ncol=2, mode="expand", borderaxespad=0.)
plt.show()












