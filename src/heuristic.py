# © Copyright IBM Corporation 2025. All Rights Reserved.
# LICENSE: Eclipse Public License - v 2.0, https://opensource.org/licenses/EPL-2.0
# SPDX-License-Identifier: EPL-2.0


import array as arr
import random


def find_edge(n, wt, elist):
    pos = 0
    e = len(wt)-1
    if e < 0:
        return e
    
    sum = 0.0
    for i in range(len(wt)):
        sum = sum + wt[i]
        
    randval = (sum-1e-6)*random.random()
    
    sum = 0.0
    for i in range(len(wt)):
        sum = sum + wt[i]
        if sum >= randval:
            e = i
            break

    te = e
    while elist[2*te] == n-1 or elist[2*te+1] == n-1:
        te = te - 1
        if te < 0:
            break
    if te < 0:
        te = e
        while elist[2*te] == n-1 or elist[2*te+1] == n-1:
            te = te + 1
            if te >= len(wt):
                break
    return te

def create_undir_edges_from_par(cnodes, cparents, cweight, elist, bdir, wt, inwt):
    edgemap = {}
    for i in range(len(inwt)):
        inwt[i] = 0.0

    for i in range(len(cnodes)):
        if cweight[i] <= 0.001:
            continue
        
        cset = cnodes[i]
        parent = cparents[i]

        if len(cset) > 1:
            for j in cset:
                bdir.append(j)
                
        for j in range(len(cset)):
            node = cset[j]
            par = parent[j]
            if node in par:
                print('mayday', cweight[i], i, cset, par)
                return []
                exit(-1)
            inwt[node] = inwt[node] + cweight[i]
            
            for k in par:
                n1 = k
                n2 = node
                
                if k > node:
                    n1 = node
                    n2 = k
                    
                nodepair = n1*1000 + n2
                if len(edgemap) == 0 or nodepair not in edgemap:
                    edgemap[nodepair] = cweight[i]
                else:
                    edgemap[nodepair] = edgemap[nodepair] + cweight[i]

    for i in edgemap:
        wt.append(edgemap[i])
        n1 = int(i / 1000)
        n2 = int(i % 1000)

        elist.append(n1)
        elist.append(n2)


def create_edges_from_par(cnodes, cparents, cweight, elist, bdir, bdirwt, wt, inwt):
    edgemap = {}
    bdirmap = {}
    
    for i in range(len(inwt)):
        inwt[i] = 0.0

    for i in range(len(cnodes)):
        if cweight[i] <= 0.001:
            continue
        
        cset = cnodes[i]
        parent = cparents[i]

        if len(cset) > 1:
            n1 = cset[0]
            n2 = cset[1]
            if n1 > n2:
                temp = n2
                n2 = n1
                n1 = temp
                    
            nodepair = n1*1000 + n2
            if nodepair not in bdirmap:
                bdirmap[nodepair] = cweight[i]
            else:
                bdirmap[nodepair] = bdirmap[nodepair] + cweight[i]
                
        for j in range(len(cset)):
            node = cset[j]
            par = parent[j]
            if node in par:
                print('mayday', cweight[i], i, cset, par)
                return []
                exit(-1)
            inwt[node] = inwt[node] + cweight[i]
            
            for k in par:

                nodepair = k*1000 + node 
                if nodepair not in edgemap:
                    edgemap[nodepair] = cweight[i]
                else:
                    edgemap[nodepair] = edgemap[nodepair] + cweight[i]

    for i in bdirmap:
        bdirwt.append(bdirmap[i])
        n1 = int(i / 1000)
        n2 = int(i % 1000)

        bdir.append(n1)
        bdir.append(n2)
        
    for i in edgemap:
        wt.append(edgemap[i])
        n1 = int(i / 1000)
        n2 = int(i % 1000)

        elist.append(n1)
        elist.append(n2)

    
def contract_edge (e, wt, elist, nodelist, cnodes, cparents, cweight, containsbdir):
    node1 = elist[2*e]
    node2 = elist[2*e+1]

    if node1 == node2:
        print('mayday2')
        exit(-1)
    
    for i in nodelist[node1]:
        nodelist[node2].append(i)
    nodelist[node1] = []
    if containsbdir[node1] == 1:
        containsbdir[node2] = 1

    for i in range(len(cnodes)):
        cset = cnodes[i]
        parent = cparents[i]

        deln = [0 for j in range(len(cset))]
        ndel = 0
        for j in range(len(cset)):
            
            node = cset[j]
            par = parent[j]
            if node == node2 and node1 in par:
                deln[j] = 1
                ndel = ndel + 1
            else:
               if node == node1:
                   if node2 in par:
                       deln[j] = 1
                   else:
                       cset[j] = node2
               else:       
                   if node1 in par:
                       ix = par.index(node1)
                       if node2 not in par:
                           par[ix] = node2
                       else:
                           par.remove(node1)

        
        if ndel > 0:
            if ndel == len(cset):
                cweight[i] = 0.0
            else:
                newcset = []
                newpar = []
                for j in range(len(cset)):
                    if deln[j] == 0:
                        newcset.append(cset[j])
                        newpar.append(parent[j])
                cnodes[i] = newcset
                cparents[i] = newpar
                del deln
    
def contract_bdir_edge (node1, node2, wt, elist, nodelist, cnodes, cparents, cweight, containsbdir):
    if node1 == node2:
        print('mayday2')
        exit(-1)

    for i in nodelist[node1]:
        nodelist[node2].append(i)
    nodelist[node1] = []
    containsbdir[node1] = 0
    containsbdir[node2] = 1

    for i in range(len(cnodes)):
        cset = cnodes[i]
        parent = cparents[i]

        if len(cset) == 1:
            if node1 in cset or node2 in cset:
                cweight[i] = 0.0
                continue
        if len(cset) > 1:
            if node1 in cset and node2 not in cset:
                cweight[i] = 0.0
                continue
            if node2 in cset and node1 not in cset:
                cweight[i] = 0.0
                continue
            if node1 in cset and node2 in cset:
                cnodes[i] = [node2]
                cparents[i] = [list(set(parent[0]).union(parent[1]))]
                continue
                                    
    for i in range(len(cnodes)):
        cset = cnodes[i]
        parent = cparents[i]

        if node1 in cset or node2 in cset:
            continue
        for j in range(len(cset)):
            node = cset[j]
            par = parent[j]
            if node1 in par:
                ix = par.index(node1)
                if node2 not in par:
                    par[ix] = node2
                else:
                    par.remove(node1)

def single_contraction_round(n, cnodes, cparents, cweight, allcyc, cutwt, dbg):
    elist = arr.array('i',[]) #array to store edges
    bdir = arr.array('i',[]) #array to store bidirected edges
    wt = arr.array('d',[]) #array to store edges
    inwt = arr.array('d', [0.0 for i in range(n)])
    nodelist = [[i] for i in range(n)]
    containsbdir = [0 for i in range(n)]

    create_undir_edges_from_par(cnodes, cparents, cweight, elist, bdir, wt, inwt)
    over = n-5
    while over >= 0:
        e = find_edge(n, wt, elist)
        
        if e < 0 or e >= len(wt):
            over = -1
            break
        if dbg > 0:
            print('e =', e, 'len =', len(wt), 'n1,n2=', elist[2*e], elist[2*e+1])
            for i in range(len(cnodes)):
                print(cnodes[i], cparents[i], cweight[i])
        contract_edge (e, wt, elist, nodelist, cnodes, cparents, cweight, containsbdir)
        
        del wt
        del elist
        del bdir
        wt = arr.array('d',[]) #array to store edges
        elist = arr.array('i',[]) #array to store edges
        bdir = arr.array('i',[]) #array to store bidirected edges

        create_undir_edges_from_par(cnodes, cparents, cweight, elist, bdir, wt, inwt)
        
        for i in range(len(inwt)):
            if inwt[i] < 0.98*cutwt and len(nodelist[i]) > 1 and nodelist[i] not in allcyc:
                allcyc.append(nodelist[i][:])

        if dbg > 0:        
            print (nodelist)
            print (inwt)
        over = over - 1

def contract_heur(n, icnodes, icparents, icweight, cutwt):
    #first locate edges with weight>= 0.99 and check if almost directed cycle exists
    wt = arr.array('d',[]) #array to store edges
    elist = arr.array('i',[]) #array to store edges
    bdir = arr.array('i',[]) #array to store bidirected edges
    bdirwt = arr.array('d',[]) #array to store bidirected edges
    allcyc = []
    inwt = arr.array('d', [0 for i in range(n)])
    cnodes = []
    cparents = []
    cweight = []
    
    # copy all input arrays as we will modify them
    for i in range(len(icnodes)):
        if icweight[i] > 0.001:
            cnodes.append(icnodes[i][:])
            cparents.append(icparents[i][:])
            cweight.append(icweight[i])
    
    nccomp = len(cnodes)

    create_edges_from_par(cnodes, cparents, cweight, elist, bdir, bdirwt, wt, inwt)
    
    
    for numrep in range(20):
        dbg = 0
        # now modify cnodes, cparents, cweight to find regular cluster inequalities
        tnodes = []
        tparents = []
        tweight = arr.array('d', [])

        # now copy all parent set info as it is used destructively inside single_contraction round
        for i in range(len(cnodes)):
            cset = []
            par = []
            for j in range(len(cnodes[i])):
                cset.append(cnodes[i][j])

            for j in range(len(cparents[i])):
                par.append(cparents[i][j][:])

            tnodes.append(cset)
            tparents.append(par)
            tweight.append(cweight[i])

        if dbg > 0:
            for i in range(len(tnodes)):
                print(cnodes[i], cparents[i], cweight[i])

        single_contraction_round (n, tnodes, tparents, tweight, allcyc, cutwt, dbg)
        
        if len(allcyc) > 10:
            break
    #end repeat
    
    return allcyc
#end of function

def contract_heur_bdir(n, icnodes, icparents, icweight):
    #first locate edges with weight>= 0.99 and check if almost directed cycle exists
    wt = arr.array('d',[]) #array to store edges
    elist = arr.array('i',[]) #array to store edges
    bdir = arr.array('i',[]) #array to store bidirected edges
    bdirwt = arr.array('d',[]) #array to store bidirected edges
    inwt = arr.array('d', [0 for i in range(n)])
    cnodes = []
    cparents = []
    cweight = []
    allcyc = []
    
    # copy all input arrays as we will modify them
    for i in range(len(icnodes)):
        if icweight[i] > 0.001:
            cnodes.append(icnodes[i][:])
            cparents.append(icparents[i][:])
            cweight.append(icweight[i])
    
    nccomp = len(cnodes)
    
    create_edges_from_par(cnodes, cparents, cweight, elist, bdir, bdirwt, wt, inwt)
    
    for b in range(int(len(bdir)/2)):

        # now modify cnodes, cparents, cweight to find regular cluster inequalities
        tnodes = []
        tparents = []
        tweight = arr.array('d', [])
        nodelist = [[i] for i in range(n)]
        containsbdir = [0 for i in range(n)]

        for i in range(len(cnodes)):
            cset = []
            par = []
            for j in range(len(cnodes[i])):
                cset.append(cnodes[i][j])

            for j in range(len(cparents[i])):
                par.append(cparents[i][j][:])

            tnodes.append(cset)
            tparents.append(par)
            tweight.append(cweight[i])


        contract_bdir_edge (bdir[2*b], bdir[2*b+1], wt, elist, nodelist, tnodes, tparents, tweight, containsbdir)

        newcyc = contract_heur (n, tnodes, tparents, tweight, bdirwt[b])
        if len(newcyc) > 0:
            for cyc in newcyc:
                if bdir[2*b+1] in cyc:
                    cyc1 = cyc
                    cyc1.remove(bdir[2*b+1])
                    cyc1.sort()
                    newc = [bdir[2*b], bdir[2*b+1]] + cyc

                    if not newc in allcyc:
                        allcyc.append(newc)
        
    return allcyc
#end of function

