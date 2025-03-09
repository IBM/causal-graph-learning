# © Copyright IBM Corporation 2025. All Rights Reserved.
# LICENSE: Eclipse Public License - v 2.0, https://opensource.org/licenses/EPL-2.0
# SPDX-License-Identifier: EPL-2.0


import array as arr

def dircyc (n , ne, elist):
    ideg = [0 for i in range(n)]
    odeg = [0 for i in range(n)] #outdegree
    oadjn = [0 for i in range(n*n)] #out adjacency list
    mark = [0 for i in range(n)] #processed nodes
    active = arr.array('i', []) #same as stack
    nchild = arr.array('i', []) #number of processed children
    isactive = [0 for i in range(n)] #idicates if node is on the stack

    nproc = 0
    allcyc = []

    for i in range(ne):
        n1 = elist[2*i]
        n2 = elist[2*i+1]
        oadjn[n1*n + odeg[n1]] = n2
        odeg[n1] += 1
        ideg[n2] += 1

    while nproc < n: 
        # Find an unmarked node to start from, preferably zero indeg node*/
        start = -1
        pstart = -1
        for i in range(n):
            if mark[i] == 0 and pstart == -1:
                pstart = i
                if ideg[i] == 0:
                    start = i
                    break
        if start == -1 and pstart >= 0:
            start = pstart
        
        #print('start =', start, ', pstart=',  pstart)
        if start == -1:
            break

        active.append(start)
        nchild.append(0)
        mark[start] = 1
        isactive[start] = 1
        
        while len(active) > 0:
            cnode = active[len(active)-1]
            nchildproc = nchild[len(nchild)-1]
            #print('cnode=',cnode,',nchildproc=',nchildproc, 'act=',len(active), len(nchild))

            if nchildproc == odeg[cnode]:
                nproc = nproc + 1
                isactive[cnode] = 0
                active.pop()
                nchild.pop()
                if len(nchild) > 0:
                    nchild[len(nchild)-1] = nchild[len(nchild)-1] + 1

            else:
                n1 = oadjn[cnode*n + nchildproc]
                #print(n1,'mark=',mark[n1])
                if mark[n1] == 0:
                    active.append(n1)
                    nchild.append(0)
                    mark[n1] = 1
                    isactive[n1] = 1
                else:
                    #if n1 in active:
                    if isactive[n1] != 0:
                        # Found previously marked node trace back to cycle
                        pos = len(active)-1
                        cur = active[pos]
                        cyc = [cnode]

                        while cur != n1 and pos > 0:
                            pos = pos -1
                            cur = active[pos]
                            cyc.append(cur)
                        if (cur == n1):
                            #print('cycle=',cyc)
                            allcyc.append(cyc[::-1])
                        else:
                            print('error in computing cycle')
                            exit(-1)

                    # at this point we have a marked child
                    nchild[len(nchild)-1] = nchild[len(nchild)-1] + 1
                    
    return allcyc

def almostdircyc (n , ne, elist, bnode, enode):
    ideg = [0 for i in range(n)]
    odeg = [0 for i in range(n)] #outdegree
    oadjn = [0 for i in range(n*n)] #out adjacency list
    mark = [0 for i in range(n)] #processed nodes
    active = arr.array('i', []) #same as stack
    nchild = arr.array('i', []) #number of processed children
    isactive = [0 for i in range(n)] #idicates if node is on the stack

    nproc = 0
    allcyc = []

    for i in range(ne):
        n1 = elist[2*i]
        n2 = elist[2*i+1]
        oadjn[n1*n + odeg[n1]] = n2
        odeg[n1] += 1
        ideg[n2] += 1


    mark[bnode] = 1
    active.append(enode)
    nchild.append(0)
    mark[enode] = 1
    isactive[enode] = 1
                
    while len(active) > 0:
        cnode = active[len(active)-1]
        nchildproc = nchild[len(nchild)-1]
        #print('cnode=',cnode,',nchildproc=',nchildproc, 'act=',len(active), len(nchild))

        if nchildproc == odeg[cnode]:
            nproc = nproc + 1
            isactive[cnode] = 0
            active.pop()
            nchild.pop()
            if len(nchild) > 0:
                nchild[len(nchild)-1] = nchild[len(nchild)-1] + 1

        else:
            n1 = oadjn[cnode*n + nchildproc]
            #print(n1,'mark=',mark[n1])
            if mark[n1] == 0:
                active.append(n1)
                nchild.append(0)
                mark[n1] = 1
                isactive[n1] = 1
            else:
                if  n1 == bnode:
                    # Found previously marked node trace back to cycle
                    pos = len(active)-1
                    cur = active[pos]
                    cyc = [cnode]

                    while cur != enode and pos > 0:
                        pos = pos -1
                        cur = active[pos]
                        cyc.append(cur)
                    if (cur == enode):
                        #print('cycle=',cyc)
                        cyc.reverse()
                        cyc.append(bnode)
                        allcyc.append(cyc)
                    else:
                        print('error in computing cycle')
                        exit(-1)
                        
                # at this point we have a marked child
                nchild[len(nchild)-1] = nchild[len(nchild)-1] + 1
                    
    return allcyc
