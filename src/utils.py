# © Copyright IBM Corporation 2025. All Rights Reserved.
# LICENSE: Eclipse Public License - v 2.0, https://opensource.org/licenses/EPL-2.0
# SPDX-License-Identifier: EPL-2.0


import numpy as np
import pickle
import math
import sys


#latent variable graph
def compute_BiC_c_comps_sets(D, c_components, parent_sets, bidirect_tuple):
    # compute the bic score
    nVars = D.shape[1]
    nSamples = D.shape[0]

    # find all c components in graph
    isParent = np.zeros(nVars)

    scores = np.zeros( (1, len(c_components)))

    covMat = np.cov(D.T)
    tol = 10e-2

    # convert to mag
    [compMag, district, parents] = find_component_mag_sets(c_components, nVars, parent_sets, bidirect_tuple)

    district = np.where(district)[0]

    scores = compute_local_BiC_comp_sets(compMag, c_components, district, len(c_components), covMat, nSamples, tol)

    score = np.sum(scores)

    return score

def compute_local_BiC_comp_sets(compMag, component, district, compSize, CovMat, nSamples,  tol):

    if district.size == 1:
        comCovMat = CovMat[district, district]
        sign, logdet = 1, np.log(2 * math.pi * comCovMat)
        sc = sign * logdet + (nSamples - 1) / nSamples
    else:
        nVar = CovMat.shape[0]
        curCovMat = CovMat[np.ix_(district, district)]
        compMag = compMag[np.ix_(district, district)]
        try:
            _, _, curHatCovMat, _, converged = RICF_fit(compMag, curCovMat, tol)
            if converged == False:
                print('RICF fit did not converge')
                return -1.0*pow(10,20)
        except np.linalg.LinAlgError as err:
            print('matrix singular error')
            return -1.0*pow(10,20)

        remParents = np.zeros(nVar)
        remParents[district] = 1
        remParents[list(component)] = 0
        parInds = np.where(remParents[district].astype(int))[0]

        sign, logdet = np.linalg.slogdet(curHatCovMat)
        if sign < 0.0:
            print('score error: sign < 0.0')
            return -1.0*pow(10,20)
        logdet_signed = sign* logdet
        if logdet_signed < -4.0:
            print('score error: logdet_signed < -4.0')
            return -1.0*pow(10,20)
        if np.sum(remParents) > 0:
            l1 = compSize * math.log(2 * math.pi)
            l2 = logdet_signed - math.log(np.prod(np.diag(curHatCovMat[np.ix_(parInds, parInds)])))
            l3 = (nSamples - 1) / nSamples * (np.trace( matrix_division_left(curHatCovMat,curCovMat)) - np.sum(remParents))
            sc = l1 + l2 +l3
        else:
            l1 = compSize * math.log(2 * math.pi)
            l2 = logdet_signed
            l3 = (nSamples - 1) / nSamples * np.trace(matrix_division_left(curHatCovMat,curCovMat))
            sc = l1 + l2 + l3

    nVars  = np.size(component)
    nEdges = np.where(compMag)[0].size/2
    bp = math.log(nSamples)/2 * (2 * nVars + nEdges)
    sc = - sc * nSamples/2 - bp

    return sc


def find_component_mag_sets(component, nVars, parent_sets, bidirect_tuple):
    compMag = np.zeros((nVars,nVars))
    if len(component) > 1:
        for edge_ind, (parent, child) in enumerate(bidirect_tuple):
            compMag[parent, child] = 2
            compMag[child, parent] = 2

    district, parents = np.zeros(( nVars,)), np.zeros((nVars,))
    district[list(component)] = True
    for index, iVar in enumerate(component):
        iParents = list(parent_sets[index])
        compMag[iParents, iVar] = 2
        compMag[iVar, iParents] = 3
        district[iParents] = True
        parents[iParents] = True

    return compMag, district, parents


def RICF_fit(smm, covMat, tol):

    converged = True
    if (smm==4).any():
        sys.exit('Graph includes bows\n')

    nVars = covMat.shape[1]
    dg = np.multiply(smm == 2, smm.T == 3)
    bg = np.multiply(smm == 2, smm.T == 2)

    # starting values
    omega = np.diag(np.diag(covMat),0)
    beta = np.eye(nVars)

    # list of parents, spouses etc
    par, sp = np.zeros((nVars, nVars)), np.zeros((nVars, nVars))
    for iVar in range(nVars):
        par[iVar, dg[:, iVar]] = True
        sp[iVar, bg[:, iVar]] = True
    iter=0

    ricf = {}
    while True:
        iter = iter+1
        ricf[iter] = {}
        ricf[iter]['omega'] = omega.copy()
        ricf[iter]['beta'] = beta.copy()

        omega = omega.copy()
        beta = beta.copy()
        # for each variable
        for iVar in range(nVars):

           vcomp = list(range(0,iVar)) + list(range(iVar+1, nVars))
           iPar= np.where(par[iVar,:])[0]
           iSp = np.where(sp[iVar, :])[0]

           if iSp.size == 0:
               if iPar.size > 0:
                   if iter==1:
                       if len(vcomp) == 1:
                           inv_mat = 1 /  (covMat[iPar, iPar])
                       else:
                           inv_mat =   np.linalg.inv(covMat[np.ix_(iPar, iPar)])
                       beta[iVar, iPar] = np.matmul(-covMat[iVar, iPar],inv_mat)
                       omega[iVar, iVar] = covMat[iVar, iVar] + np.matmul(beta[iVar, iPar],covMat[iPar, iVar])

           elif iPar.size > 0: # spouses and parents
               oInv = np.zeros((nVars, nVars))
               if len(vcomp) == 1:
                   oInv[vcomp, vcomp] = 1 / oInv[vcomp, vcomp]
               else:
                   oInv[np.ix_(vcomp, vcomp)] = np.linalg.inv(omega[np.ix_(vcomp, vcomp)])
               Z = np.matmul(oInv[np.ix_(iSp, vcomp)], beta[vcomp, :])
               if Z.ndim == 1:
                   Z = Z.reshape(1, -1)
               nPar = iPar.size
               nSp = iSp.size
               range1 = list(range(nPar))
               range2= list(range(nPar, nPar+nSp))
               # % XX
               XX = np.zeros( (nPar+nSp,(nPar+nSp)) )
               XX[:] = np.nan
               #% Upper left quadrant
               XX[np.ix_(range1, range1)] = covMat[np.ix_(iPar, iPar)]
               #% Upper right quadrant
               XX[np.ix_(range1, range2)] = np.matmul(covMat[iPar, :], Z.T)
               #% Lower left quadrant
               XX[np.ix_(range2, range1)] = XX[np.ix_(range(nPar), range(nPar,nPar+nSp))].T
               #% Lower right quadrant
               XX[np.ix_(range2, range2)] = np.matmul(np.matmul(Z,covMat),Z.T)
               YX = np.hstack((covMat[iVar, iPar].reshape(1,-1), np.array([np.matmul(covMat[iVar, :],Z.T)]).reshape(1,-1)
                                    )).reshape((-1,1))

               if YX.size == 1:
                    temp = YX / XX
                    beta[iVar, iPar] = -temp[range1]
                    omega[iVar, iSp] = temp[range2]
                    omega[iSp, iVar] = omega[iVar, iSp]
                    tempVar = covMat[iVar, iVar] - temp * YX
                    omega[iVar, iVar] = tempVar + np.matmul(omega[iVar, iSp]/ omega[np.ix_(iSp, iSp)] \
                                                            , omega[iSp, iVar])
               else:
                   temp = matrix_division_right(YX.T, XX)
                   # % update beta, omega
                   beta[iVar, iPar] = -temp[0,range1]
                   omega[iVar, iSp] = temp[0, range2]
                   omega[iSp, iVar] = omega[iVar, iSp]

                   tempVar = covMat[iVar, iVar] - np.matmul(temp, YX)

                   if iSp.size == 1:
                       omega[iVar, iVar] = tempVar + omega[iVar, iSp]/ omega[np.ix_(iSp, iSp)] * omega[iSp, iVar]
                   else:
                        omega[iVar, iVar] = tempVar + np.matmul(matrix_division_right(omega[iVar, iSp], omega[np.ix_(iSp, iSp)])\
                                       ,omega[iSp, iVar])

           else:
               oInv = np.zeros((nVars, nVars))
               if len(vcomp) == 1:
                    oInv[vcomp, vcomp] = 1/omega[vcomp, vcomp]
               else:
                    oInv[np.ix_(vcomp, vcomp)] = np.linalg.inv(omega[ np.ix_(vcomp, vcomp)])
               Z = np.matmul(oInv[np.ix_(iSp, vcomp)],beta[vcomp, :])
               XX = np.matmul( np.matmul(Z,covMat), Z.T)
               YX = np.matmul(covMat[iVar, :], Z.T).T

               if YX.size== 1:
                   omega[iVar, iSp] =  YX / XX
                   omega[iSp, iVar] = omega[iVar, iSp]
                   tempVar = covMat[iVar, iVar] -  omega[[iVar],[iSp]] * YX
                   omega[iVar, iVar] = tempVar + \
                                       np.matmul(np.matmul(omega[[iVar], [iSp]], oInv[np.ix_(iSp,iSp)]), omega[iSp, [iVar]])

               else:
                   omega[iVar, iSp] =  matrix_division_right(YX.T, XX)
                   omega[iSp, iVar] = omega[iVar, iSp]
                   tempVar = covMat[iVar, iVar] - np.matmul(omega[iVar,iSp],YX)
                   omega[iVar, iVar] = tempVar+ \
                                       np.matmul(np.matmul(omega[iVar, iSp],oInv[np.ix_(iSp,iSp)]),omega[iSp, iVar])

        if (sum(sum(abs(ricf[iter]['omega']-omega))) + sum(sum(abs(ricf[iter]['beta']-beta))) < tol):
           break

        elif iter>5000:
           print('RICF fit did not converge')
           print('ricf omega: '+str(ricf[iter]['omega']-omega))
           print('ricf beta: '+str(ricf[iter]['beta']-beta))
           converged = False
           break

    hatCovMat = np.matmul( np.matmul(np.linalg.inv(beta), omega), np.linalg.inv(beta.T))

    return beta, omega, hatCovMat, ricf, converged


def matrix_division_right(YX, XX):
    # solve for YX/XX or x XX = YX

    x = np.matmul(np.matmul(YX, XX.T), np.linalg.inv(np.matmul(XX, XX.T)))

    return x

def matrix_division_left(B, b):
    # solves for B\b or Bx = b
    x, resid, rank, s = np.linalg.lstsq(B,b)

    return x

def score_write_to_file(scores, file_name):
    # write scores to find
    with open(file_name, 'wb') as handle:
        pickle.dump(scores, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return

def check_connected_component(set, node_ids):
    # set is a list of tuples, which consists a tuple of edge tupule
    new_set = []

    # single set
    if len(set)  < 2:
        return set

    for edge_lists in set:
        # check connected
        flag = True
        for node in node_ids:
            if not (any(node in i for i in edge_lists)):
                flag = False
        if flag:
            new_set.append(edge_lists)

    return new_set

def find_bi_connected_node(iVar, comp_edges):
    # find nodes that bi-direct connect to iVar in comp edges

    nodes = []

    if len(comp_edges)  < 2:
        return nodes

    for edge_lists in comp_edges:
        if iVar in  edge_lists:
            bi_node = set( list(edge_lists)) - set([iVar])
            nodes.append( list(bi_node)[0])

    return nodes

def has_cycle_in_c_comp(parent_set_config, vars_in_comp):
    # check if parent set has a directed cycle of length 2 only inside a c-comp
    num_var_parent = len(parent_set_config)
    vars_in_comp_list = list(vars_in_comp)
    flag = False

    for index_var, iVar in enumerate(vars_in_comp_list):
        iParents = list(parent_set_config[index_var])
        for each_parent in iParents:
            if each_parent in vars_in_comp_list:
                each_parent_index = vars_in_comp_list.index(each_parent)
                par_parents = list(parent_set_config[each_parent_index])
                if iVar in par_parents:
                    flag = True
                    break

        if flag == True:
            break

    return flag
