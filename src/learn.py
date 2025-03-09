# © Copyright IBM Corporation 2025. All Rights Reserved.
# LICENSE: Eclipse Public License - v 2.0, https://opensource.org/licenses/EPL-2.0
# SPDX-License-Identifier: EPL-2.0


# To run, type:
# python3 learn.py -s score_file_name.pkl -i path_input_directory/ -o path_output_directory/ -h heuristics -c cuts -t MIP_time_limit_(seconds) -g type_of_graph_desired > log_file.txt
# For example:
# python3 learn.py -s score_sample_graph_random_bowfree_10nodes-0-cs2_sps3_ocps2.pkl -i ../data/test/ -o ../results/test/ -t 7200 -g bowfree > log_file.txt


import pickle
import copy
import time
from dircycle import dircyc,almostdircyc
from heuristic import contract_heur,contract_heur_bdir
import cplex
from utils import compute_BiC_c_comps_sets
import sys
from generate_random_graph import write_graph_to_file
from test_for_arid import test_for_arid
import random
import itertools

class BNSLlvInst:
        def __init__(self,instance,datadir,resultsdir,heuristics,cuts,milp_time_limit,bowfree,arid):
                self.instance = instance
                self.datadir = datadir
                self.resultsdir = resultsdir
                self.heuristics = heuristics
                self.cuts = cuts
                self.milp_time_limit = milp_time_limit
                self.bowfree = bowfree
                self.arid = arid
                self.originalScores = None
                self.temp1Scores = None
                self.temp2Scores = None
                self.scores = None
                self.data = None
                self.V = None
                self.cComps = None
                self.iComps = None
                self.dPars = None
                self.iPars = None
                self.biPars = None
                self.lp = None
                self.z = None
                self.ind = None
                self.indInv = None
                self.udE = None
                self.bi = None
                self.x = None
                self.clusterIP = None
                self.ConflictNodes = None
                self.ConflictEdges = None
                self.dag = None
                self.milp = None

        def readFromPkl(self):
                filename = self.datadir+self.instance
                file = open(filename, 'rb')
                [self.data,self.originalScores] = pickle.load(file)
                file.close()
                
        def readFromDag(self):
                filename = self.datadir+self.instance
                file = open(filename, 'r')
                x = file.readlines()
                file.close
                self.originalScores = {}
                Vsize = int(x[0].strip('\n'))
                lineind = 1
                for k in range(Vsize):
                        line = x[lineind]
                        z = line.strip('\n').split(' ')
                        i = int(z[0])
                        nPaSets = int(z[1])
                        self.originalScores[(i,),()] = {}
                        for p in range(nPaSets):
                                line = x[lineind+1+p]
                                z = line.strip('\n').split(' ')
                                PaSize = int(z[1])
                                PaSet = []
                                for t in range(PaSize):
                                        PaSet.append(int(z[2+t]))
                                self.originalScores[(i,),()][(tuple(PaSet),)] = float(z[0])
                        lineind = lineind+1+nPaSets
                if len(x) > lineind:
                        lineind = lineind+1
                        csize = int(x[lineind].strip('\n'))
                        for k in range(csize):
                                lineind = lineind+1
                                line = x[lineind]
                                z = line.strip('\n').split()
                                nNodes = int(z[0])
                                score = float(z[1])
                                cComp = []
                                for k in range(nNodes):
                                        cComp.append(int(z[2+k]))
                                cComp = tuple(cComp)
                                # Only consider size-2 c-components now...
                                biEdges = (cComp,)
                                if cComp not in self.originalScores.keys():
                                        self.originalScores[cComp,biEdges] = {}
                                parsets = []
                                for p in range(len(cComp)):
                                        lineind = lineind+1
                                        z = x[lineind].strip('\n').split()
                                        parset = [int(z[i]) for i in range(2,len(z))]
                                        parsets.append(tuple(parset))
                                self.originalScores[cComp,biEdges][tuple(parsets)] = score
                
        def Prune_scores(self,prune_more=False):
                t0P = time.time()
                sum1 = 0
                sum2 = 0
                counter1 = 0
                counter2 = 0
                Dpars = {}
                for D in self.originalScores.keys():
                        if len(D[0]) == 1:
                                Dpars[D] = list(self.originalScores[D].keys())
                                Dparscopy = Dpars[D].copy()
                                for ind in range(len(Dpars[D])):
                                        for i in range(ind+1,len(self.originalScores[D])):
                                                if set(Dpars[D][ind][0]).difference(set(Dpars[D][i][0])) == set() and self.originalScores[D][Dpars[D][ind]] >= self.originalScores[D][Dpars[D][i]]:
                                                        if Dpars[D][i] in Dparscopy:
                                                                Dparscopy.remove(Dpars[D][i])
                                sum1 = sum1+len(Dpars[D])
                                Dpars[D] = Dparscopy.copy()
                                sum2 = sum2+len(Dpars[D])

                        # size 2 c-component
                        elif len(D[0]) == 2:
                                Dpars[D] = list(self.originalScores[D].keys())
                                Dparscopy = Dpars[D].copy()
                                for ind in range(len(Dpars[D])):
                                        # a = D[0][0], b = D[0][1]
                                        # a<->b vs a<->b
                                        delInd = False
                                        for i in range(ind):
                                                if set(Dpars[D][i][0]).difference(set(Dpars[D][ind][0])) == set() and set(Dpars[D][i][1]).difference(set(Dpars[D][ind][1])) == set() and self.originalScores[D][Dpars[D][ind]] <= self.originalScores[D][Dpars[D][i]]:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        delInd = True
                                                        break
                                        if delInd == True:
                                                continue
                                        # a<->b vs a,b
                                        Da = ((D[0][0],),())
                                        Db = ((D[0][1],),())
                                        DaPars = list(self.originalScores[Da].keys())
                                        DbPars = list(self.originalScores[Db].keys())
                                        maxa = -float('inf')
                                        maxb = -float('inf')
                                        for ind1 in range(len(DaPars)):
                                                if set(DaPars[ind1][0]).difference(set(Dpars[D][ind][0])) == set() and self.originalScores[Da][DaPars[ind1]] > maxa:
                                                        maxa = self.originalScores[Da][DaPars[ind1]]
                                        for ind2 in range(len(DbPars)):
                                                if set(DbPars[ind2][0]).difference(set(Dpars[D][ind][1])) == set() and self.originalScores[Db][DbPars[ind2]] > maxb:
                                                        maxb = self.originalScores[Db][DbPars[ind2]]
                                        if self.originalScores[D][Dpars[D][ind]] <= maxa+maxb:
                                                Dparscopy.remove(Dpars[D][ind])
                                                counter1 += 1
                                                continue
                                        tol = 1e-10
                                        # a<->b vs a->b
                                        Da = ((D[0][0],),())
                                        Db = ((D[0][1],),())
                                        DbPars = list(self.originalScores[Db].keys())
                                        maxa = self.originalScores[Da][((),)]
                                        maxb = -float('inf')
                                        for ind2 in range(len(DbPars)):
                                                if set(DbPars[ind2][0]).difference(set(Dpars[D][ind][1])) == {D[0][0]} and self.originalScores[Db][DbPars[ind2]] > maxb:
                                                        maxb = self.originalScores[Db][DbPars[ind2]]
                                        if self.originalScores[D][Dpars[D][ind]] <= maxa+maxb+tol:
                                                Dparscopy.remove(Dpars[D][ind])
                                                counter2 += 1
                                                continue
                                        # a<->b vs a<-b
                                        Da = ((D[0][0],),())
                                        Db = ((D[0][1],),())
                                        DaPars = list(self.originalScores[Da].keys())
                                        maxa = -float('inf')
                                        maxb = self.originalScores[Db][((),)]
                                        for ind1 in range(len(DaPars)):
                                                if set(DaPars[ind1][0]).difference(set(Dpars[D][ind][0])) == {D[0][1]} and self.originalScores[Da][DaPars[ind1]] > maxa:
                                                        maxa = self.originalScores[Da][DaPars[ind1]]
                                        if self.originalScores[D][Dpars[D][ind]] <= maxa+maxb+tol:
                                                Dparscopy.remove(Dpars[D][ind])
                                                counter2 += 1
                                                continue
                                sum1 = sum1+len(Dpars[D])
                                Dpars[D] = Dparscopy.copy()
                                sum2 = sum2+len(Dpars[D])
                        elif len(D[0]) == 3:
                                Dpars[D] = list(self.originalScores[D].keys())
                                Dparscopy = Dpars[D].copy()
                                for ind in range(len(Dpars[D])):
                                        # a<->b<->c
                                        if len(D[1]) == 2:
                                                Nodeb = 0
                                                for i in D[0]:
                                                        if i in D[1][0] and i in D[1][1]:
                                                                Nodeb = i
                                                                break
                                                Nodea = 0
                                                for i in D[0]:
                                                        if i != Nodeb and i in D[1][0]:
                                                                Nodea = i
                                                                break
                                                Nodec = 0
                                                for i in D[0]:
                                                        if i != Nodeb and i in D[1][1]:
                                                                Nodec = i
                                                                break
                                                if prune_more == True and (Nodea in Dpars[D][ind][D[0].index(Nodec)] or Nodec in Dpars[D][ind][D[0].index(Nodea)]):
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                Da = ((Nodea,),())
                                                Db = ((Nodeb,),())
                                                Dc = ((Nodec,),())
                                                Dab = (tuple(sorted((Nodea,Nodeb))),(tuple(sorted((Nodea,Nodeb))),))
                                                Dbc = (tuple(sorted((Nodeb,Nodec))),(tuple(sorted((Nodeb,Nodec))),))
                                                # a<->b<->c vs a,b,c
                                                DaPars = list(self.originalScores[Da].keys())
                                                DbPars = list(self.originalScores[Db].keys())
                                                DcPars = list(self.originalScores[Dc].keys())
                                                maxa = -float('inf')
                                                maxb = -float('inf')
                                                maxc = -float('inf')
                                                for ind1 in range(len(DaPars)):
                                                        if set(DaPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and self.originalScores[Da][DaPars[ind1]] > maxa:
                                                                maxa = self.originalScores[Da][DaPars[ind1]]
                                                for ind2 in range(len(DbPars)):
                                                        if set(DbPars[ind2][0]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and self.originalScores[Db][DbPars[ind2]] > maxb:
                                                                maxb = self.originalScores[Db][DbPars[ind2]]
                                                for ind3 in range(len(DcPars)):
                                                        if set(DcPars[ind3][0]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dc][DcPars[ind3]] > maxc:
                                                                maxc = self.originalScores[Dc][DcPars[ind3]]
                                                if self.originalScores[D][Dpars[D][ind]] <= maxa+maxb+maxc:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                # a<->b<->c vs a<->b,c
                                                DabPars = list(self.originalScores[Dab].keys())
                                                DbcPars = list(self.originalScores[Dbc].keys())
                                                maxab = -float('inf')
                                                maxbc = -float('inf')
                                                abInda = 0
                                                abIndb = 1
                                                if Nodea > Nodeb:
                                                        abInda = 1
                                                        abIndb = 0
                                                for ind1 in range(len(DabPars)):
                                                        if set(DabPars[ind1][abInda]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DabPars[ind1][abIndb]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and self.originalScores[Dab][DabPars[ind1]] > maxab:
                                                                maxab = self.originalScores[Dab][DabPars[ind1]]
                                                if self.originalScores[D][Dpars[D][ind]] <= maxab+maxc:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                # a<->b<->c vs a,b<->c
                                                bcIndb = 0
                                                bcIndc = 1
                                                if Nodeb > Nodec:
                                                        bcIndb = 1
                                                        bcIndc = 0
                                                for ind1 in range(len(DbcPars)):
                                                        if set(DbcPars[ind1][bcIndb]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and set(DbcPars[ind1][bcIndc]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dbc][DbcPars[ind1]] > maxbc:
                                                                maxbc = self.originalScores[Dbc][DbcPars[ind1]]
                                                if self.originalScores[D][Dpars[D][ind]] <= maxa+maxbc:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                # a<->b<->c vs a<->b<->c
                                                delInd = False
                                                for i in range(ind):
                                                        if set(Dpars[D][i][0]).difference(set(Dpars[D][ind][0])) == set() and set(Dpars[D][i][1]).difference(set(Dpars[D][ind][1])) == set() and set(Dpars[D][i][2]).difference(set(Dpars[D][ind][2])) == set() and self.originalScores[D][Dpars[D][ind]] <= self.originalScores[D][Dpars[D][i]]:
                                                                Dparscopy.remove(Dpars[D][ind])
                                                                delInd = True
                                                                break
                                                if delInd == True:
                                                        continue
                                        # a<->b<->c<->a
                                        if len(D[1]) == 3:
                                                Nodea = D[0][0]
                                                Nodeb = D[0][1]
                                                Nodec = D[0][2]
                                                Dab = ((Nodea,Nodeb),((Nodea,Nodeb),))
                                                Dbc = ((Nodeb,Nodec),((Nodeb,Nodec),))
                                                Dac = ((Nodea,Nodec),((Nodea,Nodec),))
                                                Dabc = ((Nodea,Nodeb,Nodec),((Nodea,Nodeb),(Nodeb,Nodec)))
                                                if Dabc not in self.originalScores.keys():
                                                        Dabc = ((Nodea,Nodeb,Nodec),((Nodeb,Nodec),(Nodea,Nodeb)))
                                                Dbca = ((Nodea,Nodeb,Nodec),((Nodeb,Nodec),(Nodea,Nodec)))
                                                if Dbca not in self.originalScores.keys():
                                                        Dbca = ((Nodea,Nodeb,Nodec),((Nodea,Nodec),(Nodeb,Nodec)))
                                                Dcab = ((Nodea,Nodeb,Nodec),((Nodea,Nodec),(Nodeb,Nodec)))
                                                if Dcab not in self.originalScores.keys():
                                                        Dcab = ((Nodea,Nodeb,Nodec),((Nodeb,Nodec),(Nodea,Nodec)))
                                                DaPars = list(self.originalScores[Da].keys())
                                                DbPars = list(self.originalScores[Db].keys())
                                                DcPars = list(self.originalScores[Dc].keys())
                                                maxa = -float('inf')
                                                maxb = -float('inf')
                                                maxc = -float('inf')
                                                for ind1 in range(len(DaPars)):
                                                        if set(DaPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and self.originalScores[Da][DaPars[ind1]] > maxa:
                                                                maxa = self.originalScores[Da][DaPars[ind1]]
                                                for ind2 in range(len(DbPars)):
                                                        if set(DbPars[ind2][0]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and self.originalScores[Db][DbPars[ind2]] > maxb:
                                                                maxb = self.originalScores[Db][DbPars[ind2]]
                                                for ind3 in range(len(DcPars)):
                                                        if set(DcPars[ind3][0]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dc][DcPars[ind3]] > maxc:
                                                                maxc = self.originalScores[Dc][DcPars[ind3]]
                                                DabPars = list(self.originalScores[Dab].keys())
                                                DbcPars = list(self.originalScores[Dbc].keys())
                                                DacPars = list(self.originalScores[Dac].keys())
                                                maxab = -float('inf')
                                                maxbc = -float('inf')
                                                maxac = -float('inf')
                                                for ind1 in range(len(DabPars)):
                                                        if set(DabPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DabPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and self.originalScores[Dab][DabPars[ind1]] > maxab:
                                                                maxab = self.originalScores[Dab][DabPars[ind1]]
                                                for ind1 in range(len(DbcPars)):
                                                        if set(DbcPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and set(DbcPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dbc][DbcPars[ind1]] > maxbc:
                                                                maxbc = self.originalScores[Dbc][DbcPars[ind1]]
                                                for ind1 in range(len(DacPars)):
                                                        if set(DacPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DacPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dac][DacPars[ind1]] > maxac:
                                                                maxac = self.originalScores[Dac][DacPars[ind1]]
                                                DabcPars = list(self.originalScores[Dabc].keys())
                                                DbcaPars = list(self.originalScores[Dbca].keys())
                                                DcabPars = list(self.originalScores[Dcab].keys())
                                                maxabc = -float('inf')
                                                maxbca = -float('inf')
                                                maxcab = -float('inf')
                                                for ind1 in range(len(DabcPars)):
                                                        if set(DabcPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DabcPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and set(DabcPars[ind1][2]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dabc][DabcPars[ind1]] > maxabc:
                                                                maxabc = self.originalScores[Dabc][DabcPars[ind1]]
                                                for ind1 in range(len(DbcaPars)):
                                                        if set(DbcaPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DbcaPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and set(DbcaPars[ind1][2]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dbca][DbcaPars[ind1]] > maxbca:
                                                                maxbca = self.originalScores[Dbca][DbcaPars[ind1]]
                                                for ind1 in range(len(DcabPars)):
                                                        if set(DcabPars[ind1][0]).difference(set(Dpars[D][ind][D[0].index(Nodea)])) == set() and set(DcabPars[ind1][1]).difference(set(Dpars[D][ind][D[0].index(Nodeb)])) == set() and set(DcabPars[ind1][2]).difference(set(Dpars[D][ind][D[0].index(Nodec)])) == set() and self.originalScores[Dcab][DcabPars[ind1]] > maxcab:
                                                                maxcab = self.originalScores[Dcab][DcabPars[ind1]]
                                                # vs a,b,c
                                                if self.originalScores[D][Dpars[D][ind]] <= maxa+maxb+maxc:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                # vs a<->b,c
                                                if self.originalScores[D][Dpars[D][ind]] <= maxab+maxc:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                # vs a, b<->c
                                                if self.originalScores[D][Dpars[D][ind]] <= maxa+maxbc:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                # vs b, c<->a
                                                if self.originalScores[D][Dpars[D][ind]] <= maxac+maxb:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                # vs a<->b<->c
                                                if self.originalScores[D][Dpars[D][ind]] <= maxabc:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                # vs b<->c<->a
                                                if self.originalScores[D][Dpars[D][ind]] <= maxbca:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                # vs c<->a<->b
                                                if self.originalScores[D][Dpars[D][ind]] <= maxcab:
                                                        Dparscopy.remove(Dpars[D][ind])
                                                        continue
                                                # vs a<->b<->c<->a
                                                delInd = False
                                                for i in range(ind):
                                                        if set(Dpars[D][i][0]).difference(set(Dpars[D][ind][0])) == set() and set(Dpars[D][i][1]).difference(set(Dpars[D][ind][1])) == set() and set(Dpars[D][i][2]).difference(set(Dpars[D][ind][2])) == set() and self.originalScores[D][Dpars[D][ind]] <= self.originalScores[D][Dpars[D][i]]:
                                                                Dparscopy.remove(Dpars[D][ind])
                                                                delInd = True
                                                                break
                                                if delInd == True:
                                                        continue
                                sum1 = sum1+len(Dpars[D])
                                Dpars[D] = Dparscopy.copy()
                                sum2 = sum2+len(Dpars[D])
                        else:
                                Dpars[D] = list(self.originalScores[D].keys())
                                sum1 = sum1+len(Dpars[D])
                                sum2 = sum2+len(Dpars[D])
                print(str(sum1)+" vs "+str(sum2)+", pruning time: "+str(time.time()-t0P))
                fileName = self.resultsdir+self.instance+'_cplex.log'
                f = open(fileName,"a")
                f.write(str(sum1)+" vs "+str(sum2)+", pruning time: "+str(time.time()-t0P)+"\n")
                f.close

                nOneNodeDist = 0
                nTwoNodeDist = 0
                self.temp1Scores = {}
                for D in self.originalScores.keys():
                        for Dpar in Dpars[D]:
                                if D not in self.temp1Scores.keys():
                                        self.temp1Scores[D] = {}
                                self.temp1Scores[D][Dpar] = self.originalScores[D][Dpar]
                                # count number of c-comps with one or two node districts
                                if len(D[0]) == 1:
                                        nOneNodeDist += 1
                                elif len(D[0]) == 2:
                                        nTwoNodeDist += 1
                f = open(fileName,"a")
                f.write("number of c-comps with one node districts: "+str(nOneNodeDist)+"\n")
                f.write("number of c-comps with two node districts: "+str(nTwoNodeDist)+"\n")
                f.write("number of c_comps deleted type 1: "+str(counter1)+"\n")
                f.write("number of c_comps deleted type 2: "+str(counter2)+"\n")
                f.close

        def Prune_nonarid_ccomps(self):
                t0P = time.time()
                sum1 = 0
                sum2 = 0
                Dpars = {}
                for D in self.temp2Scores.keys():
                        Dpars[D] = list(self.temp2Scores[D].keys())
                        if len(D[0]) < 4:
                                sum1 = sum1+len(Dpars[D])
                                sum2 = sum2+len(Dpars[D])
                        else:
                                Dparscopy = Dpars[D].copy()
                                for i in range(len(Dpars[D])):
                                        is_arid = test_for_arid(D[0],Dpars[D][i],D[1])
                                        if is_arid == False:
                                                Dparscopy.remove(Dpars[D][i])
                                sum1 = sum1+len(Dpars[D])
                                Dpars[D] = Dparscopy.copy()
                                sum2 = sum2+len(Dpars[D])

                print("after prune non-arid, "+str(sum1)+" vs "+str(sum2)+", pruning time: "+str(time.time()-t0P))
                fileName = self.resultsdir+self.instance+'_cplex.log'
                f = open(fileName,"a")
                f.write("after prune non-arid, "+str(sum1)+" vs "+str(sum2)+", pruning time: "+str(time.time()-t0P)+"\n")
                f.close

                nOneNodeDist = 0
                nTwoNodeDist = 0
                self.scores = {}
                for D in self.temp2Scores.keys():
                        for Dpar in Dpars[D]:
                                if D not in self.scores.keys():
                                        self.scores[D] = {}
                                self.scores[D][Dpar] = self.temp2Scores[D][Dpar]
                                # count number of c-comps with one or two node districts
                                if len(D[0]) == 1:
                                        nOneNodeDist += 1
                                elif len(D[0]) == 2:
                                        nTwoNodeDist += 1
                f = open(fileName,"a")
                f.write("number of c-comps with one node districts: "+str(nOneNodeDist)+"\n")
                f.write("number of c-comps with two node districts: "+str(nTwoNodeDist)+"\n")
                f.close

        def Prune_bow_ccomps(self):
                t0P = time.time()
                sum1 = 0
                sum2 = 0
                Dpars = {}
                for D in self.temp1Scores.keys():
                        Dpars[D] = list(self.temp1Scores[D].keys())
                        if len(D[0]) == 1:
                                sum1 = sum1+len(Dpars[D])
                                sum2 = sum2+len(Dpars[D])
                                continue
                        Dparscopy = Dpars[D].copy()
                        for i in range(len(Dpars[D])):
                                remove = False
                                for j in range(len(Dpars[D][i])):
                                        head = D[0][j]
                                        for k in range(len(Dpars[D][i][j])):
                                                tail = Dpars[D][i][j][k]
                                                for l in range(len(D[1])):
                                                        node1 = D[1][l][0]
                                                        node2 = D[1][l][1]
                                                        if (node1==tail and node2==head) or (node1==head and node2==tail):
                                                               remove = True
                                                               break
                                                if remove == True:
                                                        break
                                        if remove == True:
                                                break
                                if remove == True:
                                        Dparscopy.remove(Dpars[D][i])
                        sum1 = sum1+len(Dpars[D])
                        Dpars[D] = Dparscopy.copy()
                        sum2 = sum2+len(Dpars[D])

                print("after prune bow, "+str(sum1)+" vs "+str(sum2)+", pruning time: "+str(time.time()-t0P))
                fileName = self.resultsdir+self.instance+'_cplex.log'
                f = open(fileName,"a")
                f.write("after prune bow, "+str(sum1)+" vs "+str(sum2)+", pruning time: "+str(time.time()-t0P)+"\n")
                f.close

                nOneNodeDist = 0
                nTwoNodeDist = 0
                self.temp2Scores = {}
                for D in self.temp1Scores.keys():
                        for Dpar in Dpars[D]:
                                if D not in self.temp2Scores.keys():
                                        self.temp2Scores[D] = {}
                                self.temp2Scores[D][Dpar] = self.temp1Scores[D][Dpar]
                                # count number of c-comps with one or two node districts
                                if len(D[0]) == 1:
                                        nOneNodeDist += 1
                                elif len(D[0]) == 2:
                                        nTwoNodeDist += 1
                f = open(fileName,"a")
                f.write("number of c-comps with one node districts: "+str(nOneNodeDist)+"\n")
                f.write("number of c-comps with two node districts: "+str(nTwoNodeDist)+"\n")
                f.close


                
        def Initialize(self,prune=True,dag=False,printsc=False,prune_parentInDistrict=False):
                                
                if prune == True:
                        self.Prune_scores(prune_more=prune_parentInDistrict)
                else:
                        self.temp1Scores = self.originalScores
                if self.bowfree == True:
                        self.Prune_bow_ccomps()
                else:
                        self.temp2Scores = self.temp1Scores
                if self.arid == True:
                        self.Prune_nonarid_ccomps()
                else:
                        self.scores = self.temp2Scores
                self.dag = dag
                self.CreateLP(printsc)

        def CreateLP(self,printsc=False):
                self.V = set()
                self.cComps = []
                for D in self.scores.keys():
                        self.cComps.append(D)
                        self.V = self.V.union(set(D[0]))
                self.iComps = {}
                for i in self.V:
                        self.iComps[i] = []
                        for D in self.cComps:
                                if i in D[0]:
                                        self.iComps[i].append(self.cComps.index(D))
                self.dPars = {}
                for d in range(len(self.cComps)):
                        self.dPars[d] = []
                        for par in self.scores[self.cComps[d]].keys():
                                self.dPars[d].append(par)
                self.iPars = {}
                for i in self.V:
                        self.iPars[i] = []
                        for d in self.iComps[i]:
                                for W in self.dPars[d]:
                                        if W[self.cComps[d][0].index(i)] not in self.iPars[i]:
                                                self.iPars[i].append(W[self.cComps[d][0].index(i)])
                if printsc == True:
                        for d in range(len(self.cComps)):
                                print(str(self.cComps[d])+':')
                                print(self.scores[self.cComps[d]])
                                print('\n')
                
                
                self.biPars = {}
                for D in self.cComps:
                        for bi in D[1]:
                                if bi not in self.biPars.keys():
                                        self.biPars[bi] = []
                                for W in self.dPars[self.cComps.index(D)]:
                                        biPar = (W[D[0].index(bi[0])],W[D[0].index(bi[1])])
                                        if biPar not in self.biPars[bi]:
                                                self.biPars[bi].append(biPar)
                self.lp = cplex.Cplex()
                self.z = {}
                for d in range(len(self.cComps)):
                        for dp in range(len(self.dPars[d])):
                                if self.dag == False or len(self.cComps[d][0]) <= 1:
                                        self.z[d,dp] = self.lp.variables.add(obj=[self.scores[self.cComps[d]][self.dPars[d][dp]]],lb=[0],ub=[1],names=['z'+str(d)+','+str(dp)])
                                else:
                                        self.z[d,dp] = self.lp.variables.add(obj=[self.scores[self.cComps[d]][self.dPars[d][dp]]],lb=[0],ub=[0],names=['z'+str(d)+','+str(dp)])

                
                self.lp.objective.set_sense(self.lp.objective.sense.maximize)
                
                
                for i in self.V:
                        self.lp.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = ['z'+str(d)+','+str(dp) for d in self.iComps[i] for dp in range(len(self.dPars[d]))], val = [1]*sum(len(self.dPars[d]) for d in self.iComps[i]))], senses=["E"], rhs=[1])

                self.indInv = []
                for i in self.V:
                        for j in range(i+1,len(self.V)):
                                self.indInv.append((i,j))

                self.udE = range(len(self.indInv))

                self.bi = {}
                for e in self.udE:
                        self.bi[e] = self.lp.variables.add(obj=[0],names=['bi'+str(e)])
                        zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if self.indInv[e] in self.cComps[d][1]]
                        self.lp.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set+['bi'+str(e)], val = [1]*len(zindex_set)+[-1])], senses=["E"], rhs=[0])
                        
                self.x = {}
                for i in self.V:
                        for ip in range(len(self.iPars[i])):
                                self.x[i,ip] = self.lp.variables.add(obj=[0],names=['x'+str(i)+','+str(ip)])
                                zindex_set = ['z'+str(d)+','+str(dp) for d in self.iComps[i] for dp in range(len(self.dPars[d])) if self.dPars[d][dp][self.cComps[d][0].index(i)] == self.iPars[i][ip]]
                                self.lp.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set+['x'+str(i)+','+str(ip)], val = [1]*len(zindex_set)+[-1])], senses=["E"], rhs=[0])

                # -------------------------------
                # new 2-node set cuts
                if 'a' in self.cuts:
                        self.cuts_2_node_set()
                # -------------------------------

        def biClusterToIneq(self,C,ii,jj):
                if jj < ii:
                        cp = jj
                        jj = ii
                        ii = cp
                ifLHS = {(d,dp):False for d in range(len(self.cComps)) for dp in range(len(self.dPars[d]))}
                for d in range(len(self.cComps)):
                        if len(set(self.cComps[d][0])&set([ii,jj]))!=1:
                                vs = [v for v in C if v in self.cComps[d][0]]
                                if len(vs) > 0:
                                        for dp in range(len(self.dPars[d])):
                                                for v in vs:
                                                        if set(self.dPars[d][dp][self.cComps[d][0].index(v)])&C.union(set([ii,jj])) == set():
                                                                ifLHS[d,dp] = True
                                                                break
                                if (ii,jj) in self.cComps[d][1]:
                                        for dp in range(len(self.dPars[d])):
                                                if set(self.dPars[d][dp][self.cComps[d][0].index(ii)])&C.union(set([ii,jj])) == set() and set(self.dPars[d][dp][self.cComps[d][0].index(jj)])&C.union(set([ii,jj])) == set():
                                                        ifLHS[d,dp] = True
                return ifLHS

        def ClusterToIneq(self,C):
                ifLHS = {(d,dp):False for d in range(len(self.cComps)) for dp in range(len(self.dPars[d]))}
                for d in range(len(self.cComps)):
                        vs = [v for v in C if v in self.cComps[d][0]]
                        if len(vs) > 0:
                                for dp in range(len(self.dPars[d])):
                                                for v in vs:
                                                        if set(self.dPars[d][dp][self.cComps[d][0].index(v)])&set(C) == set():
                                                                ifLHS[d,dp] = True
                                                                break
                return ifLHS


        def isThereBow(self, dist, pars):
                isThereBow = False
                if len(dist[0]) == 1:
                        return isThereBow
                mapNodeToPars = {}
                for i in range(len(dist[0])):
                        node = dist[0][i]
                        mapNodeToPars[node] = set(pars[i])
                for i in range(len(dist[1])):
                        node1 = dist[1][i][0]
                        node2 = dist[1][i][1]
                        if node1 in mapNodeToPars[node2] or node2 in mapNodeToPars[node1]:
                                isThereBow = True
                                break
                return isThereBow

        
        def CombineCComponents(self, z_value):
                use_our_scoring_function = True
                candDist = [] # candidate districts
                candPars = [] # candidate parents
                candScore = [] # candidate scores
                dParsFracInd = {}
                dParsOneInd = {}
                for d in range(len(self.cComps)):
                        dParsFracInd[d] = []
                        dParsOneInd[d] = []
                        for dp in range(len(self.dPars[d])):
                                if z_value[d,dp] > 0.0:
                                        dParsOneInd[d].append(dp)
                
                # increase by 1 the number of parents of c-components that has positive LP solution
                if 'a' in self.heuristics:
                        print('starting heuristic a')
                        for d in range(len(self.cComps)):                        
                                if len(dParsOneInd[d]) > 0:
                                        parents = [self.dPars[d][dParsOneInd[d][dp]] for dp in range(len(dParsOneInd[d]))]
                                        scores = []
                                        for dp in range(len(dParsOneInd[d])):
                                                scores.append(self.scores[self.cComps[d]][self.dPars[d][dParsOneInd[d][dp]]])
                                                print('positive '+'z'+str(d)+', district = '+str(self.cComps[d][0])+', parents = '+str(parents)+', scores = '+str(scores))

                                        for dp in range(len(dParsOneInd[d])):
                                                for i in range(len(self.cComps[d][0])):
                                                        dnode = self.cComps[d][0][i]
                                                        pdnode = set(self.dPars[d][dParsOneInd[d][dp]][i])
                                                        for node in self.V:
                                                                if self.cComps[d][0][i] == node or dnode in pdnode:
                                                                        continue
                                                                node = [node]
                                                                node = tuple(node)
                                                                pdnodeUnion = set(self.dPars[d][dParsOneInd[d][dp]][i]).union(node)
                                                                pdnodeUnionSorted = []
                                                                for j in self.V:
                                                                        if j in pdnodeUnion:
                                                                                pdnodeUnionSorted.append(j)

                                                                pUnionSorted = list(self.dPars[d][dParsOneInd[d][dp]])
                                                                pUnionSorted[i] = tuple(pdnodeUnionSorted)
                                                                pUnionSorted = tuple(pUnionSorted)
                                                                print('union = '+str(pUnionSorted))
                                                                if use_our_scoring_function == True:
                                                                        if len(self.cComps[d][0]) == 1:
                                                                                score = compute_BiC_c_comps_sets(self.data,self.cComps[d][0],pUnionSorted,-1)
                                                                        else:
                                                                                score = compute_BiC_c_comps_sets(self.data,self.cComps[d][0],pUnionSorted,self.cComps[d][1])

                                                                if score > 0.0:
                                                                        print('replaced score '+str(score)+' with score '+str(-1.0*pow(10,20)))
                                                                        score = -1.0*pow(10,20)
                                                                        print('score = '+str(score))
                                                                if score > self.scores[self.cComps[d]][self.dPars[d][dParsOneInd[d][dp]]]:
                                                                        print('added to LP')
                                                                        candDist.append(self.cComps[d])
                                                                        candPars.append(pUnionSorted)
                                                                        candScore.append(score)
                                                                        
                # add bidirected edge to any c-component that has positive LP solution
                # It adds a node to the district and creates a bidirected edge with an existing node in the district
                # the node added to the district cannot be a node in the parent set of the existing node
                if 'b' in self.heuristics:
                        print('starting heuristic b')
                        for d in range(len(self.cComps)):
                                if len(dParsOneInd[d]) == 0:
                                        continue # only proceed if the variable is positive
                                for dp in range(len(dParsOneInd[d])):
                                        for index in range(len(self.cComps[d][0])):
                                                n1 = self.cComps[d][0][index]
                                                for n2 in self.V:
                                                        if n1 == n2 or n2 in set(self.dPars[d][dParsOneInd[d][dp]][index]):
                                                                continue
                                                        # create a c-comp with node n2 in district and no parents
                                                        district1 = (tuple([n2]),())
                                                        parent1 = (tuple([]),)
                                                        if use_our_scoring_function == True:
                                                                score1 = compute_BiC_c_comps_sets(self.data,district1[0],parent1,district1[1])
                                                        # from this point we are going to add n2 to the district
                                                        # union of the nodes
                                                        nUnion = set(self.cComps[d][0]).union(district1[0])
                                                        nUnionSorted = []
                                                        for i in self.V:
                                                                if i in nUnion:
                                                                        nUnionSorted.append(i)
                                                        nUnionSorted = tuple(nUnionSorted)
                                                        # union of edges
                                                        edges = self.cComps[d][1]
                                                        edges1 = None
                                                        if n1 < n2:
                                                                edges1 = (n1,n2)
                                                        else:
                                                                edges1 = (n2,n1)
                                                        edges1 = (edges1,)
                                                        eUnion = []
                                                        size = len(edges)
                                                        size1 = len(edges1)
                                                        i, j = 0, 0
                                                        while i < size and j < size1:
                                                                if (edges[i][0] < edges1[j][0] or
                                                                    (edges[i][0] == edges1[j][0] and
                                                                     edges[i][1] < edges1[j][1])):
                                                                        eUnion.append(edges[i])
                                                                        i += 1
                                                                elif (edges[i][0] == edges1[j][0] and
                                                                      edges[i][1] == edges1[j][1]):
                                                                        eUnion.append(edges[i])
                                                                        i += 1
                                                                        j += 1
                                                                else:
                                                                        eUnion.append(edges1[j])
                                                                        j += 1
                                                        for iter in range(i,size):
                                                                eUnion.append(edges[iter])
                                                        for iter in range(j,size1):
                                                                eUnion.append(edges1[iter])
                                                        eUnion = tuple(eUnion)
                                                        # union of parents
                                                        parent = self.dPars[d][dParsOneInd[d][dp]]
                                                        score = self.scores[self.cComps[d]][self.dPars[d][dParsOneInd[d][dp]]]
                                                        print('positive '+'z'+str(d)+', district = '+str(self.cComps[d][0])+', parent = '+str(parent)+', score = '+str(score)+', to add: district = '+str(district1[0])+', parent = '+str(parent1)+', score = '+str(score1))
                                                        pUnion = [set() for i in range(len(nUnionSorted))]
                                                        pUnionSortedList = [[] for i in range(len(nUnionSorted))]
                                                        for i in range(len(nUnionSorted)):
                                                                node = nUnionSorted[i]
                                                                for j in range(len(self.cComps[d][0])):
                                                                        if node == self.cComps[d][0][j]:
                                                                                pUnion[i] = pUnion[i].union(self.dPars[d][dParsOneInd[d][dp]][j])
                                                                                break
                                                                for j in range(len(district1[0])):
                                                                        if node == district1[0][j]:
                                                                                pUnion[i] = pUnion[i].union(parent1[j])
                                                                                break
                                                        for i in self.V:
                                                                for j in range(len(pUnion)):
                                                                        if i in pUnion[j]:
                                                                                pUnionSortedList[j].append(i)
                                                        pUnionSorted = []
                                                        for j in range(len(pUnionSortedList)):
                                                                pUnionSorted.append(tuple(pUnionSortedList[j]))
                                                        pUnionSorted = tuple(pUnionSorted)
                                                        if use_our_scoring_function == True:
                                                                newscore = compute_BiC_c_comps_sets(self.data,nUnionSorted,pUnionSorted,eUnion)
                                                        if newscore > 0.0:
                                                                print('replaced score '+str(newscore)+' with score '+str(-1.0*pow(10,20)))
                                                                newscore = -1.0*pow(10,20)
                                                        print('new district = ('+str(nUnionSorted)+', '+str(eUnion)+'), new parents = '+str(pUnionSorted)+', score = '+str(newscore))
                                                        if newscore > score + score1:
                                                                print('added to LP')
                                                                candDist.append((nUnionSorted, eUnion))
                                                                candPars.append(pUnionSorted)
                                                                candScore.append(newscore)

                # For any c-component that has positive LP solution
                # removes a parent and makes it a spouse (basically transforming a directed edge
                # into a bidirected edge)
                if 'c' in self.heuristics:
                        print('starting heuristic c')
                        for d in range(len(self.cComps)):
                                if len(dParsOneInd[d]) == 0:
                                        continue # only proceed if the variable is positive
                                for dp in range(len(dParsOneInd[d])):
                                        for index1 in range(len(self.cComps[d][0])):
                                                n1 = self.cComps[d][0][index1]
                                                for index2 in range(len(self.dPars[d][dParsOneInd[d][dp]][index1])):
                                                        n2 = self.dPars[d][dParsOneInd[d][dp]][index1][index2]
                                                        if n2 in set(self.cComps[d][0]):
                                                                continue # n2 cannot be already part of the district.
                                                        # consider the score of the c-component formed only by node n2
                                                        district = (tuple([n2]),())
                                                        parent = (tuple([]),)
                                                        score = self.scores[district][parent]
                                                        # union of the nodes
                                                        setToAdd = set()
                                                        setToAdd.add(n2)
                                                        nUnion = set(self.cComps[d][0]).union(setToAdd)
                                                        nUnionSorted = []
                                                        for i in self.V:
                                                                if i in nUnion:
                                                                        nUnionSorted.append(i)
                                                        nUnionSorted = tuple(nUnionSorted)
                                                        # union of edges
                                                        edges = self.cComps[d][1]
                                                        edges1 = None
                                                        if n1 < n2:
                                                                edges1 = (n1,n2)
                                                        else:
                                                                edges1 = (n2,n1)
                                                        edges1 = (edges1,)
                                                        eUnion = []
                                                        size = len(edges)
                                                        size1 = len(edges1)
                                                        i, j = 0, 0
                                                        while i < size and j < size1:
                                                                if (edges[i][0] < edges1[j][0] or
                                                                    (edges[i][0] == edges1[j][0] and
                                                                     edges[i][1] < edges1[j][1])):
                                                                        eUnion.append(edges[i])
                                                                        i += 1
                                                                elif (edges[i][0] == edges1[j][0] and
                                                                      edges[i][1] == edges1[j][1]):
                                                                        eUnion.append(edges[i])
                                                                        i += 1
                                                                        j += 1
                                                                else:
                                                                        eUnion.append(edges1[j])
                                                                        j += 1
                                                        for iter in range(i,size):
                                                                eUnion.append(edges[iter])
                                                        for iter in range(j,size1):
                                                                eUnion.append(edges1[iter])
                                                        eUnion = tuple(eUnion)
                                                        # union of parents
                                                        pUnion = [set() for i in range(len(nUnionSorted))]
                                                        pUnionSortedList = [[] for i in range(len(nUnionSorted))]
                                                        for i in range(len(nUnionSorted)):
                                                                node = nUnionSorted[i]
                                                                for j in range(len(self.cComps[d][0])):
                                                                        if node == self.cComps[d][0][j]:
                                                                                if node == n1:
                                                                                        temp = set(self.dPars[d][dParsOneInd[d][dp]][j]).difference(setToAdd)
                                                                                        pUnion[i] = pUnion[i].union(temp)
                                                                                elif node == n2:
                                                                                        # do nothing because no new parent to add
                                                                                        break
                                                                                else:
                                                                                        pUnion[i] = pUnion[i].union(self.dPars[d][dParsOneInd[d][dp]][j])
                                                                                break
                                                        for i in self.V:
                                                                for j in range(len(pUnion)):
                                                                        if i in pUnion[j]:
                                                                                pUnionSortedList[j].append(i)
                                                        pUnionSorted = []
                                                        for j in range(len(pUnionSortedList)):
                                                                pUnionSorted.append(tuple(pUnionSortedList[j]))
                                                        pUnionSorted = tuple(pUnionSorted)
                                                        # score of ccomp without the parent that is going to be moved to district
                                                        if len(self.cComps[d][0]) == 1:
                                                                score1 = compute_BiC_c_comps_sets(self.data,self.cComps[d][0],pUnionSorted,-1)
                                                        else:
                                                                score1 = compute_BiC_c_comps_sets(self.data,self.cComps[d][0],pUnionSorted,self.cComps[d][1])
                                                        print('positive '+'z'+str(d)+', district = '+str(self.cComps[d][0])+', parent = '+str(self.dPars[d][dParsOneInd[d][dp]])+', score of ccomp made by parent = '+str(score)+', score of ccomp without parent = '+str(score1)+'. Moving parent '+str(n2)+' to district')
                                                        newscore = compute_BiC_c_comps_sets(self.data,nUnionSorted,pUnionSorted,eUnion)
                                                        if newscore > 0.0:
                                                                print('replaced score '+str(newscore)+' with score '+str(-1.0*pow(10,20)))
                                                                newscore = -1.0*pow(10,20)
                                                        print('new district = ('+str(nUnionSorted)+', '+str(eUnion)+'), new parents = '+str(pUnionSorted)+', score = '+str(newscore))
                                                        if newscore > score + score1:
                                                                print('added to LP')
                                                                candDist.append((nUnionSorted, eUnion))
                                                                candPars.append(pUnionSorted)
                                                                candScore.append(newscore)

                # Same district
                for d in range(len(self.cComps)):
                        if len(dParsFracInd[d]) > 1:
                                parents = [self.dPars[d][dParsFracInd[d][dp]] for dp in range(len(dParsFracInd[d]))]
                                scores = []
                                for dp in range(len(dParsFracInd[d])):
                                        scores.append(self.scores[self.cComps[d]][self.dPars[d][dParsFracInd[d][dp]]])
                                print('fractional '+'z'+str(d)+', district = '+str(self.cComps[d][0])+', parents = '+str(parents)+', scores = '+str(scores))
                                for dp in range(len(dParsFracInd[d])-1):
                                        for dp1 in range(dp+1, len(dParsFracInd[d])):
                                                if len(self.cComps[d][1]) == 0:
                                                        pUnion = set(self.dPars[d][dParsFracInd[d][dp]][0]).union(self.dPars[d][dParsFracInd[d][dp1]][0])
                                                        pUnionSorted = []
                                                        for i in self.V:
                                                                if i in pUnion:
                                                                        pUnionSorted.append(i)
                                                        pUnionSorted = (tuple(pUnionSorted),)
                                                        print('union = '+str(pUnionSorted))
                                                        score = compute_BiC_c_comps_sets(self.data,self.cComps[d][0],pUnionSorted,-1)
                                                        print('score = '+str(score))
                                                        candDist.append(self.cComps[d])
                                                        candPars.append(pUnionSorted)
                                                        candScore.append(score)
                                                else:
                                                        pUnion = []
                                                        pUnionSortedList = []
                                                        for i in range(len(self.cComps[d][0])):
                                                                pUnion.append(set(self.dPars[d][dParsFracInd[d][dp]][i]).union(self.dPars[d][dParsFracInd[d][dp1]][i]))
                                                                pUnionSortedList.append([])
                                                        for i in self.V:
                                                                for j in range(len(pUnion)):
                                                                        if i in pUnion[j]:
                                                                                pUnionSortedList[j].append(i)
                                                        pUnionSorted = []
                                                        for j in range(len(pUnionSortedList)):
                                                                pUnionSorted.append(tuple(pUnionSortedList[j]))
                                                        pUnionSorted = tuple(pUnionSorted)
                                                        print('union = '+str(pUnionSorted))
                                                        score = compute_BiC_c_comps_sets(self.data,self.cComps[d][0],pUnionSorted,self.cComps[d][1])
                                                        print('score = '+str(score))
                                                        candDist.append(self.cComps[d])
                                                        candPars.append(pUnionSorted)
                                                        candScore.append(score)

                # two different districts
                for d in range(len(self.cComps)-1):
                        if len(dParsFracInd[d]) > 0:
                                for d1 in range(d+1, len(self.cComps)):
                                        if len(dParsFracInd[d1]) > 0 and (not (len(self.cComps[d][0]) == 1 and
                                                                               len(self.cComps[d1][0]) == 1)):
                                                canMerge = False
                                                for i in range(len(self.cComps[d][0])):
                                                        for j in range(len(self.cComps[d1][0])):
                                                                if self.cComps[d][0][i] == self.cComps[d1][0][j]:
                                                                        canMerge = True
                                                                        break
                                                        if canMerge == True:
                                                                break
                                                if canMerge:
                                                        # union of the nodes
                                                        nUnion = set(self.cComps[d][0]).union(self.cComps[d1][0])
                                                        nUnionSorted = []
                                                        for i in self.V:
                                                                if i in nUnion:
                                                                        nUnionSorted.append(i)
                                                        nUnionSorted = tuple(nUnionSorted)
                                                        # union of edges
                                                        edges = self.cComps[d][1]
                                                        edges1 = self.cComps[d1][1]
                                                        eUnion = []
                                                        size = len(edges)
                                                        size1 = len(edges1)
                                                        i, j = 0, 0
                                                        while i < size and j < size1:
                                                                if (edges[i][0] < edges1[j][0] or
                                                                    (edges[i][0] == edges1[j][0] and
                                                                     edges[i][1] < edges1[j][1])):
                                                                        eUnion.append(edges[i])
                                                                        i += 1
                                                                elif (edges[i][0] == edges1[j][0] and
                                                                      edges[i][1] == edges1[j][1]):
                                                                        eUnion.append(edges[i])
                                                                        i += 1
                                                                        j += 1
                                                                else:
                                                                        eUnion.append(edges1[j])
                                                                        j += 1
                                                        for iter in range(i,size):
                                                                eUnion.append(edges[iter])
                                                        for iter in range(j,size1):
                                                                eUnion.append(edges1[iter])
                                                        eUnion = tuple(eUnion)
                                                        # union of parents
                                                        for dp in range(len(dParsFracInd[d])):
                                                                for dp1 in range(len(dParsFracInd[d1])):
                                                                        parent = self.dPars[d][dParsFracInd[d][dp]]
                                                                        score = self.scores[self.cComps[d]][self.dPars[d][dParsFracInd[d][dp]]]
                                                                        parent1 = self.dPars[d1][dParsFracInd[d1][dp1]]
                                                                        score1 = self.scores[self.cComps[d1]][self.dPars[d1][dParsFracInd[d1][dp1]]]
                                                                        print('fractional '+'z'+str(d)+', district = '+str(self.cComps[d][0])+', parent = '+str(parent)+', score = '+str(score)+', fractional '+'z'+str(d1)+', district = '+str(self.cComps[d1][0])+', parent = '+str(parent1)+', score = '+str(score1))
                                                                        pUnion = [set() for i in range(len(nUnionSorted))]
                                                                        pUnionSortedList = [[] for i in range(len(nUnionSorted))]
                                                                        for i in range(len(nUnionSorted)):
                                                                                node = nUnionSorted[i]
                                                                                for j in range(len(self.cComps[d][0])):
                                                                                        if node == self.cComps[d][0][j]:
                                                                                                pUnion[i] = pUnion[i].union(self.dPars[d][dParsFracInd[d][dp]][j])
                                                                                                break
                                                                                for j in range(len(self.cComps[d1][0])):
                                                                                        if node == self.cComps[d1][0][j]:
                                                                                                pUnion[i] = pUnion[i].union(self.dPars[d1][dParsFracInd[d1][dp1]][j])
                                                                                                break
                                                                        for i in self.V:
                                                                                for j in range(len(pUnion)):
                                                                                        if i in pUnion[j]:
                                                                                                pUnionSortedList[j].append(i)
                                                                        pUnionSorted = []
                                                                        for j in range(len(pUnionSortedList)):
                                                                                pUnionSorted.append(tuple(pUnionSortedList[j]))
                                                                        pUnionSorted = tuple(pUnionSorted)
                                                                        print('union = '+str(pUnionSorted))
                                                                        score = compute_BiC_c_comps_sets(self.data,nUnionSorted,pUnionSorted,eUnion)
                                                                        if score > 0.0:
                                                                                print('replaced score '+str(score)+' with score '+str(-1.0*pow(10,20)))
                                                                                score = -1.0*pow(10,20)
                                                                        print('new district = ('+str(nUnionSorted)+', '+str(eUnion)+'), new parents = '+str(pUnionSorted)+', score = '+str(score))
                                                                        candDist.append((nUnionSorted, eUnion))
                                                                        candPars.append(pUnionSorted)
                                                                        candScore.append(score)
                # add candidates to scores
                for i in range(len(candDist)):
                        dist = candDist[i]
                        pars = candPars[i]
                        score = candScore[i]
                        if self.isThereBow(dist,pars):
                                continue
                        if dist in self.scores.keys():
                                if pars not in self.scores[dist].keys():
                                        self.scores[dist][pars] = score
                        else:
                                self.scores[dist] = {}
                                self.scores[dist][pars] = score
                                                        
        def calculateLHS(self,z_value,ifLHS):
                lhs = 0.0
                for d in range(len(self.cComps)):
                        for dp in range(len(self.dPars[d])):
                                if ifLHS[d,dp] == True:
                                        lhs += z_value[d,dp]
                return lhs

        def cuts_2_node_set(self, z_value=None):
                print('a: cuts_2_node_set')
                if z_value != None:
                        print('cuts_2_node_set(...) works only if z_value is not provided')
                        return
                t0cuts = time.time()
                nnewcuts = 0
                ne = 0
                for e in self.udE:
                        ne += 1
                        ifLHS = {(d,dp):False for d in range(len(self.cComps)) for dp in range(len(self.dPars[d]))}
                        node1 = self.indInv[e][0]
                        node2 = self.indInv[e][1]
                        for j in range(len(self.iComps[node1])):
                                d = self.iComps[node1][j]
                                if node2 not in self.cComps[d][0]:
                                        for dp in range(len(self.dPars[d])):
                                                if node2 in self.dPars[d][dp][self.cComps[d][0].index(node1)]:
                                                        ifLHS[d,dp] = True
                        for j in range(len(self.iComps[node2])):
                                d = self.iComps[node2][j]
                                if node1 not in self.cComps[d][0]:
                                        for dp in range(len(self.dPars[d])):
                                                if node1 in self.dPars[d][dp][self.cComps[d][0].index(node2)]:
                                                        ifLHS[d,dp] = True
                        for j in range(len(self.iComps[node1])):
                                d = self.iComps[node1][j]
                                if node2 in self.cComps[d][0]:
                                        for dp in range(len(self.dPars[d])):
                                                ifLHS[d,dp] = True
                        zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
                        self.lp.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set, val = [1]*len(zindex_set))], senses=["L"], rhs=[1])
                        nnewcuts += 1

                print('Number of new cuts added: '+str(nnewcuts)+', time: '+str(time.time()-t0cuts))


        def cuts_3_node_set(self, z_value=None):
                print('b: cuts_3_node_set')
                t0cuts = time.time()
                if z_value == None:
                        print('cuts_3_node_set(...) works only if z_value provided')
                        return
                icc = []
                for d in range(len(self.cComps)):
                        for dp in range(len(self.dPars[d])):
                                if z_value[d,dp] > 1.0e-6:
                                        icc.append((d,dp))
                nnewcuts = 0
                list_nodes = [i for i in range(len(self.V))]
                all_3node_subsets = [list(i) for i in itertools.combinations(list_nodes,3)]
                few_3node_subsets = []
                # loop to select 3-node subsets where constraint is violated
                for subset in all_3node_subsets:
                        n1 = subset[0]
                        n2 = subset[1]
                        n3 = subset[2]
                        lhs = 0.0
                        for i in range(len(icc)):
                                d = icc[i][0]
                                dp = icc[i][1]
                                nodes = self.cComps[d][0]
                                pars = self.dPars[d][dp]
                                i1 = -1 # this is the index of n1 in cComps[d][0] (if n1 is not there, then i1=-1)
                                i2 = -1
                                i3 = -1
                                for j in range(len(nodes)):
                                        if n1 == nodes[j]:
                                                i1 = j
                                        elif n2 == nodes[j]:
                                                i2 = j
                                        elif n3 == nodes[j]:
                                                i3 = j
                                if i1 >= 0 and i2 == -1 and i3 == -1 and n2 in pars[i1] and n3 in pars[i1]:
                                        lhs += z_value[d,dp]
                                        continue
                                if i2 >= 0 and i1 == -1 and i3 == -1 and n1 in pars[i2] and n3 in pars[i2]:
                                        lhs += z_value[d,dp]
                                        continue
                                if i3 >= 0 and i1 == -1 and i2 == -1 and n1 in pars[i3] and n2 in pars[i3]:
                                        lhs += z_value[d,dp]
                                        continue
                                if i1 >= 0 and i2 >= 0 and i3 >= 0:
                                        lhs += z_value[d,dp]
                                        continue
                                if (i1 >= 0 and i2 >= 0 and i3 == -1) and (n3 in pars[i1] or n3 in pars[i2]):
                                        lhs += z_value[d,dp]
                                        continue
                                if (i1 >= 0 and i3 >= 0 and i2 == -1) and (n2 in pars[i1] or n2 in pars[i3]):
                                        lhs += z_value[d,dp]
                                        continue
                                if (i2 >= 0 and i3 >= 0 and i1 == -1) and (n1 in pars[i2] or n1 in pars[i3]):
                                        lhs += z_value[d,dp]
                                        continue
                        if lhs > 1.0+1.0e-6:
                                few_3node_subsets.append(subset)
                # loop to generate new cuts
                for subset in few_3node_subsets:
                        n1 = subset[0]
                        n2 = subset[1]
                        n3 = subset[2]
                        ifLHS = {(d,dp):False for d in range(len(self.cComps)) for dp in range(len(self.dPars[d]))}
                        for d in range(len(self.cComps)):
                                for dp in range(len(self.dPars[d])):
                                        nodes = self.cComps[d][0]
                                        pars = self.dPars[d][dp]
                                        i1 = -1 # this is the index of n1 in cComps[d][0] (if n1 is not there, then i1=-1)
                                        i2 = -1
                                        i3 = -1
                                        for j in range(len(nodes)):
                                                if n1 == nodes[j]:
                                                        i1 = j
                                                elif n2 == nodes[j]:
                                                        i2 = j
                                                elif n3 == nodes[j]:
                                                        i3 = j
                                        if i1 >= 0 and n2 in pars[i1] and n3 in pars[i1]:
                                                ifLHS[d,dp] = True
                                                continue
                                        if i2 >= 0 and n1 in pars[i2] and n3 in pars[i2]:
                                                ifLHS[d,dp] = True
                                                continue
                                        if i3 >= 0 and n1 in pars[i3] and n2 in pars[i3]:
                                                ifLHS[d,dp] = True
                                                continue
                                        if i1 >= 0 and i2 >= 0 and i3 >= 0:
                                                ifLHS[d,dp] = True
                                                continue
                                        if (i1 >= 0 and i2 >= 0) and (n3 in pars[i1] or n3 in pars[i2]):
                                                ifLHS[d,dp] = True
                                                continue
                                        if (i1 >= 0 and i3 >= 0) and (n2 in pars[i1] or n2 in pars[i3]):
                                                ifLHS[d,dp] = True
                                                continue
                                        if (i2 >= 0 and i3 >= 0) and (n1 in pars[i2] or n1 in pars[i3]):
                                                ifLHS[d,dp] = True
                                                continue
                        addcut = True
                        if addcut == True:
                                debug = False
                                if debug == True:
                                        print('-------- 3-node set cut -------')
                                        print('set: '+str(n1)+', '+str(n2)+', '+str(n3))
                                        for d in range(len(self.cComps)):
                                                for dp in range(len(self.dPars[d])):
                                                        if ifLHS[d,dp] == True:
                                                                print('dist: '+str(self.cComps[d])+', pars: '+str(self.dPars[d][dp]))
                                        print('-------------------------------')
                                zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
                                self.lp.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set, val = [1]*len(zindex_set))], senses=["L"], rhs=[1])
                                nnewcuts += 1
                print('Number of new 3-node set cuts added: '+str(nnewcuts)+', time: '+str(time.time()-t0cuts))

        def cuts_facets_3_node_set(self, z_value=None):
                print('c: cuts_facets_3_node_set')
                t0cuts = time.time()
                if z_value == None:
                        print('cuts_facets_3_node_set(...) works only if z_value provided')
                        return
                icc = []
                for d in range(len(self.cComps)):
                        for dp in range(len(self.dPars[d])):
                                if z_value[d,dp] > 1.0e-6:
                                        icc.append((d,dp))
                nnewcuts = 0
                list_nodes = [i for i in range(len(self.V))]
                all_3node_subsets = [list(i) for i in itertools.combinations(list_nodes,3)]
                few_3node_subsets = []
                # loop to select 3-node subsets where constraint is violated
                for subset in all_3node_subsets:
                        n1 = subset[0]
                        n2 = subset[1]
                        n3 = subset[2]
                        lhs = 0.0
                        for i in range(len(icc)):
                                d = icc[i][0]
                                dp = icc[i][1]
                                nodes = self.cComps[d][0]
                                pars = self.dPars[d][dp]
                                i1 = -1 # this is the index of n1 in cComps[d][0] (if n1 is not there, then i1=-1)
                                i2 = -1
                                i3 = -1
                                for j in range(len(nodes)):
                                        if n1 == nodes[j]:
                                                i1 = j
                                        elif n2 == nodes[j]:
                                                i2 = j
                                        elif n3 == nodes[j]:
                                                i3 = j
                                if (i1 >= 0 and i2 >= 0) or (i1 >= 0 and i3 >= 0) or (i2 >= 0 and i3 >= 0): 
                                        lhs += z_value[d,dp]
                                        continue
                        if lhs > 1.0+1.0e-6:
                                few_3node_subsets.append(subset)
                # loop to generate new cuts
                for subset in few_3node_subsets:
                        n1 = subset[0]
                        n2 = subset[1]
                        n3 = subset[2]
                        ifLHS = {(d,dp):False for d in range(len(self.cComps)) for dp in range(len(self.dPars[d]))}
                        for d in range(len(self.cComps)):
                                for dp in range(len(self.dPars[d])):
                                        nodes = self.cComps[d][0]
                                        pars = self.dPars[d][dp]
                                        i1 = -1 # this is the index of n1 in cComps[d][0] (if n1 is not there, then i1=-1)
                                        i2 = -1
                                        i3 = -1
                                        for j in range(len(nodes)):
                                                if n1 == nodes[j]:
                                                        i1 = j
                                                elif n2 == nodes[j]:
                                                        i2 = j
                                                elif n3 == nodes[j]:
                                                        i3 = j
                                        if (i1 >= 0 and i2 >= 0) or (i1 >= 0 and i3 >= 0) or (i2 >= 0 and i3 >= 0): 
                                                ifLHS[d,dp] = True
                                                continue
                        addcut = True
                        if addcut == True:
                                debug = False
                                if debug == True:
                                        print('-------- 3-node set facets cut -------')
                                        print('set: '+str(n1)+', '+str(n2)+', '+str(n3))
                                        for d in range(len(self.cComps)):
                                                for dp in range(len(self.dPars[d])):
                                                        if ifLHS[d,dp] == True:
                                                                print('dist: '+str(self.cComps[d])+', pars: '+str(self.dPars[d][dp]))
                                        print('-------------------------------')
                                zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
                                self.lp.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set, val = [1]*len(zindex_set))], senses=["L"], rhs=[1])
                                nnewcuts += 1
                print('Number of new 3-node set facets cuts added: '+str(nnewcuts)+', time: '+str(time.time()-t0cuts))


        def RunLP(self,t0):
                ContinueCondt = True
                LPiter = 0
                Objvalue = float('inf')
                nvar = sum(len(self.dPars[d]) for d in range(len(self.cComps)))
                nbi = len(self.udE)
                ncluster = 0
                nbicluster = 0
                
                z_value = {} # this stores the values of the LP variables
                
                MaxLPiter = 100
                while ContinueCondt == True and LPiter < MaxLPiter:
                        nz_bi = 0
                        nz_z = 0
                        ContinueCondt = False
                        LPiter = LPiter+1
                        model_soln = self.lp.solve()
                        PrevObjvalue = Objvalue
                        Objvalue = self.lp.solution.get_objective_value()
                        
                        for d in range(len(self.cComps)):
                                for dp in range(len(self.dPars[d])):
                                        if self.lp.solution.get_values(['z'+str(d)+','+str(dp)])[0] > 0:
                                                nz_z = nz_z+1
                                                if self.lp.solution.get_values(['z'+str(d)+','+str(dp)])[0] < 1.0:
                                                        print('noninteger '+'z'+str(d)+','+str(dp)+' = '+str(self.lp.solution.get_values(['z'+str(d)+','+str(dp)])[0])+'; '+str(self.cComps[d])+': '+str(self.dPars[d][dp])+', score: '+str(self.scores[self.cComps[d]][self.dPars[d][dp]]))
                                                else:
                                                        print('equal to one '+'z'+str(d)+','+str(dp)+' = '+str(self.lp.solution.get_values(['z'+str(d)+','+str(dp)])[0])+'; '+str(self.cComps[d])+': '+str(self.dPars[d][dp])+', score: '+str(self.scores[self.cComps[d]][self.dPars[d][dp]]))

                        for e in self.udE:
                                if self.lp.solution.get_values(['bi'+str(e)])[0] > 0:
                                        nz_bi = nz_bi+1
                        
                        print('LP iter '+str(LPiter)+', ObjVal: '+str(Objvalue)+', time: '+str(time.time()-t0)+', frac. of nonzero variables: '+\
                                str(nz_z)+'/'+str(nvar)+', frac. of nonzero bidirected edges: '+str(nz_bi)+'/'+str(nbi)+', # cluster: '+str(ncluster)+', # bi-cluster: '+str(nbicluster))

                
                        ncluster = 0
                        nbicluster = 0
                        bi_value = {}
                        x_value = {}
                        for e in self.udE:
                                bi_value[e] = self.lp.solution.get_values(['bi'+str(e)])[0]
                        for i in self.V:
                                for ip in range(len(self.iPars[i])):
                                        x_value[i,ip] = self.lp.solution.get_values(['x'+str(i)+','+str(ip)])[0]
                        for d in range(len(self.cComps)):
                                for dp in range(len(self.dPars[d])):
                                        z_value[d,dp] = self.lp.solution.get_values(['z'+str(d)+','+str(dp)])[0]

                        # -------------------------------
                        # new 3-node set cuts
                        if 'b' in self.cuts:
                                self.cuts_3_node_set(z_value)
                        # -------------------------------
                        # new 3-node set facets cuts
                        if 'c' in self.cuts:
                                self.cuts_facets_3_node_set(z_value)
                        # -------------------------------
                        
                        wt = {}
                        for i in self.V:
                                for ip in range(len(self.iPars[i])):
                                        if x_value[i,ip] > 1e-6:
                                                for par in self.iPars[i][ip]:
                                                        if (par,i) not in wt.keys():
                                                                wt[(par,i)] = x_value[i,ip]
                                                        else:
                                                                wt[(par,i)] = wt[(par,i)]+x_value[i,ip]
                        telist = []
                        for (i,j) in wt.keys():
                                if wt[(i,j)] >= 1-1e-6:
                                        telist.append(i)
                                        telist.append(j)
                        allcyc = []
                        
                        allcyc = dircyc(len(self.V),int(len(telist)/2),telist)
                        for Cluster in allcyc:
                                ifLHS = self.ClusterToIneq(Cluster)
                                lhs = 0.0
                                for d in range(len(self.cComps)):
                                        for dp in range(len(self.dPars[d])):
                                                if ifLHS[d,dp] == True:
                                                        lhs += z_value[d,dp]
                                                        if z_value[d,dp] > 0.0:
                                                                print('z['+str(self.cComps[d][0])+','+str(self.dPars[d][dp])+']= '+str(z_value[d,dp]))
                                zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
                                self.lp.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set, val = [1]*len(zindex_set))], senses=["G"], rhs=[1])
                                
                        if len(allcyc) > 0:
                                ncluster = len(allcyc)
                                ContinueCondt = True
                        else:
                                nnode = len(self.V)
                                gcnodes = []
                                gcweight = []
                                gcparents = []
                                for d in range(len(self.cComps)):
                                        for dp in range(len(self.dPars[d])):
                                                if z_value[d,dp] > 0:
                                                        gcnodes.append(list(self.cComps[d][0]))
                                                        gcweight.append(z_value[d,dp])
                                                        gcparents.append(list(list(pars) for pars in self.dPars[d][dp]))
                                for i in range(len(gcparents)):
                                        for j in range(len(gcparents[i])):
                                                if len(gcparents[i][j]) == 0:
                                                        gcparents[i][j].append(nnode)
                                dircycs = contract_heur(nnode+1, gcnodes, gcparents, gcweight,1.0)
                                for Cluster in dircycs:
                                        ifLHS = self.ClusterToIneq(Cluster)
                                        lhs = 0.0
                                        for d in range(len(self.cComps)):
                                                for dp in range(len(self.dPars[d])):
                                                        if ifLHS[d,dp] == True:
                                                                lhs += z_value[d,dp]
                                                                if z_value[d,dp] > 0.0:
                                                                        print('z['+str(self.cComps[d][0])+','+str(self.dPars[d][dp])+']= '+str(z_value[d,dp]))
                                        zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
                                        self.lp.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set, val = [1]*len(zindex_set))], senses=["G"], rhs=[1])

                                ncluster = ncluster+len(dircycs)

                                if self.bowfree == False:
                                        aldircycs = contract_heur_bdir(nnode+1, gcnodes, gcparents, gcweight)
                                        for Cluster in aldircycs:
                                                ii = Cluster[0]
                                                jj = Cluster[1]
                                                e = self.indInv.index((ii,jj))
                                                C = set(Cluster[2:])
                                                ifLHS = self.biClusterToIneq(C,ii,jj)
                                                zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
                                                self.lp.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = zindex_set+['bi'+str(e)], val = [1]*len(zindex_set)+[-1])], senses=["G"], rhs=[0])
                                        nbicluster = nbicluster+len(aldircycs)

                                if ncluster+nbicluster > 0 and Objvalue < PrevObjvalue:
                                        ContinueCondt = True

                        if ContinueCondt == False:
                                print('Last LP iter '+str(LPiter)+', ObjVal: '+str(Objvalue)+', time: '+str(time.time()-t0)+', frac. of nonzero variables: '+str(nz_z)+'/'+str(nvar)+', frac. of nonzero bidirected edges: '+str(nz_bi)+'/'+str(nbi)+', # cluster: '+str(ncluster)+', # bi-cluster: '+str(nbicluster))
                                        
                return z_value


        def RunMILP(self,t0):
                # starts MILP
                self.milp = cplex.Cplex(self.lp)                
                self.milp.set_problem_type(self.milp.problem_type.MILP)
                # make variables z binary
                for d in range(len(self.cComps)):
                        for dp in range(len(self.dPars[d])):
                                if self.dag == False or len(self.cComps[d][0]) <= 1:
                                        self.milp.variables.set_types('z'+str(d)+','+str(dp),self.milp.variables.type.binary)
                                else:
                                        self.milp.variables.set_types('z'+str(d)+','+str(dp),self.milp.variables.type.binary)
                self.milp.parameters.mip.limits.nodes = 9223372036800000000
                tRoot = time.time()-t0

                class LazyCallback(cplex.callbacks.LazyConstraintCallback):
                        def __init__(self, env):
                                super().__init__(env)
                        
                        def __call__(self):
                                NoCluster = False
                                bi_value = {}
                                x_value = {}
                                for e in self.udE:
                                        bi_value[e] =  self.get_values(['bi'+str(e)])[0]
                                for i in self.V:
                                        for ip in range(len(self.iPars[i])):
                                                x_value[i,ip] = self.get_values(['x'+str(i)+','+str(ip)])[0]
                                
                                ActiveEdgeList = []
                                
                                for i in self.V:
                                        for ip in range(len(self.iPars[i])):
                                                if x_value[i,ip] > 0.5:
                                                        for par in self.iPars[i][ip]:
                                                                ActiveEdgeList.append(par)
                                                                ActiveEdgeList.append(i)
                                ne = int(len(ActiveEdgeList)/2)

                                # Detecting directed cycles and adding cluster inequalities
                                cycList = dircyc(len(self.V),ne,ActiveEdgeList)
                                for Cluster in cycList:
                                        ifLHS = self.ClusterToIneq(Cluster)
                                        zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
                                        self.add(constraint=cplex.SparsePair(ind = zindex_set, val = [1]*len(zindex_set)), sense="G", rhs=1)

                                # Detecting almost directed cycles and adding bi-cluster inequalities
                                if self.bowfree == False:
                                        for e in self.udE:
                                                if bi_value[e] > 0.5:
                                                        ii = self.indInv[e][0]
                                                        jj = self.indInv[e][1]
                                                        adcycList = almostdircyc(len(self.V),ne,ActiveEdgeList,ii,jj)+almostdircyc(len(self.V),ne,ActiveEdgeList,jj,ii)
                                                        for Cluster in adcycList:
                                                                C = set(Cluster[1:-1])
                                                                ifLHS = self.biClusterToIneq(C,ii,jj)
                                                                zindex_set = ['z'+str(d)+','+str(dp) for d in range(len(self.cComps)) for dp in range(len(self.dPars[d])) if ifLHS[d,dp]==True]
                                                                self.add(constraint=cplex.SparsePair(ind = zindex_set+['bi'+str(e)], val = [1]*len(zindex_set)+[-1]), sense="G", rhs=0)

                
                self.milp.parameters.preprocessing.presolve.set(0)
                self.milp.parameters.mip.strategy.search.set(1)
                self.milp.parameters.threads.set(32)
                self.milp.parameters.timelimit.set(self.milp_time_limit)
                self.milp.parameters.mip.tolerances.absmipgap.set(0.0) # absolute mip gap
                self.milp.parameters.mip.tolerances.mipgap.set(0.0) # relative mip gap
                lazyModel = self.milp.register_callback(LazyCallback)
                lazyModel.udE = self.udE
                lazyModel.bi = self.bi
                lazyModel.iPars = self.iPars
                lazyModel.x = self.x
                lazyModel.z = self.z
                lazyModel.V = self.V
                lazyModel.dPars = self.dPars
                lazyModel.cComps = self.cComps
                lazyModel.indInv = self.indInv
                lazyModel.ClusterToIneq = self.ClusterToIneq
                lazyModel.biClusterToIneq = self.biClusterToIneq
                lazyModel.bowfree = self.bowfree
                print("start last solve")
                #fileName = self.resultsdir+self.instance+'_cplex_milp.lp'
                #self.milp.write(fileName)
                self.milp.solve()
                print('Score: '+str(self.milp.solution.get_objective_value()))
                return tRoot
        
        def Solve_with_cb(self,CG=False):
                t0 = time.time()

                if self.heuristics == "":
                        z_value = self.RunLP(t0)
                else:
                        iter = 0
                        maxiter = 5
                        while iter < maxiter:
                        
                                z_value = self.RunLP(t0)
                                iter += 1
                        
                                if iter < maxiter:
                                        # adds combinations of c-components and creates the corresponding LP
                                        self.CombineCComponents(z_value)
                                        self.CreateLP()
                
                #self.write_ccomps_to_file()
                tRoot = self.RunMILP(t0)

                fileName = self.resultsdir+self.instance+'_cplex.log'
                wrtStr = 'Time at root: '+str(tRoot)+'\n'
                wrtStr = wrtStr+'Total solution time (LPs + MIP): '+str(time.time()-t0)+'\n'
                wrtStr = wrtStr+'Score: '+str(self.milp.solution.get_objective_value())+'\n'
                print('Score: '+str(self.milp.solution.get_objective_value()))
                print('Bidirected edges: ')
                wrtStr = wrtStr+'Bidirected edges: \n'
                for e in self.udE:
                        if self.milp.solution.get_values(['bi'+str(e)])[0] > 0.5:
                                print(self.indInv[e])
                                wrtStr = wrtStr+str(self.indInv[e])+'\n'
                print('Parent sets: ')
                wrtStr = wrtStr+'Parent sets: \n'
                for i in self.V:
                        for ip in range(len(self.iPars[i])):
                                if self.milp.solution.get_values(['x'+str(i)+','+str(ip)])[0] > 0.5:
                                        print(str(i)+': '+str(self.iPars[i][ip]))
                                        wrtStr = wrtStr+str(i)+': '+str(self.iPars[i][ip])+'\n'
                print('z solution: ')
                wrtStr = wrtStr+'z solution: '
                for d in range(len(self.cComps)):
                        for dp in range(len(self.dPars[d])):
                                if self.milp.solution.get_values(['z'+str(d)+','+str(dp)])[0] > 0.5:
                                        print(str(self.cComps[d])+': '+str(self.dPars[d][dp])+', score: '+str(self.scores[self.cComps[d]][self.dPars[d][dp]]))
                                        wrtStr = wrtStr+str(self.cComps[d])+': '+str(self.dPars[d][dp])+', score: '+str(self.scores[self.cComps[d]][self.dPars[d][dp]])+'\n'
                f = open(fileName,"a")
                f.write(wrtStr)
                f.close

        def write_graph_to_file(self):
                nnodes = len(self.V)
                D = [[0]*nnodes for i in range(nnodes)]
                B = [[0]*nnodes for i in range(nnodes)]

                # get directed edges from parent sets
                for i in self.V:
                        for ip in range(len(self.iPars[i])):
                                if self.milp.solution.get_values(['x'+str(i)+','+str(ip)])[0] > 0.5:
                                        for j in range(len(self.iPars[i][ip])):
                                                D[self.iPars[i][ip][j]][i] = 1

                # get bidirected edges
                for e in self.udE:
                        if self.milp.solution.get_values(['bi'+str(e)])[0] > 0.5:
                                i = self.indInv[e][0]
                                j = self.indInv[e][1]
                                B[i][j] = 1
                                B[j][i] = 1

                graph = {}
                graph[0] = D
                graph[1] = B

                graph_file_name = self.resultsdir+'output_graph_'+self.instance+'_cplex.log'
                write_graph_to_file(graph, graph_file_name)

        def write_ccomps_to_file(self):
                wrtStr = ''
                for d in range(len(self.cComps)):
                        for dp in range(len(self.dPars[d])):
                                wrtStr = wrtStr+str(self.cComps[d])+': '+str(self.dPars[d][dp])+', score: '+str(self.scores[self.cComps[d]][self.dPars[d][dp]])+'\n'
                fileName = self.resultsdir+'all_ccomps_'+self.instance+'.txt'
                f = open(fileName,"w")
                f.write(wrtStr)
                f.close


        def write_samples_and_scores_to_files(self):
                scores_file_name = self.datadir+self.instance+'_scoresread.txt'
                f = open(scores_file_name,"w")
                Dpars = {}
                for D in self.originalScores.keys():
                        Dpars[D] = list(self.originalScores[D].keys())
                        for Dpar in Dpars[D]:
                                f.write(str(D)+' '+str(Dpar)+' '+str(self.originalScores[D][Dpar])+'\n')
                f.close()
                sample_file_name = self.datadir+self.instance+'_sampleread.txt'
                nrows,ncols = self.data.shape
                wrt = ''
                for i in range(nrows):
                        for j in range(ncols-1):
                                wrt = wrt + str(self.data[i][j]) + ' '
                        wrt = wrt + str(self.data[i][ncols-1]) + '\n'
                f = open(sample_file_name,"w")
                f.write(wrt)
                f.close()



if __name__ == '__main__':
        random.seed(1)
        t0 = time.time()
        datadir = ''
        resultsdir = ''
        # heuristic a: add a node to the parent set of the c-component
        # heuristic b: add a node to the district of the c-component
        # heuristic c: Move a node that is currently a parent to the district of the c-component
        heuristics = '' #'abc'
        # cuts a: cuts_2_node_set
        # cuts b: cuts_3_node_set
        # cuts c: cuts_facets_3_node_set
        cuts = '' #'abc'
        bowfree = True
        arid = False
        milp_time_limit = 3600
        scoresets = ['']

        if len(sys.argv) <= 1:
                print('To run, type:')
                print('python3 learn.py -s score_file_name.pkl -i path_input_directory/ -o path_output_directory/ -h heuristics -c cuts -t MIP_time_limit_(seconds) -g type_of_graph_desired')
                print('where:')
                print('score_file_name.pkl is a file with input c-components and corresponding scores,')
                print('path_input_directory/ is the path to the directory where score_file_name.pkl resides,')
                print('path_output_directory/ is the path to the directory where the output is written,')
                print('heuristics is a string that can be any combination of the letters a, b, c (e.g., abc or ac),')
                print('cuts is a string that can be any combination of the letters a, b, c (e.g., abc or ac),')
                print('MIP_time_limit_(seconds) is an integer with the time limit for the MIP,')
                print('type_of_graph_desired is either aadmg, arid, or bowfree.')
                exit()

        for i in range(1,len(sys.argv)):
                temp = sys.argv[i]
                if temp == '-s':
                        scoresets = [sys.argv[i+1]]
                elif temp == '-i':
                        datadir = sys.argv[i+1]
                elif temp == '-o':
                        resultsdir = sys.argv[i+1]
                elif temp == '-h':
                        heuristics = sys.argv[i+1]
                elif temp == '-c':
                        cuts = sys.argv[i+1]
                elif temp == '-t':
                        milp_time_limit = int(sys.argv[i+1])
                elif temp == '-g':
                        if sys.argv[i+1] == 'aadmg':
                                bowfree = False
                                arid = False
                        if sys.argv[i+1] == 'arid':
                                bowfree = True
                                arid = True
                        if sys.argv[i+1] == 'bowfree':
                                bowfree = True
                                arid = False
                        

        for instName in scoresets:
                inst = BNSLlvInst(instName, datadir, resultsdir, heuristics, cuts, milp_time_limit, bowfree, arid)
                
                inst.readFromPkl()
                #inst.readFromDag()
                inst.Initialize(prune=True,printsc=False)
                #inst.Initialize(prune=True,dag=True)
                inst.Solve_with_cb()
                inst.write_graph_to_file()

        print('Total time including pruning: '+str(time.time()-t0)+'\n')
