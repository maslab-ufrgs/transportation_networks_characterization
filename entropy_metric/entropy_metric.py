# -*- coding: utf-8 -*-

import math

import random
from random_allocation import random_allocation

import numpy as np
import matplotlib.pyplot as plt

from KSP3 import KShortestPaths

from successive_averages import dijkstra
from successive_averages import generateGraph
from successive_averages import getPathAsEdges
from successive_averages import pathToStr

import argparse

class SmartTable:
    
    def __init__(self):
        
        self._table = dict()        
        self._preference = None
        
    def updateTable(self, key, value):
        self._table[key] = value
        
    def updatePreference(self, key):
        self._preference = key
        
    def get(self, key):
        
        if key not in self._table:
            return 1000000000
        elif key in self._table:
            return self._table[key]
    
    def preference(self):
        
        if (self._preference == None):
            
            if self._table == dict():
                
                return None
            
            elif self._table != dict():
                
                for x in self._table:
                    # Return first element
                    # of self._table
                    return x

        elif self._preference != None:
            
            return self._preference
        
    def __repr__(self):
        return str(self._table)


class Dtmc:
    
    def __init__(self):
        
        self.__dtmc = dict()
        self.__visiting = dict()
        self.__counter  = 0
    
    def addOne(self, i, j=None):
        
        self.__counter += 1
        
        if j == None:
            j = i
            
        # Atualizar as visitações de cada estado
        if j not in self.__visiting:
            self.__visiting[j] = 1
        elif j in self.__visiting:
            self.__visiting[j] += 1
            
        if i not in self.__dtmc:
            
            self.__dtmc[i] = dict()
            self.__dtmc[i][j] = 1
            
        elif i in self.__dtmc:
            
            if j in self.__dtmc[i]:
                
                self.__dtmc[i][j] += 1
            
            elif j not in self.__dtmc[i]:
                
                self.__dtmc[i][j] = 1
                
    def get(self, i, j=None):
        
        if j == None:
            
            j = i
            
        if i not in self.__dtmc:
            
            return None
            
        elif i in self.__dtmc:
            
            if j in self.__dtmc[i]:
                
                return self.__dtmc[i][j]
            
            elif j not in self.__dtmc[i]:
                
                return None
            
    def __repr__(self):
        
        return str(self.__dtmc)
    
    
    def getVisiting(self):
        return self.__visiting
    
    def entropy(self):
        
        if self.__counter == 0:
            return 0.0
        
        entropy = 0.0
        
        route_choice = self.getVisiting()
        
        for i in route_choice:
            
            mu = route_choice[i] / float(self.__counter)
                
            if i not in self.__dtmc:
                
                continue
            
            elif i in self.__dtmc:
                
                if self.__dtmc[i] == dict():
                    
                    continue
                
                #sum_Pij = float(sum([self.__dtmc[i][j] for j in self.__dtmc[i]]))
                
                for j in self.__dtmc[i]:
                    
                    #Pij = self.__dtmc[i][j] / sum_Pij
                    Pij = self.__dtmc[i][j] / route_choice[i]
                    
                    entropy -= mu * Pij * math.log(Pij, 2.0)

        return entropy
    

    def asymptoty(self):

        if self.__visiting == dict():
            
            return dict()
        
        mu = dict()
        
        for j in self.__visiting:
            
            mu[j] = self.__visiting[j] / float(self.__counter)
            
        return mu
    
    def probability(self):
        
        my_dict = dict()
        
        for i in self.__dtmc:
            
            sum_Pij = float(sum([self.__dtmc[i][j] for j in self.__dtmc[i]]))
            
            if sum_Pij == 0:
                
                my_dict[i] = dict()
                
            elif sum_Pij > 0:
                
                if i not in my_dict:
                    
                    my_dict[i] = dict()
                
                for j in self.__dtmc[i]:
                    
                    Pij = self.__dtmc[i][j] / sum_Pij
                    
                    my_dict[i][j] = Pij
                    
        return my_dict
            

class entropy_metric:
    
    def __init__(self, V, E, OD, k, Ep, ep, ksp=True):
        
        # Set of nodes (class Nodes imported from successive_averages.py)
        self.V  = V
        
        # dict of edges (class Nodes imported from successive_averages.py)        
        self.E  = dict()
        for edge in E:
            self.E[edge.name] = edge
        
        # dict of od pairs of the form
        '''
        {
            key1 : demand of key1,
            key2 : demand of key2,
            ...,
            keyn : demand of keyn
        }
        '''
        # where keys are of the form "origin|destination"
        self.OD = OD
        
        # dictionary of routes with od labels as keys
        self.routes  = dict()
        for od in self.OD:
            # each self.routes[od] is a dictionary of routes
            self.routes[od] = dict()
        
        self.new_routes = dict()
        self.new_routes_str = dict()
        
        # It may be the number of k shortest paths used
        # if ksp flag is True, but it can also be the
        # maximum number of routes if ksp argument is False
        self.k           = int(k)
        
        # Episodes for calculating the entropy        
        self.Episodes    = int(Ep)
        
        # Episodes for discovery of new routes
        self.episodes    = int(ep)
        
        # Episode counter variable
        self.t = 0        
        
        # Initialize smart tables, one for each OD pair
        self.smartTables = dict()
        for od in self.OD:
            self.smartTables[od] = SmartTable()
        
        # Initialize entropies dictionary and Dtmc,
        # one for each OD pair.
        self.entropy      = dict()
        self.dtmc         = dict()
        for od in self.OD:
            # For each OD pair, keep a list of the entropy values
            self.entropy[od]      = list()
            # For each OD pair, keep last Dtmc state
            self.dtmc[od]         = Dtmc()
            
        # If ksp == True, calculate the k shortest paths
        self.ksp = ksp
        if self.ksp == True:
            
            for od in self.OD:
                
                [o,d] = od.split('|')
                
                k_shortest_paths = KShortestPaths(self.V, E,
                                                  o,d,self.k)
                
                for sp in k_shortest_paths:
                    
                    key = pathToStr(sp, self.V, E)
                    route_edge = getPathAsEdges(sp, E)
                
                    self.routes[od][key] = route_edge

    def search_iteration(self):
                
        for od in self.OD:
            
            if len(self.routes[od]) == self.k:
                # If maximum number of routes reached
                continue
            
            [o, d] = od.split("|")
            
            edge_list = list(self.E.values())
            
            route = getPathAsEdges(dijkstra(self.V,edge_list,o,d,[]),edge_list)
            route_str = pathToStr(route,self.V,edge_list)
            
            if route_str in self.routes[od]:
                continue
            elif route_str not in self.routes[od]:
                self.routes[od][route_str] = route
        
    def assignment_iteration(self):
        
        # Reset flow to 0 on each edge
        for edge in self.E:
            self.E[edge].flow = 0
        
        for od in self.OD:
            
            # New random assingment for OD pair od         
            random_flow = random_allocation(int(self.OD[od]),
                                            len(self.routes[od]))
            
            for z, route_key in enumerate(self.routes[od]):
                flow = random_flow[ z ]
                for edge in self.routes[od][route_key]:
                    self.E[edge.name].flow += flow
                
        # Update cost of each edge
        for edge in self.E:
            self.E[edge].update_cost()
        
        key = dict()
        
        for od in self.OD:

            # Sorteia uma rota para atualizar o custo
            key[od] = random.choice( list(self.routes[od].keys()) )
                
            # Calcula o custo da rota sorteada
            cost = 0.0
            for edge in self.routes[od][key[od]]:
                cost += self.E[edge.name].cost
            # -----------------------
                
            # Atualiza o custo dessa rota
            self.smartTables[od].updateTable(key[od], cost)
                
        for od in self.OD:
                
            if self.smartTables[od].get(key[od]) < \
                   self.smartTables[od].get(self.smartTables[od].preference()):
                    
                # Atualizar o dtmc
                self.dtmc[od].addOne(self.smartTables[od].preference(), key[od])
                    
                # Atualizar preferencia
                self.smartTables[od].updatePreference(key[od])
                    
            else:
                pref = self.smartTables[od].preference()
                # Atualizar o dtmc
                self.dtmc[od].addOne(pref, pref) 

    def entropy_iteration(self):
        # Cálculo das entropias
        for od in self.OD:
            self.entropy[od].append(self.dtmc[od].entropy())

    def run(self):
        
        if (self.ksp == True):
            
            while (self.t <= self.Episodes):
            
                self.t += 1
                self.assignment_iteration()
                self.entropy_iteration()  
        
        elif (self.ksp == False):
            
            while (self.t <= self.episodes):
            
                self.t += 1
                self.search_iteration()
                self.assignment_iteration()
                
            self.t = 0
            for od in self.OD:
                self.dtmc[od] = Dtmc()
            
            while (self.t <= self.Episodes):
            
                self.t += 1                
                self.assignment_iteration()
                self.entropy_iteration()                        
            
    def iteration(self):
        
        self.search_iteration()
        self.assignment_iteration()
        self.entropy_iteration() 
        
    def plotEntropy(self, od=None):
        
        if od == None:
            x = np.arange(self.t)
            for od in self.OD:
                plt.plot(x, np.array(self.entropy[od]))
            plt.legend([od for od in self.OD], loc='upper right')
            plt.show()
            
        elif od in self.OD:
            x = np.arange(self.t)
            plt.plot(x, np.array(self.entropy[od]))
            plt.legend([od,], loc='upper right')
            plt.show()            
        
    def route_choice(self, od):
        
        return self.dtmc[od].getVisiting()
        

if __name__ == '__main__':
    
    prs = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description="""The entropy metric.""")

    prs.add_argument("-f", dest="network_file",
                     required=True, help="The network file.\n")

    prs.add_argument("-K", dest="maxroutes", required=False,
                     default=0, type=int,
                      help="Max number of routes for each OD pair.\n")
    
    prs.add_argument("-k", dest="ksp", required=False, type=int,
                     default=0, help="Use or not ksp.\n")    
    
    prs.add_argument("-E", dest="Episodes", type=int,
                     required=True, default=1000,
                      help="Number of episodes calculating entropy.\n")
    
    prs.add_argument("-e", dest="episodes", type=int,
                     required=False, default=0,
                      help="Number of episodes estimating new route.\n")
    
    args = prs.parse_args()

    V, E, OD = generateGraph(args.network_file)
    
    
    if args.maxroutes == 0 and args.ksp == 0:
        # Parâmetros de chamada inconsistentes
        print('Must set ksp or kmax flag!')
    
    elif args.maxroutes > 0 and args.ksp > 0:
        # Parâmetros de chamada inconsistentes
        print("Can't set ksp and kmax simultaneosly!")
    
    elif args.maxroutes == 0 and args.ksp > 0:
            
        # Chamada com o ksp
        H = entropy_metric(V, E, OD,
                      args.ksp,
                      args.Episodes,
                      args.episodes, ksp=True)
        H.run()

        for od in OD:
            print(od, round(H.entropy[od][-1],3))
                    
    elif args.maxroutes > 0 and args.ksp == 0:
        
        # Chamada com o descubrimento de rotas
        H = entropy_metric(V, E, OD,
                       args.maxroutes,
                       args.Episodes,
                       args.episodes, ksp=False)
        H.run()
        
        for od in OD:
            print(od, H.entropy[od][-1])
                
