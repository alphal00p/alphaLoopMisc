import numpy as np
import itertools
import time
import os
import sys
import copy
from functools import wraps 


def adjacentedge(e1,e2):
    if (e1[0]==e2[0] and e1[1]!=e2[1]) or (e1[1]==e2[1] and e1[0]!=e2[0]) or (e1[1]==e2[0] and e1[0]!=e2[1]) or (e1[0]==e2[1] and e1[1]!=e2[0]):
        return True
    
    return False

def adjacentVertex(v1,v2,eset):
    for e in list(eset):
        if (e[0]==v1 and e[1]==v2) or (e[1]==v1 and e[0]==v2):
            return 1
    
    return 0

def adjacentset(vset,eset,vertex):
    adjacent=set()

    for e in list(eset):
        if e[0]==vertex or e[1]==vertex:
            adjacent.add(e)

    return adjacent

def spanning_tree_finder(vset,eset):

    spanning_tree_v={list(vset)[0]}
    spanning_tree_e=set()

    for i in range(0,len(vset)):
        
        for v in list(vset):
            adj=adjacentset(vset,eset,v)
            
            for e_adj in list(adj):

                if ((e_adj[0] in spanning_tree_v) and (e_adj[1] not in spanning_tree_v)) or ((e_adj[1] in spanning_tree_v) and (e_adj[0] not in spanning_tree_v)):
                    spanning_tree_e.add(e_adj)
                    spanning_tree_v.add(e_adj[1])
                    spanning_tree_v.add(e_adj[0])

    return [spanning_tree_v,spanning_tree_e]


    return spanning_tree_v

                    
def generators_finder(vset,eset):

    spanning_tree=spanning_tree_finder(vset,eset)   

    cycle_sets=[]

    for e in list(eset):

        if e not in spanning_tree[1]:
            
            cycle=set()
            st_cycle=set(spanning_tree[1])
            st_cycle.add(e)

            for e1 in list(st_cycle):

                new_st_cycle=set(st_cycle)
                new_st_cycle.remove(e1)
                st_spanning_tree=spanning_tree_finder(vset, new_st_cycle)[1]

                if len(st_spanning_tree)==len(new_st_cycle):
                    cycle.add(e1)

            cycle.add(e)
            cycle_sets.append(cycle)

    return cycle_sets

def cycle_merger(cycle1,cycle2):

    return cycle1.union(cycle2).difference(cycle1.intersection(cycle2))
        

def powerset(s):
    x = len(s)
    powa=[]
    for i in range(1 << x):
        powa.append([s[j] for j in range(x) if (i & (1 << j))])
    return powa

def winding_cycle_finder(vset,eset,cycle_set,w_cycle):

    zero_w_cycle=[]
    one_w_cycle=[]
    tot_one_w_cycle=[]

    for i in range(0,len(cycle_set)):
        if w_cycle[i]==0:
            zero_w_cycle.append(cycle_set[i])
        else:
            one_w_cycle.append(cycle_set[i])

    powa_zero_w_cycle=powerset(zero_w_cycle)

    for cycle in one_w_cycle:

        for cycle_set_p in powa_zero_w_cycle:

            f_cycle=cycle
            for cycle2 in cycle_set_p:
                f_cycle=cycle_merger(f_cycle,cycle2)

            tot_one_w_cycle.append(f_cycle)


    return tot_one_w_cycle




def cutter(vset,eset,cycle_set,w_cycle):
    winding_cycles=winding_cycle_finder(vset,eset,cycle_set,w_cycle)

    #cycle_set_list=[list(cycle) for cycle in cycle_set]
    #print(cycle_set_list)
    prod=itertools.product(*winding_cycles)
    prod_set=[set(prods) for prods in prod]
    print("prod_set")
    print(prod_set)

    cuts=[]
    #print("-----")
    for st in prod_set:    
#        print(st)
        graph_del_e=eset.difference(st)
        counter1=0
        

#        print("deleted graph")
#        print(graph_del_e)
        sp_tree=spanning_tree_finder(vset,graph_del_e)
        st_v=set([e[0] for e in graph_del_e]).union(set([e[1] for e in graph_del_e]))
#        print("spanning_tree")
#        print(sp_tree)
#        print(st_v)
        if len(sp_tree[0])==len(st_v):
#            print("yes")
            counter1+=1

        else:
           print("no")

        counter2=0
        #in order to check minimality, add back the edges in the cut once at a time, and if adding back the edge does not restore at least one of the winding cycles, then it is not minimal
#        print("-------")
#        print("cut")
#        print(st)
        
        for e in st:
#            print("*")
#            print("winding_cycles")
#            print(winding_cycles)
#            print("edge")
#            print(e)
            gg=graph_del_e.difference(st).union(set([e]))
#            print(gg)
            counter3=0
            for cyc in winding_cycles:
#                print("cycle")
#                print(cyc)
                if len(cyc.difference(gg))==0:
#                    print("yes")
                    counter3+=1

            if counter3==0:
                counter2=1




        if counter1>0 and counter2==0:
            print("yes the cut is approved")
            cuts.append(st)


    return cuts

  











        
        


if __name__ == '__main__':
    g=set([1])
    g=set([(1,2)])
    
    vset={1,2,3,4,5,6}
    eset=set([(1,2),(2,3),(3,1),(4,3),(4,5),(5,6),(6,2)])

    vset={1,2,3,4}
    eset=set([(1,2),(2,3),(3,1),(4,3),(1,4),(2,4)])

    listy=[1,2,3]
    listy[1]
    adjacentset(vset,eset,2)
    print(list(adjacentset(vset,eset,2)))
    print(spanning_tree_finder(vset,eset))
    print(generators_finder(vset,eset))
    #print(spanning_tree_finder(vset,{(6, 2), (2, 3), (4, 3), (4, 5), (3, 1)}))
    print(winding_cycle_finder(vset,eset,generators_finder(vset,eset),[1,0,1]))
    lst = list(itertools.product([0, 1], repeat=3))
    print(lst)
    
    cuts=cutter(vset,eset,generators_finder(vset,eset),[1,0,1])

    print(cuts)
