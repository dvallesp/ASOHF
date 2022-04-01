##################################################
### EXAMPLE SCRIPT
#### Going up the main branch of the merger tree
##################################################

import json, numpy as np, sys

arg1 = int(sys.argv[1])
arg2 = int(sys.argv[2])

######### Parameters
itfin=arg1 # last iteration
itini=100 # first iteration
every=50 # spacing of the snapshots

max_iterations_back=4 # same as in mtree.f

halo=arg2 # ID of the halo where to start to navigate, in iteration itfin
#########

itpost=itfin
it_back=0
print(itpost, halo)
while itpost > itini:
    itprev=itpost-(it_back+1)*every
    if itprev < itini:
        break

    with open('mtree_{:05d}_{:05d}.json'.format(itprev,itpost), 'r') as f:
        mtree=json.load(f)
    if str(halo) in mtree.keys():
        mtree=mtree[str(halo)]
    else:
        mtree=[]
    flag=False
    for prog in mtree:
        if prog['containsMostBound'] is True:
            halotry=prog['id']
            flag=True
            break
    if flag:
        it_back=0
        halo=halotry
        itpost=itprev
        print(itprev, halo, 'M={:.3e} Msun, R={:.3f} cMpc'.format(prog['mass'], prog['radius']))
    if not flag:
        print(itprev, '--')
        it_back+=1
        if it_back > max_iterations_back:
            break


