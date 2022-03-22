##################################################
### EXAMPLE SCRIPT
#### Going up the main branch of the merger tree
##################################################

import json, numpy as np

itfin=1550 # last iteration
itini=500 # first iteration
every=50 # spacing of the snapshots

halo=50 # ID of the halo where to start to navigate, in iteration itfin

for itpost in range(itfin,itini-every,-every):
    itprev=itpost-every
    with open('mtree_{:05d}_{:05d}.json'.format(itprev,itpost), 'r') as f:
        mtree=json.load(f)
    print(itpost, halo)
    mtree=mtree[str(halo)]
    flag=False
    for prog in mtree:
        if prog['containsMostBound'] is True:
            halo=prog['id']
            flag=True
            break
    if not flag:
        print('Broken!')
        break

