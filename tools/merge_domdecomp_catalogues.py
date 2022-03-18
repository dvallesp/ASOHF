import os, sys, numpy as np
from cython_fortran_file import FortranFile as FF

######### DOMAIN DECOMPOSITION PARAMETERS ########
# Everything is given in input coordinate units
#sidelength_x=25000.
#sidelength_y=25000.
#sidelength_z=25000.
#domain_xleft=0.
#domain_yleft=0.
#domain_zleft=0.

## units 
h=0.6711
length_unit_in_cMpc=0.001/h
##

itini=0
itfin=0
every=50

div_x=2
div_y=2
div_z=2
sec_boundary=2500.

dd_foldername='domain_{:}_{:}_{:}'
domains_filename='domains.out'

are_stars=True
are_particles=True
##################################################

# Read domain file
with open(domains_filename, 'r') as f:
    domains=np.zeros((div_x,div_y,div_z,6))
    for l in f:
        i,j,k,x1,x2,y1,y2,z1,z2=l.split()
        i,j,k=(int(scr) for scr in (i,j,k))
        x1,x2,y1,y2,z1,z2=(float(scr) for scr in (x1,x2,y1,y2,z1,z2))
        domains[i,j,k,:]=np.array([x1,x2,y1,y2,z1,z2])

# Loop over the iterations
for it in range(itini, itfin+every, every):
    print('*****************************')
    print('* Iteration {:05d}           *'.format(it))
    print('*****************************')
    ##### FAMILIES FILE ##################
    nclus=0
    newcatalogue=[]

    old_to_new_id_global = {i: {j: {k: {} for k in range(div_z)} for j in range(div_y)} for i in range(div_x)}
    new_to_old_id_global = {}

    # Loop over the domains
    for i in range(div_x):
        for j in range(div_y):
            for k in range(div_z):
                x1,x2,y1,y2,z1,z2=domains[i,j,k,:]*length_unit_in_cMpc
            
                # Read the catalogue
                folder=dd_foldername.format(i,j,k)
                with open(os.path.join(folder,'output_files/families{:05d}'.format(it)), 'r') as f:
                    families_header=[]
                    for iline in range(7):
                        families_header.append(f.readline())
                        
                    splitlines=[]
                    ilist=0
                    id_to_list={}
                    list_to_id={}

                    for line in f:
                        halo=line.split()
                        splitlines.append(halo)
                            
                        idhalo,realclus,xc,yc,zc=int(halo[0]),int(halo[1]),float(halo[2]),float(halo[3]),float(halo[4])

                        id_to_list[idhalo]=ilist
                        list_to_id[ilist]={'id':idhalo, 'realclus':realclus, 'xc':xc, 'yc':yc, 'zc':zc}

                        ilist+=1
                
                # Keep only haloes which were in the 'good' part of each domain
                keep=[] # list index of the ones to keep
                caution_substructures=[]
                for ilist in range(len(splitlines)):
                    thisclus=list_to_id[ilist]
                    if thisclus['realclus']==-1:
                        xc,yc,zc=thisclus['xc'],thisclus['yc'],thisclus['zc']
                        if x1<=xc<x2 and y1<=yc<y2 and z1<=zc<z2:
                            keep.append(ilist)
                    else: # substructure
                        pare_ilist=ilist
                        while list_to_id[pare_ilist]['realclus']!=-1:
                            thatclus=list_to_id[pare_ilist]
                            pare=thatclus['realclus']
                            pare_ilist=id_to_list[pare]
                        #print(ilist, pare_ilist)
                        thatclus=list_to_id[pare_ilist] 
                        xcp,ycp,zcp=thatclus['xc'],thatclus['yc'],thatclus['zc']
                        if x1<=xcp<x2 and y1<=ycp<y2 and z1<=zcp<z2:
                            keep.append(ilist)
                            xc,yc,zc=thisclus['xc'],thisclus['yc'],thisclus['zc']
                            if not (x1<=xc<x2 and y1<=yc<y2 and z1<=zc<z2):
                                caution_substructures.append(ilist)

                # Build the new catalogue
                list_to_newid={}
                oldid_to_newid={}
                for ilist in keep:
                    nclus+=1
                    line=splitlines[ilist]
                    oldid=int(line[0])

                    list_to_newid[ilist]=nclus
                    oldid_to_newid[oldid]=nclus
                    line[0]=str(nclus)

                    realclus=int(line[1])
                    if realclus==-1:
                        line[1]=str(-1)
                    else:
                        line[1]=str(oldid_to_newid[realclus])
                        #print(line[1])
                    newcatalogue.append(line)
                
                # Save the relation between input and output catalogues
                for k_old, v_new in oldid_to_newid.items():
                    old_to_new_id_global[i][j][k][k_old] = v_new
                    new_to_old_id_global[v_new] = [i,j,k,k_old]

    ## check subs
    #subs=[]
    #for iline in range(len(newcatalogue)):
    #    if int(newcatalogue[iline][1]) > 0:
    #        subs.append(iline)
    #
    #for idx_sub in range(len(subs)):
    #    x1,y1,z1,r1=(float(newcatalogue[idx_sub][a]) for a in [2,3,4,8])
    #    for jdx_sub in range(idx_sub+1,len(subs)):
    #        x2,y2,z1,r2=(float(newcatalogue[jdx_sub][a]) for a in [2,3,4,8])
    #        if (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 < (r1+r2)**2:
    #            print('Warning! Duplicated substructures',idx_sub,jdx_sub)
    ## end check subs

    divisions_x=domains[1:,0,0,0]
    divisions_y=domains[0,1:,0,2]
    divisions_z=domains[0,0,1:,4]

    divisions_x1=(divisions_x-sec_boundary)*length_unit_in_cMpc
    divisions_x2=(divisions_x+sec_boundary)*length_unit_in_cMpc
    divisions_y1=(divisions_y-sec_boundary)*length_unit_in_cMpc
    divisions_y2=(divisions_y+sec_boundary)*length_unit_in_cMpc
    divisions_z1=(divisions_z-sec_boundary)*length_unit_in_cMpc
    divisions_z2=(divisions_z+sec_boundary)*length_unit_in_cMpc
    divisions_x=divisions_x.size
    divisions_y=divisions_y.size
    divisions_z=divisions_z.size

    in_the_intersection=[]
    #in_the_intersection_subs=[]
    for iline in range(len(newcatalogue)):
        x1,y1,z1,realclus1=*(float(newcatalogue[iline][ii]) for ii in [2,3,4]), int(newcatalogue[iline][1])
        
        if realclus1>0:
            continue
        
        flag=False
        if not flag:
            for ii in range(divisions_x):
                if divisions_x1[ii]<x1<divisions_x2[ii]:
                    flag=True
                    break
        if not flag:
            for ii in range(divisions_y):
                if divisions_y1[ii]<y1<divisions_y2[ii]:
                    flag=True
                    break
        if not flag:
            for ii in range(divisions_z):
                if divisions_z1[ii]<z1<divisions_z2[ii]:
                    flag=True
                    break
        if flag:
            in_the_intersection.append(iline)
    
    #print(len(in_the_intersection),len(in_the_intersection_subs))

    for isub in caution_substructures:
        print(isub)
        ilist=newcatalogue[isub]
        x1,y1,z1,r1,realclus1 = *(float(newcatalogue[isub][a]) for a in [2,3,4,8]),int(newcatalogue[isub][1])
          
        for jnotsub in in_the_intersection:
            x2,y2,z2,r2 = (float(newcatalogue[jnotsub][a]) for a in [2,3,4,6])
            dista2=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2

            if dista2<r1**2:
                print(isub,'--',jnotsub,dista2**0.5,r1,r2)

    with open('output_files/families{:05d}'.format(it), 'w') as f:
        scr = families_header[1].split()
        scr[1] = str(len(newcatalogue))
        scr[2] = str(len(newcatalogue))
        families_header[1] = ('{:>10}\t'*3+'{:>10}\n').format(*scr)
        for l in families_header:
            f.write(l)
        for l in newcatalogue:
            f.write(('{:>14}'*50+'{:>14}\n').format(*l))

    ##### STELLAR_HALOES FILE ##################
    if are_stars:
        nclusst=0
        newcataloguest=[]

        old_to_new_id_global_st = {i: {j: {k: {} for k in range(div_z)} for j in range(div_y)} for i in range(div_x)}
        new_to_old_id_global_st = {}

        # Loop over the domains
        for i in range(div_x):
            for j in range(div_y):
                for k in range(div_z):
                    x1,x2,y1,y2,z1,z2=domains[i,j,k,:]*length_unit_in_cMpc
            
                    # Read the catalogue
                    folder=dd_foldername.format(i,j,k)
                    with open(os.path.join(folder,'output_files/stellar_haloes{:05d}'.format(it))) as f:
                        stars_header=[]
                        for iline in range(7):
                            stars_header.append(f.readline())

                        the_keys = old_to_new_id_global[i][j][k].keys()
                        for line in f:
                            halo = line.split()
                            stid = int(halo[0])
                            dmid = int(halo[1])
                            if dmid in the_keys:
                                newdmid = old_to_new_id_global[i][j][k][dmid]
                                nclusst += 1
                                halo[0] = int(nclusst)
                                halo[1] = int(newdmid)
                                
                                newcataloguest.append(halo)

                                old_to_new_id_global_st[i][j][k][stid] = nclusst
                                new_to_old_id_global_st[nclusst] = [i,j,k,stid]

        with open('output_files/stellar_haloes{:05d}'.format(it), 'w') as f:
            scr = stars_header[1].split()
            scr[1] = str(len(newcataloguest))
            scr[2] = str(len(newcataloguest))
            stars_header[1] = ('{:>10}\t'*3+'{:>10}\n').format(*scr)
            for l in stars_header:
                f.write(l)
            for l in newcataloguest:
                f.write(('{:>14}'*29+'{:>14}\n').format(*l))
    
    if are_particles:
        particles_old_lut = {i: {j: {k: {} for k in range(div_z)} for j in range(div_y)} for i in range(div_x)}
        particles_lists = {i: {j: {k: None for k in range(div_z)} for j in range(div_y)} for i in range(div_x)}
        if are_stars:
            particles_old_lut_st = {i: {j: {k: {} for k in range(div_z)} for j in range(div_y)} for i in range(div_x)}
            particles_lists_st = {i: {j: {k: None for k in range(div_z)} for j in range(div_y)} for i in range(div_x)}

        for i in range(div_x):
            for j in range(div_y):
                for k in range(div_z): 
                    folder=dd_foldername.format(i,j,k)
                    with FF(os.path.join(folder, 'output_files/particles{:05d}'.format(it))) as f:
                        hm = f.read_vector('i4')[0]
                        for ii in range(hm):
                            nclusold,i1,i2=f.read_vector('i4')
                            i1-=1
                            i2-=1
                            particles_old_lut[i][j][k][nclusold]=[i1,i2]
                        f.read_vector('i4')
                        particles_lists[i][j][k]=f.read_vector('i4')

                    if are_stars:
                        with FF(os.path.join(folder, 'output_files/stellar_particles{:05d}'.format(it))) as f:
                            hm = f.read_vector('i4')[0]
                            for ii in range(hm):
                                nclusold,i1,i2=f.read_vector('i4')
                                i1-=1
                                i2-=1
                                particles_old_lut_st[i][j][k][nclusold]=[i1,i2]
                            f.read_vector('i4')
                            particles_lists_st[i][j][k]=f.read_vector('i4')
        particles=[]
        new_lut=[]
        inew1=0
        inew2=0
        for halo in newcatalogue:
            nclus=int(halo[0])
            i,j,k,nclusold = new_to_old_id_global[nclus]
            i1,i2 = particles_old_lut[i][j][k][nclusold]
            particles.append(particles_lists[i][j][k][i1:i2+1])
            numpart = i2-i1+1
            inew1 = inew2+1
            inew2 = inew1+numpart-1
            new_lut.append([nclus,inew1,inew2])

        with FF('output_files/particles{:05d}'.format(it), 'w') as f:
            f.write_vector(np.array([len(newcatalogue)], dtype='i4'))
            for lut_entry in new_lut:
                f.write_vector(np.array(lut_entry, dtype='i4'))
            particles=np.concatenate(particles)
            print('Particles check:', inew2, particles.size)
            f.write_vector(np.array([particles.size], dtype='i4'))
            f.write_vector(particles.astype('i4'))
    
        if are_stars: 
            particles=[]
            new_lut=[]
            inew1=0
            inew2=0
            for halo in newcataloguest:
                nclus=int(halo[0])
                i,j,k,nclusold = new_to_old_id_global_st[nclus]
                i1,i2 = particles_old_lut_st[i][j][k][nclusold]
                particles.append(particles_lists_st[i][j][k][i1:i2+1])
                numpart = i2-i1+1
                inew1 = inew2+1
                inew2 = inew1+numpart-1
                new_lut.append([nclus,inew1,inew2])

            with FF('output_files/particles_stellar{:05d}'.format(it), 'w') as f:
                f.write_vector(np.array([len(newcataloguest)], dtype='i4'))
                for lut_entry in new_lut:
                    f.write_vector(np.array(lut_entry, dtype='i4'))
                particles=np.concatenate(particles)
                print('Particles check:', inew2, particles.size)
                f.write_vector(np.array([particles.size], dtype='i4'))
                f.write_vector(particles.astype('i4'))

    print('Merged!')
    nclus=len(newcatalogue)
    print('{:} haloes'.format(nclus))
    if are_stars:
        nclusst=len(newcataloguest)
        print('{:} stellar haloes'.format(nclusst))
    nsubs=len([1 for a in newcatalogue if int(a[1])>0])
    print('{:} substructures'.format(nsubs))
    #newcatalogue = sorted(newcatalogue, key=lambda a: (int(a[1])==-1)*float(a[5]) + (int(a[1])>0)*float(a[7])*1.e-30, reverse=True)
