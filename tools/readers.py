import os, sys, numpy as np
from cython_fortran_file import FortranFile as FF

def read_families(it, path='', output_format='dictionaries', output_redshift=False,
                  min_mass=None, exclude_subhaloes=False):
    '''
    Reads the families files, containing the halo catalogue.
    Can be outputted as a list of dictionaries, one per halo
     (output_format='dictionaries') or as a dictionary of 
     arrays (output_format='arrays').
    '''
    with open(os.path.join(path, 'families{:05d}'.format(it)), 'r') as f:
        _=f.readline()
        _,_,_,zeta = f.readline().split()
        zeta=float(zeta)
        for i in range(5):
            _=f.readline()
        haloes=[]
        for l in f:
            l=l.split()
            halo={}
            halo['id']=int(l[0])
            halo['substructureOf']=int(l[1])
            halo['x']=float(l[2])
            halo['y']=float(l[3])
            halo['z']=float(l[4])
            halo['Mvir']=float(l[5])
            halo['Rvir']=float(l[6])
            if halo['substructureOf']==-1:
                halo['M']=halo['Mvir']
                halo['R']=halo['Rvir']
            else:
                halo['M']=float(l[7])
                halo['R']=float(l[8])
            halo['partNum']=int(l[9])
            halo['mostBoundPart']=int(l[10])
            halo['xcm']=float(l[11])
            halo['ycm']=float(l[12])
            halo['zcm']=float(l[13])
            halo['majorSemiaxis']=float(l[14])
            halo['intermediateSemiaxis']=float(l[15])
            halo['minorSemiaxis']=float(l[16])
            halo['Ixx']=float(l[17])
            halo['Ixy']=float(l[18])
            halo['Ixz']=float(l[19])
            halo['Iyy']=float(l[20])
            halo['Iyz']=float(l[21])
            halo['Izz']=float(l[22])
            halo['Lx']=float(l[23])
            halo['Ly']=float(l[24])
            halo['Lz']=float(l[25])
            halo['sigma_v']=float(l[26])
            halo['vx']=float(l[27])
            halo['vy']=float(l[28])
            halo['vz']=float(l[29]) 
            halo['max_v']=float(l[30])
            halo['mean_vr']=float(l[31])
            halo['Ekin']=float(l[32])
            halo['Epot']=float(l[33])
            halo['Vcmax']=float(l[34])
            halo['Mcmax']=float(l[35])
            halo['Rcmax']=float(l[36])
            halo['R200m']=float(l[37])
            halo['M200m']=float(l[38])
            halo['R200c']=float(l[39])
            halo['M200c']=float(l[40])
            halo['R500m']=float(l[41])
            halo['M500m']=float(l[42])
            halo['R500c']=float(l[43])
            halo['M500c']=float(l[44])
            halo['R2500m']=float(l[45])
            halo['M2500m']=float(l[46])
            halo['R2500c']=float(l[47])
            halo['M2500c']=float(l[48])
            halo['fsub']=float(l[49])
            halo['Nsubs']=int(l[50])

            haloes.append(halo)
    
    if exclude_subhaloes:
        haloes=[halo for halo in haloes if halo['substructureOf']==-1]
    if min_mass is not None:
        haloes=[halo for halo in haloes if halo['M']>min_mass]

    if output_format=='dictionaries':
        if output_redshift:
            return haloes, zeta
        else:
            return haloes
    elif output_format=='arrays':
        if output_redshift:
            return {k: np.array([h[k] for h in haloes]) for k in haloes[0].keys()}, zeta
        else:
            return {k: np.array([h[k] for h in haloes]) for k in haloes[0].keys()}


def read_stellar_haloes(it, path='', output_format='dictionaries'):
    '''
    Reads the stellar haloes files, containing the halo catalogue.
    Can be outputted as a list of dictionaries, one per halo
     (output_format='dictionaries') or as a dictionary of 
     arrays (output_format='arrays').
    '''
    with open(os.path.join(path, 'stellar_haloes{:05d}'.format(it)), 'r') as f:
        for i in range(7):
            f.readline()
        haloes=[]
        for l in f:
            l=l.split()
            halo={}
            halo['id']=int(l[0])
            halo['DMid']=int(l[1])
            halo['xDM']=float(l[2])
            halo['yDM']=float(l[3])
            halo['zDM']=float(l[4])
            halo['x']=float(l[5])
            halo['y']=float(l[6])
            halo['z']=float(l[7]) 
            halo['Mhalf']=float(l[8])
            halo['Rhalf']=float(l[9])
            halo['numpart']=int(l[10])
            halo['xcm']=float(l[11])
            halo['ycm']=float(l[12])
            halo['zcm']=float(l[13])
            halo['majorSemiaxis']=float(l[14])
            halo['intermediateSemiaxis']=float(l[15])
            halo['minorSemiaxis']=float(l[16])
            halo['Ixx']=float(l[17])
            halo['Ixy']=float(l[18])
            halo['Ixz']=float(l[19])
            halo['Iyy']=float(l[20])
            halo['Iyz']=float(l[21])
            halo['Izz']=float(l[22])
            halo['Lx']=float(l[23])
            halo['Ly']=float(l[24])
            halo['Lz']=float(l[25])
            halo['sigma_v']=float(l[26])
            halo['vx']=float(l[27])
            halo['vy']=float(l[28])
            halo['vz']=float(l[29]) 

            haloes.append(halo)
    
    if output_format=='dictionaries':
        return haloes
    elif output_format=='arrays':
        return {k: np.array([h[k] for h in haloes]) for k in haloes[0].keys()}


def read_particles(it, path='', parttype='DM', sort='oripa'):
    '''
    Reads the particles list. Set parttype='DM' for DM particles,
    'stellar' for stellar particles.
    - sort: 'oripa' (sort particles by increasing ID) or 'r' (sort by 
       increasing radius).
    '''
    if parttype=='DM':
        filename='particles'
    elif parttype=='stellar':
        filename='stellar_particles'

    particles_oripa={}
    particles_lut={}
    with FF(os.path.join(path, filename+'{:05d}'.format(it)), 'r') as f:
        nhaloes=f.read_vector('i4')[0]
        for i in range(nhaloes):
            nclus,i1,i2=f.read_vector('i4')
            particles_lut[nclus]=[i1-1,i2-1]
        f.read_vector('i4')
        particles=f.read_vector('i4')
    
    if sort=='oripa':
        for k,(i1,i2) in particles_lut.items():
            particles_oripa[k]=np.sort(particles[i1:i2+1])
    elif sort=='r':
        for k,(i1,i2) in particles_lut.items():
            particles_oripa[k]=particles[i1:i2+1]
    else:
        print('Error! sort should be either oripa or r')

    return particles_oripa

