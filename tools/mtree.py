import os, sys, numpy as np, sortednp as snp, json, datetime
from scipy import integrate
from readers import *
from multiprocessing import Pool
import gc

from tqdm import tqdm
## if tqdm not installed, comment the above line and uncomment the following
#def tqdm(x): return x

######### Parameters
outputs_ASOHF='output_files/'
simulation_results='simu_masclet/'
itini=100
itfin=1550
every=50

h=0.6711
omegam=0.3026
omegalambda=1-omegam

ncores_distance=20
ncores_intersect=4

min_given_mass = 0.001 # fraction of mass of the progenitor given to the descendant to report this progenitor

look_further_iterations = False # if True, looks at previous iterations for the, until it finds a matching halo which contains the mostBoundParticle
max_iterations_back = 2 # set it to a very large number to go always until the first iteration
#########

### function to read output a dictionary mapping particle IDs to particles masses
### implement it to suit your simulation code.
def IDs_to_masses(it, folder_simu):
    from masclet_framework.read_masclet import read_cldm
    from masclet_framework.units import mass_to_sun as munit
    mdm,oripa = read_cldm(it, path=folder_simu, output_deltadm=False, output_position=False, output_velocity=False,
                          output_mass=True, output_id=True)
    # fix IDs from MASCLET
    oripa[mdm>mdm.max()/2]=-np.abs(oripa[mdm>mdm.max()/2])
    return {partid: mi for mi,partid in zip(mdm*munit,oripa)}
###
###

def LCDM_time(z1,z2,omega_m,omega_lambda,h):
    """
    Computes the time between z1 and z2, in years, using a quadrature method
    
    TAKEN FROM masclet_framework: https://github.com/dvallesp/masclet_framework

    Args:
        z1: initial redshift
        z2: final redshift
        omega_m: matter density parameter, at z=0
        omega_lambda: dark energy density parameter, at z=0
        h: dimensionless Hubble constant

    Returns:
        Time, in years, from z1 to z2
    """
    
    def E(z,omega_m, omega_lambda, omega_r=0):
        """
        Computes the evolution function of the Hubble parameter (H(z) = H_0 E(z)).

        TAKEN FROM masclet_framework: https://github.com/dvallesp/masclet_framework
    
        Args:
            z: redshift where to evaluate the function
            omega_m: matter density parameter, at z=0
            omega_lambda: dark energy density parameter, at z=0
            omega_r: radiation density parameter, at z=0. Don't specify if negligible

        Returns:
            Value of the E(z) function for this cosmology
        """
        omega_k = 1 - omega_m - omega_lambda - omega_r
        a = 1/(1+z)
        E = np.sqrt(omega_lambda + omega_k/a**2 + omega_m/a**3 + omega_r/a**4)
        return E

    if z1<z2:
        raise ValueError('Initial redshift must be larger than final redshift! Negative time will be got.')

    tH = 9784597488

    # nested function for the integrand
    def integrand(z, omega_m, omega_lambda):
        return 1/(E(z, omega_m, omega_lambda)*(1+z))

    t = integrate.quad(integrand, z2, z1, (omega_m, omega_lambda))

    if t[1] > 0.001*t[0]:
        raise ValueError("Error greater than 1 per mil")

    return tH * t[0] / h

def p_consider(j_post):
    xj,yj,zj=(haloes_post[aaa][j_post] for aaa in ['x','y','z'])
    consider_j=[]
    for i_prev in range(n_prev):
        xi,yi,zi,ri=(haloes_prev[aaa][i_prev] for aaa in ['x','y','z','R'])
        dista2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
        if dista2 < (max_dista+ri)**2:
            consider_j.append(i_prev)
    return consider_j

def p_intersect(j_post):
    id_j_post = int(haloes_post['id'][j_post])
    merger_tree_j_post = []

    oripas_post = particles_haloes_post[id_j_post]
    mpost = sum([particledict[oripa_i] for oripa_i in oripas_post])
    mostBoundPost = haloes_post['mostBoundPart'][j_post]

    for i_prev in consider[j_post]:
        id_i_prev = int(haloes_prev['id'][i_prev])

        oripas_prev = particles_haloes_prev[id_i_prev]
        mprev = sum([particledict[oripa_i] for oripa_i in oripas_prev])
        intersection = sum([particledict[oripa_i] for oripa_i in snp.intersect(oripas_post, oripas_prev)])
        mostBoundPrev = haloes_prev['mostBoundPart'][i_prev]
        if intersection > min_given_mass*mpost:
            #print(j_post, 'from', i_prev, intersection)
            merger_tree_j_post.append({'id': id_i_prev,
                                       'givenMassFraction': intersection/mpost,
                                       'retainedMassFraction': intersection/mprev,
                                       'normalisedIntersection': intersection/np.sqrt(mprev*mpost),
                                       'containsMostBound': (mostBoundPost in oripas_prev) and (mostBoundPrev in oripas_post),
                                       'position': [haloes_prev[aaa][i_prev] for aaa in ['x','y','z']],
                                       'radius': haloes_prev['R'][i_prev], 'mass': haloes_prev['M'][i_prev]})
        del oripas_prev

    del oripas_post 
    return merger_tree_j_post

particledict = IDs_to_masses(itini, simulation_results)

print('Building merger tree... From iter {:} to iter {:}'.format(itini, itfin))
for it_post in range(itini+every, itfin+every, every):
    it_prev = it_post-every
    print('******************************')
    print('NEW STEP: {:} --> {:}'.format(it_prev, it_post))
    print(datetime.datetime.now())
    print('******************************')

    haloes_prev, zeta_prev = read_families(it_prev, path=outputs_ASOHF, 
                                           output_format='arrays', output_redshift=True)
    particles_haloes_prev = read_particles(it_prev, path=outputs_ASOHF)
    print('Iteration {:}, redshift {:.2f}: read! {:} haloes, {:} substructure'.format(it_prev, zeta_prev,
                                                                                      haloes_prev['substructureOf'].size,
                                                                                      (haloes_prev['substructureOf']>0).sum()))

    haloes_post, zeta_post = read_families(it_post, path=outputs_ASOHF, 
                                           output_format='arrays', output_redshift=True)
    particles_haloes_post = read_particles(it_post, path=outputs_ASOHF)
    print('Iteration {:}, redshift {:.2f}: read! {:} haloes, {:} substructure'.format(it_post, zeta_post,
                                                                                      haloes_post['substructureOf'].size,
                                                                                      (haloes_post['substructureOf']>0).sum()))

    n_post = haloes_post['id'].size
    n_prev = haloes_prev['id'].size
    
    maxv = max([haloes_prev['max_v'].max()*(1+zeta_prev), haloes_post['max_v'].max()*(1+zeta_post)])
    maxv2 = max([max([haloes_prev['vx'].max(), haloes_prev['vy'].max(), haloes_prev['vz'].max()])*(1+zeta_prev),
                 max([haloes_post['vx'].max(), haloes_post['vy'].max(), haloes_post['vz'].max()])*(1+zeta_post)])
    maxv = max([maxv, maxv2])/299792.458
    Dt = LCDM_time(zeta_prev, zeta_post, omegam, omegalambda, h) 
    max_dista = Dt * (maxv * 3.06392e-7)
    print('Time span, max velocity, max distance to consider: {:.3f} Gyr, {:.3e} c, {:.2f} Mpc'.format(Dt/1e9, maxv, max_dista))
    print('Finding candidates...')

    with Pool(ncores_distance) as p:
        consider=list(tqdm(p.imap(p_consider, range(n_post)), total=n_post))
    consider = {i: ci for i,ci in enumerate(consider)}

    print('Possible relations to check: ', sum([len(v) for v in consider.values()]))

    with Pool(ncores_intersect) as p:
        merger_tree=list(tqdm(p.imap(p_intersect, range(n_post)), total=n_post))
    merger_tree = {int(haloes_post['id'][j_post]): v for j_post, v in enumerate(merger_tree)}
    gc.collect()

    with open(os.path.join(outputs_ASOHF, 'mtree_{:05d}_{:05d}.json'.format(it_prev, it_post)), 'w') as f:
        json.dump(merger_tree, f, indent=4)

    # find 'Lost haloes'
    if look_further_iterations:
        lost = []
        for j_post in range(n_post):
            found=False
            for prog in merger_tree[int(haloes_post['id'][j_post])]:
                if prog['containsMostBound'] is True:
                    found=True
                    break
            if not found:
                lost.append(j_post)

        print('** Number of lost haloes:', len(lost))

        for it_prev in range(it_post-2*every, max([itini-every, it_post-(2+max_iterations_back)*every]), -every):
            haloes_prev, zeta_prev = read_families(it_prev, path=outputs_ASOHF,
                                                   output_format='arrays', output_redshift=True)
            particles_haloes_prev = read_particles(it_prev, path=outputs_ASOHF)
            print('\n==> Iteration {:}, redshift {:.2f}: read! {:} haloes, {:} substructure <=='.format(it_prev, zeta_prev,
                                                                                              haloes_prev['substructureOf'].size,
                                                                                              (haloes_prev['substructureOf'] > 0).sum()))
            n_prev = haloes_prev['id'].size

            maxv = max([haloes_prev['max_v'].max() * (1 + zeta_prev), maxv*299792.458])
            maxv2 = max([haloes_prev['vx'].max(), haloes_prev['vy'].max(), haloes_prev['vz'].max()]) * (1 + zeta_prev)
            maxv = max([maxv, maxv2]) / 299792.458
            Dt = LCDM_time(zeta_prev, zeta_post, omegam, omegalambda, h)
            max_dista = Dt * (maxv * 3.06392e-7)
            print('Time span, max velocity, max distance to consider: {:.3f} Gyr, {:.3e} c, {:.2f} Mpc'.format(Dt / 1e9, maxv, max_dista))
            print('Finding candidates...')

            with Pool(ncores_distance) as p:
                consider = list(tqdm(p.imap(p_consider, lost), total=len(lost)))
            consider = {i: ci for i, ci in zip(lost, consider)}

            print('Possible relations to check: ', sum([len(v) for v in consider.values()]))

            with Pool(ncores_intersect) as p:
                merger_tree = list(tqdm(p.imap(p_intersect, consider), total=len(consider)))
            merger_tree = {int(haloes_post['id'][j_post]): v for j_post, v in zip(lost, merger_tree)}

            lost0 = [i for i in lost]
            lost = []
            for j_post in lost0:
                found = False
                for prog in merger_tree[int(haloes_post['id'][j_post])]:
                    if prog['containsMostBound'] is True:
                        found = True
                        break
                if not found:
                    lost.append(j_post)

            print('Number of lost haloes: {:}. Recovered: {:}'.format(len(lost), len(lost0)-len(lost)))

            with open(os.path.join(outputs_ASOHF, 'mtree_{:05d}_{:05d}.json'.format(it_prev, it_post)), 'w') as f:
                json.dump(merger_tree, f, indent=4)

            gc.collect()
