import numpy as np

# Cosmology
h = 0.6766
omega_m = 0.3111
omega_b = 0.0490

# Box size
size = 100. / h
nx = 256
zini = 100.

#### DEFINITIONS
class Units:
    pass
units = Units()
units.mpc_to_m = 3.085678e22 # m/Mpc
units.m_to_mpc = 1/units.mpc_to_m
units.G_isu = 6.67430e-11
units.sun_to_kg = 1.98848e30
units.kg_to_sun = 1/units.sun_to_kg

def critical_density(h, z=0, omega_m=0, omega_lambda=0):
    """
    Computes the value of the critical density of the universe (for making k=0, Lambda=0) at a given redshift, z.
    Units: solar masses per Mpc^3.

    Args:
        h: dimensionless Hubble parameter, H_0 = 100 h km/s/Mpc
        z: redshift (don't specify if values at the present time required)
        omega_m: matter density parameter, at z=0
        omega_lambda: dark energy density parameter, at z=0

    Returns:
        Critical density, in Msun / Mpc^3
    """
    isu_0 = 3 * (100000 * h / units.mpc_to_m) ** 2 / (8*np.pi*units.G_isu)
    if z == 0:
        isu = isu_0
    else:
        isu = isu_0 * E(z, omega_m, omega_lambda) ** 2

    return isu * units.kg_to_sun / units.m_to_mpc**3

def background_density(h, omega_m, z=0):
    """
    Computes the value of the background (matter) density of the universe, at a given redshift z.
    Units: solar masses per Mpc^3.

    Args:
        h: dimensionless Hubble parameter, H_0 = 100 h km/s/Mpc
        z: redshift (don't specify if values at the present time required)
        omega_m: matter density parameter, at z=0

    Returns:
        Bkg density, in Msun / Mpc^3
    """
    # rhoB = omega_m * rhoC0 / (1+z)^3
    return omega_m * critical_density(h) * (1+z)**3
####

## MASS UNITS (to Msun)
mass_box = background_density(h, omega_m) * size**3
fb = (omega_m - omega_b) / omega_m
mass_particle_combined = mass_box / nx**3
# masses in the file are (1-fb) * (Dx_inigrid / Dx_grid)^3
# for now I am doing unigrids so they are just 1-fb; 
# so that the actual mass in Msun will be that * mass_particle_combined
# And therefore the unit mass to specify will be mass_particle_combined
um = mass_particle_combined

## LENGTH UNITS (to Mpc)
ul = size # i put it positive and therefore don't have to correct it by h (it already is when defining it)

## Velocity untis (to km/s)
length_unit_isu = ul * units.mpc_to_m / (1 + zini)
rhoB0 = background_density(h, omega_m) * units.sun_to_kg / units.mpc_to_m**3
time_unit_isu = 1 / np.sqrt(4 * np.pi * units.G_isu * rhoB0 * (1 + zini)**3)
velocity_unit_isu = length_unit_isu / time_unit_isu
velocity_unit_kms = velocity_unit_isu / 1000.

print('Hubble constant (h):', h)
print('Omega matter:', omega_m)
print('Fraction of dark matter', (omega_m-omega_b)/omega_m)


print('Mass unit (Msun):', um)
print('Length unit (Mpc):', ul)
print('Velocity unit (km/s):', velocity_unit_kms)
