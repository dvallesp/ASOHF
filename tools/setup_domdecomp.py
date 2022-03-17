import os, sys, shutil

######### DOMAIN DECOMPOSITION PARAMETERS ########
# Everything is given in input coordinate units
sidelength_x=25000.
sidelength_y=25000.
sidelength_z=25000.
domain_xleft=0.
domain_yleft=0.
domain_zleft=0.

div_x=2
div_y=2
div_z=2
sec_boundary=2500.

dd_foldername='domain_{:}_{:}_{:}'
executable='asohf_sa.x'
other_files=['run.sh']
simulation_foldername='simulation'
##################################################

domains_spec=[]

for i in range(div_x):
    for j in range(div_y):
        for k in range(div_z):
            x1 = domain_xleft + i*sidelength_x/div_x
            y1 = domain_yleft + j*sidelength_y/div_y
            z1 = domain_zleft + k*sidelength_z/div_z
            x2 = domain_xleft + (i+1)*sidelength_x/div_x
            y2 = domain_yleft + (j+1)*sidelength_y/div_y
            z2 = domain_zleft + (k+1)*sidelength_z/div_z

            domains_spec.append(('{:}\t'*3+'{:.4f}\t'*5+'{:.4f}\n').format(i,j,k,x1,x2,y1,y2,z1,z2))

            if i!=0:
                x1=x1-sec_boundary
            if j!=0:
                y1=y1-sec_boundary
            if k!=0:
                z1=z1-sec_boundary
            if i!=div_x-1:
                x2=x2+sec_boundary
            if j!=div_y-1:
                y2=y2+sec_boundary
            if k!=div_z-1:
                z2=z2+sec_boundary

            print(i,j,k,x1,x2,y1,y2,z1,z2)
            
            folder=dd_foldername.format(i,j,k)
            if not os.path.exists(folder):
                os.mkdir(folder)
            shutil.copy(executable, folder)
            for otherfile in other_files:
                shutil.copy(otherfile, folder)

            folder2=os.path.join(folder, 'output_files')
            if not os.path.exists(folder2):
                os.mkdir(folder2)
                
            folder2=os.path.join(folder, 'input_files')
            if not os.path.exists(folder2):
                os.mkdir(folder2)
            with open('input_files/asohf.dat', 'r') as f:
                lines=f.readlines()
            for iline,line in enumerate(lines):
                if 'Domain to keep particles' in line:
                    break
            lines[iline+1]='{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n'.format(x1,x2,y1,y2,z1,z2)
            with open(os.path.join(folder2, 'asohf.dat'), 'w') as f:
                f.writelines(lines)
            #shutil.copy('input_files/asohf.dat', folder2)

            folder2=os.path.join(folder, simulation_foldername)
            if os.path.exists(folder2):
                os.unlink(os.path.join(folder2))
            
            os.symlink(os.path.join('..',simulation_foldername), folder2)

            
with open('domains.out', 'w') as f:
    for line in domains_spec:
        f.write(line)
