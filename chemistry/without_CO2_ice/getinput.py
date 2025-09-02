import sys
sys.path.append('../../')
# sys.path.append('../../wrapper_scripts/')
import numpy as np
from wrapper_parameters import *
from wrapper_scripts.constants import *
import matplotlib.pyplot as plt

save_name = chemical_save_name
print(save_name)

# NEEDS TO BE IN THE CSV FILE
# nr = 100
# nz = 200
with open('../../COinitgrid-' + save_name + '-cyl.dat', 'r') as f:
    lines = f.readlines()
    nr = int(lines[1][3:-1])
    nz = int(lines[2][3:-1])

###############################################################
# Compute few quantities needed by the chemical model
###############################################################
#
# 1. Compute 0th, 1st and 2nd moment of the dust distribution
#
def dust_fun(nn):
    return (amax_all**(nn-pla_dustsize) - amin_all**(nn-pla_dustsize))/(nn-pla_dustsize)
eta =  1.4*mp/rhobulk/ratio_g2d_total
nd   = 3.0*eta*dust_fun(1)/4.0/np.pi/dust_fun(4)
nda  = 3.0*eta*dust_fun(2)/4.0/np.pi/dust_fun(4)
nda2 = 3.0*eta*dust_fun(3)/4.0/np.pi/dust_fun(4)
#
# 2. Compute G0 from spectrum (in Habing units)
#
f   = np.loadtxt('../../'+fmodel_filename)
wl  = f[:,0]
Fnu = f[:,1]*(3.08572e18/rstar)**2
nu  = clum/wl/1.0e-4
wll = hh*clum/13.6/eV/1.0e-4
wlu = hh*clum/6.0/eV/1.0e-4
Fnu[wl<wll] = 0.0
Fnu[wl>wlu] = 0.0
G0Hab = np.trapz(Fnu[::-1],nu[::-1])/1.6e-3
###############################################################


###############################################################
# Read structure computed with RADMC-3D and write it in a     #
# format compatible with the chemical model                   #
###############################################################

f = np.loadtxt('../../COinitgrid-' + save_name + '-cyl.dat',skiprows=4,delimiter=',')
rdisk  = np.reshape(f[:,0],(nr,nz))
zrdisk = np.reshape(f[:,1],(nr,nz))
zdisk  = np.reshape(f[:,2],(nr,nz))
ndisk  = np.reshape(f[:,3],(nr,nz))
Tgdisk = np.reshape(f[:,4],(nr,nz))
Avdisk = np.reshape(1.086*f[:,5],(nr,nz))
Avup   = np.reshape(1.086*f[:,6],(nr,nz))
Avdn   = np.reshape(1.086*f[:,7],(nr,nz))
Tddisk = Tgdisk.copy() # the area weighted Tdust = Tgas.

#
nrtmp = 0
nztmp = np.zeros(nr,dtype=int)
for i in range(0,nr):
    for j in range(nz-1,-1,-1):
        if(ndisk[i,j]>1.0):
            nztmp[i]+=1
    if (nztmp[i]>0): nrtmp += 1
# Write parameters files
with open('data/struc/diskparameters','w+') as f:
    f.write('%s \n'%rstar)
    f.write('%s \n'%G0Hab)
    f.write('%s \n'%nd)
    f.write('%s \n'%nda)
    f.write('%s \n'%nda2)
with open('data/struc/diskstructure','w+') as f:
    f.write('%s \n'%nrtmp)
    for i in range(0,nr):
        f.write('%s \n'%nztmp[i])
        for j in range(0,nz):
            if(ndisk[i,j]>1.0):
                f.write('%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e \n'%(np.log10(rdisk[i,j]/au),
                np.log10(zdisk[i,j]/au),np.log10(ndisk[i,j]),np.log10(Tgdisk[i,j]),np.log10(Tddisk[i,j]),
                Avdisk[i,j],Avup[i,j],Avdn[i,j]))
with open('data/struc/uvfield','w+') as f:
    for i in range(0,nr):
        for j in range(0,nz):
            if(ndisk[i,j]>1.0):
                f.write('%13.6e %13.6e %13.6e \n'%(Avdisk[i,j],Avup[i,j],Avdn[i,j]))