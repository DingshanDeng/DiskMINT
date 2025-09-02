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

f = np.loadtxt('../../COinitgrid-' + save_name + '-cyl.dat',skiprows=4,delimiter=',')
rdisk  = f[:,0]
zdisk  = f[:,2]
ndisk  = f[:,3]
abh2 = np.zeros(nr*nz)
abh2[:] = -99
abc18o = np.zeros(nr*nz)
abc18o[:] = -99

f = np.loadtxt('out/ab/H2.txt',skiprows=1)
abh2tmp = f[:,2]
f = np.loadtxt('out/ab/CO@.txt',skiprows=1)
abc18otmp = f[:,2]

abh2[ndisk>1.0]   = abh2tmp
abc18o[ndisk>1.0] = abc18otmp
        
with open('COendgrid-withCO2-' + save_name + '-cyl.chem','w+') as f:
    k=0
    for i in range(nr):
        for j in range(nz):
            f.write('%13.5e %13.5e %13.5e %13.5e\n'%(rdisk[k], zdisk[k], abh2[k], abc18o[k]))
            k+=1
