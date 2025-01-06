# mesh summary
'''
LES and RANS
'''

import numpy as np

def compute_mesh(nxOuter, nQuarter, nxOuter2, nxestension, nyOuter):
    x_mesh = nxOuter + 2*nQuarter + nxOuter2 + nxestension
    y_mesh = 2*nyOuter + 2*nQuarter
    return x_mesh, y_mesh


'''
RANS
'''


# RANS - sim
nQuarter    =  14
nxOuter     =  20
nxOuter2    = 30
nyOuter     = 28
nxestension = 150
nz = 1
ncells = 23296

x_mesh, y_mesh = compute_mesh(nxOuter, nQuarter, nxOuter2, nxestension, nyOuter)
print(f'mesh for RANS - sim   : x = {x_mesh}, y = {y_mesh}, z = {nz}, ncells = {ncells}')



# RANS - coarse
nQuarter    = np.int64(np.round( 14*0.8,0))
nxOuter     = np.int64(np.round( 20*0.8,0))
nxOuter2    = np.int64(np.round( 30*0.8,0))
nyOuter     = np.int64(np.round( 28*0.8,0))
nxestension = np.int64(np.round(150*0.8,0))
nz = 1

x_mesh, y_mesh = compute_mesh(nxOuter, nQuarter, nxOuter2, nxestension, nyOuter)
print(f'mesh for RANS - coarse: x = {x_mesh}, y = {y_mesh}')


# RANS - finer
nQuarter    = np.int64(np.round( 14*1.2,0))
nxOuter     = np.int64(np.round( 20*1.2,0))
nxOuter2    = np.int64(np.round( 30*1.2,0))
nyOuter     = np.int64(np.round( 28*1.2,0))
nxestension = np.int64(np.round(150*1.2,0))
nz = 1

x_mesh, y_mesh = compute_mesh(nxOuter, nQuarter, nxOuter2, nxestension, nyOuter)
print(f'mesh for RANS - fine  : x = {x_mesh}, y = {y_mesh}')





'''
LES
'''

# LES - keq
nQuarter    =  14
nxOuter     =  20
nxOuter2    = 100
nyOuter     = 28
nxestension = 500 
nz = 14 
ncells = 2529800

x_mesh, y_mesh = compute_mesh(nxOuter, nQuarter, nxOuter2, nxestension, nyOuter)
print(f'mesh for LES - keq    : x = {x_mesh}, y = {y_mesh}, z = {nz}, ncells = {ncells}')

# LES - keq turb
nQuarter    =  14
nxOuter     =  20
nxOuter2    = 100
nyOuter     = 28
nxestension = 500 
nz = 14 
ncells = 2529800

x_mesh, y_mesh = compute_mesh(nxOuter, nQuarter, nxOuter2, nxestension, nyOuter)
print(f'mesh for LES - keqTurb: x = {x_mesh}, y = {y_mesh}, z = {nz}, ncells = {ncells}')


# LES - smago
nQuarter    =  14
nxOuter     =  20
nxOuter2    = 100
nyOuter     = 28
nxestension = 500 
nz = 14 
ncells = 2529800

x_mesh, y_mesh = compute_mesh(nxOuter, nQuarter, nxOuter2, nxestension, nyOuter)
print(f'mesh for LES - sma    : x = {x_mesh}, y = {y_mesh}, z = {nz}, ncells = {ncells}')



# LES - smago turb
nQuarter    =  14
nxOuter     =  20
nxOuter2    = 100
nyOuter     = 28
nxestension = 500
nz = 14 
ncells = 2529800

x_mesh, y_mesh = compute_mesh(nxOuter, nQuarter, nxOuter2, nxestension, nyOuter)
print(f'mesh for LES - smaTurb: x = {x_mesh}, y = {y_mesh}, z = {nz}, ncells = {ncells}')




les_smago = 287163
rans_komega = 1201
print('time les smago [h]', les_smago/60/60)
print('time rans komega [min]', rans_komega/60)
