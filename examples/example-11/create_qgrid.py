import numpy as np
qgrid = np.array([3,3,3])

# we canot create only half of the k-point pairs, since the in the nnkpts block all kpoints must have the same number of neighbours
with open('./tmds-kpoints-bse.dat','w') as f:
  f.write(f'{int(qgrid[0]*qgrid[1]*qgrid[2])}\n')
  f.write(f'2\n')
  for i in range(0,qgrid[0]):
    for j in range(0,qgrid[1]):
        for k in range(0,qgrid[2]):
            f.write(f"{1/qgrid[0]*i:.8f} {1/qgrid[1]*j:.8f} {1/qgrid[2]*k:.8f}\n")
