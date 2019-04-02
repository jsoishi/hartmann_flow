import sys
import pathlib
import h5py
import matplotlib.pyplot as plt

basedir = pathlib.Path(sys.argv[-1])
dfile = basedir / 'integrals/integrals.h5'

data = h5py.File(str(dfile), "r")

y = data['scales/y/1'][:]
vx = data['tasks/<vx>_x'][-1,0,:]

plt.semilogy(y, vx)
plt.xlabel("y")
plt.ylabel(r"$<v_x>$")

plt.savefig('vx_yprof.png')

