'''
Generates .in problem definition file for a 3D torus.
Usage:
python hollow_sphere.py [--N N] [--R R] [--r r] [--o OUT_DIR]

> --N N : Indicates size of volume. Volume is assumed to be a rectangular prism with dimension (2N x 2N x N)
> --R R : Indicates the radius of the overall torus
> --r r : Indicates the radius of the ring girth of the torus
> --o OUT_DIR : Indicates the directory to write output files
'''

import numpy as np
import matplotlib as mpl
mpl.use('Agg')  # suppress Xterm invocation
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(
  description='Heatmap generator for a 3D torus'
)
parser.add_argument(
  '--N',
  type=int,
  default=50
)
parser.add_argument(
  '--R',
  type=float,
  default=0.39
)
parser.add_argument(
  '--r',
  type=float,
  default=0.1
)
parser.add_argument(
  '--o',
  type=str,
  default='.'
)

args = parser.parse_args()
N = args.N
N2 = 2*N
R = args.R
r = args.r

out_dir = Path(args.o)

NN2 = np.arange(N2+1)
x = -0.5 + NN2/N2
y = -0.5 + NN2/N2
NN = np.arange(N+1)
z = -0.25 + NN/N2
d = 1/N2

h3 = np.empty((N2+1,N2+1,N+1))
for i in range(N2+1):
  for j in range(N2+1):
    for k in range(N+1):
      rr = np.sqrt(x[i]*x[i] + y[j]*y[j]) #xy projection 
      if rr == 0: #projects to origin
        x1 = R    #arbitrary point on R-radius
        y1 = 0 
      else:
        x1 = x[i]/rr*R  #corresponding point on the R-circle
        y1 = y[j]/rr*R
      r0 = np.sqrt((x[i]-x1)*(x[i]-x1) + (y[j]-y1)*(y[j]-y1) + z[k]*z[k])
      if r0<r:
        h3[i,j,k] = 1.0
      elif r<=r0 and r0<(r+d):
        h3[i,j,k] = (r+d-r0)/d
      else:
        h3[i,j,k]=0.0

fig = plt.figure()
extent = np.min(x), np.max(x), np.min(z), np.max(z)
fig.suptitle('Torus, X-Z plane')
for i in range(1,16):
  ax = fig.add_subplot(3, 5, i)
  ax.set_aspect('equal', adjustable='box')
  ax.imshow(1-h3[:,i*N2//16,:].transpose(1,0),interpolation='nearest', extent=extent,cmap='gray')
  ax.set_xticks([], []) 
  ax.set_yticks([], []) 
fig.savefig('torus_xz.png',dpi=200)
plt.close()

fig = plt.figure()
extent = np.min(x), np.max(x), np.min(y), np.max(y)
fig.suptitle('Torus, X-Y plane')
for i in range(13,28):
  ax = fig.add_subplot(3, 5, i-12)
  ax.set_aspect('equal', adjustable='box')
  ax.imshow(1-h3[:,:,i*N//39].transpose(1,0),interpolation='nearest', extent=extent,cmap='gray')
  ax.set_xticks([], [])
  ax.set_yticks([], [])
fig.savefig(out_dir / 'torus_xy.png',dpi=200)
plt.close()

with open(out_dir / F"torus-{N}.in", mode='w') as f:
  f.write('3\n')
  f.write(F'{N2} {N2} {N}\n')
  f.write('800 1600\n')
  f.write('1 1\n')
  for i in range(N+1):
    for j in range(N2+1):
      for k in range(N2+1):
        f.write(F'{h3[k,j,i]}\n')
