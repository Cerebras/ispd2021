'''
Generates .in problem definition files for a 2D ring and 3D hollow spheres.
Usage:
python hollow_sphere.py [--N N] [--r1 R1] [--r2 R2] [--o OUT_DIR]

> --N N : Indicates size of volume. Volume is assumed to be a cube
> --r1 R1 : Indicates the inner radius of the ring/hollow sphere
> --r2 R2 : Indicates the outer radius of the ring/hollow sphere
> --o OUT_DIR : Indicates the directory to write output files
'''

import numpy as np
import matplotlib as mpl
mpl.use('Agg')  # suppress Xterm invocation
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(
  description='Heatmap generator for 2D ring and 3D hollow sphere'
)
parser.add_argument(
  '--N',
  type=int,
  default=100
)
parser.add_argument(
  '--r1',
  type=float,
  default=0.45
)
parser.add_argument(
  '--r2',
  type=float,
  default=0.48
)
parser.add_argument(
  '--o',
  type=str,
  default='.'
)

args = parser.parse_args()
N = args.N
r1 = args.r1
r2 = args.r2

out_dir = Path(args.o)

NN = np.arange(N+1)
x = -0.5 + NN/N
y = -0.5 + NN/N
d = 1/N

h2 = np.empty((N+1,N+1))
for i in range(N+1):
  for j in range(N+1):
    r = np.sqrt(x[i]*x[i] + y[j]*y[j])
    if r<(r1-d):
      h2[i,j]=0.0
    elif (r1-d)<=r and r<r1:
      h2[i,j] = (d + r - r1)/d
    elif r1<=r and r<r2:
      h2[i,j] = 1.0
    elif r2<=r and r<(r2+d):
      h2[i,j] = (r2+d-r)/d
    else:
      h2[i,j]=0.0

extent = np.min(x), np.max(x), np.min(y), np.max(y)
fig = plt.figure()
fig.suptitle('2D Ring')
ax1 = fig.add_subplot(1, 1, 1)
ax1.set_aspect('equal', adjustable='box')
ax1.imshow(1-h2,interpolation='nearest', extent=extent,cmap='gray')
fig.savefig(out_dir / 'ring.png',dpi=200)
plt.close()

with open(F"ring-{N}.in", mode='w') as f:
  f.write('2\n')
  f.write(F'{N} {N}\n')
  f.write('800 1600\n')
  f.write('1 1\n')
  for i in range(N+1):
    for j in range(N+1):
      f.write(F'{h2[i,j]}\n')

z = -0.5 + NN/N

h3 = np.empty((N+1,N+1,N+1))
for i in range(N+1):
  for j in range(N+1):
    for k in range(N+1):
      r = np.sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
      if r<(r1-d):
        h3[i,j,k]=0.0
      elif (r1-d)<=r and r<r1:
        h3[i,j,k] = (d + r - r1)/d
      elif r1<=r and r<r2:
        h3[i,j,k] = 1.0
      elif r2<=r and r<(r2+d):
        h3[i,j,k] = (r2+d-r)/d
      else:
        h3[i,j,k]=0.0

fig = plt.figure()
extent = np.min(x), np.max(x), np.min(y), np.max(y)
fig.suptitle('3D Sphere, X-Y plane')
for i in range(1,16):
  ax = fig.add_subplot(3, 5, i)
  ax.set_aspect('equal', adjustable='box')
  ax.imshow(1-h3[:,i*N//16,:].transpose(1,0),interpolation='nearest', extent=extent,cmap='gray')
  ax.set_xticks([], [])
  ax.set_yticks([], [])
fig.savefig(out_dir / 'hollow_sphere_xz.png',dpi=200)
plt.close()


with open(out_dir / F"hollow_sphere-{N}.in", mode='w') as f:
  f.write('3\n')
  f.write(F'{N} {N} {N}\n')
  f.write('800 1600\n')
  f.write('1 1\n')
  for i in range(N+1):
    for j in range(N+1):
      for k in range(N+1):
        f.write(F'{h3[k,j,i]}\n')
