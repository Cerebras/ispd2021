import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')  # suppress Xterm invocation
import matplotlib.pyplot as plt
import re

def get_next_line(file: str) -> str:
    with open(file) as f:
        line = f.readline()
        while line != '':
            if not line.startswith('#'):
                yield line

            line = f.readline()

parser = argparse.ArgumentParser(description="Visualize a .in problem file")
parser.add_argument('filename', metavar='filename', type=str, help='.in file to parse')


args = parser.parse_args()

N = 0

iterator = get_next_line(args.filename)

dim = int(next(iterator))
volume = [int(x) for x in re.findall(r'\d+', next(iterator))]
fabric = [int(x) for x in re.findall(r'\d+', next(iterator))]
cost = [float(x) for x in re.findall(r'\d+\.?\d*', next(iterator))]

N = volume[0]

heatmap = np.empty((N+1,N+1,N+1))

for i in range(N+1):
    for j in range(N+1):
        for k in range(N+1):
            heatmap[k, j, i] = float(next(iterator).strip())

NN = np.arange(N+1)
x = NN/N
y = NN/N
z = NN/N
d = 1/N

fig = plt.figure()
extent = np.min(x), np.max(x), np.min(z), np.max(z)
fig.suptitle('X-Z plane')
for i in range(1,16):
  ax = fig.add_subplot(3, 5, i)
  ax.set_aspect('equal', adjustable='box')
  ax.imshow(1-heatmap[:,(i-1)*N//16,:].transpose(1,0),interpolation='nearest', extent=extent,cmap='gray')
  ax.set_xticks([], []) 
  ax.set_yticks([], []) 
fig.savefig('out_xz.png',dpi=200)
plt.close()

fig = plt.figure()
extent = np.min(x), np.max(x), np.min(y), np.max(y)
fig.suptitle('X-Y plane')
for i in range(1,16):
  ax = fig.add_subplot(3, 5, i)
  ax.set_aspect('equal', adjustable='box')
  ax.imshow(1-heatmap[:,:,(i-1)*N//16].transpose(1,0),interpolation='nearest', extent=extent,cmap='gray')
  ax.set_xticks([], [])
  ax.set_yticks([], [])
fig.savefig('out_xy.png',dpi=200)
plt.close()