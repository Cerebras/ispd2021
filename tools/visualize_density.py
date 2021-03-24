import argparse
from typing import List
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import re

def get_next_line(file: str) -> str:
    with open(file) as f:
        line = f.readline()
        while line != '':
            if not line.startswith('#'):
                yield line

            line = f.readline()


parser = argparse.ArgumentParser(description="Visualize a .in problem file")
parser.add_argument('input', metavar='input', type=str, help='.in file to parse')
parser.add_argument('output', metavar='output', type=str, help='.out file to parse')
parser.add_argument('-f', '--filename', default='out.png', type=str, help='output directory')

args = parser.parse_args()

in_iterator = get_next_line(args.input)

dim = int(next(in_iterator))
volume = [int(x) for x in re.findall(r'\d+', next(in_iterator))]
fabric = [int(x) for x in re.findall(r'\d+', next(in_iterator))]
cost = [float(x) for x in re.findall(r'\d+\.?\d*', next(in_iterator))]

out_iterator = get_next_line(args.output)

inv_sampling_step = int(next(out_iterator))

prisms = []

line = str(next(out_iterator)).strip()

max_res = 0

while not line.startswith("compute_map:"):
    prism = [int(x) for x in re.findall(r'-?\d+', line)]

    if prism[0] >= 0:
        prisms.append(prism)
        if max_res < prism[0]:
            max_res = prism[0]

    line = str(next(out_iterator)).strip()

z = int(volume[2] * inv_sampling_step / 2)

fig, ax = plt.subplots()

rects = []
resols = []
# Write adapters
for prism in prisms:

    z_low = prism[3]
    z_high = z_low + (prism[6] * 10) * 2** prism[0] - 1
    if z_low <= z and z <= z_high:
        rect = patches.Rectangle((prism[1], prism[2]), (prism[4] * 10) * 2** prism[0] - 1, (prism[5] * 10) * 2** prism[0] - 1)
        rects.append(rect)
        resols.append(1.0 / (2**prism[0]))

r = PatchCollection(rects, cmap=matplotlib.cm.jet)
r.set_array(np.array(resols))

ax.add_collection(r)
ax.set_title(F'{args.output} PE Placement')
ax.axis([0, volume[0] * inv_sampling_step, 0, volume[1] * inv_sampling_step])
fig.colorbar(r)

fig.savefig(args.filename, dpi=500)