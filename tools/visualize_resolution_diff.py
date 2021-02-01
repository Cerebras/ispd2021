import argparse
from typing import List
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

def visualize_resolution_diff(filename: str, out: str):
    # Find resolution list
    with open(filename) as f:
        line = f.readline().strip()
        while line != '' and not line.startswith('Resolution Difference'):
            line = f.readline().strip()
        
        target_res = []
        tile_res = []

        line = f.readline().strip()

        scale = 1

        while not line.endswith('End Resolution Difference'):
            res = [float(x) for x in re.findall(r'\d+\.?\d*', line)]

            target_res.append(res[1])
            tile_res.append(res[0])

            scale = max(scale, 1 + (res[1] - res[0]) / res[0])
            line = f.readline().strip()

        scaled_res = [x * scale for x in tile_res]

        x = np.arange(0, len(tile_res), 1, dtype=float)

        width = 0.3
        fig, ax = plt.subplots()

        ax.set_title(F'{filename} Tile Resolution Difference')
        ax.set_xlabel('Tile ID')
        ax.set_ylabel('Resolution')

        ax.bar(x - width, tile_res, width=width, color='r', align='center', label='Tile Resolution')
        ax.bar(x, scaled_res, width=width, color='g', align='center', label='Scaled Tile Resolution')
        ax.bar(x + width, target_res, width=width, color='b', align='center', label='Target Resolution')
        ax.xaxis.set_ticks(np.arange(0, len(x), 1))

        ax.legend()

        fig.savefig(out, dpi=500)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize PE placements of solution for a given problem")
    parser.add_argument('data', metavar='data', type=str, help='validator output data')

    parser.add_argument('-o', '--output', default='out.png', type=str, help='output directory')

    args = parser.parse_args()

    visualize_resolution_diff(args.data, args.output)
