import argparse
from typing import List
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re

def get_next_line(file: str) -> str:
    with open(file) as f:
        line = f.readline()
        while line != '':
            if not line.startswith('#'):
                yield line

            line = f.readline()

    while True:
        yield ''


def visualize_placements(problem: str, solution: str, out: str):
    prob_iter = get_next_line(problem)

    dim = int(next(prob_iter))
    volume = [int(x) for x in re.findall(r'\d+', next(prob_iter))]
    gpte = int(next(prob_iter))
    fabric = [int(x) for x in re.findall(r'\d+', next(prob_iter))]
    cost = [float(x) for x in re.findall(r'\d+', next(prob_iter))]

    print(dim)
    print(volume)
    print(gpte)
    print(fabric)
    print(cost)
    
    sol_iter = get_next_line(solution)

    sample_step = int(next(sol_iter))

    line = next(sol_iter)

    prism_list = []
    tile_list = []

    while not line.startswith('compute_map:'):
        prism = [int (x) for x in re.findall(r'\d+', line)]

        prism_list.append(prism)
        if dim == 2:
            for i in range(0, prism[3]):
                for j in range(0, prism[4]):
                    tile_list.append(1.0 / (2 ** prism[0]))
        else:
            for i in range(0, prism[4]):
                for j in range(0, prism[5]):
                    for k in range(0, prism[6]):
                        tile_list.append(1.0 / (2 ** prism[0]))

        line = next(sol_iter)

    x = []
    y = []

    compute_tiles = len(tile_list)

    line = next(sol_iter)
    while not line.startswith('adapter_map:'):
        point = [int(x) for x in re.findall(r'\d+', line)]
        x.append(point[0])
        y.append(point[1])

        line = next(sol_iter)

    line = next(sol_iter)

    x_adapter = []
    y_adapter = []
    while line != '':
        point = [int(x) for x in re.findall(r'\d+', line)]
        x_adapter.append(point[0])
        y_adapter.append(point[1])

        line = next(sol_iter)

    pe_map = np.zeros((fabric[0], fabric[1]), dtype=float)
    tile_id_map = -np.ones((fabric[0], fabric[1]), dtype=int)
    for i in range(0, len(tile_list)):
        pe_map[x[i]][y[i]] = tile_list[i]
        if tile_list[i] >= 0:
            tile_id_map[x[i]][y[i]] = i

    fig, ax = plt.subplots()
    c_mesh = ax.pcolormesh(pe_map, cmap='viridis',  vmin=0, vmax=1)

    # Write adapters
    for i in range(0, len(x_adapter)):
        rect = patches.Rectangle((y_adapter[i], x_adapter[i]), 1, 1, edgecolor='r',facecolor='none',linewidth=1)
        ax.add_patch(rect)

    ax.set_title(F'{solution} PE Placement')
    ax.axis([0, fabric[0], 0, fabric[1]])

    if compute_tiles < 400:
        ax.xaxis.set_ticks(np.arange(0, fabric[0] + 1, 1))
        ax.yaxis.set_ticks(np.arange(0, fabric[1] + 1, 1))
        for i in range(0, len(tile_list)):
            tile_id = tile_id_map[x[i]][y[i]]
            if tile_id >= 0:
                ax.text(y[i] + 0.5, x[i] + 0.5, tile_id, ha='center', va='center', color='black')

    fig.colorbar(c_mesh)
    fig.savefig(out, dpi=500)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize PE placements of solution for a given problem")
    parser.add_argument('problem', metavar='problem', type=str, help='.in file to parse')
    parser.add_argument('solution', metavar='solution', type=str, help='.out file to parse')

    parser.add_argument('-o', '--output', default='out.png', type=str, help='output directory')

    args = parser.parse_args()

    assert args.problem.endswith('.in')
    assert args.solution.endswith('.out')
    visualize_placements(args.problem, args.solution, args.output)
