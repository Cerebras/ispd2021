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
    fabric = [int(x) for x in re.findall(r'\d+', next(prob_iter))]
    cost = [float(x) for x in re.findall(r'\d+', next(prob_iter))]

    print(dim)
    print(volume)
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

    compute_tiles = len(tile_list)

    line = next(sol_iter)
    while not line.startswith('adapter_map:'):
        line = next(sol_iter)

    line = next(sol_iter)

    adapter_tiles = 0

    while line != '':
        adapter_tiles += 1
        line = next(sol_iter)

    tile_by_res = {}
    for tile in tile_list:
        if tile in tile_by_res:
            tile_by_res[tile] += 1
        else:
            tile_by_res[tile] = 1

    tile_labels = []
    tile_counts = []

    for key, value in tile_by_res.items():
        tile_labels.append(str(key))
        tile_counts.append(value)

    tile_labels.append("Adapters")
    tile_counts.append(adapter_tiles)

    tile_labels.append("Empty")
    tile_counts.append(fabric[0] * fabric[1] - compute_tiles - adapter_tiles)

    fig, ax = plt.subplots()
    
    ax.set_title(F'{out} PE Breakdown')
    ax.pie(tile_counts, labels=tile_labels, autopct='%1.2f%%')
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
