import argparse
from typing import List
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

def get_next_line(file: str) -> str:
    with open(file) as f:
        line = f.readline()
        while line != '':
            if not line.startswith('#'):
                yield line

            line = f.readline()


def visualize_heatmap(filename: str, out: str, interpolate: bool):

    iterator = get_next_line(filename)

    dim = int(next(iterator))
    volume = [int(x) for x in re.findall(r'\d+', next(iterator))]
    fabric = [int(x) for x in re.findall(r'\d+', next(iterator))]
    cost = [float(x) for x in re.findall(r'\d+\.?\d*', next(iterator))]

    print(dim)
    print(volume)
    print(fabric)
    print(cost)

    out_file = Path(out)

    if dim == 2:
        
        data = []
        data_row = []
        count = 0
        for line in iterator:
            data_row.append(float(line))

            if count == volume[0]:
                data.append(data_row)
                data_row = []
                count = 0
            else:
                count +=1
        
        heatmap = np.array(data)

        fig, ax = plt.subplots()

        c_mesh = ax.pcolormesh(heatmap, cmap='gist_heat_r',  shading='gouraud' if interpolate else 'flat', vmin=np.min(heatmap), vmax=np.max(heatmap))

        ax.set_title(F'{filename} Heatmap')
        ax.axis([0, volume[0], 0, volume[1]])
        fig.colorbar(c_mesh)

        fig.savefig(str(out_file), dpi=500)

    else:

        data = []
        data_page = []
        data_row = []
        x_count = 0
        y_count = 0
        for line in iterator:
            data_row.append(float(line))

            if x_count == volume[0]:
                data_page.append(data_row)
                data_row = []
                x_count = 0

                if y_count == volume[1]:
                    data.append(data_page)
                    data_page = []
                    y_count = 0
                else:
                    y_count += 1
            else:
                x_count +=1
        
        heatmap = np.array(data)

        x = []
        y = []
        z = []
        c = []

        for i in range(0, volume[2] + 1):
            for j in range(0, volume[1] + 1):
                for k in range(0, volume[0] + 1):
                    x.append(k)
                    y.append(j)
                    z.append(i)
                    c.append(heatmap[i, j, k])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        c_map = plt.get_cmap('gist_heat_r')
        c_norm = matplotlib.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
        scalar_map = matplotlib.cm.ScalarMappable(norm=c_norm, cmap=c_map)

        ax.set_title(F'{filename} Heatmap')
        ax.scatter(x, y, z, c=scalar_map.to_rgba(c), marker='|')
        fig.colorbar(scalar_map)
        fig.savefig(str(out_file), dpi=500)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize a .in problem file")
    parser.add_argument('filename', metavar='filename', type=str, help='.in file to parse')
    parser.add_argument('-o', '--output', default='out.png', type=str, help='output directory')
    parser.add_argument('-i', '--interpolate', action='store_true', help='interpolate the heatmap values')

    args = parser.parse_args()

    assert args.filename.endswith('.in')
    visualize_heatmap(args.filename, args.output, args.interpolate)