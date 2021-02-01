dim = 3
vol_x = 20
vol_y = 40
vol_z = 200
gridpoints_per_tile_edge = 10
fab_x = 400
fab_y = 400
cost_a = 1
cost_b = 1

inv_sampling_step = gridpoints_per_tile_edge

with open("ispd/tests/3d/flat/flat_big.in", "w+") as f:
    f.write(F"# Dimension\n{dim}\n")
    f.write(F"# Volume X, Y, Z\n{vol_x} {vol_y} {vol_z}\n")
    f.write(F"# Fabric X, Y\n{fab_x} {fab_y}\n")
    f.write(F"# Cost Params\n{cost_a} {cost_b}\n")

    f.write(F"# Heatmap\n")
    for i in range(0, vol_z + 1):
        for j in range(0, vol_y + 1):
            f.write(F"# z, y = {i}, {j}\n")
            for k in range(0, vol_x + 1):
                f.write("1\n")


with open("ispd/tests/3d/flat/flat_big.out", "w+") as f:
    f.write(F"# Inv Sampling Step\n{inv_sampling_step}\n")

    f.write(
        F"0 0 0 0 {vol_x} {vol_y} {vol_z}\n")

    f.write(F"compute_map:\n")

    for i in range(0, fab_x):
        for j in range(0, fab_y):
            f.write(F"{i} {j}\n")

    f.write(F"adapter_map:\n")
