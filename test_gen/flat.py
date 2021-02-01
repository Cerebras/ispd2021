dim = 3
vol_x = 20
vol_y = 40
vol_z = 200
gridpoints_per_tile_edge = 1
fab_x = 400
fab_y = 400
cost_a = 1
cost_b = 1

inv_sampling_step = 1

with open("ispd/tests/flat/flat_big.in", "w+") as f:
    f.write(F"{dim}\n")
    f.write(F"{vol_x} {vol_y} {vol_z}\n")
    f.write(F"{gridpoints_per_tile_edge}\n")
    f.write(F"{fab_x} {fab_y}\n")
    f.write(F"{cost_a} {cost_b}\n")

    for i in range(0, vol_z + 1):
        for j in range(0, vol_y + 1):
            for k in range(0, vol_x + 1):
                f.write("1\n")


with open("ispd/tests/flat/flat_big.out", "w+") as f:
    f.write(F"{inv_sampling_step}\n")

    f.write(
        F"0 0 0 0 {vol_x * inv_sampling_step} {vol_y * inv_sampling_step} {vol_z * inv_sampling_step}\n")

    f.write(F"compute_map:\n")

    for i in range(0, fab_x):
        for j in range(0, fab_y):
            f.write(F"{i} {j}\n")

    f.write(F"adapter_map:\n")
