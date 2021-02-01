dim = 2
vol_x = 32
vol_y = 32
gridpoints_per_tile_edge = 1
fab_x = 256
fab_y = 256
cost_a = 1
cost_b = 1

inv_sampling_step = 8

with open("ispd/tests/gradient/gradient_big.in", "w+") as f:
    f.write(F"{dim}\n")
    f.write(F"{vol_x} {vol_y}\n")
    f.write(F"{gridpoints_per_tile_edge}\n")
    f.write(F"{fab_x} {fab_y}\n")
    f.write(F"{cost_a} {cost_b}\n")

    for i in range(0, vol_y + 1):
        for j in range(0, vol_x + 1):
            f.write(F"{float(i + 1) / (vol_y + 1)}\n")

with open("ispd/tests/gradient/gradient_big.out", "w+") as f:
    f.write(F"{inv_sampling_step}\n")
    
    y_sum = 0
    y_height = vol_y / 4
    compute_num = 0
    for i in reversed(range(0, 4)):
        x_dim = int(vol_x * inv_sampling_step / (2 ** i * gridpoints_per_tile_edge))
        y_dim = int(vol_y * inv_sampling_step / (2 ** i * gridpoints_per_tile_edge * 4))
        f.write(F"{i} 0 {y_sum} {x_dim} {y_dim}\n")
        y_sum += int(vol_y * inv_sampling_step / 4)
        compute_num += int(x_dim * y_dim)

    f.write("compute_map:\n")
    curr_x = 0
    curr_y = 0
    for i in range(0, compute_num):
        f.write(F"{curr_x} {curr_y}\n")
        curr_x += 1
        if curr_x >= fab_x:
            curr_y += 1
            curr_x = 0
    f.write("adapter_map:\n")
