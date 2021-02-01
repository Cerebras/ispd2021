dim = 2
vol_x = 9
vol_y = 9
gridpoints_per_tile_edge = 1
fab_x = 200
fab_y = 200
cost_a = 1
cost_b = 1

inv_sampling_step = 4

with open("ispd/tests/circle/circle.in", "w+") as f:
    f.write(F"{dim}\n")
    f.write(F"{vol_x} {vol_y}\n")
    f.write(F"{gridpoints_per_tile_edge}\n")
    f.write(F"{fab_x} {fab_y}\n")
    f.write(F"{cost_a} {cost_b}\n")

    half_y = int((vol_y + 1) / 2)
    half_x = int((vol_x + 1) / 2)

    radius = min(vol_x + 1, vol_y + 1)

    for i in range(0, vol_y + 1):
        for j in range(0, vol_x + 1):
            value = max(1 - ((i - half_y) ** 2 + (j - half_x) ** 2) / radius, 0.01)
            f.write(F"{value}\n")

with open("ispd/tests/circle/circle.out", "w+") as f:
    f.write(F"{inv_sampling_step}\n")

    compute_num = 0

    for i in range(0, 3):
        for j in range(0, 3):
            x_point = int(vol_x * inv_sampling_step / 3 * j)
            y_point = int(vol_y * inv_sampling_step / 3 * i)

            x_shape = int(vol_x * inv_sampling_step / (3 * 4))
            y_shape = int(vol_y * inv_sampling_step / (3 * 4))
            if i == 1 and j == 1:
                x_shape *= 4
                y_shape *= 4
                f.write(
                    F"0 {x_point} {y_point} {x_shape} {y_shape}\n")
            elif i == 1 or j == 1:
                x_shape *= 2
                y_shape *= 2
                f.write(
                    F"1 {x_point} {y_point} {x_shape} {y_shape}\n")
            else:
                f.write(
                    F"2 {x_point} {y_point} {x_shape} {y_shape}\n")
            compute_num += x_shape * y_shape

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
