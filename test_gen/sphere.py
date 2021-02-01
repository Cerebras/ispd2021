dim = 3
vol_x = 15
vol_y = 15
vol_z = 15
gridpoints_per_tile_edge = 10
fab_x = 400
fab_y = 400
cost_a = 1
cost_b = 1

inv_sampling_step = 80

with open("ispd/tests/3d/sphere/sphere.in", "w+") as f:
    f.write(F"# Dimension\n{dim}\n")
    f.write(F"# Volume X, Y, Z\n{vol_x} {vol_y} {vol_z}\n")
    f.write(F"# Fabric X, Y\n{fab_x} {fab_y}\n")
    f.write(F"# Cost Params\n{cost_a} {cost_b}\n")

    half_z = int((vol_z + 1) / 2)
    half_y = int((vol_y + 1) / 2)
    half_x = int((vol_x + 1) / 2)

    radius = min(vol_x + 1, vol_y + 1, vol_z + 1)

    f.write(F"# Heatmap\n")

    for i in range(0, vol_z + 1):
        for j in range(0, vol_y + 1):
            f.write(F"# z, y = {i}, {j}\n")
            for k in range(0, vol_x + 1):
                value = max(1 - ((i - half_z) ** 2 + (j - half_y) ** 2 + (k - half_x) ** 2) / radius, 0.01)
                f.write(F"{value}\n")

with open("ispd/tests/3d/sphere/sphere.out", "w+") as f:
    f.write(F"#Inverse sampling step\n{inv_sampling_step}\n")

    compute_num = 0
    f.write("# Prisms\n")
    for i in range(0, 3):
        for j in range(0, 3):
            for k in range(0, 3):
                x_point = int(vol_x * inv_sampling_step / 3 * k)
                y_point = int(vol_y * inv_sampling_step / 3 * j)
                z_point = int(vol_z * inv_sampling_step / 3 * i)

                x_shape = int(vol_x * inv_sampling_step / gridpoints_per_tile_edge / (3 * 8))
                y_shape = int(vol_y * inv_sampling_step / gridpoints_per_tile_edge / (3 * 8))
                z_shape = int(vol_z * inv_sampling_step / gridpoints_per_tile_edge / (3 * 8))

                if i == 1 and j == 1 and k == 1:
                    x_shape *= 8
                    y_shape *= 8
                    z_shape *= 8
                    f.write(
                        F"0 {x_point} {y_point} {z_point} {x_shape} {y_shape} {z_shape}\n")
                elif i == 1 or j == 1 or k == 1:
                    if (i == 1 and j == 1) or (i == 1 and k == 1) or (k == 1 and j == 1):
                        x_shape *= 4
                        y_shape *= 4
                        z_shape *= 4
                        f.write(
                            F"1 {x_point} {y_point} {z_point} {x_shape} {y_shape} {z_shape}\n")
                    else:
                        x_shape *= 2
                        y_shape *= 2
                        z_shape *= 2
                        f.write(
                            F"2 {x_point} {y_point} {z_point} {x_shape} {y_shape} {z_shape}\n")
                else:
                    f.write(
                        F"3 {x_point} {y_point} {z_point} {x_shape} {y_shape} {z_shape}\n")
                compute_num += x_shape * y_shape * z_shape

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
