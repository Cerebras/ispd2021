dim = 3
vol_x = 15
vol_y = 15
vol_z = 15
gridpoints_per_tile_edge = 10
fab_x = 200
fab_y = 200
cost_a = 1
cost_b = 1

inv_sampling_step = 20

with open("ispd/tests/3d/box/box_big.in", "w+") as f:
    f.write(F"# Dimension\n{dim}\n")
    f.write(F"# Volume X Y Z\n{vol_x} {vol_y} {vol_z}\n")
    f.write(F"# Fabric X Y\n{fab_x} {fab_y}\n")
    f.write(F"# Cost Params\n{cost_a} {cost_b}\n")

    third_z = int(vol_z / 3)
    third_y = int(vol_y / 3)
    third_x = int(vol_x / 3)
    f.write(F"# Heatmap\n")

    for i in range(0, vol_z + 1):
        for j in range(0, vol_y + 1):
            f.write(F"# z, y = {i}, {j}\n")
            for k in range(0, vol_x + 1):
                if i >= third_z and i < 2 * third_z and j >= third_y and j < 2 * third_y and j >= third_x and j < 2 * third_y:
                    f.write("1\n")
                else:
                    f.write("0.5\n")

with open("ispd/tests/3d/box/box_big.out", "w+") as f:
    f.write(F"# Inverse sampling step\n{inv_sampling_step}\n")

    compute_num = 0

    f.write(F"# Prisms\n")

    for i in range(0, 3):
        for j in range(0, 3):
            for k in range(0, 3):
                x_point = int(vol_x * inv_sampling_step / 3 * k)
                y_point = int(vol_y * inv_sampling_step / 3 * j)
                z_point = int(vol_z * inv_sampling_step / 3 * i)

                x_shape = int(vol_x * inv_sampling_step / gridpoints_per_tile_edge / (3 * 2))
                y_shape = int(vol_y * inv_sampling_step / gridpoints_per_tile_edge / (3 * 2))
                z_shape = int(vol_z * inv_sampling_step / gridpoints_per_tile_edge / (3 * 2))

                if i == 1 and j == 1 and k == 1:
                    x_shape *= 2
                    y_shape *= 2
                    z_shape *= 2
                    f.write(F"0 {x_point} {y_point} {z_point} {x_shape} {y_shape} {z_shape}\n")
                else:    
                    f.write(F"1 {x_point} {y_point} {z_point} {x_shape} {y_shape} {z_shape}\n")
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
    