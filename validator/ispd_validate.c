#include "ispd_validate.h"

#define DOUBLE_EPSILON 0.000000001 // Epsilon value for approximate double equality

#define WIRE_POW 1.5 // Power to which each wire score is raised (L1.5 Norm)

// Global variables
struct option opt;                        // Initialize a global option struct
int dimension = NDIM;                     // Dimension of problem/solution
struct problem *prob = 0;                 // Parsed problem
struct solution *sol = 0;                 // Parsed solution
node_bounds volume_bounds;                // Bounds of volume given by problem
int ncompute = 0;                         // Number of compute tiles given by solution
int nadaptor = 0;                         // Number of adapter tiles given by solution
tile_tree *tree = 0;                      // R-tree which holds the tiles of a solution
tile **tile_list = 0;                     // List of tiles in solution
struct connectivity *need_adapters = 0;   // Holds connectivity of all tiles that need adapters
struct connectivity *all_connections = 0; // Holds connectivity of all tiles

// Creates a tile
tile *get_tile(struct prism *p, int pid, int prev_tiles, int *coord)
{
    tile *prism_tile = calloc(1, sizeof(tile));

    for (int i = 0; i < dimension; i += 1)
    {
        prism_tile->tile_id[i + 1] = coord[i];
        prism_tile->bounds.low_bounds[i] = p->origin[i] + ((coord[i] * prob->gridpoints_per_tile_edge) << p->resolution);
        prism_tile->bounds.up_bounds[i] = p->origin[i] + (((coord[i] + 1) * prob->gridpoints_per_tile_edge) << p->resolution) - 1;
    }
    prism_tile->parent = p;
    prism_tile->prev_tiles = prev_tiles;
    prism_tile->tile_id[0] = pid;

    return prism_tile;
}

// Gets the tile id relating to a tile (prev_tile plus the position of the tile in its prism)
int get_tile_id(tile *t)
{
    int tid = t->prev_tiles;

    struct prism *p = t->parent;

    if (dimension == 3)
    {
        tid += t->tile_id[3] * p->shape[1] * p->shape[0]; // Z position
    }
    tid += t->tile_id[2] * p->shape[0]; // Y position
    tid += t->tile_id[1];               // X position
    return tid;
}

// Get upper bound of a prism in sampling space at the finest resolution
int get_close(struct prism *prism, int axis)
{
    return prism->origin[axis] + ((prism->shape[axis] * prob->gridpoints_per_tile_edge) << prism->resolution) - 1;
}

void print_tabs(int tabs)
{
    printf("%*s", 2 * tabs, "");
}

void print_problem(struct problem *prob)
{
    printf("Problem:\n");
    printf("Volume Shape:");
    for (int i = 0; i < dimension; i += 1)
        printf(" %d", prob->volume_shape[i]);
    printf("\n");

    printf("Gridpoints Per Tile Edge: %d\n", prob->gridpoints_per_tile_edge);

    printf("Fabric Shape: Length: %d Width: %d\n", prob->fabric_shape[0], prob->fabric_shape[1]);

    printf("Cost Params: %f %f\n", prob->cost.alpha, prob->cost.beta);

    printf("Heatmap:\n");
    int heatmap_num = 0;
    if (dimension == 2)
    {
        heatmap_num = prob->volume_shape[0] * prob->volume_shape[1];
    }
    else
    {
        heatmap_num = prob->volume_shape[0] * prob->volume_shape[1] * prob->volume_shape[2];
    }
    for (int i = 0; i < heatmap_num; i += 1)
    {
        printf("Point %d\n", i);
        print_tabs(1);
        printf("%f\n", prob->heatmap[i]);
    }
    printf("\n");
}

void print_point(struct pt *point, int tabs)
{
    print_tabs(tabs);
    printf("X: %d Y: %d\n", point->x, point->y);
}

// Prints a prism with given tabs. If print_range is not 0, also print the range of the prism
void print_prism(struct prism *p, int tabs, bool print_range)
{
    print_tabs(tabs);
    printf("Resolution: %d\n", p->resolution);

    print_tabs(tabs);
    printf("Origin:");
    for (int i = 0; i < dimension; i += 1)
    {
        printf(" %d", p->origin[i]);
    }
    printf("\n");

    print_tabs(tabs);
    printf("Shape:");
    for (int i = 0; i < dimension; i += 1)
    {
        printf(" %d", p->shape[i]);
    }
    printf("\n");

    // Optional range output
    if (print_range)
    {
        print_tabs(tabs);
        printf("Range:");
        for (int j = 0; j < dimension; j += 1)
        {
            printf(" %c: [%d,%d]", j == 0 ? 'x' : j == 1 ? 'y' : 'z',
                   p->origin[j],
                   get_close(p, j));
        }
        printf("\n");
    }
}

void print_solution(struct solution *sol, int ncompute, int nadapter)
{
    printf("Solution:\n");
    printf("Sampling Step: %d\n", sol->inv_sampling_step);
    printf("Number of Prisms: %d\n", sol->nprims);
    for (int i = 0; i < sol->nprims; i += 1)
    {
        printf("Prism %d:\n", i);
        print_prism(&sol->prism[i], 1, 0);
    }
    printf("Number of Compute Tiles: %d\n", ncompute);
    for (int i = 0; i < ncompute; i += 1)
    {
        printf("Point %d:\n", i);
        print_point(&sol->compute_map[i], 1);
    }
    printf("Number of Adapter Tiles: %d\n", nadapter);
    for (int i = 0; i < nadapter; i += 1)
    {
        printf("Point %d:\n", i);
        print_point(&sol->adapter_map[i], 1);
    }
    printf("\n");
}
// Check if line is compute map header
int check_compute_map(char *line)
{
    char *cmp = "compute_map:";
    for (int i = 0; i < 12; i += 1)
    {
        if (line[i] == '\0' || line[i] != cmp[i])
        {
            return 0;
        }
    }
    return 1;
}
// Check if line is adapter map header
int check_adapter_map(char *line)
{
    char *cmp = "adapter_map:";
    for (int i = 0; i < 12; i += 1)
    {
        if (line[i] == '\0' || line[i] != cmp[i])
        {
            return 0;
        }
    }
    return 1;
}

// Gets the next line of input to parse, will ignore empty lines and comments (lines starting with #)
ssize_t get_next_line(char **line_ptr, size_t *len, FILE *input)
{
    ssize_t read;
    while ((read = getline(line_ptr, len, input)) != -1)
    {
        char *line = *line_ptr;
        // Check for comment/empty line
        bool is_comment = 1;
        for (int i = 0; i < *len; i += 1)
        {
            if (!isspace(line[i]))
            {
                is_comment = line[i] == '#';

                break;
            }
        }

        if (is_comment)
        {
            continue;
        }

        // Trim trailing whitespace before returning
        for (int i = *len - 1; i >= 0; i -= 1)
        {
            if (isspace(line[i]))
            {
                line[i] = '\0';
            }
            else if (line[i] != '\0')
            {
                break;
            }
        }
        break;
    }
    return read;
}

int max(int a, int b)
{
    return a > b ? a : b;
}
int min(int a, int b)
{
    return a < b ? a : b;
}

double max_d(double a, double b)
{
    return a > b ? a : b;
}
double min_d(double a, double b)
{
    return a < b ? a : b;
}

double abs_d(double n)
{
    return n < 0 ? -n : n;
}

// Initialize node to starting values
void init_node(tile_node *node, node_type type, tile_node *parent)
{
    node->type = type;
    node->children = 0;
    node->num_children = 0;
    node->parent = parent;
    for (int i = 0; i < dimension; i += 1)
    {
        node->bounds.low_bounds[i] = 0;
        node->bounds.up_bounds[i] = 0;
    }
}

void print_bounds(const node_bounds *bounds, int tabs)
{
    print_tabs(tabs);
    printf("Bounds:");
    for (int j = 0; j < dimension; j += 1)
    {
        printf(" %c: [%d,%d]", j == 0 ? 'x' : j == 1 ? 'y' : 'z',
               bounds->low_bounds[j],
               bounds->up_bounds[j]);
    }
    printf("\n");
}

void print_tile(tile *t, int tabs)
{
    print_tabs(tabs);
    if (dimension == 2)
    {
        printf("Tile ID: %d i: %d x: %d y: %d\n", get_tile_id(t), t->tile_id[0], t->tile_id[1], t->tile_id[2]);
    }
    else
    {
        printf("Tile ID: %d i: %d x: %d y: %d z: %d\n", get_tile_id(t), t->tile_id[0], t->tile_id[1], t->tile_id[2], t->tile_id[3]);
    }

    print_bounds(&t->bounds, tabs);
}

// Test if two bounds intersect
bool intersects(const node_bounds *a, const node_bounds *b)
{
    for (int i = 0; i < dimension; i += 1)
    {
        int lower = a->low_bounds[i];
        int upper = a->up_bounds[i];
        if (a->low_bounds[i] > a->up_bounds[i])
        {
            lower ^= upper;
            upper ^= lower;
            lower ^= upper;
        }
        if (lower > b->low_bounds[i] && lower > b->up_bounds[i])
        {
            return 0;
        }
        if (upper < b->low_bounds[i] && upper < b->up_bounds[i])
        {
            return 0;
        }
    }

    return 1;
}

// Combines two bounds together and return the new bounds
node_bounds add_bounds(const node_bounds *a, const node_bounds *b)
{
    node_bounds new_bounds;
    for (int i = 0; i < dimension; i += 1)
    {
        new_bounds.low_bounds[i] = min(a->low_bounds[i], b->low_bounds[i]);
        new_bounds.up_bounds[i] = max(a->up_bounds[i], b->up_bounds[i]);
    }

    return new_bounds;
}

// Get volume of space denoted by bounds
int get_bounds_volume(node_bounds bounds)
{
    int value = 1;
    for (int i = 0; i < dimension; i += 1)
    {
        int dis = bounds.up_bounds[i] - bounds.low_bounds[i];
        fatalif(dis < 0, "Bounds distance is negative");
        value *= (dis + 1);
    }

    return value;
}

void free_node(tile_node *node)
{
    if (!node)
    {
        return;
    }

    if (node->type == INTERIOR)
    {
        tile_node **children = (tile_node **)node->children;
        if (children != 0)
        {
            for (int i = 0; i < node->num_children; i += 1)
            {
                free_node(children[i]);
            }
            v_free(children);
        }
    }
    free(node);
}

void free_tree(tile_tree *tree)
{
    if (!tree)
    {
        return;
    }

    free_node(tree->head);

    free(tree);
}

void print_node(tile_node *node, int depth)
{
    switch (node->type)
    {
    case INTERIOR:
    {
        print_tabs(depth);
        printf("Interior Node:\n");
        print_tabs(depth + 1);
        printf("Children: %d\n", node->num_children);
        print_bounds(&node->bounds, depth + 1);
        for (int i = 0; i < node->num_children; i += 1)
        {
            print_node(((tile_node **)node->children)[i], depth + 1);
        }
        break;
    }
    case LEAF:
    {
        print_tabs(depth);
        printf("Leaf Node:\n");
        print_tile((tile *)node->children, depth + 1);
        break;
    }
    default:
    {
        fatalif(1, "Unknown type");
    }
    }
}

// Inserts list of tile_node into a tile_tree based of Hilbert R-tree
tile_node **insert_tile_node(tile_tree *tree, tile_node **tiles, int ntiles)
{
    if (ntiles == 0)
    {
        return 0;
    }

    tile_node **internal_nodes = 0;

    tile_node *node = 0;
    int curr_children = 0;

    tile_node **children = 0;

    for (int i = 0; i < ntiles; i += 1)
    {
        if (curr_children == 0)
        {
            node = calloc(1, sizeof(tile_node));
            init_node(node, INTERIOR, 0);
            v_adjoin(internal_nodes, node);

            node->bounds = tiles[i]->bounds;
            children = (tile_node **)node->children;
        }
        else
        {
            node->bounds = add_bounds(&node->bounds, &tiles[i]->bounds);
        }

        tiles[i]->parent = node;
        v_adjoin(children, tiles[i]);
        node->children = children;
        node->num_children += 1;

        curr_children += 1;

        // Finish this level
        if (curr_children == tree->max_children)
        {
            curr_children = 0;
        }
    }

    return internal_nodes;
}

// Inserts list of tile into a tile_tree based of Hilbert R-tree
tile_node **insert_tile(tile_tree *tree, tile **tiles, int ntiles)
{
    if (ntiles == 0)
    {
        return 0;
    }

    tile_node **internal_nodes = 0;

    tile_node *node = 0;
    int curr_children = 0;

    tile_node **children = 0;

    for (int i = 0; i < ntiles; i += 1)
    {
        if (curr_children == 0)
        {
            node = calloc(1, sizeof(tile_node));
            init_node(node, INTERIOR, 0);
            v_adjoin(internal_nodes, node);

            node->bounds = tiles[i]->bounds;
            children = (tile_node **)node->children;
        }
        else
        {
            node->bounds = add_bounds(&node->bounds, &tiles[i]->bounds);
        }

        tile_node *leaf = calloc(1, sizeof(tile_node));
        init_node(leaf, LEAF, 0);
        leaf->children = tiles[i];
        leaf->num_children = 1;
        leaf->bounds = tiles[i]->bounds;

        leaf->parent = node;
        v_adjoin(children, leaf);
        node->children = children;
        node->num_children += 1;

        curr_children += 1;

        // Finish this level
        if (curr_children == tree->max_children)
        {
            curr_children = 0;
        }
    }

    return internal_nodes;
}

// Compare two tiles by spatial position, order goes x, then y then z
int cmp_tiles(const void *a, const void *b)
{
    tile **a_tile = (tile **)a;
    tile **b_tile = (tile **)b;
    int x_diff = (*a_tile)->bounds.low_bounds[0] - (*b_tile)->bounds.low_bounds[0];
    if (x_diff == 0)
    {
        int y_diff = (*a_tile)->bounds.low_bounds[1] - (*b_tile)->bounds.low_bounds[1];
        if (y_diff == 0 && dimension == 3)
        {
            return (*a_tile)->bounds.low_bounds[2] - (*b_tile)->bounds.low_bounds[2];
        }
        else
        {
            return y_diff;
        }
    }
    else
    {
        return x_diff;
    }
}

// Compare two tiles by tile id
int cmp_tiles_id(const void *a, const void *b)
{
    tile **a_tile = (tile **)a;
    tile **b_tile = (tile **)b;
    int diff = get_tile_id(*a_tile) - get_tile_id(*b_tile);
    return diff;
}

// Build tree using prism sized tiles (For overlap testing)
tile_tree *build_prism_tree(tile ***tiles, struct prism *prisms, int nprims, int max_children)
{
    tile_tree *tree = calloc(1, sizeof(tile_tree));

    tile **added_tiles = 0;

    tree->max_children = max_children;
    for (int i = 0; i < nprims; i += 1)
    {
        struct prism *p = &prisms[i];

        tile *prism_tile = calloc(1, sizeof(tile));

        for (int i = 0; i < dimension; i += 1)
        {
            prism_tile->tile_id[i + 1] = 0;
            prism_tile->bounds.low_bounds[i] = p->origin[i];
            prism_tile->bounds.up_bounds[i] = p->origin[i] + ((p->shape[i] * prob->gridpoints_per_tile_edge) << p->resolution) - 1;
        }
        prism_tile->parent = p;
        prism_tile->prev_tiles = i;
        prism_tile->tile_id[0] = i;

        v_adjoin(added_tiles, prism_tile);
    }

    // Sort tiles spatially for node construction
    qsort(added_tiles, nprims, sizeof(tile *), cmp_tiles);

    // Create lowest level interior nodes
    tile_node **interior_nodes = insert_tile(tree, added_tiles, nprims);

    // Build rest of interior nodes
    for (; v_count(interior_nodes) > 1;)
    {
        tile_node **new_level = insert_tile_node(tree, interior_nodes, v_count(interior_nodes));
        v_free(interior_nodes);
        interior_nodes = new_level;
    }

    tree->head = interior_nodes[0];

    v_free(interior_nodes);

    // Return tiles added to tree in tiles pointer
    *tiles = added_tiles;
    return tree;
}

// Build tree using PE sized tiles (Adjacency checking)
tile_tree *build_tile_tree(tile ***tiles, struct prism *prisms, int nprims, int max_children)
{
    tile_tree *tree = calloc(1, sizeof(tile_tree));

    tile **added_tiles = 0;

    tree->max_children = max_children;

    // Create tiles by iterating through prisms and their shape
    int ntiles = 0;
    for (int i = 0; i < nprims; i += 1)
    {
        int prev_prism_tiles = ntiles;
        struct prism *prism = &prisms[i];
        if (dimension == 2)
        {
            for (int x = 0; x < prism->shape[0]; x += 1)
            {
                for (int y = 0; y < prism->shape[1]; y += 1)
                {
                    int coords[2] = {x, y};
                    tile *curr_tile = get_tile(prism, i, prev_prism_tiles, coords);
                    v_adjoin(added_tiles, curr_tile);
                    ntiles += 1;
                }
            }
        }
        else
        {
            for (int x = 0; x < prism->shape[0]; x += 1)
            {
                for (int y = 0; y < prism->shape[1]; y += 1)
                {
                    for (int z = 0; z < prism->shape[2]; z += 1)
                    {
                        int coords[3] = {x, y, z};
                        tile *curr_tile = get_tile(prism, i, prev_prism_tiles, coords);
                        v_adjoin(added_tiles, curr_tile);
                        ntiles += 1;
                    }
                }
            }
        }
    }

    // Sort tiles spatially for building tree
    qsort(added_tiles, ntiles, sizeof(tile *), cmp_tiles);

    // Create lowest level interior nodes
    tile_node **interior_nodes = insert_tile(tree, added_tiles, ntiles);

    // Build rest of interior nodes
    for (; v_count(interior_nodes) > 1;)
    {
        tile_node **new_level = insert_tile_node(tree, interior_nodes, v_count(interior_nodes));
        v_free(interior_nodes);
        interior_nodes = new_level;
    }

    tree->head = interior_nodes[0];

    v_free(interior_nodes);

    // Sort tiles by id for later lookup
    qsort(added_tiles, ntiles, sizeof(tile *), cmp_tiles_id);

    // Return added tiles in tiles pointer
    *tiles = added_tiles;
    return tree;
}

// Range search a tile_node and return tiles found in tiles_found vector
void search_node(tile ***tiles_found, tile_node *node, const node_bounds *range)
{
    switch (node->type)
    {
    case INTERIOR:
    {
        tile_node **children = (tile_node **)node->children;
        for (int i = 0; i < node->num_children; i += 1)
        {
            if (intersects(&(children[i]->bounds), range))
            {
                search_node(tiles_found, children[i], range);
            }
        }
        break;
    }
    case LEAF:
    {
        if (intersects(&node->bounds, range))
        {
            tile **temp_list = *tiles_found;
            v_adjoin(temp_list, (tile *)node->children);
            *tiles_found = temp_list;
        }
        break;
    }
    default:
        break;
    }
}

// Range search tree for tiles that intersect a specified range
tile **search_tree(tile_tree *tree, node_bounds range)
{
    tile **tiles_found = 0;
    search_node(&tiles_found, tree->head, &range);
    return tiles_found;
}

void print_tree(tile_tree *tree)
{
    printf("Prism Tree:\nMax Children: %d\n", tree->max_children);
    print_node(tree->head, 0);
}

/* Connectivity in tile space, not prism space */
/* derives the connectivity between tiles 
 * based on the problem and solution */
void build_connectivity()
{
    need_adapters = calloc(1, sizeof(struct connectivity));
    all_connections = calloc(1, sizeof(struct connectivity));

    int ntiles = v_count(tile_list);

    for (int i = 0; i < ntiles; i += 1)
    {
        tile *t = tile_list[i];
        int curr_tid = get_tile_id(t);
        int has_res_diff = 0;
        v_adjoin(all_connections->current_tile, curr_tid);
        all_connections->ntiles += 1;
        for (int j = 0; j < dimension; j += 1)
        {
            if (t->bounds.low_bounds[j] > volume_bounds.low_bounds[j])
            {
                node_bounds adjacent_range;
                for (int k = 0; k < dimension; k += 1)
                {
                    adjacent_range.low_bounds[k] = t->bounds.low_bounds[k];
                    adjacent_range.up_bounds[k] = t->bounds.up_bounds[k];
                }

                // Change current dimension to 1 unit slice to search adjacent
                adjacent_range.low_bounds[j] -= 1;
                adjacent_range.up_bounds[j] = adjacent_range.low_bounds[j];

                // Range search for adjacent tiles
                tile **curr_adj = search_tree(tree, adjacent_range);

                if (v_count(curr_adj) > (dimension == 2 ? 2 : 4))
                {
                    print_tile(t, 0);
                    fatalif(1, "Tile has more than 4 neighbours for this face");
                }
                // Assume curr_adj has at most 4 adjacent
                for (int k = 0; k < 4; k += 1)
                {
                    // Continue if no more adjacent tiles
                    if (k >= v_count(curr_adj))
                    {
                        v_adjoin(all_connections->adjacent_tiles[j * 2][k], -1);
                        continue;
                    }
                    // Only go from low res to high res (e.g. 1 to 0)
                    int res_diff = t->parent->resolution - curr_adj[k]->parent->resolution;
                    int adj_tid = get_tile_id(curr_adj[k]);

                    v_adjoin(all_connections->adjacent_tiles[j * 2][k], adj_tid);

                    // Ignore other types of connections, (only adapter connections)
                    if (res_diff == 1)
                    {
                        fatalif(v_count(curr_adj) != (dimension == 2 ? 2 : 4),
                                "Tile %d alignment is off, adjacent to more than %d (%d) tile of low to high resolution",
                                curr_tid, (dimension == 2 ? 2 : 4), v_count(curr_adj));

                        v_adjoin(need_adapters->adjacent_tiles[j * 2][k], adj_tid);
                        if (!has_res_diff)
                        {
                            v_adjoin(need_adapters->current_tile, curr_tid);
                            need_adapters->ntiles += 1;
                        }
                        has_res_diff = 1;
                    }
                    else if (res_diff == 0)
                    {
                        fatalif(v_count(curr_adj) != 1, "Tile %d alignment is off, adjacent to more than 1 (%d) tile of same resolution", curr_tid, v_count(curr_adj));
                    }
                    else if (res_diff == -1)
                    {
                        fatalif(v_count(curr_adj) != 1,
                                "Tile %d alignment is off, adjacent to more than 1 (%d) tile of high to low resolution",
                                curr_tid, v_count(curr_adj));
                    }
                    else if (res_diff != 0 && res_diff != -1)
                    {
                        fatalif(1, "Resolution can only jump by power of 2 between tile %d and %d", curr_tid, adj_tid);
                    }
                }

                v_free(curr_adj);
            }
            else
            {
                // If edge, just add no connections to all_connection list
                for (int k = 0; k < 4; k += 1)
                {
                    v_adjoin(all_connections->adjacent_tiles[j * 2][k], -1);
                }
            }

            if (t->bounds.low_bounds[j] < volume_bounds.up_bounds[j])
            {
                node_bounds adjacent_range;
                for (int k = 0; k < dimension; k += 1)
                {
                    adjacent_range.low_bounds[k] = t->bounds.low_bounds[k];
                    adjacent_range.up_bounds[k] = t->bounds.up_bounds[k];
                }

                // Change current dimension to 1 unit slice
                adjacent_range.up_bounds[j] += 1;
                adjacent_range.low_bounds[j] = adjacent_range.up_bounds[j];

                tile **curr_adj = search_tree(tree, adjacent_range);

                if (v_count(curr_adj) > (dimension == 2 ? 2 : 4))
                {
                    print_tile(t, 0);
                    fatalif(1, "Tile has more than 4 neighbours for this face");
                }
                // Assume curr_adj has at most 4 adjacent
                for (int k = 0; k < 4; k += 1)
                {
                    // Continue if no more adjacent tiles
                    if (k >= v_count(curr_adj))
                    {
                        v_adjoin(all_connections->adjacent_tiles[j * 2 + 1][k], -1);
                        continue;
                    }
                    // Only go from low res to high res (e.g. 1 to 0)
                    int res_diff = t->parent->resolution - curr_adj[k]->parent->resolution;
                    int adj_tid = get_tile_id(curr_adj[k]);

                    v_adjoin(all_connections->adjacent_tiles[j * 2 + 1][k], adj_tid);

                    // Ignore other types of connections, (only adapter connections)
                    if (res_diff == 1)
                    {
                        fatalif(v_count(curr_adj) != (dimension == 2 ? 2 : 4),
                                "Tile %d alignment is off, adjacent to more than %d (%d) tile of low to high resolution",
                                curr_tid, (dimension == 2 ? 2 : 4), v_count(curr_adj));

                        v_adjoin(need_adapters->adjacent_tiles[j * 2 + 1][k], adj_tid);
                        if (!has_res_diff)
                        {
                            v_adjoin(need_adapters->current_tile, curr_tid);
                            need_adapters->ntiles += 1;
                        }
                        has_res_diff = 1;
                    }
                    else if (res_diff == 0)
                    {
                        fatalif(v_count(curr_adj) != 1, "Tile %d alignment is off, adjacent to more than 1 (%d) tile of same resolution", curr_tid, v_count(curr_adj));
                    }
                    else if (res_diff == -1)
                    {
                        fatalif(v_count(curr_adj) != 1,
                                "Tile %d alignment is off, adjacent to more than 1 (%d) tile of low to high resolution",
                                curr_tid, v_count(curr_adj));
                    }
                    else if (res_diff != 0 && res_diff != -1)
                    {
                        fatalif(1, "Resolution can only jump by power of 2 between tile %d and %d", curr_tid, adj_tid);
                    }
                }
                v_free(curr_adj);
            }
            else
            {
                // If edge, just add no connections to all_connection list
                for (int k = 0; k < 4; k += 1)
                {
                    v_adjoin(all_connections->adjacent_tiles[j * 2 + 1][k], -1);
                }
            }
        }

        // Fill rest of adjacent tiles with -1 if there was a res_diff
        if (has_res_diff)
        {
            for (int j = 0; j < dimension; j += 1)
            {
                for (int k = 0; k < 4; k += 1)
                {
                    if (v_count(need_adapters->adjacent_tiles[j * 2][k]) != need_adapters->ntiles)
                    {
                        v_adjoin(need_adapters->adjacent_tiles[j * 2][k], -1);
                    }
                    if (v_count(need_adapters->adjacent_tiles[j * 2 + 1][k]) != need_adapters->ntiles)
                    {
                        v_adjoin(need_adapters->adjacent_tiles[j * 2 + 1][k], -1);
                    }
                }
            }
        }
    }
}

// Compare two neighbours, order is curr tile id, face number
// Faces are numbered 0(low x), 1(high x), 2(low y), 3(high y), 4(low z), 5 (high z)
int cmp_neighbours(const void *a, const void *b)
{
    neighbour *a_neighbour = (neighbour *)a;
    neighbour *b_neighbour = (neighbour *)b;

    int curr_diff = a_neighbour->curr_tile - b_neighbour->curr_tile;
    if (curr_diff == 0)
    {
        int face_diff = a_neighbour->face - b_neighbour->face;
        if (face_diff == 0)
        {
            int adj_diff = a_neighbour->adjacent_tile - b_neighbour->adjacent_tile;
            return adj_diff;
        }
        else
        {
            return face_diff;
        }
    }
    else
    {
        return curr_diff;
    }
}
// Parses a connectivity struct and returns only adjacent connections between tiles, tiles are sorted using cmp_neighbours
neighbour *get_adjacent_list(struct connectivity *connection)
{
    neighbour *list = 0;
    for (int i = 0; i < connection->ntiles; i += 1)
    {
        for (int j = 0; j < dimension; j += 1)
        {
            for (int k = 0; k < 4; k += 1)
            {

                if (connection->adjacent_tiles[j * 2][k][i] != -1)
                {
                    neighbour entry;
                    entry.curr_tile = connection->current_tile[i];
                    entry.face = j * 2;
                    entry.adjacent_tile = connection->adjacent_tiles[j * 2][k][i];
                    v_adjoin(list, entry);
                }

                if (connection->adjacent_tiles[j * 2 + 1][k][i] != -1)
                {
                    neighbour entry;
                    entry.curr_tile = connection->current_tile[i];
                    entry.face = j * 2 + 1;
                    entry.adjacent_tile = connection->adjacent_tiles[j * 2 + 1][k][i];
                    v_adjoin(list, entry);
                }
            }
        }
    }

    qsort(list, v_count(list), sizeof(neighbour), cmp_neighbours);

    return list;
}

void free_connectivity(struct connectivity *connection)
{
    if (!connection)
    {
        return;
    }
    for (int i = 0; i < dimension; i += 1)
    {
        for (int j = 0; j < 4; j += 1)
        {
            v_free(connection->adjacent_tiles[i * 2][j]);
            v_free(connection->adjacent_tiles[i * 2 + 1][j]);
        }
    }

    v_free(connection->current_tile);

    free(connection);
}

void print_connectivity(struct connectivity *connection, int tabs)
{

    neighbour *connection_list = get_adjacent_list(connection);
    print_tabs(tabs);
    printf("Curr\tFace\tAdj\n");
    for (int i = 0; i < v_count(connection_list); i += 1)
    {
        print_tabs(tabs);
        printf("%d\t%d\t%d\n", connection_list[i].curr_tile, connection_list[i].face, connection_list[i].adjacent_tile);
    }
    v_free(connection_list);
}

void validate_pe_fit()
{
    // Check that ncompute is equal to prism definition
    int expected_compute = 0;
    for (int i = 0; i < sol->nprims; i += 1)
    {
        int curr_pes = 1;
        for (int j = 0; j < dimension; j += 1)
        {
            curr_pes *= sol->prism[i].shape[j];
        }
        expected_compute += curr_pes;
    }

    fatalif(expected_compute != ncompute,
            "Expected %d compute PEs, but got %d", expected_compute, ncompute);

    // Cache to check used PEs
    int *pe_cache = calloc(prob->fabric_shape[0] * prob->fabric_shape[1], sizeof(int));
    fatalif(pe_cache == 0,
            "Unable to allocate PE cache for validation, fabric shape %d X %d",
            prob->fabric_shape[0], prob->fabric_shape[1]);

    // Cache to check used adpaters PEs
    int *adapter_cache = calloc(prob->fabric_shape[0] * prob->fabric_shape[1], sizeof(int));
    fatalif(adapter_cache == 0,
            "Unable to allocate adapter_cache cache for validation, fabric shape %d X %d",
            prob->fabric_shape[0], prob->fabric_shape[1]);

    // Validate fit for compute_map
    for (int i = 0; i < ncompute; i += 1)
    {
        fatalif(sol->compute_map[i].x < 0 ||
                    sol->compute_map[i].x >= prob->fabric_shape[0],
                "Compute map point %d [%d, %d]: Invalid x coordinate, does not fit in fabric [0, %d]",
                i, sol->compute_map[i].x, sol->compute_map[i].y, prob->fabric_shape[0] - 1);
        fatalif(sol->compute_map[i].y < 0 ||
                    sol->compute_map[i].y >= prob->fabric_shape[1],
                "Compute map point %d [%d, %d]: Invalid y coordinate, does not fit in fabric [0, %d]",
                i, sol->compute_map[i].x, sol->compute_map[i].y, prob->fabric_shape[1] - 1);

        int pe_index = sol->compute_map[i].x + sol->compute_map[i].y * prob->fabric_shape[0];
        fatalif(pe_cache[pe_index] > 0,
                "Compute map point %d [%d, %d]: Invalid as compute map point %d already exists here",
                i, sol->compute_map[i].x, sol->compute_map[i].y, pe_cache[pe_index] - 1);
        pe_cache[pe_index] = i + 1;
    }

    // Validate fit for adapter_map
    for (int i = 0; i < nadaptor; i += 1)
    {
        fatalif(sol->adapter_map[i].x < 0 ||
                    sol->adapter_map[i].x >= prob->fabric_shape[0],
                "Adapter map point %d [%d, %d]: Invalid x coordinate, does not fit in fabric [0, %d]",
                i, sol->adapter_map[i].x, sol->adapter_map[i].y, prob->fabric_shape[0] - 1);
        fatalif(sol->adapter_map[i].y < 0 ||
                    sol->adapter_map[i].y >= prob->fabric_shape[1],
                "Adapter map point %d [%d, %d]: Invalid y coordinate, does not fit in fabric [0, %d]",
                i, sol->adapter_map[i].x, sol->adapter_map[i].y, prob->fabric_shape[1] - 1);

        int pe_index = sol->adapter_map[i].x + sol->adapter_map[i].y * prob->fabric_shape[0];
        fatalif(pe_cache[pe_index] > 0,
                "Adapter map point %d [%d, %d]: Invalid as compute map point %d already exists here",
                i, sol->adapter_map[i].x, sol->adapter_map[i].y, pe_cache[pe_index] - 1);
        fatalif(adapter_cache[pe_index] == (dimension == 2 ? 2 : 4),
                "Adapter map point %d [%d, %d]: Invalid as adapter map point already used for maximum adapter connections",
                i, sol->adapter_map[i].x, sol->adapter_map[i].y);
        adapter_cache[pe_index] += 1;
    }

    free(adapter_cache);
    free(pe_cache);
}

void validate_prism_overlap()
{
    // Build R tree of prisms to test overlap
    tile **prism_list = 0;
    tile_tree *prism_tree = build_prism_tree(&prism_list, sol->prism, sol->nprims, 2);

    struct prism **overlaps = 0;

    // Check overlap
    for (int i = 0; i < sol->nprims; i += 1)
    {
        tile **curr_prism = search_tree(prism_tree, prism_list[i]->bounds);

        fatalif(v_count(curr_prism) == 0, "Unable to find prism %d in tree during overlap validation", i);

        // Look for prisms in range search that are not the current prism
        for (int j = 0; j < v_count(curr_prism); j += 1)
        {
            if (curr_prism[j]->parent != prism_list[i]->parent)
            {
                v_adjoin(overlaps, prism_list[i]->parent);
            }
        }
        v_free(curr_prism);
    }

    int count = v_count(overlaps);

    if (count > 0)
    {
        printf("Overlapping prisms exist:\n");
        for (int i = 0; i < count; i += 1)
        {
            print_prism(overlaps[i], 1, 1);
        }
    }

    // Cleanup

    v_free(overlaps);
    free_tree(prism_tree);
    for (int i = 0; i < v_count(prism_list); i += 1)
    {
        free(prism_list[i]);
    }
    v_free(prism_list);

    if (count > 0)
    {
        exit(1); // Exit with failure
    }
}

void validate_coverage()
{
    // Assuming no overlap can just check if prism volume = sample space volume
    int expected_volume = 1;
    for (int i = 0; i < dimension; i += 1)
    {
        expected_volume *= prob->volume_shape[i] * sol->inv_sampling_step;
    }

    int real_volume = 0;

    for (int i = 0; i < sol->nprims; i += 1)
    {
        int prism_volume = 1;
        for (int j = 0; j < dimension; j += 1)
        {
            prism_volume *= (get_close(&(sol->prism[i]), j) -
                             sol->prism[i].origin[j] + 1);
        }
        real_volume += prism_volume;
    }

    fatalif(real_volume != expected_volume, "Coverage failed, expected: %d, real %d", expected_volume, real_volume);
}

int validate_adapters()
{
    int success = 1;

    // Build r tree of tiles
    tile_list = 0;
    tree = build_tile_tree(&tile_list, sol->prism, sol->nprims, 4);

    // Print tree build
    if (opt.flags & DEBUG)
    {
        print_tree(tree);
    }

    // Build need_adapters and all_connections structs
    build_connectivity();

    neighbour *adapter_list = get_adjacent_list(need_adapters);

    if (v_count(adapter_list) != nadaptor)
    {
        printf("Invalid number of adapters, see below for required connections:\n");
        print_connectivity(need_adapters, 1);
        success = 0;
    }

    v_free(adapter_list);

    return success;
}

/* feasibility check:
 * nonoverlap 
 * full coverage 
 * required adapters are 1:2 or 2:1 only */
int validate()
{
    validate_prism_overlap();
    printf("[OK] - Overlap Validated\n");
    validate_coverage();
    printf("[OK] - Coverage Validated\n");
    int adapters_good = validate_adapters();
    if (adapters_good)
    {
        printf("[OK] - Connectivity Validated\n");
        validate_pe_fit();
        printf("[OK] - PE Fit Validated\n");
    }

    return adapters_good;
}

// Get the exact heatmap value in heatmap array at inputted coords (in heatmap space)
double get_heatmap_value(int *coords)
{
    if (dimension == 2)
    {
        return prob->heatmap[coords[0] + coords[1] * (prob->volume_shape[0] + 1)];
    }
    else
    {
        return prob->heatmap[coords[0] + coords[1] * (prob->volume_shape[0] + 1) + coords[2] * (prob->volume_shape[0] + 1) * (prob->volume_shape[1] + 1)];
    }
}

// Interpolate heat map value at inputted coords (in heatmap space)
double interpolate_heatmap_value(double *scaled_coords)
{

    int floored_coords[NDIM];
    int ceil_coords[NDIM];

    for (int i = 0; i < dimension; i += 1)
    {
        floored_coords[i] = (int)floor(scaled_coords[i]);
        ceil_coords[i] = (int)ceil(scaled_coords[i]);
    }

    if (dimension == 2)
    {
        // Bilinear interpolation for 2D
        double x_diff = min_d(scaled_coords[0] - floored_coords[0], ceil_coords[0] - scaled_coords[0]);
        double y_diff = min_d(scaled_coords[1] - floored_coords[1], ceil_coords[1] - scaled_coords[1]);

        int p_11[2] = {floored_coords[0], floored_coords[1]};
        int p_21[2] = {ceil_coords[0], floored_coords[1]};

        int p_12[2] = {floored_coords[0], ceil_coords[1]};
        int p_22[2] = {ceil_coords[0], ceil_coords[1]};

        double f_11 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_11) : (ceil_coords[0] - scaled_coords[0]) * get_heatmap_value(p_11);
        double f_21 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_21) : (scaled_coords[0] - floored_coords[0]) * get_heatmap_value(p_21);

        double f_12 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_12) : (ceil_coords[0] - scaled_coords[0]) * get_heatmap_value(p_12);
        double f_22 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_22) : (scaled_coords[0] - floored_coords[0]) * get_heatmap_value(p_22);

        double f_1 = y_diff < DOUBLE_EPSILON ? (f_11 + f_21) : (ceil_coords[1] - scaled_coords[1]) * (f_11 + f_21);
        double f_2 = y_diff < DOUBLE_EPSILON ? (f_12 + f_22) : (scaled_coords[1] - floored_coords[1]) * (f_12 + f_22);

        double inv_x = x_diff < DOUBLE_EPSILON ? 0.5 : 1.0 / (ceil_coords[0] - floored_coords[0]);
        double inv_y = y_diff < DOUBLE_EPSILON ? 0.5 : 1.0 / (ceil_coords[1] - floored_coords[1]);

        double value = inv_x * inv_y * (f_1 + f_2);

        return value;
    }
    else
    {
        // Trilinear interpolation for 3D

        double x_diff = min_d(scaled_coords[0] - floored_coords[0], ceil_coords[0] - scaled_coords[0]);
        double y_diff = min_d(scaled_coords[1] - floored_coords[1], ceil_coords[1] - scaled_coords[1]);
        double z_diff = min_d(scaled_coords[2] - floored_coords[2], ceil_coords[2] - scaled_coords[2]);

        // Z Low
        int p_111[3] = {floored_coords[0], floored_coords[1], floored_coords[2]};
        int p_211[3] = {ceil_coords[0], floored_coords[1], floored_coords[2]};

        int p_121[3] = {floored_coords[0], ceil_coords[1], floored_coords[2]};
        int p_221[3] = {ceil_coords[0], ceil_coords[1], floored_coords[2]};

        // Z High
        int p_112[3] = {floored_coords[0], floored_coords[1], ceil_coords[2]};
        int p_212[3] = {ceil_coords[0], floored_coords[1], ceil_coords[2]};

        int p_122[3] = {floored_coords[0], ceil_coords[1], ceil_coords[2]};
        int p_222[3] = {ceil_coords[0], ceil_coords[1], ceil_coords[2]};

        // Z Low
        double f_111 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_111) : (ceil_coords[0] - scaled_coords[0]) * get_heatmap_value(p_111);
        double f_211 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_211) : (scaled_coords[0] - floored_coords[0]) * get_heatmap_value(p_211);

        double f_121 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_121) : (ceil_coords[0] - scaled_coords[0]) * get_heatmap_value(p_121);
        double f_221 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_221) : (scaled_coords[0] - floored_coords[0]) * get_heatmap_value(p_221);

        // Z High
        double f_112 = y_diff < DOUBLE_EPSILON ? get_heatmap_value(p_112) : (ceil_coords[0] - scaled_coords[0]) * get_heatmap_value(p_112);
        double f_212 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_212) : (scaled_coords[0] - floored_coords[0]) * get_heatmap_value(p_212);

        double f_122 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_122) : (ceil_coords[0] - scaled_coords[0]) * get_heatmap_value(p_122);
        double f_222 = x_diff < DOUBLE_EPSILON ? get_heatmap_value(p_222) : (scaled_coords[0] - floored_coords[0]) * get_heatmap_value(p_222);

        double f_11 = y_diff < DOUBLE_EPSILON ? (f_111 + f_211) : (ceil_coords[1] - scaled_coords[1]) * (f_111 + f_211);
        double f_21 = y_diff < DOUBLE_EPSILON ? (f_121 + f_221) : (scaled_coords[1] - floored_coords[1]) * (f_121 + f_221);

        double f_12 = y_diff < DOUBLE_EPSILON ? (f_112 + f_212) : (ceil_coords[1] - scaled_coords[1]) * (f_112 + f_212);
        double f_22 = y_diff < DOUBLE_EPSILON ? (f_122 + f_222) : (scaled_coords[1] - floored_coords[1]) * (f_122 + f_222);

        double f_1 = z_diff < DOUBLE_EPSILON ? (f_11 + f_21) : (ceil_coords[2] - scaled_coords[2]) * (f_11 + f_21);
        double f_2 = z_diff < DOUBLE_EPSILON ? (f_12 + f_22) : (scaled_coords[2] - floored_coords[2]) * (f_12 + f_22);

        double inv_x = x_diff < DOUBLE_EPSILON ? 0.5 : 1.0 / (ceil_coords[0] - floored_coords[0]);
        double inv_y = y_diff < DOUBLE_EPSILON ? 0.5 : 1.0 / (ceil_coords[1] - floored_coords[1]);
        double inv_z = z_diff < DOUBLE_EPSILON ? 0.5 : 1.0 / (ceil_coords[2] - floored_coords[2]);

        double value = inv_x * inv_y * inv_z * (f_1 + f_2);

        return value;
    }
}

// Interpolate heat map value at inputted coords (in sampling space)
double scale_interpolate_heatmap_value(int *coords)
{
    double scaled_coords[NDIM];
    for (int i = 0; i < dimension; i += 1)
    {
        scaled_coords[i] = (double)coords[i] / sol->inv_sampling_step;
        scaled_coords[i] = min_d(scaled_coords[i], prob->volume_shape[i]);
        scaled_coords[i] = max_d(scaled_coords[i], 0);
    }

    return interpolate_heatmap_value(scaled_coords);
}

// Integrate using trapezoid rule over heatmap at given y, z over given x range
double integrate_x(double x_low, double x_high, double y, double z)
{
    fatalif(x_low > x_high, "Integrating x wrong way!");

    int range_low = (int)ceil(x_low);
    int range_high = (int)floor(x_high);

    // Occurs when both are in same grid space
    if (range_low > range_high)
    {
        double low[NDIM] = {x_low, y, z};
        double high[NDIM] = {x_high, y, z};

        double low_val = interpolate_heatmap_value(low);
        double high_val = interpolate_heatmap_value(high);

        return 0.5 * (high_val + low_val) * (x_high - x_low);
    }

    double sum = 0.0;

    for (int i = range_low; i < range_high; i += 1)
    {
        double low[NDIM] = {(double)i, y, z};
        double high[NDIM] = {(double)i + 1.0, y, z};

        double low_val = interpolate_heatmap_value(low);
        double high_val = interpolate_heatmap_value(high);

        // Since units of 1, can just integrate directly without delta
        sum += 0.5 * (high_val + low_val);
    }

    double low_diff = (double)range_low - x_low;
    double high_diff = x_high - (double)range_high;

    // Add endpoints if they were not exactly on heatmap
    if (low_diff > DOUBLE_EPSILON)
    {
        double low[NDIM] = {x_low, y, z};
        double high[NDIM] = {(double)range_low, y, z};

        double low_val = interpolate_heatmap_value(low);
        double high_val = interpolate_heatmap_value(high);

        sum += 0.5 * (high_val + low_val) * low_diff;
    }

    if (high_diff > DOUBLE_EPSILON)
    {
        double low[NDIM] = {(double)range_high, y, z};
        double high[NDIM] = {x_high, y, z};

        double low_val = interpolate_heatmap_value(low);
        double high_val = interpolate_heatmap_value(high);

        sum += 0.5 * (high_val + low_val) * high_diff;
    }

    return sum;
}

// Double integrate heatmap using trapezoid rule at given z over given x, y range
double integrate_y(double x_low, double x_high, double y_low, double y_high, double z)
{

    fatalif(y_low > y_high, "Integrating y wrong way!");

    int range_low = (int)ceil(y_low);
    int range_high = (int)floor(y_high);

    // Occurs when both are in same grid space
    if (range_low > range_high)
    {
        double low_val = integrate_x(x_low, x_high, y_low, z);
        double high_val = integrate_x(x_low, x_high, y_high, z);

        return 0.5 * (high_val + low_val) * (y_high - y_low);
    }

    double sum = 0.0;

    for (int i = range_low; i < range_high; i += 1)
    {
        double low_val = integrate_x(x_low, x_high, i, z);
        double high_val = integrate_x(x_low, x_high, i + 1, z);

        // Trapezoidal rule (unit step)
        sum += 0.5 * (high_val + low_val);
    }

    double low_diff = (double)range_low - y_low;
    double high_diff = y_high - (double)range_high;

    // Add endpoints if they were not exactly on heatmap
    if (low_diff > DOUBLE_EPSILON)
    {

        double low_val = integrate_x(x_low, x_high, y_low, z);
        double high_val = integrate_x(x_low, x_high, range_low, z);

        sum += 0.5 * (high_val + low_val) * low_diff;
    }

    if (high_diff > DOUBLE_EPSILON)
    {

        double low_val = integrate_x(x_low, x_high, range_high, z);
        double high_val = integrate_x(x_low, x_high, y_high, z);

        sum += 0.5 * (high_val + low_val) * high_diff;
    }

    return sum;
}

// Triple integrate heatmap using trapezoid rule over given x, y, z range
double integrate_z(double x_low, double x_high, double y_low, double y_high, double z_low, double z_high)
{

    fatalif(z_low > z_high, "Integrating z wrong way!");

    int range_low = (int)ceil(z_low);
    int range_high = (int)floor(z_high);

    // Occurs when both are in same grid space
    if (range_low > range_high)
    {
        double low_val = integrate_y(x_low, x_high, y_low, y_high, z_low);
        double high_val = integrate_y(x_low, x_high, y_low, y_high, z_high);

        return 0.5 * (high_val + low_val) * (z_high - z_low);
    }

    double sum = 0.0;

    for (int i = range_low; i < range_high; i += 1)
    {
        double low_val = integrate_y(x_low, x_high, y_low, y_high, i);
        double high_val = integrate_y(x_low, x_high, y_low, y_high, i + 1);

        // Trapezoidal rule (unit step)
        sum += 0.5 * (high_val + low_val);
    }

    double low_diff = (double)range_low - z_low;
    double high_diff = z_high - (double)range_high;

    // Add endpoints if they were not exactly on heatmap
    if (low_diff > DOUBLE_EPSILON)
    {

        double low_val = integrate_y(x_low, x_high, y_low, y_high, z_low);
        double high_val = integrate_y(x_low, x_high, y_low, y_high, range_low);

        sum += 0.5 * (high_val + low_val) * low_diff;
    }

    if (high_diff > DOUBLE_EPSILON)
    {

        double low_val = integrate_y(x_low, x_high, y_low, y_high, range_high);
        double high_val = integrate_y(x_low, x_high, y_low, y_high, z_high);

        sum += 0.5 * (high_val + low_val) * high_diff;
    }

    return sum;
}

// Find max overshoot between target heatmap and solution resolution at given y, z over given x range
double find_max_overshoot_x(double x_low, double x_high, double y, double z, double sol_val)
{
    fatalif(x_low > x_high, "x wrong way!");

    int range_low = (int)ceil(x_low);
    int range_high = (int)floor(x_high);

    // Occurs when both are in same grid space
    if (range_low > range_high)
    {
        double low[NDIM] = {x_low, y, z};
        double high[NDIM] = {x_high, y, z};

        double low_val = interpolate_heatmap_value(low);
        double high_val = interpolate_heatmap_value(high);

        return max_d(0, max_d(low_val - sol_val, high_val - sol_val));
    }

    double max_overshoot = 0;

    for (int i = range_low; i <= range_high; i += 1)
    {
        double curr_coord[NDIM] = {(double)i, y, z};
        double val = interpolate_heatmap_value(curr_coord);

        max_overshoot = max_d(max_overshoot, val - sol_val);
    }

    double low_diff = (double)range_low - x_low;
    double high_diff = x_high - (double)range_high;

    // Add endpoints if they were not exactly on heatmap
    if (low_diff > DOUBLE_EPSILON)
    {
        double curr_coord[NDIM] = {(double)x_low, y, z};

        double val = interpolate_heatmap_value(curr_coord);

        max_overshoot = max_d(max_overshoot, val - sol_val);
    }

    if (high_diff > DOUBLE_EPSILON)
    {
        double curr_coord[NDIM] = {(double)x_high, y, z};

        double val = interpolate_heatmap_value(curr_coord);

        max_overshoot = max_d(max_overshoot, val - sol_val);
    }

    return max_overshoot;
}

// Find max overshoot between target heatmap and solution resolution at given z over given x, y range
double find_max_overshoot_y(double x_low, double x_high, double y_low, double y_high, double z, double sol_val)
{

    fatalif(y_low > y_high, "y wrong way!");

    int range_low = (int)ceil(y_low);
    int range_high = (int)floor(y_high);

    // Occurs when both are in same grid space
    if (range_low > range_high)
    {
        double low_val = find_max_overshoot_x(x_low, x_high, y_low, z, sol_val);
        double high_val = find_max_overshoot_x(x_low, x_high, y_high, z, sol_val);

        return max_d(0, max_d(low_val, high_val));
    }

    double max_overshoot = 0.0;

    for (int i = range_low; i <= range_high; i += 1)
    {
        double val = find_max_overshoot_x(x_low, x_high, i, z, sol_val);

        max_overshoot = max_d(max_overshoot, val);
    }

    double low_diff = (double)range_low - y_low;
    double high_diff = y_high - (double)range_high;

    // Add endpoints if they were not exactly on heatmap
    if (low_diff > DOUBLE_EPSILON)
    {
        double val = find_max_overshoot_x(x_low, x_high, y_low, z, sol_val);

        max_overshoot = max_d(max_overshoot, val);
    }

    if (high_diff > DOUBLE_EPSILON)
    {
        double value = find_max_overshoot_x(x_low, x_high, y_high, z, sol_val);

        max_overshoot = max_d(max_overshoot, value);
    }

    return max_overshoot;
}

// Find max overshoot between target heatmap and solution resolution over given x, y, z range
double find_max_overshoot_z(double x_low, double x_high, double y_low, double y_high, double z_low, double z_high, double sol_val)
{

    fatalif(z_low > z_high, "z wrong way!");

    int range_low = (int)ceil(z_low);
    int range_high = (int)floor(z_high);

    // Occurs when both are in same grid space
    if (range_low > range_high)
    {
        double low_val = find_max_overshoot_y(x_low, x_high, y_low, y_high, z_low, sol_val);
        double high_val = find_max_overshoot_y(x_low, x_high, y_low, y_high, z_high, sol_val);

        return max_d(0, max_d(low_val, high_val));
    }

    double max_overshoot = 0.0;

    for (int i = range_low; i < range_high; i += 1)
    {
        double val = find_max_overshoot_y(x_low, x_high, y_low, y_high, i, sol_val);

        max_overshoot = max_d(max_overshoot, val);
    }

    double low_diff = (double)range_low - z_low;
    double high_diff = z_high - (double)range_high;

    // Add endpoints if they were not exactly on heatmap
    if (low_diff > DOUBLE_EPSILON)
    {
        double val = find_max_overshoot_y(x_low, x_high, y_low, y_high, z_low, sol_val);

        max_overshoot = max_d(max_overshoot, val);
    }

    if (high_diff > DOUBLE_EPSILON)
    {
        double val = find_max_overshoot_y(x_low, x_high, y_low, y_high, z_high, sol_val);

        max_overshoot = max_d(max_overshoot, val);
    }

    return max_overshoot;
}

// Convert tile node bounds from sampling space to heatmap space
node_bounds_d scale_tile_to_heatmap(struct tile *t)
{
    node_bounds_d scaled;

    for (int i = 0; i < dimension; i += 1)
    {
        double low = (double)t->bounds.low_bounds[i] / sol->inv_sampling_step;
        low = min_d(low, prob->volume_shape[i]);
        low = max_d(low, 0.0);

        double high = ((double)t->bounds.up_bounds[i] + 1.0) / sol->inv_sampling_step;
        high = min_d(high, prob->volume_shape[i]);
        high = max_d(high, 0);

        scaled.low_bounds[i] = low;
        scaled.up_bounds[i] = high;
    }

    return scaled;
}

// Finds the max overshoot related to a single tile
double find_max_overshoot(struct tile *t)
{
    // Get scaled volume of tile
    node_bounds_d scaled = scale_tile_to_heatmap(t);

    double max_overshoot = 0.0;

    // Compute resolution of tile based on parent prism
    double solution_res = 1.0 / pow(2, t->parent->resolution);

    if (dimension == 2)
    {
        max_overshoot = find_max_overshoot_y(scaled.low_bounds[0], scaled.up_bounds[0], scaled.low_bounds[1], scaled.up_bounds[1], 0, solution_res);
    }
    else
    {
        max_overshoot = find_max_overshoot_z(scaled.low_bounds[0], scaled.up_bounds[0], scaled.low_bounds[1], scaled.up_bounds[1], scaled.low_bounds[2], scaled.up_bounds[2], solution_res);
    }

    return max_overshoot;
}

// Finds the maximum overshoot and use this to compute the scaling factor to scale each tile by
double find_max_scale()
{
    int ntiles = v_count(tile_list);
    double max_overshoot = 0;
    double max_scale = 1.0;

    for (int i = 0; i < ntiles; i += 1)
    {
        double curr_overshoot = find_max_overshoot(tile_list[i]);
        if (curr_overshoot > max_overshoot)
        {
            // curr_overshoot = target_res - sol_res => scale = curr_overshoot/sol_res + 1
            double sol_res = 1.0 / pow(2, tile_list[i]->parent->resolution);
            max_scale = 1.0 + curr_overshoot / sol_res;
            max_overshoot = curr_overshoot;
        }
    }

    return max_scale;
}

// Integrate heatmap and tile in tile volume and returns the scaled ratio
double integrate_tile(struct tile *t, double scale)
{

    node_bounds_d scaled = scale_tile_to_heatmap(t);

    double i1 = 1.0 / pow(2, t->parent->resolution);

    for (int i = 0; i < dimension; i += 1)
    {
        i1 *= (scaled.up_bounds[i] - scaled.low_bounds[i]);
    }

    double i2 = 0.0;

    if (dimension == 2)
    {
        i2 = integrate_y(scaled.low_bounds[0], scaled.up_bounds[0], scaled.low_bounds[1], scaled.up_bounds[1], 0);
    }
    else
    {
        i2 = integrate_z(scaled.low_bounds[0], scaled.up_bounds[0], scaled.low_bounds[1], scaled.up_bounds[1], scaled.low_bounds[2], scaled.up_bounds[2]);
    }

    fatalif(i1 == 0, "Tile volume zero!");
    if (opt.flags & RESOLUTION)
    {
        printf("%f %f\n", i1, i2);
    }

    return i2 / (i1 * scale);
}

// Scores a tile accuracy with a given scale
double score_tile_accuracy(struct tile *t, double scale)
{
    return integrate_tile(t, scale);
}

// Score accuracy of solution
double score_accuracy()
{
    double scale = find_max_scale();

    fatalif(scale < 1.0, "Scale cannot be less than 1");

    if (opt.flags & RESOLUTION)
    {
        printf("Resolution Difference:\n");
    }

    double score = 0.0;

    // Iterate through all tiles and add their scores to cumulative
    int ntiles = v_count(tile_list);
    for (int i = 0; i < ntiles; i += 1)
    {
        double tile_score = score_tile_accuracy(tile_list[i], scale);
        score += tile_score;
    }

    if (opt.flags & RESOLUTION)
    {
        printf("End Resolution Difference\n");
    }
    return score;
}

// Score accuracy of wires
double score_wires()
{
    // Uses the global connectivity struct
    fatalif(all_connections == 0, "Did not create connectivity struct");

    neighbour *connection_list = get_adjacent_list(all_connections);

    int nconnection = v_count(connection_list);
    double sum_1_5_norm = 0;
    double sum_wires = 0;
    // Sum all non-adapter connections first
    for (int i = 0; i < nconnection; i += 1)
    {
        int cid = connection_list[i].curr_tile;
        int aid = connection_list[i].adjacent_tile;

        int curr_res = tile_list[cid]->parent->resolution;
        int adj_res = tile_list[cid]->parent->resolution;

        if (curr_res == adj_res)
        {
            // Only do direct wire calculation for connections without res diff (i.e. adapters)
            double wire_x = sol->compute_map[cid].x - sol->compute_map[aid].x;
            double wire_y = sol->compute_map[cid].y - sol->compute_map[aid].y;

            sum_wires += abs_d(wire_x) + abs_d(wire_y);
            sum_1_5_norm += pow(abs_d(wire_x) + abs_d(wire_y), WIRE_POW);
        }
    }

    v_free(connection_list);

    // Sum all adapter connections
    connection_list = get_adjacent_list(need_adapters);
    nconnection = v_count(connection_list);
    int curr_cid = -1;
    bool curr_added = 0;
    for (int i = 0; i < nconnection; i += 1)
    {
        int cid = connection_list[i].curr_tile;
        int aid = connection_list[i].adjacent_tile;

        // this works since connection list is sorted by tile id
        if (curr_cid != cid)
        {
            curr_cid = cid;
            curr_added = 0;
        }

        // Add score of current to adapter if not already added
        if (!curr_added)
        {
            double c_x = sol->compute_map[cid].x - sol->adapter_map[i].x;
            double c_y = sol->compute_map[cid].y - sol->adapter_map[i].y;

            // Multiply by 2 as bidirectional
            sum_wires += 2 * (abs_d(c_x) + abs_d(c_y));
            sum_1_5_norm += 2 * pow(abs_d(c_x) + abs_d(c_y), WIRE_POW);

            curr_added = 1;
        }

        // Add score of adapter to adjacent
        double a_x = sol->compute_map[aid].x - sol->adapter_map[i].x;
        double a_y = sol->compute_map[aid].y - sol->adapter_map[i].y;

        // Multiply by 2 as bidrectional
        sum_wires += 2 * (abs_d(a_x) + abs_d(a_y));
        sum_1_5_norm += 2 * pow(abs_d(a_x) + abs_d(a_y), WIRE_POW);
    }
    if (opt.flags & SCORE)
    {
        printf("Total Wire Length: %f\n", sum_wires);
    }

    return pow(sum_1_5_norm, 1.0 / WIRE_POW);
}

// Compute the score of a solution
double score_solution()
{
    double accuracy = score_accuracy();
    double wires = score_wires();

    int max_tiles = prob->fabric_shape[0] * prob->fabric_shape[1];
    int num_tiles = ncompute + nadaptor;

    double max_accuracy = max_tiles; // TODO figure out this metric

    double max_wires_tot = pow(4 * max_tiles, 1.0 / WIRE_POW);
    double max_wires_curr = pow(4 * num_tiles, 1.0 / WIRE_POW);

    if (dimension == 3)
    {
        max_wires_tot = pow(4 * max_tiles + 2 * pow(pow(max_tiles, 1.0 / 3.0), WIRE_POW), 1.0 / WIRE_POW);
        max_wires_curr = pow(4 * num_tiles + 2 * pow(pow(num_tiles, 1.0 / 3.0), WIRE_POW), 1.0 / WIRE_POW);
    }
    double accuracy_norm = accuracy / max_accuracy;
    double wires_norm = min_d(max_wires_tot / wires, max_wires_curr / wires);

    if (opt.flags & SCORE)
    {
        printf("Pure Accuracy: %f / %f\n", accuracy, max_accuracy);
        printf("Pure Connectivity: %f / (%f or %f)\n", wires, max_wires_tot, max_wires_curr);
        printf("Normalized Accuracy: %f\n", accuracy_norm);
        printf("Normalized Connectivity: %f\n", wires_norm);
    }

    return prob->cost.alpha * accuracy_norm + prob->cost.beta * wires_norm;
}

int main(int argc, char **argv)
{
    // The answer input is ignored since we do not need it to validate the solution
    fatalif(argc < 2, "Invalid number of arguments");
    opt.input = argv[1];
    opt.feedbackdir = "";

    opt.flags = 0;
    // Check for flags
    if (argc > 2)
    {
        for (int i = 2; i < argc; i += 1)
        {
            if (strcmp(argv[i], "-d") == 0)
            {
                opt.flags |= DEBUG;
            }
            if (strcmp(argv[i], "-c") == 0)
            {
                opt.flags |= CONNECTIVITY;
            }
            if (strcmp(argv[i], "-r") == 0)
            {
                opt.flags |= RESOLUTION;
            }
            if (strcmp(argv[i], "-s") == 0)
            {
                opt.flags |= SCORE;
            }
            if (strcmp(argv[i], "-o") == 0)
            {
                i += 1;
                opt.feedbackdir = argv[i];
            }
        }
    }

    char *line = 0;
    size_t len = 0;
    ssize_t read;
    int res = 0;

    // Read in problem from file
    prob = calloc(1, sizeof(struct problem));

    FILE *input = fopen(opt.input, "r");

    fatalif(!input, "%s: not found", opt.input);

    if ((read = get_next_line(&line, &len, input)) != -1)
    {
        res = sscanf(line, "%d", &dimension);
        fatalif(res != 1 || dimension > 3 || dimension < 2,
                "%s: Invalid input for dimension, must be 2 or 3", line);
    }
    else
    {
        fatalif(1, "Dimension missing");
    }

    if ((read = get_next_line(&line, &len, input)) != -1)
    {
        if (dimension == 2)
        {
            res = sscanf(line, "%d%d", &prob->volume_shape[0], &prob->volume_shape[1]);
            fatalif(res != 2 ||
                        prob->volume_shape[0] <= 0 ||
                        prob->volume_shape[1] <= 0,
                    "%s: Invalid input for volume shape, must be 2 positive integers",
                    line);
        }
        else if (dimension == 3)
        {
            res = sscanf(line, "%d%d%d", &prob->volume_shape[0], &prob->volume_shape[1], &prob->volume_shape[2]);
            fatalif(res != 3 ||
                        prob->volume_shape[0] <= 0 ||
                        prob->volume_shape[1] <= 0 ||
                        prob->volume_shape[2] <= 0,
                    "%s: Invalid input for volume shape, must be 3 positive integers",
                    line);
        }
    }
    else
    {
        fatalif(1, "Volume Shape missing");
    }

    // Gridpoints per tile edge is always set to 10

    prob->gridpoints_per_tile_edge = 10;

    if ((read = get_next_line(&line, &len, input)) != -1)
    {
        res = sscanf(line, "%d%d",
                     &prob->fabric_shape[0], &prob->fabric_shape[1]);
        fatalif(res != 2 ||
                    prob->fabric_shape[0] <= 0 ||
                    prob->fabric_shape[1] <= 0,
                "%s: Invalid input for fabric shape, must be 2 positive integers",
                line);
    }
    else
    {
        fatalif(1, "Fabric Shape missing");
    }

    if ((read = get_next_line(&line, &len, input)) != -1)
    {
        res = sscanf(line, "%f%f",
                     &prob->cost.alpha, &prob->cost.beta);
        fatalif(res != 2,
                "%s: Invalid input for cost paramters, must be 2 numbers",
                line);
    }
    else
    {
        fatalif(1, "Cost parameters missing");
    }

    int heatmap_num = 0;
    if (dimension == 2)
    {
        heatmap_num = (prob->volume_shape[0] + 1) * (prob->volume_shape[1] + 1);
    }
    else
    {
        heatmap_num = (prob->volume_shape[0] + 1) * (prob->volume_shape[1] + 1) * (prob->volume_shape[2] + 1);
    }

    prob->heatmap = (float *)calloc(heatmap_num, sizeof(float));

    for (int i = 0; i < heatmap_num; i += 1)
    {
        if ((read = get_next_line(&line, &len, input)) != -1)
        {
            res = sscanf(line, "%f",
                         &prob->heatmap[i]);
            fatalif(res != 1 ||
                        prob->heatmap[i] < 0.0f ||
                        prob->heatmap[i] > 1.0f,
                    "%s: Invalid input for heatmap, must be number between [0,1]",
                    line);
        }
        else
        {
            fatalif(1, "Expected %d heatmap entries. Got %d", heatmap_num, i);
        }
    }

    // Read in solution from stdin
    sol = calloc(1, sizeof(struct solution));
    if ((read = get_next_line(&line, &len, stdin) != -1))
    {
        res = sscanf(line, "%d",
                     &sol->inv_sampling_step);
        fatalif(res != 1 ||
                    sol->inv_sampling_step <= 0,
                "%s: Invalid input for inv_sampling_step, must be positive number",
                line);
    }
    else
    {
        fatalif(1, "Sampling Step required");
    }

    for (int i = 0; i < dimension; i += 1)
    {
        volume_bounds.low_bounds[i] = 0;
        volume_bounds.up_bounds[i] = prob->volume_shape[i] * sol->inv_sampling_step - 1;
    }

    struct prism *prisms = 0;

    // Prism parsing

    while ((read = get_next_line(&line, &len, stdin)) != -1 && !check_compute_map(line))
    {
        struct prism p;
        if (dimension == 2)
        {
            res = sscanf(line, "%d%d%d%d%d",
                         &p.resolution,
                         &p.origin[0], &p.origin[1],
                         &p.shape[0], &p.shape[1]);
            fatalif(res != 5, "%s: Invalid paramters for prism, need 5 integers", line);
        }
        else
        {
            res = sscanf(line, "%d%d%d%d%d%d%d",
                         &p.resolution,
                         &p.origin[0], &p.origin[1], &p.origin[2],
                         &p.shape[0], &p.shape[1], &p.shape[2]);
            fatalif(res != 7, "%s: Invalid paramters for prism, need 7 integers", line);
        }

        fatalif(p.resolution < 0, "%s: Invalid resolution, should be non negative integer", line);
        for (int i = 0; i < dimension; i += 1)
        {
            fatalif(p.origin[i] < volume_bounds.low_bounds[i] ||
                        p.origin[i] > volume_bounds.up_bounds[i],
                    "%s: Invalid origin coordinate %d in dimension %d, should be within volume bounds [%d,%d]",
                    line, p.origin[i],
                    i, volume_bounds.low_bounds[i], volume_bounds.up_bounds[i]);

            fatalif(p.shape[i] <= 0, "%s: Invalid shape dimension %d, should be positive integer", line, i);

            int close = get_close(&p, i);
            fatalif(close > volume_bounds.up_bounds[i] || close < volume_bounds.low_bounds[i],
                    "%s: Shape length %d too large in dimension %d, should be within volume bounds [%d,%d]",
                    line, close,
                    i, volume_bounds.low_bounds[i], volume_bounds.up_bounds[i]);
        }

        v_adjoin(prisms, p);
    }

    sol->prism = prisms;
    sol->nprims = v_count(prisms);

    fatalif(!check_compute_map(line), "'compute_map:' header missing");

    struct pt *compute_map = 0;

    while ((read = get_next_line(&line, &len, stdin)) != -1 && !check_adapter_map(line))
    {
        struct pt compute_point;
        res = sscanf(line, "%d%d", &compute_point.x, &compute_point.y);
        fatalif(res != 2 ||
                    compute_point.x < 0 ||
                    compute_point.y < 0,
                "%s: Invalid compute map point, should be non negative integer", line);

        v_adjoin(compute_map, compute_point);
    }

    sol->compute_map = compute_map;

    fatalif(!check_adapter_map(line), "'adapter_map:' header missing");

    struct pt *adapter_map = 0;

    while ((read = get_next_line(&line, &len, stdin)) != -1)
    {
        struct pt adapter_point;
        res = sscanf(line, "%d%d", &adapter_point.x, &adapter_point.y);
        fatalif(res != 2 ||
                    adapter_point.x < 0 ||
                    adapter_point.y < 0,
                "%s: Invalid adapter map point, should be non negative integer", line);

        v_adjoin(adapter_map, adapter_point);
    }

    sol->adapter_map = adapter_map;

    ncompute = v_count(compute_map);
    nadaptor = v_count(adapter_map);

    // If debug flag is specified, print parsed problem and solution
    if (opt.flags & DEBUG)
    {
        print_problem(prob);
        print_solution(sol, ncompute, nadaptor);
    }

    // Only score solution or extra info if validated
    if (validate())
    {
        // Score
        double score = score_solution();

        // Print score
        int feedback_len = strlen(opt.feedbackdir);
        char *score_file = calloc(feedback_len + 11, sizeof(char));

        strcpy(score_file, opt.feedbackdir);
        strcat(score_file, "score.txt");

        FILE *score_output = fopen(score_file, "w");
        fprintf(score_output, "%f", score);
        fclose(score_output);
        free(score_file);
        score_file = NULL;

        // Output connectivity if desired
        if (opt.flags & CONNECTIVITY)
        {
            printf("Tile Connection List:\n");
            print_connectivity(all_connections, 0);
            printf("End Connection List\n");
        }
    }

    // Clean up
    free_connectivity(all_connections);
    free_connectivity(need_adapters);

    free_tree(tree);

    for (int i = 0; i < v_count(tile_list); i += 1)
    {
        free(tile_list[i]);
    }
    v_free(tile_list);

    if (prob->heatmap)
    {
        free(prob->heatmap);
    }
    free(prob);

    v_free(sol->prism);
    v_free(sol->compute_map);
    v_free(sol->adapter_map);

    free(sol);

    if (line)
    {
        free(line);
        line = NULL;
    }
    fclose(input);
    return 0;
}
