#define NDIM 3 /* 2 or 3 */
#define NUM_PRISM_FACES (NDIM == 3 ? 6 : 4)

#include <stdbool.h>

/* [ ] problem specific parameters, the defining cost landscape */
struct costparams
{
    float alpha;
    float beta;
};

/* problem definition */
struct problem
{
    /* heatmap space parameters */
    int volume_shape[NDIM]; /* heatmap space shape */
    float *heatmap;         /* heatmap */

    /* solution space parameters */
    int gridpoints_per_tile_edge; /* edge of the cube mapped to a tile */
    int fabric_shape[2];          /* fabric shape */

    struct costparams cost; /* coeffs that parameterize the cost */
};

/* coordinate of a tile on the fabric */
struct pt
{
    int x; /* x coordinate of point */
    int y; /* y coordinate of point */
};

/* partition prism parameters */
struct prism
{
    int resolution;   /* log2 of spacing between adjecent gridpoints */
    bool is_empty;     /* is this prism considered empty */
    int origin[NDIM]; /* prism's origin designated using the finest resolution coordinates  (sampling space)    */
    int shape[NDIM];  /* prism's dimensions in units of tiles */
};

/* partition of the heatmap space */
struct solution
{
    int inv_sampling_step;  /* number of gridpoints per volume unit */
    int nprims;             /* number of prisms    */
    struct prism *prism;    /* array of prisims    */
    struct pt *compute_map; /* compute tile array  */
    struct pt *adapter_map; /* adapter tile array  */
};

/* undirectional connectivity [derived] */
struct connectivity
{
    int ntiles;                              /* total number of tiles */
    int *adjacent_tiles[NUM_PRISM_FACES][4]; /* Id of adjacent tiles, 6 faces, max of 4 neighbors per face */
    int *current_tile;                       /* Current tile id */
};

/* option flags */
typedef enum flag_type
{
    DEBUG = 1,        // Debug flag
    CONNECTIVITY = 2, // Print connectivity flag
    RESOLUTION = 4,   // Print resolution flag
    SCORE = 8,        // Print score info
} flag_type;

/* Options taken from command line */
struct option
{
    char *input;       // Location of problem file
    char *feedbackdir; // Location of feedback directory
    int flags;         // Flags [-d | -c | -r | -s]
};

/* Defines the bounds of a space in NDIM dimensions */
typedef struct node_bounds
{
    int low_bounds[NDIM]; // Lower bounds, always less than or equal to upper bounds
    int up_bounds[NDIM];  // Upper bounds, always greater than or equal to lower bounds
} node_bounds;

/* Defines bounds of a space but using doubles instead of integers */
typedef struct node_bounds_d
{
    double low_bounds[NDIM]; // Lower bounds, always less than or equal to upper bounds
    double up_bounds[NDIM];  // Upper bounds, always greater than or equal to lower bounds
} node_bounds_d;

/* Represents a single tile(PE) in the solution */
typedef struct tile
{
    int prev_tiles;        // Tiles that come before it in order of tile id
    struct prism *parent;  // Prism that tile is a part of
    int tile_id[NDIM + 1]; // Stores values used to compute tile id
    node_bounds bounds;    // Bounds of tile
} tile;

/* Data structure to hold information regarding a tile and its adjacent tiles */
typedef struct neighbour
{
    int curr_tile;     // Current tile id
    int face;          // Which face of current tile this refers to
    int adjacent_tile; // Adjacent tile id
} neighbour;

/* Used to define nodes in R tree */
typedef enum node_type
{
    INTERIOR = 0, // Only containt tile_nodes
    LEAF,         // Leaf nodes contain tile pointers
} node_type;

/* Node of a R tree */
typedef struct tile_node
{
    node_type type;           // Type of node
    void *children;           // Pointer to children, (tile_node** for INTERIOR nodes, tile* for LEAF nodes)
    struct tile_node *parent; // Parent node in tree
    int num_children;         // Number of immediate children
    node_bounds bounds;       // Bounding box surronding bounds of each child
} tile_node;

/* Hilbert variation of R Trees: https://en.wikipedia.org/wiki/Hilbert_R-tree */
typedef struct tile_tree
{
    tile_node *head;  // Head of tile tree
    int max_children; // Max children for leaf nodes
} tile_tree;
