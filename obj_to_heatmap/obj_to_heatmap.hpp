#include <vector>
#include <string>
#include <iostream>

struct Vec3
{
    double x;
    double y;
    double z;

    Vec3 &operator+=(const Vec3 &rhs);
    Vec3 &operator-=(const Vec3 &rhs);

    Vec3 &operator*=(double scalar);
    Vec3 &operator/=(double scalar);

    friend Vec3 operator+(Vec3 lhs, const Vec3 &rhs);
    friend Vec3 operator-(Vec3 lhs, const Vec3 &rhs);

    friend Vec3 operator/(Vec3 lhs, double rhs);
    friend Vec3 operator*(double lhs, Vec3 rhs);

    friend std::ostream &operator<<(std::ostream &out, const Vec3 &vec);
};

struct AABB
{
    Vec3 min;
    Vec3 max;

    double distance2To(const Vec3 &point);

    AABB &operator+=(const AABB &rhs);
    friend AABB operator+(AABB lhs, const AABB &rhs);
};

struct Face
{
    int v0;
    int v1;
    int v2;
    bool has_norms;
};

class Triangle
{
    Vec3 v_0_1;
    Vec3 v_0_2;

public:
    Vec3 verts[3];
    Vec3 norms[3];
    AABB bounds;

    Triangle(Vec3 *verts, Vec3 *norms = nullptr);
    double distance2To(const Vec3 &point);
    bool inFront(const Vec3 &point);

    friend std::ostream &operator<<(std::ostream &out, const Triangle &mesh);
};

struct Node
{
    int height;
    std::vector<Node *> children;
    Triangle *data;
    AABB bounds;
    Node();
    ~Node();
    double distance2To(const Vec3 &point, double best);

    friend std::ostream &operator<<(std::ostream &out, const Node &node);
};

struct Tree
{
    Node *head;
    Tree();
    ~Tree();
    void build(const std::vector<Triangle*> &data, int max_children);
    double distance2To(const Vec3 &point);

    friend std::ostream &operator<<(std::ostream &out, const Tree &tree);
};

struct Mesh
{
    std::vector<Face> faces;
    std::vector<Vec3> vertices;
    std::vector<Vec3> vert_norms;
    std::vector<int> facecount;

    std::vector<Vec3> norms;

    std::vector<Triangle*> triangles;

    Mesh(const std::string &filename);
    ~Mesh();

    friend std::ostream &operator<<(std::ostream &out, const Mesh &mesh);
};
