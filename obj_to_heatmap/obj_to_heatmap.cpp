#include "obj_to_heatmap.hpp"
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>

Vec3 &Vec3::operator+=(const Vec3 &rhs)
{
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
}

Vec3 &Vec3::operator-=(const Vec3 &rhs)
{
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
}

Vec3 &Vec3::operator*=(double scalar)
{
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}

Vec3 &Vec3::operator/=(double scalar)
{
    x /= scalar;
    y /= scalar;
    z /= scalar;
    return *this;
}

Vec3 operator+(Vec3 lhs, const Vec3 &rhs)
{
    lhs += rhs;
    return lhs;
}

Vec3 operator-(Vec3 lhs, const Vec3 &rhs)
{
    lhs -= rhs;
    return lhs;
}

Vec3 operator/(Vec3 lhs, double rhs)
{
    lhs /= rhs;
    return lhs;
}

Vec3 operator*(double lhs, Vec3 rhs)
{
    rhs *= lhs;
    return rhs;
}

std::ostream &operator<<(std::ostream &out, const Vec3 &vec)
{
    return out << "[" << vec.x << " " << vec.y << " " << vec.y << "]";
}

Vec3 normalize(const Vec3 &v)
{
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return Vec3{v.x / len, v.y / len, v.z / len};
}

Vec3 cross(const Vec3 &lhs, const Vec3 &rhs)
{
    return Vec3{lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x};
}

double dot(const Vec3 &lhs, const Vec3 &rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Vec3 min(const Vec3 &lhs, const Vec3 &rhs)
{
    return Vec3{std::min(lhs.x, rhs.x), std::min(lhs.y, rhs.y), std::min(lhs.z, rhs.z)};
}

Vec3 max(const Vec3 &lhs, const Vec3 &rhs)
{
    return Vec3{std::max(lhs.x, rhs.x), std::max(lhs.y, rhs.y), std::max(lhs.z, rhs.z)};
}

double length2(const Vec3 &vec)
{
    return dot(vec, vec);
}

AABB &AABB::operator+=(const AABB &rhs)
{
    this->min = ::min(this->min, rhs.min);
    this->max = ::max(this->max, rhs.max);
    return *this;
}

AABB operator+(AABB lhs, const AABB &rhs)
{
    lhs += rhs;
    return lhs;
}

double distance2ToRectangle(const Vec3 &origin, const Vec3 &dim, const Vec3 &point)
{
    // Adapted from https://www.geometrictools.com/Documentation/DistanceLine3Rectangle3.pdf

    double e_0, e_1;
    double x_0, x_1, x_2;

    if (dim.x < 0.00001)
    {
        // rectangle is y, z
        e_0 = dim.y;
        e_1 = dim.z;

        x_0 = point.y - origin.y;
        x_1 = point.z - origin.z;
        x_2 = point.x - origin.x;
    }
    else if (dim.y < 0.00001)
    {
        // rectangle is x, z
        e_0 = dim.x;
        e_1 = dim.z;

        x_0 = point.x - origin.x;
        x_1 = point.z - origin.z;
        x_2 = point.y - origin.y;
    }
    else
    {
        // rectangle is x, y
        e_0 = dim.x;
        e_1 = dim.y;

        x_0 = point.x - origin.x;
        x_1 = point.y - origin.y;
        x_2 = point.z - origin.z;
    }

    double dist;

    if (x_1 > e_1)
    {
        if (x_0 < -e_0)
        {
            dist = (x_0 + e_0) * (x_0 + e_0) + (x_1 - e_1) * (x_1 - e_1) + x_2 * x_2;
        }
        else if (x_0 > e_0)
        {
            dist = (x_0 - e_0) * (x_0 - e_0) + (x_1 - e_1) * (x_1 - e_1) + x_2 * x_2;
        }
        else
        {
            dist = (x_1 - e_1) * (x_1 - e_1) + x_2 * x_2;
        }
    }
    else if (x_1 < -e_1)
    {
        if (x_0 < -e_0)
        {
            dist = (x_0 + e_0) * (x_0 + e_0) + (x_1 + e_1) * (x_1 + e_1) + x_2 * x_2;
        }
        else if (x_0 > e_0)
        {
            dist = (x_0 - e_0) * (x_0 - e_0) + (x_1 + e_1) * (x_1 + e_1) + x_2 * x_2;
        }
        else
        {
            dist = (x_1 + e_1) * (x_1 + e_1) + x_2 * x_2;
        }
    }
    else
    {
        if (x_0 < -e_0)
        {
            dist = (x_0 + e_0) * (x_0 + e_0) + x_2 * x_2;
        }
        else if (x_0 > e_0)
        {
            dist = (x_0 - e_0) * (x_0 - e_0) + x_2 * x_2;
        }
        else
        {
            dist = x_2 * x_2;
        }
    }

    return dist;
}

double AABB::distance2To(const Vec3 &point)
{
    double dist = -1;
    Vec3 origin;
    Vec3 dim;

    {
        // +x face
        if (point.x >= max.x)
        {
            origin.x = max.x;
            origin.y = min.y;
            origin.z = min.z;

            dim = max - min;
            dim.x = 0.0;

            double curr_dist = distance2ToRectangle(origin, dim, point);

            if (dist < 0 || curr_dist < dist)
            {
                dist = curr_dist;
            }
        }
    }

    {
        // -x face
        if (point.x <= min.x)
        {
            origin.x = min.x;
            origin.y = min.y;
            origin.z = min.z;

            dim = max - min;
            dim.x = 0.0;

            double curr_dist = distance2ToRectangle(origin, dim, point);

            if (dist < 0 || curr_dist < dist)
            {
                dist = curr_dist;
            }
        }
    }

    {
        // +y face
        if (point.y >= max.y)
        {
            origin.x = min.x;
            origin.y = max.y;
            origin.z = min.z;

            dim = max - min;
            dim.y = 0.0;

            double curr_dist = distance2ToRectangle(origin, dim, point);

            if (dist < 0 || curr_dist < dist)
            {
                dist = curr_dist;
            }
        }
    }

    {
        // -y face
        if (point.y <= min.y)
        {
            origin.x = min.x;
            origin.y = min.y;
            origin.z = min.z;

            dim = max - min;
            dim.y = 0.0;

            double curr_dist = distance2ToRectangle(origin, dim, point);

            if (dist < 0 || curr_dist < dist)
            {
                dist = curr_dist;
            }
        }
    }

    {
        // +z face
        if (point.z >= max.z)
        {
            origin.x = min.x;
            origin.y = min.y;
            origin.z = max.z;

            dim = max - min;
            dim.z = 0.0;

            double curr_dist = distance2ToRectangle(origin, dim, point);

            if (dist < 0 || curr_dist < dist)
            {
                dist = curr_dist;
            }
        }
    }

    {
        // -z face
        if (point.z <= min.z)
        {
            origin.x = min.x;
            origin.y = min.y;
            origin.z = min.z;

            dim = max - min;
            dim.z = 0.0;

            double curr_dist = distance2ToRectangle(origin, dim, point);

            if (dist < 0 || curr_dist < dist)
            {
                dist = curr_dist;
            }
        }
    }

    return dist;
}

Triangle::Triangle(Vec3 *verts, Vec3 *norms)
{
    for (int i = 0; i < 3; i += 1)
    {
        this->verts[i] = verts[i];
        if (norms != nullptr)
        {
            this->norms[i] = norms[i];
        }
    }

    bounds.min = verts[0];
    bounds.max = verts[0];
    for (int i = 1; i < 3; i += 1)
    {
        bounds.min = min(bounds.min, verts[i]);
        bounds.max = max(bounds.max, verts[i]);
    }
    v_0_2 = verts[2] - verts[0];
    v_0_1 = verts[1] - verts[0];
        std::cout << v_0_1 << v_0_2 << std::endl;

}

double Triangle::distance2To(const Vec3 &point)
{
    // Adapted from https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

    Vec3 D = verts[0] - point;
    double a = dot(v_0_1, v_0_1);
    double b = dot(v_0_1, v_0_2);
    double c = dot(v_0_2, v_0_2);
    double d = dot(v_0_1, D);
    double e = dot(v_0_2, D);

    double delta = a * c - b * b;
    if (delta < 0)
    {
        delta = -delta;
    }

    double s = b * e - c * d;
    double t = b * d - a * e;

    if (s + t <= delta)
    {
        if (s < 0)
        {
            if (t < 0)
            {
                // Region 4
                if (d < 0)
                {
                    // On edge t = 0 with s > 0
                    t = 0;
                    if (-d >= a)
                    {
                        s = 1;
                    }
                    else
                    {
                        s = -d / a;
                    }
                }
                else
                {
                    // On edge s = 0
                    s = 0;
                    if (e >= 0)
                    {
                        t = 0;
                    }
                    else if (-e >= c)
                    {
                        t = 1;
                    }
                    else
                    {
                        t = -e / c;
                    }
                }
            }
            else
            {
                // Region 3

                // Line at s = 0, F(t) = 0 => t = -e/c, if t <= 0 or t >= 1 use endpoints
                s = 0;
                if (e >= 0)
                {
                    t = 0;
                }
                else if (-e >= c)
                {
                    t = 1;
                }
                else
                {
                    t = -e / c;
                }
            }
        }
        else if (t < 0)
        {
            // Region 5

            // Line at t = 0, F(s) = 0 -> s = -d/a, if s <=0 or s >=1 use endpoints
            t = 0;
            if (d >= 0)
            {
                s = 0;
            }
            else if (-d >= a)
            {
                s = 1;
            }
            else
            {
                s = -d / a;
            }
        }
        else
        {
            // Region 0

            // In triangle
            s /= delta;
            t /= delta;
        }
    }
    else
    {
        if (s < 0)
        {
            // Region 2

            double tmp0 = b + d;
            double tmp1 = c + e;
            if (tmp1 > tmp0)
            {
                // On edge s+t = 1 with s > 0
                double numer = tmp1 - tmp0;
                double denom = a - 2 * b + c;
                if (numer >= denom)
                {
                    s = 1;
                }
                else
                {
                    s = numer / denom;
                }
                t = 1 - s;
            }
            else
            {
                // On edge s = 0 with t <= 1
                s = 0;
                if (tmp1 <= 0)
                {
                    t = 1;
                }
                else if (e >= 0)
                {
                    t = 0;
                }
                else
                {
                    t = -e / c;
                }
            }
        }
        else if (t < 0)
        {
            // Region 6

            double tmp0 = b + e;
            double tmp1 = a + d;
            if (tmp1 > tmp0)
            {
                // Minimum on edge s + 1 = 1 with t > 0
                double numer = tmp1 - tmp0;
                double denom = a - 2 * b + c;
                if (numer >= denom)
                {
                    t = 1;
                }
                else
                {
                    t = numer / denom;
                }
                s = 1 - t;
            }
            else
            {
                // On edge t = 0 with s <= 1
                t = 0;
                if (tmp1 <= 0)
                {
                    s = 1;
                }
                else if (d >= 0)
                {
                    s = 0;
                }
                else
                {
                    s = -d / a;
                }
            }
        }
        else
        {
            // Region 1

            // This is a line (s+t = 1),  F(s) = 0 => s = ((c+e)-(b+d))/(a-2b+c). If s <= 0, or s >= 1 use endpoints instead
            double numer = (c + e) - (b + d);
            if (numer <= 0)
            {
                s = 0;
            }
            else
            {
                double denom = a - 2 * b + c;
                if (numer >= denom)
                {
                    s = 1;
                }
                else
                {
                    s = numer / denom;
                }
            }
            t = 1 - s;
        }
    }

    Vec3 closest_point = verts[0] + s * v_0_1 + t * v_0_2;
    double dist = length2(closest_point - point);
    std::cout << "Point: " << point << "\nClosest: " << closest_point << "\nDist: " << dist << std::endl;
    std::cout << "S: " << s << "\nT: " << t << std::endl;
    std::cout << "E0: " << v_0_1 << "\nE1: " << v_0_2 << std::endl;

    return dist;
}

bool Triangle::inFront(const Vec3 &point)
{
    Vec3 dir = point - verts[0];

    return dot(dir, norms[0]) >= 0;
}

Node::Node() : children{}, data{nullptr}, bounds{} {}
Node::~Node()
{
    if (children.size() > 0)
    {
        for (auto child : children)
        {
            delete child;
        }
    }
}

double Node::distance2To(const Vec3 &point, double best)
{
    if (data == nullptr)
    {
        for (auto child : children)
        {
            if (best < 0 || child->bounds.distance2To(point) < best)
            {
                double curr_dist = child->distance2To(point, best);
                if (best < 0 || (curr_dist >= 0 && curr_dist < best))
                {
                    best = curr_dist;
                }
            }
        }
    }
    else
    {
        if (data->inFront(point))
        {
            double curr_dist = data->distance2To(point);
            if (best < 0 || (curr_dist >= 0 && curr_dist < best))
            {
                best = curr_dist;
            }
        }
    }

    return best;
}

Tree::Tree() : head{nullptr} {}
Tree::~Tree()
{
    if (head)
        delete head;

    head = nullptr;
}

void updateHeight(Node *node, int height)
{
    node->height = height - node->height;
    for (auto child : node->children)
    {
        updateHeight(child, height);
    }
}

void Tree::build(const std::vector<Triangle*> &data, int max_children)
{
    std::vector<Node *> level;

    for (auto triangle : data)
    {
        Node *node = new Node();
        node->data = triangle;
        node->height = 0;
        node->bounds = triangle->bounds;
        level.push_back(node);
    }

    std::sort(level.begin(), level.end(), [](const Node *l, const Node *r) -> bool {
        int x_diff = l->data->verts[0].x - r->data->verts[0].x;
        if (x_diff == 0)
        {
            int y_diff = l->data->verts[0].y - r->data->verts[0].y;
            if (y_diff == 0)
            {
                int z_diff = l->data->verts[0].z - r->data->verts[0].z;
                return z_diff > 0;
            }
            else
            {
                return y_diff > 0;
            }
        }
        else
        {
            return x_diff > 0;
        }
    });

    int height = 1;
    for (; level.size() > 1;)
    {
        std::vector<Node *> next_level;
        int curr_children = 0;
        Node *curr_node = nullptr;
        for (int i = 0; i < level.size(); i += 1)
        {
            if (curr_children == 0)
            {
                curr_node = new Node();
                curr_node->data = nullptr;
                curr_node->height = height;
                curr_node->bounds = level[i]->bounds;
                next_level.push_back(curr_node);
            }
            else
            {
                curr_node->bounds += level[i]->bounds;
            }
            curr_node->children.push_back(level[i]);
            curr_children += 1;

            if (curr_children == max_children)
            {
                curr_children = 0;
            }
        }

        std::swap(level, next_level);
        next_level.clear();
        height += 1;
    }

    head = level.at(0);

    height -= 1;
    updateHeight(head, height);
}

double Tree::distance2To(const Vec3 &point)
{
    if (head == nullptr)
    {
        return -1;
    }
    double best = head->distance2To(point, -1);
    return best;
}

Mesh::Mesh(const std::string &filename)
{
    std::ifstream file_in(filename.c_str());

    if (file_in.good())
    {
        std::string line;
        for (; std::getline(file_in, line);)
        {
            std::stringstream ss(line);

            std::string code;

            ss >> code;
            if (code == "v")
            {
                double vx, vy, vz;
                ss >> vx >> vy >> vz;

                vertices.push_back(Vec3{vx, vy, vz});
                vert_norms.push_back(Vec3{0, 0, 0});
                facecount.push_back(0);
            }
            else if (code == "vn")
            {
                double nx, ny, nz;
                ss >> nx >> ny >> nz;

                norms.push_back(Vec3{nx, ny, nz});
            }
            else if (code == "f")
            {
                int v0, v1, v2;
                v0 = v1 = v2 = -1;

                bool vert_only = false;

                std::string index;
                ss >> index;

                std::stringstream index_ss(index);

                std::string token;
                getline(index_ss, token, '/');
                v0 = std::stoi(token) - 1;

                // Skip texture index
                if (getline(index_ss, token, '/'))
                {
                    if (getline(index_ss, token, '/'))
                    {
                        int n = std::stoi(token) - 1;
                        vert_norms[v0] += norms[n];
                        facecount[v0] += 1;
                    }
                }
                else
                {
                    // Vertex only
                    vert_only = true;
                }

                ss >> index;
                index_ss = std::stringstream(index);

                getline(index_ss, token, '/');
                v1 = std::stoi(token) - 1;

                // Skip texture index
                if (!vert_only && getline(index_ss, token, '/'))
                {
                    if (getline(index_ss, token, '/'))
                    {
                        int n = std::stoi(token) - 1;
                        vert_norms[v1] += norms[n];
                        facecount[v1] += 1;
                    }
                }

                ss >> index;
                index_ss = std::stringstream(index);

                getline(index_ss, token, '/');
                v2 = std::stoi(token) - 1;

                // Skip texture index
                if (!vert_only && getline(index_ss, token, '/'))
                {
                    if (getline(index_ss, token, '/'))
                    {
                        int n = std::stoi(token) - 1;
                        vert_norms[v2] += norms[n];
                        facecount[v2] += 1;
                    }
                }

                faces.push_back(Face{v0, v1, v2, !vert_only});
            }
        }

        for (auto face : faces)
        {
            if (face.has_norms)
            {
                continue;
            }

            // Needs to calculate normals
            Vec3 v_0_2 = vertices[face.v2] - vertices[face.v0];
            Vec3 v_0_1 = vertices[face.v1] - vertices[face.v0];

            Vec3 normal = normalize(cross(v_0_1, v_0_2));

            vert_norms[face.v0] += normal;
            facecount[face.v0] += 1;

            vert_norms[face.v1] += normal;
            facecount[face.v1] += 1;

            vert_norms[face.v2] += normal;
            facecount[face.v2] += 1;
        }

        for (int i = 0; i < vert_norms.size(); i += 1)
        {
            vert_norms[i] = normalize(vert_norms[i] / facecount[i]);
        }

        // Create triangles
        for (auto face : faces)
        {
            Vec3 tri_verts[3];
            Vec3 tri_norms[3];
            tri_verts[0] = vertices[face.v0];
            tri_verts[1] = vertices[face.v1];
            tri_verts[2] = vertices[face.v2];

            tri_norms[0] = vert_norms[face.v0];
            tri_norms[1] = vert_norms[face.v1];
            tri_norms[2] = vert_norms[face.v2];

            triangles.push_back(new Triangle{tri_verts, tri_norms});
        }
    }
}

Mesh::~Mesh() {
    for (auto tri : triangles) {
        delete tri;
    }
}

std::ostream &operator<<(std::ostream &out, const Triangle &tri)
{
    return out << "Vertices:\n"
               << tri.verts[0] << tri.verts[1] << tri.verts[2]
               << "Normals:\n"
               << tri.norms[0] << tri.norms[1] << tri.norms[2];
}

std::ostream &operator<<(std::ostream &out, const Mesh &mesh)
{
    out << "Triangles:\n";
    for (auto tri : mesh.triangles)
    {
        out << *tri;
    }
    return out;
}

std::ostream &operator<<(std::ostream &out, const Node &node)
{
    out << std::string(node.height * 2, ' ') << "Node\n";
    for (auto child : node.children)
    {
        out << *child;
    }
    return out;
}

std::ostream &operator<<(std::ostream &out, const Tree &tree)
{
    return out << *(tree.head);
}

int main(int argc, char **argv)
{
    if (argc >= 2)
    {
        std::string file_name(argv[1]);

        Mesh mesh(file_name);

        Tree tree;
        tree.build(mesh.triangles, 5);

        std::cout << "Min: " << tree.head->bounds.min << "\nMax: " << tree.head->bounds.max << std::endl;
        double dist = tree.distance2To(Vec3{100, 100, 100});
        std::cout << dist << std::endl;
    }

    // AABB bounds;
    // bounds.min = Vec3{0, 0, 0};
    // bounds.max = Vec3{1, 1, 1};

    // std::cout << bounds.distance2To(Vec3{0, 2, 0}) << std::endl;
    // std::cout << bounds.distance2To(Vec3{1, 1, 1}) << std::endl;
    // std::cout << bounds.distance2To(Vec3{1, 3, 1}) << std::endl;
}