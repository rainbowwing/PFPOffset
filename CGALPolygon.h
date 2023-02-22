//
// Created by te1t0ch1phead on 2022/6/6.
//

#ifndef THICKEN2_LIB_IMPL_H
#define THICKEN2_LIB_IMPL_H

#define myeps 1e-6

#define merge_eps 1e-3


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <queue>
#include <functional>
#include <set>
#include <memory>
#include <algorithm>
#include <bitset>
#include "MeshKernel/Mesh.h"



#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/enum.h>
#include <vector>
#include <fstream>
#include <limits>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;
typedef K::Plane_3 Plane_3;
typedef K::Intersect_3 Intersect_3;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef CGAL::Exact_predicates_exact_constructions_kernel K2;
typedef std::list< K2::Triangle_3>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K2, Iterator> Primitive;
typedef CGAL::AABB_traits<K2, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
namespace PMP = CGAL::Polygon_mesh_processing;



using namespace std;


class CGALPolygon {
    CGAL::Polyhedron_3<K2> * poly;
    CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2 > * inside;

public:
    CGALPolygon() {}
    CGALPolygon(const shared_ptr <MeshKernel::SurfaceMesh> &mesh);
    bool inMesh(K2::Point_3 v);
    bool outMesh(K2::Point_3 v);
};


double cgal_vertex_triangle_dist(MeshKernel::iGameFace f, MeshKernel::iGameVertex v, std::shared_ptr<MeshKernel::SurfaceMesh>mesh) ;


#endif //THICKEN2_LIB_IMPL_H
