//
// Created by te1t0ch1phead on 2022/6/6.
//

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


using namespace std;

#include "CGALPolygon.h"

CGALPolygon::CGALPolygon(const shared_ptr<MeshKernel::SurfaceMesh> &mesh) {
    std::vector<K2::Point_3> ps;
    std::vector<std::vector<std::size_t> > fs;
    for (int i = 0; i < mesh->VertexSize(); i++) {
        ps.emplace_back(mesh->vertices(MeshKernel::iGameVertexHandle(i)).x(),
                        mesh->vertices(MeshKernel::iGameVertexHandle(i)).y(),
                        mesh->vertices(MeshKernel::iGameVertexHandle(i)).z());

    }
    for (int i = 0; i < mesh->FaceSize(); i++) {
        fs.push_back({static_cast<unsigned long long>(mesh->faces(MeshKernel::iGameFaceHandle(i)).vh(0).idx()),
                      static_cast<unsigned long long>(mesh->faces(MeshKernel::iGameFaceHandle(i)).vh(1).idx()),
                      static_cast<unsigned long long>(mesh->faces(MeshKernel::iGameFaceHandle(i)).vh(2).idx())});
    }
    poly = new CGAL::Polyhedron_3<K2>();
    PMP::polygon_soup_to_polygon_mesh(ps, fs, *poly, CGAL::parameters::all_default());
    inside = new CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2>(*poly);

};

bool CGALPolygon::inMesh(K2::Point_3 v) {
    // CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2 > inside(poly);
    CGAL::Bounded_side res = (*inside)(v);
    if (res == CGAL::ON_BOUNDED_SIDE)
        return true;
    return false;
}

bool CGALPolygon::outMesh(K2::Point_3 v) {
    // CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2 > inside(poly);
    CGAL::Bounded_side res = (*inside)(v);
    if (res == CGAL::ON_UNBOUNDED_SIDE)
        return true;
    return false;
}


double cgal_vertex_triangle_dist(MeshKernel::iGameFace f, MeshKernel::iGameVertex v,
                                 std::shared_ptr<MeshKernel::SurfaceMesh> mesh) {
    Point a(mesh->fast_iGameVertex.at(f.vh(0)).x(), mesh->fast_iGameVertex.at(f.vh(0)).y(),
            mesh->fast_iGameVertex.at(f.vh(0)).z());
    Point b(mesh->fast_iGameVertex.at(f.vh(1)).x(), mesh->fast_iGameVertex.at(f.vh(1)).y(),
            mesh->fast_iGameVertex.at(f.vh(1)).z());
    Point c(mesh->fast_iGameVertex.at(f.vh(2)).x(), mesh->fast_iGameVertex.at(f.vh(2)).y(),
            mesh->fast_iGameVertex.at(f.vh(2)).z());
    Point point_query(v.x(), v.y(), v.z());
    K::FT sqd = squared_distance(K::Triangle_3(a, b, c), point_query);
    double res = sqrt(sqd);
    return res;

}



