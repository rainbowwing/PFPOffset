#define CGAL_HAS_THREADS
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
#include <thread>
#include "MeshKernel/Mesh.h"
#include "CGALPolygon.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "OsqpEigen/OsqpEigen.h"
#include <vector>
#include <fstream>
#include <limits>
#include <unordered_set>
#include <sstream>
#include "DSU.h"
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/version.h>
typedef CGAL::Constrained_Delaunay_triangulation_2<K2> CDT;
//#include "BVH.h"

//#define MODEIN

using namespace std;

shared_ptr <MeshKernel::SurfaceMesh> mesh;

shared_ptr<CGALPolygon>cgal_polygon;

double default_move = 0.1;
int thread_num = 16;
//Thicken2 ../data/
MeshKernel::SurfaceMesh ReadObjFile(const std::string &_InputFile) {
    //std::ifstream inputfile(_InputFile, std::ios::in);
    std::vector<MeshKernel::iGameVertex> vertices;
    std::vector<std::vector<MeshKernel::iGameVertexHandle> > faces;
    std::vector<double> move_dist;
    std::vector<std::vector<double>> normals;
    std::vector<std::vector<double>> uvs;
    std::unordered_map<int, int> V2N;// vertex to normal
    std::unordered_map<int, int> V2T;// vertex to uv
    std::cout << "Reading " << _InputFile << " File" << std::endl;
    FILE *inputfile = fopen(_InputFile.c_str(), "r");
    char inputs[100];
    while (fscanf(inputfile, "%[^\n]\n", inputs) != EOF) {
        string line(inputs);
        if (line[0] == '#') {
            continue;
        }
        std::stringstream linestream;
        linestream.str(line);
        std::string flag;
        linestream >> flag;
        if (flag == "v") {
            double x, y, z;
            linestream >> x >> y >> z;
            vertices.push_back(MeshKernel::iGameVertex(x, y, z));
        } else if (flag == "f") {
            std::vector<std::string> vex;
            std::string tmp;
            while (linestream >> tmp) vex.push_back(tmp);
            std::vector<MeshKernel::iGameVertexHandle> face(3);
            for (size_t i = 0; i < 3; i++) {
                size_t idx = 0;
                while (idx < vex[i].length() && std::isdigit(vex[i][idx])) idx++;
                int vh = std::stoi(vex[i].substr(0, idx)) - 1;//  obj start v from 1
                face[i] = (MeshKernel::iGameVertexHandle)(vh);
                if (vh >= vertices.size()) {
                    std::cerr << vh << " " << vertices.size() << std::endl;
                }
            }
            if (vex.size() >= 4)
                move_dist.push_back(std::stod(vex[3]));
            else
                move_dist.push_back(default_move);
            faces.push_back(face);
        }
    }
    if (!normals.empty()) {
        int ncnt = normals.size();
        for (int i = 0; i < vertices.size(); ++i) {
            int nidx = V2N[i];
            //if (nidx < 0 || nidx >= ncnt) printf("error: nidx = %d\n", nidx);// debug 用
            assert(nidx >= 0 && nidx < ncnt);
            vertices[i].setNormal(normals[nidx]);
        }
    }
    auto mesh = MeshKernel::SurfaceMesh(vertices, faces, move_dist);
    fclose(inputfile);
    return mesh;
}

double stx, sty, stz;
double grid_len = 1e10;

struct grid {
    int x;
    int y;
    int z;
    grid(){
        x=-1;
        y=-1;
        z=-1;
    }
    grid(int x,int y,int z){
        this->x=x;
        this->y=y;
        this->z=z;
    }
    friend bool operator<(const grid &a, const grid &b) {
        if (a.x != b.x)
            return a.x < b.x;
        if (a.y != b.y)
            return a.y < b.y;
        return a.z < b.z;
    }

    friend bool operator==(const grid &a, const grid &b) {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    }

    bool valid(){
        return x>=0 && y>=0 && z>=0;
    }
};

std::hash<int>int_hash;
struct grid_hash{
    size_t operator () (grid x) const {
        return int_hash(x.x) ^ int_hash(x.y<<6) ^ int_hash(x.z<<12);
    }};
struct grid_equal
{
    bool operator() (grid a,  grid b) const {
        return a.x == b.x  &&  a.y == b.y &&  a.z == b.z;
    }
};



vector <vector<int>> GridVertexDir = {{0, 0, 0}, //0
                                      {1, 0, 0}, //1
                                      {0, 1, 0}, //2
                                      {0, 0, 1}, //3
                                      {1, 1, 0}, //4
                                      {1, 0, 1}, //5
                                      {0, 1, 1}, //6
                                      {1, 1, 1}}; //7


vector <vector<int> > container_grid_face = {{0,3,5,1},
                                             {0,2,6,3},
                                             {0,1,4,2},
                                             {2,4,7,6},
                                             {3,6,7,5},
                                             {1,5,7,4}
};


MeshKernel::iGameVertex getGridVertex(grid g, int k) {
    double gsx = stx + g.x * grid_len;
    double gsy = sty + g.y * grid_len;
    double gsz = stz + g.z * grid_len;
    return MeshKernel::iGameVertex(gsx + GridVertexDir[k][0] * grid_len,
                                   gsy + GridVertexDir[k][1] * grid_len,
                                   gsz + GridVertexDir[k][2] * grid_len);
}

vector <vector<int>> DirectedGridEdge = {{1, 2, 3},
                                         { 4, 5},
                                         { 4, 6},
                                         {5, 6},
                                         { 7},
                                         { 7},
                                         { 7}};

MeshKernel::iGameVertex getGridiGameVertex(const MeshKernel::iGameVertex& small, const MeshKernel::iGameVertex& big, int k) {

    double x = small.x();
    double y = small.y();
    double z = small.z();
    if(GridVertexDir[k][0] > 0 )
        x = big.x();
    if(GridVertexDir[k][1] > 0 )
        y = big.y();
    if(GridVertexDir[k][2] > 0 )
        z = big.z();
    return MeshKernel::iGameVertex(x,y,z);
}

K2::Point_3 getGridK2Vertex(const K2::Point_3& small, const K2::Point_3& big, int k) {

    CGAL::Epeck::FT x = small.x();
    CGAL::Epeck::FT y = small.y();
    CGAL::Epeck::FT z = small.z();
    if(GridVertexDir[k][0] > 0 )
        x = big.x();
    if(GridVertexDir[k][1] > 0 )
        y = big.y();
    if(GridVertexDir[k][2] > 0 )
        z = big.z();
    return K2::Point_3(x,y,z);
}


grid getGridFrameVertex(grid g, int k) {
    grid v = g;
    v.x += GridVertexDir[k][0];
    v.y += GridVertexDir[k][1];
    v.z += GridVertexDir[k][2];
    return v;
}

grid vertex_to_grid(MeshKernel::iGameVertex v) {
    int x = int((v.x() - stx) / grid_len + myeps); // 先不管
    int y = int((v.y() - sty) / grid_len + myeps); // 先不管
    int z = int((v.z() - stz) / grid_len + myeps); // 先不管
    return grid{x, y, z};
}



struct GridVertex {
    int grid_type;
    vector <MeshKernel::iGameFaceHandle> face_list;
    vector<MeshKernel::iGameVertex>build_v;
    vector<int>build_v_global;
    vector<vector<size_t> >build_face;
    vector<K2::Triangle_3> generate_face_list;
    vector<K2::Triangle_3> generate_face_list_depth[5];
    GridVertex() {
        grid_type = -1;
    }
};



vector<MeshKernel::iGameVertex> field_move_vertex;
vector<vector<MeshKernel::iGameVertex> > field_move_face;
vector<K2::Triangle_3> field_move_K2_triangle;

//vector<vector<int> >approximate_field_face_table = {{0,1,2},{3,4,5},{0,2,4},{4,3,0},{1,5,4},{4,2,1},{0,3,5},{5,1,0}};

inline K::Point_3 iGameVertex_to_Point(const MeshKernel::iGameVertex& v){
    return K::Point_3(v.x(),v.y(),v.z());
}

inline K2::Point_3 iGameVertex_to_Point_K2(const MeshKernel::iGameVertex& v){
    return K2::Point_3(v.x(),v.y(),v.z());
}


inline  MeshKernel::iGameVertex Point_K2_to_iGameVertex(const K2::Point_3& v){
    return MeshKernel::iGameVertex(CGAL::to_double(v.x()),CGAL::to_double(v.y()),CGAL::to_double(v.z()));
}


bool segment_in_line(K2::Segment_3 a,K2::Segment_3  b){
    CGAL::Epeck::FT d0 = CGAL::squared_distance(a.supporting_line(),b.vertex(0));
    CGAL::Epeck::FT d1 = CGAL::squared_distance(a.supporting_line(),b.vertex(1));
    if(d0 <= CGAL::Epeck::FT(0) &&
       d1 <= CGAL::Epeck::FT(0))
        return true;
    return false;
}


bool segment_in_line2(K2::Segment_3 a,K2::Segment_3  b){
    CGAL::Epeck::FT d0 = CGAL::squared_distance(a.supporting_line(),b.vertex(0));
    CGAL::Epeck::FT d1 = CGAL::squared_distance(a.supporting_line(),b.vertex(1));
    if(d0 <= CGAL::Epeck::FT(myeps/10) &&
       d1 <= CGAL::Epeck::FT(myeps/10))
        return true;
    return false;
}

bool segment_in_line(K2::Segment_2 a,K2::Segment_2  b){
    K2::FT d0 = CGAL::squared_distance(a.supporting_line(),b.vertex(0));
    K2::FT d1 = CGAL::squared_distance(a.supporting_line(),b.vertex(1));
    if(d0 <= CGAL::Epeck::FT(0) &&
       d1 <= CGAL::Epeck::FT(0))
        return true;
    return false;
}


CGAL::Epeck::FT triangle_squared_aspect_ratio(vector<K2::Point_2>v){
    vector<CGAL::Epeck::FT> length;
    for(int i=0;i<3;i++){
        length.push_back((v[i] - v[(i+1)%3]).squared_length());
    }
    return *max_element(length.begin(),length.end()) / *min_element(length.begin(),length.end());
}

double triangle_squared_aspect_ratio(vector<MeshKernel::iGameVertex > v){
    vector<double> length;
    for(int i=0;i<3;i++){
        length.push_back((v[i] - v[(i+1)%3]).norm2());
    }
    return *max_element(length.begin(),length.end()) / *min_element(length.begin(),length.end());
}

//


vector<vector<K2::Point_3> > CGAL_CDT(vector<K2::Point_3> sorted_bound_vertex, vector<K2::Segment_3> cs,K2::Triangle_3 origin_face) {

    //CGAL::make_conforming_Delaunay_2();
    // cout << sorted_bound_vertex.size() <<" "<< cs.size() << endl;
    K2::Point_3 base_point_3d = origin_face.vertex(0);


    K2::Vector_3 origin_face_v0 = origin_face.vertex(1) - origin_face.vertex(0);
    K2::Vector_3 origin_face_v1 = origin_face.vertex(2) - origin_face.vertex(0);
    K2::Vector_3 new_direct = CGAL::cross_product(origin_face_v0,origin_face_v1);
    new_direct = new_direct /  CGAL::Epeck::FT(Point_K2_to_iGameVertex(K2::Point_3(0,0,0) +  new_direct).norm());
    K2::Vector_3 X_axis = origin_face_v0 / CGAL::Epeck::FT(Point_K2_to_iGameVertex(K2::Point_3(0,0,0) + origin_face_v0).norm());
    K2::Vector_3 Y_axis = CGAL::cross_product(new_direct,X_axis);


    vector<CDT::Vertex_handle> sorted_bound_vertex_in_cdt(sorted_bound_vertex.size());


    CDT cdt;
    map<CDT::Vertex_handle,K2::Point_3> mp;
    for(int i=0;i<sorted_bound_vertex.size();i++){
        CGAL::Epeck::FT x = (sorted_bound_vertex[i]-base_point_3d) * X_axis;
        CGAL::Epeck::FT y = (sorted_bound_vertex[i]-base_point_3d) * Y_axis;
        auto t = Point_K2_to_iGameVertex(sorted_bound_vertex[i]);
        // cout <<"X Y : " <<CGAL::to_double(x) <<" "<<CGAL::to_double(y) << endl;
        //cout <<"t: "<< t.x() <<" "<<t.y()<<" "<<t.z() << endl;
        auto tt = base_point_3d + x*X_axis + y*Y_axis;
        //  cout <<"tt: "<< tt.x() <<" "<<tt.y()<<" "<<tt.z() << endl;

        sorted_bound_vertex_in_cdt[i] = cdt.insert(K2::Point_2(x,y));
        mp[ sorted_bound_vertex_in_cdt[i] ] = sorted_bound_vertex[i];
    }


    CGAL::Polygon_2<K2>poly;

    for(int i=0;i<sorted_bound_vertex.size();i++){
        poly.push_back(cdt.point(sorted_bound_vertex_in_cdt[i]));
    }

    CGAL::make_conforming_Gabriel_2(cdt);

    //TODO 约束裁剪 ，，新增德劳内点怎么计算 ！！！！！！这里还有问题
    vector<vector<K2::Point_2> > faces;
    for(auto fit = cdt.finite_faces_begin();fit != cdt.finite_faces_end();fit++){
        faces.push_back({cdt.point(fit->vertex(0)),cdt.point(fit->vertex(1)),cdt.point(fit->vertex(2))});
    }
    for(K2::Segment_3 seg: cs){
        vector<vector<K2::Point_2> > faces_new;
        CGAL::Epeck::FT x0 = (seg.vertex(0) - base_point_3d) * X_axis;
        CGAL::Epeck::FT y0 = (seg.vertex(0) - base_point_3d) * Y_axis;
        CGAL::Epeck::FT x1 = (seg.vertex(1) - base_point_3d) * X_axis;
        CGAL::Epeck::FT y1 = (seg.vertex(1) - base_point_3d) * Y_axis;

        K2::Segment_2 seg2(K2::Point_2(x0,y0),K2::Point_2(x1,y1));

        for(auto face : faces){
            K2::Triangle_2 tri(face[0],face[1],face[2]);
            CGAL::cpp11::result_of<K2::Intersect_2(K2::Segment_2 , K2::Triangle_2)>::type
                    res_st = intersection(seg2,tri);

            if (res_st) {
                if (const K2::Segment_2 *s = boost::get<K2::Segment_2>(&*res_st)) {
                    bool is_same_edge = false;
                    for(int i=0;i<3;i++){
                        K2::Segment_2 edge(face[i],face[(i+1)%3]);
                        if(segment_in_line(edge,*s))
                            is_same_edge = true;
                    }
                    if(!is_same_edge){
                        CGAL::cpp11::result_of<K2::Intersect_2(K2::Line_2 , K2::Triangle_2)>::type
                                res_lt = intersection(seg2.supporting_line(),tri);
                        if (const K2::Segment_2 *ss = boost::get<K2::Segment_2>(&*res_lt)) {
                            K2::Point_2 v0 = ss->vertex(0);
                            K2::Point_2 v1 = ss->vertex(1);

                            vector <K2::Point_2> positive_side;
                            vector <K2::Point_2> negative_side;
                            for (int i = 0; i < 3; i++) {
                                if (seg2.supporting_line().has_on_positive_side(face[i])) {
                                    positive_side.push_back(face[i]);
                                } else if (seg2.supporting_line().has_on_negative_side(face[i])) {
                                    negative_side.push_back(face[i]);
                                }
                            }
                            for (vector <K2::Point_2> vs: {positive_side, negative_side}) {
                                if (vs.size() == 1) {
                                    faces_new.push_back({v0,v1,vs[0]});
                                } else if (vs.size() == 2) {
                                    vector<K2::Point_2> ans0v0;
                                    vector<K2::Point_2> ans0v1;
                                    vector<K2::Point_2> ans1v0;
                                    vector<K2::Point_2> ans1v1;
                                    K2::Segment_2 v0s0(v0,vs[0]);
                                    K2::Segment_2 v1s1(v1,vs[1]);
                                    CGAL::cpp11::result_of<K2::Intersect_2(K2::Segment_2 , K2::Segment_2)>::type
                                            res_ss2 = intersection(v0s0,v1s1);
                                    if(!res_ss2){
                                        ans0v0 = {v0,vs[0],v1};
                                        ans1v0 = {v0,v1,vs[1]};
                                    }
                                    else{
                                        ans0v0 = {v0,vs[1],v1};
                                        ans1v0 = {v0,v1,vs[0]};
                                    }
                                    ans0v1 = {vs[0],vs[1],v1};
                                    ans1v1 = {vs[0],vs[1],v0};
                                    if(triangle_squared_aspect_ratio(ans0v0) + triangle_squared_aspect_ratio(ans0v1) <
                                       triangle_squared_aspect_ratio(ans1v0) + triangle_squared_aspect_ratio(ans1v1)
                                            ){
                                        faces_new.push_back(ans0v0);
                                        faces_new.push_back(ans0v1);
                                    }
                                    else{
                                        faces_new.push_back(ans1v0);
                                        faces_new.push_back(ans1v1);
                                    }

                                }

                            }
                            continue;
                        }
                    }
                }
            }
            faces_new.push_back(face);
        }
        swap(faces_new,faces);
    }

    vector<vector<K2::Point_3> > ret;

    for(auto i : faces){
        vector<K2::Point_3>tmp(3);
        for(int j=0;j<3;j++){
            tmp[j] = base_point_3d + i[j].x()*X_axis + i[j].y()*Y_axis;
        }
        ret.push_back(tmp);
    }

    return ret;
}


struct ApproximateField {
    vector<MeshKernel::iGameVertex> origin_vertices;
    vector<MeshKernel::iGameVertex> extend_vertices;
    vector<K2::Tetrahedron_3>tet_list;
    vector<vector<MeshKernel::iGameVertex> > outer_face;
    vector<vector<MeshKernel::iGameVertex> > inner_face;
    vector<vector<MeshKernel::iGameVertex> > side_face;
    vector<MeshKernel::iGameVertex > bound_face_vertex;
    vector<vector<int> > bound_face_id;
    vector<bool>bound_face_useful;
    K2::Point_3 center;

    // vector<vector<MeshKernel::iGameVertex> > bounded_face;
    // vector<K2::Triangle_3 > bounded_face_k2;
    MeshKernel::iGameFaceHandle fh;
    ApproximateField(){}
    ~ApproximateField(){}
    ApproximateField(MeshKernel::iGameFaceHandle fh) {

        this->fh=fh;
        //this->mu = make_shared<std::mutex>();
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)]);
        //TODO 注意分解
        MeshKernel::iGameVertex normal = ((origin_vertices[1] - origin_vertices[0]) % (origin_vertices[2] - origin_vertices[0])).normalize();
        for(int i=0;i<3;i++){
            //extend_vertices.push_back(field_move_vertex[mesh->fast_iGameFace[fh].vh(i)]);
            //if(!in_triangle_positive_side(this->fh,field_move_vertex[mesh->fast_iGameFace[fh].vh(i)])){
            extend_vertices.push_back(field_move_vertex[mesh->fast_iGameFace[fh].vh(i)]);
            //}
            /* else{
                 extend_vertices.push_back(origin_vertices[i] + normal * mesh->fast_iGameFace[fh].move_dist);
             }*/
        }
        MeshKernel::iGameVertex new_normal = ((extend_vertices[1] - extend_vertices[0]) % (extend_vertices[2] - extend_vertices[0])).normalize();
        if(new_normal * normal <0)
            swap(extend_vertices[0],extend_vertices[1]);

        for(int i=0;i<3;i++) {
            tet_list.emplace_back(iGameVertex_to_Point_K2(origin_vertices[0]),
                                  iGameVertex_to_Point_K2(origin_vertices[1]),
                                  iGameVertex_to_Point_K2(origin_vertices[2]),
                                  iGameVertex_to_Point_K2(extend_vertices[i]));

        }
        for(int i=0;i<3;i++) {
            tet_list.emplace_back(iGameVertex_to_Point_K2(extend_vertices[0]),
                                  iGameVertex_to_Point_K2(extend_vertices[1]),
                                  iGameVertex_to_Point_K2(extend_vertices[2]),
                                  iGameVertex_to_Point_K2(origin_vertices[i]));
        }
        std::vector<K2::Triangle_3 >face_list;

        vector<vector<int> >side_face_generate;
        side_face_generate.push_back({0,1,2});
        side_face_generate.push_back({0+3,2+3,1+3});

        for(int i=0;i<3;i++)
            bound_face_vertex.push_back(origin_vertices[i]);
        for(int i=0;i<3;i++)
            bound_face_vertex.push_back(extend_vertices[i]);
        // MeshKernel::iGameVertex c(0,0,0);
        K2::Point_3 zero(0,0,0);
        K2::Vector_3 vv(0,0,0);
        for(int i=0;i<3;i++){
            vv += (iGameVertex_to_Point_K2(origin_vertices[i]) -zero);
            vv += (iGameVertex_to_Point_K2(extend_vertices[i]) -zero);
        }
        //c/=6;
        center = zero + (vv / 6);

        int cnt = 0;
        for(int i=0;i<3;i++) {
            auto iv0 = origin_vertices[(i+1)%3];
            auto iv1 = origin_vertices[i];
            for(int j=0;j<3;j++){
                auto ov = extend_vertices[j];
                vector<MeshKernel::iGameVertex> new_face{iv0,ov,iv1};
                K2::Triangle_3 this_tri(iGameVertex_to_Point_K2(iv0),iGameVertex_to_Point_K2(ov),
                                        iGameVertex_to_Point_K2(iv1));

                K2::Plane_3 this_plane_K1(iGameVertex_to_Point_K2(iv0),iGameVertex_to_Point_K2(ov),
                                          iGameVertex_to_Point_K2(iv1));
                bool positive_side = false;
                bool negative_side = false;
                for(int k=0;k<tet_list.size();k++){
                    positive_side |= this_tri.supporting_plane().has_on_positive_side(centroid(tet_list[k]));
                    negative_side |= this_tri.supporting_plane().has_on_negative_side(centroid(tet_list[k]));
                }
                if( positive_side ^ negative_side){
                    bool flag = true;
                    for(auto k : side_face_generate)
                    {

                        K2::Point_3 other_center = centroid(K2::Triangle_3 (iGameVertex_to_Point_K2(bound_face_vertex[k[0]]),
                                                                            iGameVertex_to_Point_K2(bound_face_vertex[k[1]]),
                                                                            iGameVertex_to_Point_K2(bound_face_vertex[k[2]])
                        ));
                        if(CGAL::squared_distance(other_center,this_plane_K1) == CGAL::Epeck::FT(0)){
                            K2::Point_3 other_v0 = iGameVertex_to_Point_K2(bound_face_vertex[k[0]]);
                            K2::Point_3 other_v1 = iGameVertex_to_Point_K2(bound_face_vertex[k[1]]);
                            K2::Point_3 other_v2 = iGameVertex_to_Point_K2(bound_face_vertex[k[2]]);
                            K2::Triangle_3 other_tri(other_v0,other_v1,other_v2);

                            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
                                    res_tt = intersection(this_tri,other_tri);
                            if (res_tt) {
                                if (const K2::Point_3 *p = boost::get<K2::Point_3>(&*res_tt)) {
                                    continue;
                                }
                                if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
                                    continue;
                                }
                                flag = false;
                                break;
                            }
                        }
                    }
                    if(flag){
                        side_face_generate.push_back({(i+1)%3,i,j+3});
                    }
                    //side_face_generate.push_back(new_face);
                }
            }
        }



        for(int i=0;i<3;i++){
            auto iv0 = extend_vertices[i];
            auto iv1 = extend_vertices[(i+1)%3];
            for(int j=0;j<3;j++){
                auto ov = origin_vertices[j];
                vector<MeshKernel::iGameVertex> new_face{iv0,ov,iv1};
                K2::Triangle_3 this_tri(iGameVertex_to_Point_K2(iv0),iGameVertex_to_Point_K2(ov),
                                        iGameVertex_to_Point_K2(iv1));
                K2::Plane_3 this_plane_K1(iGameVertex_to_Point_K2(iv0),iGameVertex_to_Point_K2(ov),
                                          iGameVertex_to_Point_K2(iv1));
                bool positive_side = false;
                bool negative_side = false;
                for(int k=0;k<tet_list.size();k++){
                    positive_side |= this_tri.supporting_plane().has_on_positive_side(centroid(tet_list[k]));
                    negative_side |= this_tri.supporting_plane().has_on_negative_side(centroid(tet_list[k]));
                }
                if( positive_side ^ negative_side){
                    bool flag = true;
                    for(auto k : side_face_generate)
                    {
                        K2::Point_3 other_center = centroid(K2::Triangle_3 (iGameVertex_to_Point_K2(bound_face_vertex[k[0]]),
                                                                            iGameVertex_to_Point_K2(bound_face_vertex[k[1]]),
                                                                            iGameVertex_to_Point_K2(bound_face_vertex[k[2]])
                        ));

                        if(CGAL::squared_distance(other_center,this_plane_K1) == CGAL::Epeck::FT(0) ){
                            K2::Point_3 other_v0 = iGameVertex_to_Point_K2(bound_face_vertex[k[0]]);
                            K2::Point_3 other_v1 = iGameVertex_to_Point_K2(bound_face_vertex[k[1]]);
                            K2::Point_3 other_v2 = iGameVertex_to_Point_K2(bound_face_vertex[k[2]]);
                            K2::Triangle_3 other_tri(other_v0,other_v1,other_v2);

                            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
                                    res_tt = intersection(this_tri,other_tri);
                            if (res_tt) {
                                if (const K2::Point_3 *p = boost::get<K2::Point_3>(&*res_tt)) {
                                    continue;
                                }
                                if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
                                    continue;
                                }
                                flag = false;
                                break;
                            }
                        }
                    }
                    if(flag){
                        side_face_generate.push_back({i+3,(i+1)%3+3,j});
                    }
                    //side_face_generate.push_back(new_face);
                }
            }
        }


        for(int i=0;i<side_face_generate.size();i++) {

            K2::Point_3 this_v0 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[i][0]]);
            K2::Point_3 this_v1 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[i][1]]);
            K2::Point_3 this_v2 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[i][2]]);
            K2::Triangle_3 this_tri(this_v0,this_v1,this_v2);
            K2::Ray_3 ray(CGAL::centroid(this_tri),CGAL::centroid(this_tri) + this_tri.supporting_plane().orthogonal_vector());
            K2::Ray_3 ray2(CGAL::centroid(this_tri),CGAL::centroid(this_tri) - this_tri.supporting_plane().orthogonal_vector());

            bool flag1 = true;
            bool flag2 = true;

            int yy= 0;
            for(int j=0;j<side_face_generate.size();j++) {
                if(i==j)continue;
                K2::Point_3 other_v0 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[j][0]]);
                K2::Point_3 other_v1 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[j][1]]);
                K2::Point_3 other_v2 = iGameVertex_to_Point_K2(bound_face_vertex[side_face_generate[j][2]]);
                K2::Triangle_3 other_tri(other_v0,other_v1,other_v2);
                CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Ray_3)>::type
                        res_tt = intersection(other_tri,ray);
                if(res_tt) {
                    flag1 = false;
                }
                CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Ray_3)>::type
                        res_tt2 = intersection(other_tri,ray2);
                if(res_tt2) {
                    flag2 = false;
                }
//                if(xx == 10){
//                    auto vvv0 = bound_face_vertex[side_face_generate[j][0]];
//                    auto vvv1 = bound_face_vertex[side_face_generate[j][1]];
//                    auto vvv2 = bound_face_vertex[side_face_generate[j][2]];
//                    printf("v %lf %lf %lf\n",vvv0.x(),vvv0.y(),vvv0.z());
//                    printf("v %lf %lf %lf\n",vvv1.x(),vvv1.y(),vvv1.z());
//                    printf("v %lf %lf %lf\n",vvv2.x(),vvv2.y(),vvv2.z());
//                    printf("f %d %d %d\n",yy+1,yy+2,yy+3);
//                    yy+=3;
//                }
            }
            if(flag1 || flag2 ){
                bound_face_id.push_back({side_face_generate[i][0],side_face_generate[i][1],side_face_generate[i][2]});
            }
            //cout << xx <<" "<< flag1 <<" "<< flag2 << endl; //10 坏了
            //xx+=3;
        }
        set<int >se;
        queue<int>q;
        q.push(0);
        while(!q.empty()){
            int id = q.front();
            q.pop();
            if(se.count(id))continue;
            se.insert(id);
            map<int,int> mp_id_vertex;
            for(int i=0;i<3;i++)
                mp_id_vertex[bound_face_id[id][i]] = i;

            for(int j=0;j<bound_face_id.size();j++){ // 找其他所有
                if(!se.count(j) && mp_id_vertex.count(bound_face_id[j][0])
                                   + mp_id_vertex.count(bound_face_id[j][1]) + mp_id_vertex.count(bound_face_id[j][2]) >= 2  ){
                    for(int k=0;k<3;k++){
                        if(mp_id_vertex.count(bound_face_id[j][k]) && mp_id_vertex.count(bound_face_id[j][(k + 1) % 3])) {
                            int in_id_k = mp_id_vertex[bound_face_id[j][k]];
                            int in_id_kp1 = mp_id_vertex[bound_face_id[j][(k+1)%3]];
                            if( (in_id_kp1 + 1) % 3 != in_id_k){
                                swap(bound_face_id[j][k],bound_face_id[j][(k+1)%3]);
                                break;
                            }
                        }
                    }
                    q.push(j);
                }
            }
        }





        // bounded_face.push_back({extend_vertices[0],extend_vertices[2],extend_vertices[1]});

        outer_face.push_back({extend_vertices[0],extend_vertices[2],extend_vertices[1]});
        inner_face.push_back({origin_vertices[0],origin_vertices[2],origin_vertices[1]});
        // bounded_face.push_back(outer_face[0]);


        vector<K2::Point_3> vertices_list;
        std::vector<std::vector<std::size_t> > faces_list;

        for (int i = 0; i < 6; i++) {
            vertices_list.push_back(iGameVertex_to_Point_K2(bound_face_vertex[i]));
        }
        for (auto i: bound_face_id) {
            faces_list.push_back({std::size_t(i[0]), std::size_t(i[1]), std::size_t(i[2])});
            bound_face_useful.push_back(true);
        }
        poly = new CGAL::Polyhedron_3<K2>();

        PMP::polygon_soup_to_polygon_mesh(vertices_list, faces_list, *poly, CGAL::parameters::all_default());
        inside_ptr = new CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2>(*poly);

        // 6个四面体法，防止相交处理麻烦 并且用tet 来判断内外这样每一个面的偏移就是6个tet，；
    }

    bool in_field(K2::Point_3 v) {
        //CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> inside(*poly);
        if ((*inside_ptr)(v) == CGAL::ON_BOUNDED_SIDE)
            return true;
        return false;
    }

    bool in_or_on_field(K2::Point_3 v) {
        //CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> inside(*poly);
        auto side = (*inside_ptr)(v);
        if (side== CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY)
            return true;
        return false;
    }

public:
    CGAL::Polyhedron_3<K2> * poly;

    CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> * inside_ptr;
    // std::shared_ptr<std::mutex> mu;

};
vector<ApproximateField>faces_approximate_field;


vector<MeshKernel::iGameVertex> min_move_g;
vector<MeshKernel::iGameVertex> max_move_g;
MeshKernel::iGameVertex do_quadratic_error_metric(MeshKernel::iGameVertexHandle vh,bool &is_succ,int depth=0){

    MeshKernel::iGameVertex v = mesh->fast_iGameVertex[vh];
    int m = mesh->FastNeighborFhOfVertex_[vh].size();
    Eigen::SparseMatrix<double> hessian(3, 3);     //P: n*n正定矩阵,必须为稀疏矩阵SparseMatrix
    hessian.setZero();
    Eigen::VectorXd gradient(3);                  //Q: n*1向量
    gradient.setZero();
    Eigen::SparseMatrix<double> linearMatrix(m, 3); //A: m*n矩阵,必须为稀疏矩阵SparseMatrix
    linearMatrix.setZero();
    Eigen::VectorXd lowerBound(m);                  //L: m*1下限向量
    lowerBound.setZero();
    Eigen::VectorXd upperBound(m);                  //U: m*1上限向量
    upperBound.setZero();
    int cnt = 0;

    double avg_move_dist=0;
    double max_move_dist = 0;
    MeshKernel::iGameVertex avg_move_vertex(0,0,0);
    for (auto f : mesh->FastNeighborFhOfVertex_[vh]){
        avg_move_dist += mesh->fast_iGameFace[f].move_dist;
        max_move_dist = max(max_move_dist, mesh->fast_iGameFace[f].move_dist);

    }
    avg_move_dist /= (1.0*mesh->FastNeighborFhOfVertex_[vh].size());

    for (auto f : mesh->FastNeighborFhOfVertex_[vh]) {
        MeshKernel::iGameVertex normal
                = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]
                    - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]) %
                   (mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]
                    - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])).normalize();


//        MeshKernel::iGameVertex new_v = v + normal * mesh->fast_iGameFace[f].move_dist;
//
//        MeshKernel::iGameVertex move_max_v = v + normal * mesh->fast_iGameFace[f].move_dist*1.35;
//        MeshKernel::iGameVertex move_min_v = v + normal * mesh->fast_iGameFace[f].move_dist*0.8;
#ifdef MODEIN
        normal = normal * -1;
#endif

        MeshKernel::iGameVertex new_v = v + normal * avg_move_dist;
        avg_move_vertex += normal * avg_move_dist;

        MeshKernel::iGameVertex move_max_v = v + normal * avg_move_dist * (1.35+0.15*depth);
        MeshKernel::iGameVertex move_min_v = v + normal * avg_move_dist * (0.80-0.15*depth);

        double d = -(normal.x() * new_v.x() + normal.y() * new_v.y() +  normal.z() * new_v.z());

        Eigen::Vector3d p(normal.x(),  normal.y(), normal.z());
        Eigen::Matrix3d A = p * p.transpose();
        Eigen::Vector3d D2(2*d*normal.x(),2*d*normal.y(),2*d*normal.z());
        Eigen::Vector3d LimitV(normal.x(),normal.y(),normal.z());

        double lower = normal.x()*move_min_v.x() + normal.y()*move_min_v.y() + normal.z()*move_min_v.z();
        double upper = normal.x()*move_max_v.x() + normal.y()*move_max_v.y() + normal.z()*move_max_v.z();



//        for(int i=0;i<3;i++)
//            for(int j=0;j<3;j++)
//                hessian.coeffRef(i,j)+=A.coeff(i,j)*2;
//        for(int i=0;i<3;i++)
//            gradient.coeffRef(i) +=  D2.coeffRef(i);


        for(int i=0;i<3;i++)
            linearMatrix.coeffRef(cnt,i) = LimitV[i];
        lowerBound.coeffRef(cnt) = lower;
        upperBound.coeffRef(cnt) = upper;

        cnt++;
    }
    avg_move_vertex/= mesh->FastNeighborFhOfVertex_[vh].size();

    min_move_g[vh] =  v + avg_move_vertex;

    avg_move_vertex = v + avg_move_vertex;;
    hessian.coeffRef(0,0) += (2.0);
    hessian.coeffRef(1,1) += (2.0);
    hessian.coeffRef(2,2) += (2.0);

    gradient.coeffRef(0) -= (2.0) * v.x();
    gradient.coeffRef(1) -= (2.0) * v.y();
    gradient.coeffRef(2) -= (2.0) * v.z();

    OsqpEigen::Solver solver;
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables(3);   //变量数n
    solver.data()->setNumberOfConstraints(m); //约束数m
    solver.data()->setHessianMatrix(hessian);
    solver.data()->setGradient(gradient);
    solver.data()->setLinearConstraintsMatrix(linearMatrix);
    solver.data()->setLowerBound(lowerBound);
    solver.data()->setUpperBound(upperBound);
    solver.initSolver();
    // return avg_move_vertex;
    if(solver.solve()) {
        is_succ = true;
        Eigen::VectorXd QPSolution = solver.getSolution();
        MeshKernel::iGameVertex res(QPSolution.coeffRef(0), QPSolution.coeffRef(1), QPSolution.coeffRef(2));
        if( (res-v).norm() > 2*avg_move_dist){
            return v+ (res-v)/(res-v).norm()*(2*avg_move_dist);
            // return avg_move_vertex;
        }

        return {QPSolution.coeffRef(0), QPSolution.coeffRef(1), QPSolution.coeffRef(2)};
    }
    else{
        is_succ = false;
        return do_quadratic_error_metric(vh,is_succ,depth+1);
        /* cout <<"v "<< avg_move_vertex.x()<<" "<<avg_move_vertex.y()<<" "<<avg_move_vertex.z() << endl;
         puts("use avg_move_vertex instead");
         return avg_move_vertex;*/
    }
}
double mix_factor = 0.5;
void mix(int T){
    for (int times = 0; times < T; times++) {
        std::vector<double> fix_move_dist;
        fix_move_dist.resize(mesh->FaceSize() + 1);
        for (auto i: mesh->allfaces()) {
            double dist = 0;
            for (auto j: mesh->NeighborFh(i.first)) {
                dist += mesh->faces(j).move_dist;
            }
            dist /= mesh->NeighborFh(i.first).size();
            fix_move_dist[i.first] = (1.0 - mix_factor) * i.second.move_dist + mix_factor * dist;
        }
        for (auto i: mesh->allfaces()) {
            mesh->faces(i.first).move_dist = fix_move_dist[i.first];
        }
    }
}


/*
 *  N*N*N 来确定点
 *   然后保证后，通过格点表面采样法，然后连接边还原
 *    正负不管了
 */

vector<vector<int> > container_grid_dir{{-1,-1,-1},{-1,-1,0},{-1,-1,1},
                                        {-1,0,-1},{-1,0,0},{-1,0,1},
                                        {-1,1,-1},{-1,1,0},{-1,1,1},
                                        {0,-1,-1},{0,-1,0},{0,-1,1},
                                        {0,0,-1},{0,0,1},
                                        {0,1,-1},{0,1,0},{0,1,1},
                                        {1,-1,-1},{1,-1,0},{1,-1,1},
                                        {1,0,-1},{1,0,0},{1,0,1},
                                        {1,1,-1},{1,1,0},{1,1,1},};


void sort_by_polar_order(vector<K2::Point_3>& v,K2::Vector_3 orthogonal_direction){
    K2::Vector_3 center_v(0,0,0);
    for(auto j: v){
        center_v += j - K2::Point_3 (0,0,0);
    }
    center_v /= v.size();
    if(v.size() <=1)
        return ;
    K2::Point_3 center =  K2::Point_3 (0,0,0) + center_v;

    function<int(CGAL::Epeck::FT,CGAL::Epeck::FT)> quadrant = [](CGAL::Epeck::FT x,CGAL::Epeck::FT y){
        auto zero = CGAL::Epeck::FT(0);
        if(x>  zero && y > zero)return 1;
        else if(x<= zero && y >  zero)return 2;
        else if(x<= zero && y <= zero)return 3;
        return 4;
    };
    Plane_3 p;
    K2::Vector_3 x_axis = (v[0] - center) / sqrt(CGAL::to_double(CGAL::squared_distance(v[0],center)));
    K2::Vector_3 y_axis = CGAL::cross_product(orthogonal_direction,x_axis);
    y_axis = y_axis / sqrt(CGAL::to_double(y_axis.squared_length()));

    sort(v.begin(),v.end(),[&](K2::Point_3 a, K2::Point_3 b){

        CGAL::Epeck::FT x1 = (a - center) * x_axis;
        CGAL::Epeck::FT y1 = (a - center) * y_axis;
        CGAL::Epeck::FT x2 = (b - center) * x_axis;
        CGAL::Epeck::FT y2 = (b - center) * y_axis;
        int q1 = quadrant(x1,y1);
        int q2 = quadrant(x2,y2);
        if(q1!=q2)return q1<q2;
        else
            return x1*y2 - x2*y1 > 0;

    });
    return ;
}

class MeshBuilder{
    vector<K2::Triangle_3>face_list;
    vector<K2::Point_3> v;
    vector<set<int> >face_near_point;
    vector<vector<set<int>>>edge_add;
    unordered_map<size_t,int>vmp;
public:
    MeshBuilder(){};

    MeshBuilder(const vector<K2::Triangle_3>&generate_face_final){
        this->face_list = generate_face_final;
        for(auto i : this->face_list){
            int sz = vmp.size();
            vmp[i.id()] = sz;
            v.push_back(i.vertex(0));
            v.push_back(i.vertex(1));
            v.push_back(i.vertex(2));
        }
        face_near_point.resize(generate_face_final.size());
    }
    void build(double eps = myeps*5){

        //  return ;
        DSU dsu(face_list.size() * 3);
        edge_add.resize(face_list.size());
        for(auto &i : edge_add){
            i.resize(3);
        }
        std::list<K2::Triangle_3> tri_list;
        for(auto i : face_list){
            tri_list.push_back(i);
        }

        Tree aabb_tree_final(tri_list.begin(), tri_list.end());
        FILE *file2 = fopen( "ansbud22.obj", "w");
        for(int face_id = 0; face_id < face_list.size(); face_id++){
            for(int id=0;id<3;id++){
                K2::Point_3 this_v = face_list[face_id].vertex(id);
                K2::Tetrahedron_3 tet(K2::Point_3(this_v.x(),this_v.y(),this_v.z()+eps),
                                      K2::Point_3(this_v.x(),this_v.y()+eps,this_v.z()-eps/2),
                                      K2::Point_3(this_v.x()-eps/2,this_v.y()-eps/2,this_v.z()-eps/2),
                                      K2::Point_3(this_v.x()+eps/2,this_v.y()-eps/2,this_v.z()-eps/2));
                std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
                aabb_tree_final.all_intersections(tet,std::back_inserter(intersections));
                for(auto i : intersections){
                    if(face_list[face_id].id() != i.second->id()){
                        int other_field_belong_id = vmp[i.second->id()];
                        face_near_point[other_field_belong_id].insert(face_id*3+id);
                    }
                }
            }
        }
        for(int i=0;i<face_list.size();i++){
            for(auto j : face_near_point[i]){
                bool flag = false;
                for(int k=0;k<3;k++) {
                    if (CGAL::squared_distance(v[j], v[i * 3 + k]) <= CGAL::Epeck::FT(eps)*CGAL::Epeck::FT(eps)) {
//                       static int xxxx=1;
//                       fprintf(file2,"v %lf %lf %lf\n",CGAL::to_double(v[j].x()),CGAL::to_double(v[j].y()),CGAL::to_double(v[j].z()));
//                       fprintf(file2,"v %lf %lf %lf\n",CGAL::to_double(v[i * 3 + k].x()),CGAL::to_double(v[i * 3 + k].y()),CGAL::to_double(v[i * 3 + k].z()));
//                       fprintf(file2,"l %d %d\n",xxxx,xxxx+1);
//                       xxxx+=2;
                        dsu.join(j,i * 3 + k);
                        flag = true;
                    }
                }
                if(!flag){

                    for(int k=0;k<3;k++) {
                        K2::Segment_3 se(v[i * 3 + k],v[ i * 3 + (k+1)%3]);
                        if (CGAL::squared_distance(v[j], se) <= CGAL::Epeck::FT(eps)*CGAL::Epeck::FT(eps)) {
                            edge_add[i][k].insert(j);
                            cout << i<<" "<<k<<" "<<j << endl;
                            flag = true;
                        }
                    }
                }
                //cout << i <<" "<<v[j].x() <<" "<< v[j].y()<<" "<< v[j].z() << endl;
            }
        }
        vector<int>obj_id(face_list.size() * 3);
        int cnt=0;
        FILE *file = fopen( "ansbud.obj", "w");
        for(int i=0;i<v.size();i++){
            if(dsu.find_root(i) == i){
                obj_id[i] = cnt;
                fprintf(file,"v %lf %lf %lf\n",CGAL::to_double(v[i].x()),CGAL::to_double(v[i].y()),CGAL::to_double(v[i].z()));
                cnt++;
            }
        }
        for(int i=0;i<face_list.size();i++){

        }

        for(int i=0;i<v.size();i+=3){
            if(set<int>{dsu.find_root(i)+1,
                        dsu.find_root(i+1)+1,
                        dsu.find_root(i+2)+1}.size()==3)
                fprintf(file,"f %d %d %d\n",obj_id[dsu.find_root(i)]+1,
                        obj_id[dsu.find_root(i+1)]+1,
                        obj_id[dsu.find_root(i+2)]+1);
        }
    }

};


using namespace std;

namespace CH3d{
    const int MAXN=305000;
    const double eps=1e-8;
    struct Point{
        double x,y,z;
        Point(){}
        Point(double xx,double yy,double zz):x(xx),y(yy),z(zz){}
        //两向量之差
        Point operator -(const Point p1){
            return Point(x-p1.x,y-p1.y,z-p1.z);
        }
        //两向量之和
        Point operator +(const Point p1){
            return Point(x+p1.x,y+p1.y,z+p1.z);
        }
        //叉乘
        Point operator *(const Point p){
            return Point(y*p.z-z*p.y,z*p.x-x*p.z,x*p.y-y*p.x);
        }
        Point operator *(double d){
            return Point(x*d,y*d,z*d);
        }
        Point operator / (double d){
            return Point(x/d,y/d,z/d);
        }
        //点乘
        double  operator ^(Point p){
            return (x*p.x+y*p.y+z*p.z);
        }
    };
    struct CH3D{
        struct face{
            //表示凸包一个面上的三个点的编号
            int a,b,c;
            //表示该面是否属于最终凸包上的面
            bool ok;
        };
        int n;//初始顶点数
        Point P[MAXN];//初始顶点
        int num; //凸包表面的三角形数
        face F[8*MAXN]; //凸包表面的三角形
        map<int,int>g[MAXN];//凸包表面的三角形
        double vlen(Point a){//向量长度
            return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        }
        //叉乘
        Point cross(const Point &a,const Point &b,const Point &c){
            return Point((b.y-a.y)*(c.z-a.z)-(b.z-a.z)*(c.y-a.y),(b.z-a.z)*(c.x-a.x)-(b.x-a.x)*(c.z-a.z),(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x));
        }
        //三角形面积*2
        double area(Point a,Point b,Point c){
            return vlen((b-a)*(c-a));
        }
        //四面体有向体积*6
        double volume(Point a,Point b,Point c,Point d){
            return (b-a)*(c-a)^(d-a);
        }
        //正：点在面同向
        double dblcmp(Point &p,face &f){
            Point m=P[f.b]-P[f.a];
            Point n=P[f.c]-P[f.a];
            Point t=p-P[f.a];
            return (m*n)^t;
        }
        void deal(int p,int a,int b){
            int f=g[a][b];//搜索与该边相邻的另一个平面
            face add;
            if(F[f].ok){
                if(dblcmp(P[p],F[f])>eps)dfs(p,f);
                else{
                    add.a=b;
                    add.b=a;
                    add.c=p;//这里注意顺序，要成右手系
                    add.ok=true;
                    g[p][b]=g[a][p]=g[b][a]=num;
                    F[num++]=add;
                }
            }
        }
        void dfs(int p,int now)//递归搜索所有应该从凸包内删除的面
        {
            F[now].ok=false;//删除可以看到的面
            deal(p,F[now].b,F[now].a);
            deal(p,F[now].c,F[now].b);
            deal(p,F[now].a,F[now].c);
        }
        bool same(int s,int t){
            Point &a=P[F[s].a];
            Point &b=P[F[s].b];
            Point &c=P[F[s].c];
            return fabs(volume(a,b,c,P[F[t].a]))<eps &&
                   fabs(volume(a,b,c,P[F[t].b]))<eps &&
                   fabs(volume(a,b,c,P[F[t].c]))<eps;
        }
        //构建三维凸包
        void create(){
            int i,j,tmp;
            face add;

            num=0;
            if(n<4)return;
            //**********************************************
            //此段是为了保证前四个点不共面
            bool flag=true;
            for(i=1;i<n;i++){
                if(vlen(P[0]-P[i])>eps){
                    swap(P[1],P[i]);
                    flag=false;
                    break;
                }
            }
            if(flag)return;
            flag=true;
            //使前三个点不共线
            for(i=2;i<n;i++){
                if(vlen((P[0]-P[1])*(P[1]-P[i]))>eps){
                    swap(P[2],P[i]);
                    flag=false;
                    break;
                }
            }
            if(flag)return;
            flag=true;
            //使前四个点不共面
            for(int i=3;i<n;i++){
                if(fabs((P[0]-P[1])*(P[1]-P[2])^(P[0]-P[i]))>eps){
                    swap(P[3],P[i]);
                    flag=false;
                    break;
                }
            }
            if(flag)return;
            //*****************************************
            for(i=0;i<4;i++){
                add.a=(i+1)%4;
                add.b=(i+2)%4;
                add.c=(i+3)%4;
                add.ok=true;
                if(dblcmp(P[i],add)>0)swap(add.b,add.c);
                g[add.a][add.b]=g[add.b][add.c]=g[add.c][add.a]=num;
                F[num++]=add;
            }
            for(i=4;i<n;i++){
                for(j=0;j<num;j++){
                    if(F[j].ok&&dblcmp(P[i],F[j])>eps){
                        dfs(i,j);
                        break;
                    }
                }
            }
            tmp=num;
            for(i=num=0;i<tmp;i++)if(F[i].ok)F[num++]=F[i];
        }
        set<int> get_v(){
            set<int>se;
            for(int i=0;i<num;i++){
                se.insert(F[i].a);
                se.insert(F[i].b);
                se.insert(F[i].c);
            }
            return se;
        }
        //表面积
        double area(){
            double res=0;
            if(n==3){
                Point p=cross(P[0],P[1],P[2]);
                res=vlen(p)/2.0;
                return res;
            }
            for(int i=0;i<num;i++)res+=area(P[F[i].a],P[F[i].b],P[F[i].c]);
            return res/2.0;
        }
        //体积
        double volume(){
            double res=0;
            Point tmp(0,0,0);
            for(int i=0;i<num;i++)res+=volume(tmp,P[F[i].a],P[F[i].b],P[F[i].c]);
            return fabs(res/6.0);
        }
        //表面三角形个数
        int triangle(){
            return num;
        }
        //表面多边形个数
        int polygon(){
            int i,j,res,flag;
            for(i=res=0;i<num;i++){
                flag=1;
                for(j=0;j<i;j++)
                    if(same(i,j)){
                        flag=0;
                        break;
                    }
                res+=flag;
            }
            return res;
        }
        //三维凸包重心
        Point barycenter(){
            Point ans(0,0,0),o(0,0,0);
            double all=0;
            for(int i=0;i<num;i++){
                double vol=volume(o,P[F[i].a],P[F[i].b],P[F[i].c]);
                ans=ans+(o+P[F[i].a]+P[F[i].b]+P[F[i].c])/4.0*vol;
                all+=vol;
            }
            ans=ans/all;
            return ans;
        }
        //点到面的距离
        double ptoface(Point p,int i){
            return fabs(volume(P[F[i].a],P[F[i].b],P[F[i].c],p)/vlen((P[F[i].b]-P[F[i].a])*(P[F[i].c]-P[F[i].a])));
        }
    };
}

CH3d::CH3D chg;

int main(int argc, char* argv[]) {


    cout <<"CGAL_RELEASE_DATE:" << CGAL_RELEASE_DATE << endl;
    string input_filename(argv[1]);
    vector<chrono::time_point<chrono::system_clock> > time_path;// start = std::chrono::system_clock::now();
    // FILE *file9 = fopen( (input_filename + "_9.obj").c_str(), "w");
    //FILE *file10 = fopen( (input_filename + "_10.obj").c_str(), "w");
//    FILE *file13 = fopen( (input_filename + "_13.off").c_str(), "w");
//
//
//    FILE *file5_1 = fopen( (input_filename + "_5_1.obj").c_str(), "w");
//    FILE *file5_3 = fopen( (input_filename + "_5_3.obj").c_str(), "w");
//    FILE *file5_2 = fopen( (input_filename + "_5_2.obj").c_str(), "w");
//    FILE *file5_4 = fopen( (input_filename + "_5_4.obj").c_str(), "w");
//    FILE *file5_5 = fopen( (input_filename + "_5_5.obj").c_str(), "w");
//    FILE *file5_6 = fopen( (input_filename + "_5_6.obj").c_str(), "w");
//    FILE *file5_7 = fopen( (input_filename + "_5_7.obj").c_str(), "w");
//    FILE *file5_8 = fopen( (input_filename + "_5_8.obj").c_str(), "w");
//    FILE *file5 = fopen( (input_filename + "_5.obj").c_str(), "w");
    FILE *file4 = fopen( (input_filename + "_4.obj").c_str(), "w");
//    FILE *file7 = fopen( (input_filename + "_7.obj").c_str(), "w");
    // freopen("../debugoutput.txt","w",stdout);
    default_move = 0.01;
    grid_len = 2.5;
    // 2.5 552531870234
    // 1.7 550170141470
    cout << grid_len <<endl;
    mix_factor = 0.5;
    //3.59
    mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile(input_filename)); grid_len = 0.1;
    mesh->initBBox();
    mesh->build_fast();
    cout <<"mesh->build_fast() succ" << endl;
    double default_move_dist = 0.05;
    if(argc > 2 )
        grid_len = stod(string(argv[2]));
    else
    {
        double sum = 0;
        for(int i=0;i<mesh->FaceSize();i++){
            double minx = *set<double>{(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)]).norm(),
                                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)]).norm(),
                                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)]).norm()
            }.begin();
            double maxx = *set<double>{(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)]).norm(),
                                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)]).norm(),
                                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)]).norm()
            }.rbegin();

            sum += sqrt(minx*maxx);

        }
        //grid_len = sum/mesh->FaceSize()*1.5;
        grid_len = sum/mesh->FaceSize()*2.5;
//        for(int i = 0;i<mesh->EdgeSize();i++ ){
//            sum += (mesh->fast_iGameVertex[mesh->fast_iGameEdge[i].vh(0)] - mesh->fast_iGameVertex[mesh->fast_iGameEdge[i].vh(1)]).norm();
//        }
//        sum /= mesh->EdgeSize();
        cout<<"GL:" <<grid_len << endl;
    }
    //mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/test_orgv2.obj2")); grid_len = 12.5; double default_move_dist = 0.8;
    if(argc > 3 ) {
        cout <<"default_move_dist : "<< stod(string(argv[3])) << endl;
        for (int i = 0; i < mesh->FaceSize(); i++) {
            mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist =  stod(string(argv[3]));
        }
    }
    else if(*input_filename.rbegin() != '2') {
        auto tmp = mesh->BBoxMax -  mesh->BBoxMin;
        default_move_dist = grid_len/2.5;//abs(*set<double>{tmp.x(),tmp.y(),tmp.z()}.begin());
        double minx = mesh->BBoxMin.z();
        double maxx = mesh->BBoxMax.z();
        cout << "default_move_dist: "<< default_move_dist << endl;
        for(int i=0;i<mesh->FaceSize();i++){
            double this_z = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(0)]+
                              mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(1)]+
                              mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(2)])/3).z();
            double ratio = (maxx-this_z)/(maxx-minx);

            mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist =default_move_dist; //max(default_move_dist*1.5*ratio*ratio,default_move_dist*0.5);
            //default_move_dist;
        }
    }
    else{
        cout <<"mix start" << endl;
        for(int times = 0; times <50;times++) {
            for (int i = 0; i < mesh->FaceSize(); i++) {
                double avg = mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist;
                int cnt = 1;
                for (auto j: mesh->NeighborFh(MeshKernel::iGameFaceHandle(i))) {
                    avg += mesh->fast_iGameFace[j].move_dist;
                    cnt++;
                }
                mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist = avg / cnt;
            }
        }
        cout <<"mix end" << endl;
    }


    // grid_len = max(grid_len,default_move_dist/2);
//    cout <<"mix start" << endl;
// //   mesh->build_fast();
//
//    for(int times = 0; times <50;times++) {
//        for (int i = 0; i < mesh->FaceSize(); i++) {
//            double avg = mesh->faces(MeshKernel::iGameFaceHandle(i)).move_dist;
//            int cnt = 1;
//            for (auto j: mesh->NeighborFh(MeshKernel::iGameFaceHandle(i))) {
//                avg += mesh->faces(j).move_dist;
//                cnt++;
//            }
//            mesh->faces(MeshKernel::iGameFaceHandle(i)).move_dist = avg / cnt;
//        }
//    }
//    cout <<"mix end" << endl;

    //只动xy
    faces_approximate_field.resize(mesh->FaceSize());
    field_move_vertex.resize(mesh->VertexSize());
    min_move_g.resize(mesh->VertexSize());
    max_move_g.resize(mesh->VertexSize());

    field_move_face.resize(mesh->FaceSize());
    field_move_K2_triangle.resize(mesh->FaceSize());
//    int f4id = 1;
//    cout <<"st do_quadratic_error_metric" << endl;
    time_path.push_back(std::chrono::system_clock::now());
    for(int i=0;i<mesh->VertexSize();i++){
        bool is_succ = true;

        field_move_vertex[i] = do_quadratic_error_metric(MeshKernel::iGameVertexHandle(i),is_succ);
//        fprintf(file4, "v %lf %lf %lf\n", mesh->fast_iGameVertex[i].x(), mesh->fast_iGameVertex[i].y(),
//                mesh->fast_iGameVertex[i].z());
//        fprintf(file4, "v %lf %lf %lf\n", field_move_vertex[i].x(), field_move_vertex[i].y(),
//                field_move_vertex[i].z());
////        fprintf(file5, "v %lf %lf %lf\n", field_move_vertex[i].x(), field_move_vertex[i].y(),
////                field_move_vertex[i].z());
//        fprintf(file4, "l %d %d\n", f4id, f4id + 1);
//        f4id+=2;
    }
    time_path.push_back(std::chrono::system_clock::now());

//    chg.n = mesh->VertexSize();
//    for(int i=0;i<mesh->VertexSize();i++){
//        chg.P[i].x = field_move_vertex[i].x();
//        chg.P[i].y = field_move_vertex[i].y();
//        chg.P[i].z = field_move_vertex[i].z();
//    }
//    chg.create();
//    set<int>se = chg.get_v();

    // cout <<"seseize" << se.size() << endl;
//    list<K2::Triangle_3> cutalg_aabbtree;
//    for(int i=0;i<mesh->FaceSize();i++){
//        K2::Triangle_3 tri(iGameVertex_to_Point_K2( field_move_vertex[mesh->fast_iGameFace[i].vh(0)]),
//                       iGameVertex_to_Point_K2( field_move_vertex[mesh->fast_iGameFace[i].vh(1)]),
//                       iGameVertex_to_Point_K2( field_move_vertex[mesh->fast_iGameFace[i].vh(2)])
//                       );
//        cutalg_aabbtree.push_back(tri);
//    }
//    Tree cut_aabb(cutalg_aabbtree.begin(),cutalg_aabbtree.end());
//
//    vector<vector<K2::Triangle_3 > >st_list;
//    vector<bool>is_cut;
//    vector< K2::Triangle_3>this_tri_list;
//    vector<vector<K2::Segment_3> >cs_list;
//    for(int i=0;i<mesh->FaceSize();i++){
//        std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_r;
//        K2::Triangle_3 tri(iGameVertex_to_Point_K2( field_move_vertex[mesh->fast_iGameFace[i].vh(0)]),
//                           iGameVertex_to_Point_K2( field_move_vertex[mesh->fast_iGameFace[i].vh(1)]),
//                           iGameVertex_to_Point_K2( field_move_vertex[mesh->fast_iGameFace[i].vh(2)])
//        );
//        cut_aabb.all_intersections(tri,std::back_inserter(intersections_r));
//        this_tri_list.push_back(tri);
//        //vector<bool>cutting_field_id(field_through_list.size(),false);
//        vector<K2::Segment_3>cs;
//        for(auto item : intersections_r) {
//            if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
//
//            }
//            else if(K2::Segment_3 * se = boost::get<K2::Segment_3>(&(item.first))){
//                if(!segment_in_line(*se,K2::Segment_3(tri.vertex(0),tri.vertex(1))) &&
//                    !segment_in_line(*se,K2::Segment_3(tri.vertex(2),tri.vertex(1))) &&
//                    !segment_in_line(*se,K2::Segment_3(tri.vertex(0),tri.vertex(2))))
//                cs.push_back(*se);
//            }
//        }
//        if(cs.size()){
//            is_cut.push_back(true);
//            cs_list.push_back(cs);
//            auto res = CGAL_CDT({tri.vertex(0),tri.vertex(1),tri.vertex(2)},cs,tri);
//            vector<K2::Triangle_3> ts;
//            st_list.push_back(ts);
//        }
//        else{
//            is_cut.push_back(false);
//            st_list.emplace_back();
//            cs_list.emplace_back();
//        }
//    }
//    set<int>se;
//    for(int i=0;i<mesh->FaceSize();i++){
//        int cnt=0;
//        for(int k=0;k<3;k++) {
//            std::list<Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections_r;
//            std::list<Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections_r2;
//            auto this_tri = this_tri_list[i];
//            K2::Ray_3 ray(this_tri.vertex(k),
//                          this_tri.vertex(k) + this_tri.supporting_plane().orthogonal_vector());
//            K2::Ray_3 ray2(this_tri.vertex(k),
//                           this_tri.vertex(k) - this_tri.supporting_plane().orthogonal_vector());
//
//            bool flag1 = true;
//            bool flag2 = true;
//
//            cut_aabb.all_intersections(ray, std::back_inserter(intersections_r));
//            cut_aabb.all_intersections(ray2, std::back_inserter(intersections_r2));
//            for (auto j: intersections_r) {
//
//                if (const K2::Point_3 *p = boost::get<K2::Point_3>(&j.first)) {
//                    if(CGAL::squared_distance(*p,this_tri.vertex(k))>CGAL::Epeck::FT(0))
//                        if (j.second->id() != this_tri.id())
//                            flag1 = false;
//                }
//
//            }
//            for (auto j: intersections_r2) {
//                if (const K2::Point_3 *p = boost::get<K2::Point_3>(&j.first)) {
//                    if(CGAL::squared_distance(*p,this_tri.vertex(k))>CGAL::Epeck::FT(0))
//                        if (j.second->id() != this_tri.id())
//                            flag2 = false;
//                }
//            }
//            if(flag1 || flag2 ){
//                cnt++;
//            }
//        }
//        if(cnt!=3){
//            se.insert(mesh->fast_iGameFace[i].vh(0));
//            se.insert(mesh->fast_iGameFace[i].vh(1));
//            se.insert(mesh->fast_iGameFace[i].vh(2));
//        }
//        //cout << xx <<" "<< flag1 <<" "<< flag2 << endl; //10 坏了
//        //xx+=3;
//    }
//    cout <<"ses "<< se.size() << endl;

////
//
//
//    vector<int>bfs_flag(mesh->FaceSize());
//    fill(bfs_flag.begin(),bfs_flag.end(),0);
//    vector<vector<MeshKernel::iGameVertex> > answ;
//    vector<pair<K2::Triangle_3,int> >vg;
//    queue<int>qq ;
//    for(int i=0;i<mesh->FaceSize();i++){
//        if(is_cut[i]){
//            for(int j=0;j<st_list[i].size();j++){
//                int state = 0;
//                for(int ii=0;ii<3;ii++){
//                    K2::Segment_3 sei( st_list[i][j].vertex(ii),
//                                       st_list[i][j].vertex((ii+1)%3));
//                    for(int jj=0;jj<cs_list[i].size();jj++) {
//                        if(segment_in_line2(sei,cs_list[i][jj])){
//                            state |= (1<<ii);
//                        }
//                    }
//                }
//                if(cs_list[i].size())cout <<"state" <<" "<< state << endl;
//                vg.emplace_back(st_list[i][j],state);
//            }
//        }
//        else  {
//            qq.push(vg.size());
//            bfs_flag[vg.size()] = 1;
//            K2::Triangle_3 tri(iGameVertex_to_Point_K2( field_move_vertex[mesh->fast_iGameFace[i].vh(0)]),
//                               iGameVertex_to_Point_K2( field_move_vertex[mesh->fast_iGameFace[i].vh(1)]),
//                               iGameVertex_to_Point_K2( field_move_vertex[mesh->fast_iGameFace[i].vh(2)]));
//            vg.push_back({tri,0});
//        }
//    }
//
//
//    list<K2::Triangle_3>li;
//    map<size_t,int>mp;
//    for(auto i : vg) {
//        mp[i.first.id()] = li.size();
//        li.push_back(i.first);
//    }
//    Tree aabb2(li.begin(),li.end());
//
//    while(!qq.empty()){
//        int nowid = qq.front();
//        qq.pop();
//        answ.push_back(vector<MeshKernel::iGameVertex>{
//                MeshKernel::iGameVertex(CGAL::to_double(vg[nowid].first.vertex(0).x()),
//                                        CGAL::to_double(vg[nowid].first.vertex(0).y()),
//                                        CGAL::to_double(vg[nowid].first.vertex(0).z())),
//                MeshKernel::iGameVertex(CGAL::to_double(vg[nowid].first.vertex(1).x()),
//                                        CGAL::to_double(vg[nowid].first.vertex(1).y()),
//                                        CGAL::to_double(vg[nowid].first.vertex(1).z())),
//                MeshKernel::iGameVertex(CGAL::to_double(vg[nowid].first.vertex(2).x()),
//                                        CGAL::to_double(vg[nowid].first.vertex(2).y()),
//                                        CGAL::to_double(vg[nowid].first.vertex(2).z()))
//
//        });
//
//        for(int ii=0;ii<3;ii++){
//           // if((vg[nowid].second) == 0){
//            if(((1<<ii)&(vg[nowid].second)) == 0){
//                K2::Segment_3 sei( vg[nowid].first.vertex(ii),
//                                   vg[nowid].first.vertex((ii+1)%3));
//                std::list< Tree::Intersection_and_primitive_id<K2::Segment_3>::Type> intersections_r;
//                aabb2.all_intersections(sei,std::back_inserter(intersections_r));
//                for(auto item : intersections_r){
//                    if(K2::Segment_3 * it = boost::get<K2::Segment_3>(&(item.first))) {
//                        int j = mp[item.second->id()];
//                        if(bfs_flag[j]==0){
//                            bfs_flag[j]=1;
//                            qq.push(j);
//                        }
//                    }
//                }
//            }
//
//        }
//    }
//
//
//
///*
//for(int k=0;k< st_list[j].size();k++){
//                            bool uuu = false;
//                            for(int ii=0;ii<3;ii++){
//                                K2::Segment_3 sei( tri.vertex(ii),
//                                                   tri.vertex((ii+1)%3));
//                                for(int jj=0;jj<3;jj++) {
//                                    K2::Segment_3 sej(st_list[j][k].vertex(jj),st_list[j][k].vertex((jj+1)%3));
//                                    if(segment_in_line2(sei,sej)){
//                                        uuu = true;
//                                    }
//                                }
//                            }
//                            if(uuu){
//                                answ.push_back({Point_K2_to_iGameVertex(st_list[j][k].vertex(0)),
//                                                Point_K2_to_iGameVertex(st_list[j][k].vertex(1)),
//                                                Point_K2_to_iGameVertex(st_list[j][k].vertex(2))});
//                                set<int>selo;
//                                selo.insert(k);
//                                queue<int>qq;
//                                qq.push(k);
//                                while(!qq.empty()){
//                                    int now = qq.front();
//                                    qq.pop();
//                                    for(int l=0;l<3;l++){
//                                        K2::Segment_3  ls(st_list[j][now].vertex(l), st_list[j][now].vertex((l+1)%3));
//                                        bool lf = false;
//                                        for(int m=0;m<cs_list[j].size();m++){
//                                            lf |= segment_in_line(cs_list[j][m],ls);
//                                        }
//                                        if(!lf){
//                                            for(int kk=0; kk < st_list[j].size(); kk++){
//                                                bool aaa = false;
//                                                if(!selo.count(kk)){
//                                                    for(int ii=0;ii<3;ii++){
//                                                        K2::Segment_3 sei( st_list[j][kk].vertex(ii),
//                                                                           st_list[j][kk].vertex((ii + 1) % 3));
//                                                            if(segment_in_line(sei,ls)){
//                                                                aaa = true;
//                                                            }
//                                                    }
//                                                }
//                                                if(aaa)
//                                                {
//                                                    answ.push_back({Point_K2_to_iGameVertex(st_list[j][kk].vertex(0)),
//                                                                    Point_K2_to_iGameVertex(st_list[j][kk].vertex(1)),
//                                                                    Point_K2_to_iGameVertex(st_list[j][kk].vertex(2))});
//                                                    selo.insert(kk);
//                                                    qq.push(kk);
//                                                }
//                                            }
//                                        }
//                                    }
//
//
//                                }
//
//                                break;
//                            }
//                        }
// */
//
//
//    FILE * ffff = fopen("../ress.obj","w");
//    FILE * gggg = fopen("../regg.obj","w");
//    FILE * llll = fopen("../rell.obj","w");
//    int f67id = 1;
//    for(int i=0;i<mesh->VertexSize();i++){
//        if(!se.count(i))
//        fprintf(gggg, "v %lf %lf %lf \n",field_move_vertex[i].x(),
//                field_move_vertex[i].y(),
//                field_move_vertex[i].z());
//    }
////    for(int i=0;i<mesh->FaceSize();i++){
////        if(is_cut[i]){
//////            if(se.count(mesh->fast_iGameFace[i].vh(0)) ||
//////               se.count(mesh->fast_iGameFace[i].vh(1))||
//////               se.count(mesh->fast_iGameFace[i].vh(2))
//////                    )continue;
////            for(auto kk:st_list[i]) {
////
////                fprintf(llll, "v %lf %lf %lf \n", field_move_vertex[mesh->fast_iGameFace[i].vh(0)].x(),
////                        field_move_vertex[mesh->fast_iGameFace[i].vh(0)].y(),
////                        field_move_vertex[mesh->fast_iGameFace[i].vh(0)].z());
////                fprintf(llll, "v %lf %lf %lf \n", field_move_vertex[mesh->fast_iGameFace[i].vh(1)].x(),
////                        field_move_vertex[mesh->fast_iGameFace[i].vh(1)].y(),
////                        field_move_vertex[mesh->fast_iGameFace[i].vh(1)].z());
////                fprintf(llll, "v %lf %lf %lf \n", field_move_vertex[mesh->fast_iGameFace[i].vh(2)].x(),
////                        field_move_vertex[mesh->fast_iGameFace[i].vh(2)].y(),
////                        field_move_vertex[mesh->fast_iGameFace[i].vh(2)].z());
////                fprintf(llll, "f %d %d %d\n", f67id, f67id + 1, f67id + 2);
////                f67id += 3;
////            }
////        }
////    }
//    for(int i=0;i<mesh->FaceSize();i++){
//        if(is_cut[i]){
////            if(se.count(mesh->fast_iGameFace[i].vh(0)) ||
////               se.count(mesh->fast_iGameFace[i].vh(1))||
////               se.count(mesh->fast_iGameFace[i].vh(2))
////                    )continue;
//            for(auto kk:st_list[i]) {
//
//                fprintf(llll, "v %lf %lf %lf \n",CGAL::to_double(kk.vertex(0).x()),
//                        CGAL::to_double(kk.vertex(0).y()),
//                        CGAL::to_double(kk.vertex(0).z()));
//                fprintf(llll, "v %lf %lf %lf \n",CGAL::to_double(kk.vertex(1).x()),
//                        CGAL::to_double(kk.vertex(1).y()),
//                        CGAL::to_double(kk.vertex(1).z()));
//                fprintf(llll, "v %lf %lf %lf \n",CGAL::to_double(kk.vertex(2).x()),
//                        CGAL::to_double(kk.vertex(2).y()),
//                        CGAL::to_double(kk.vertex(2).z()));
//                fprintf(llll, "f %d %d %d\n", f67id, f67id + 1, f67id + 2);
//                f67id += 3;
//            }
//        }
//    }
//
//
//    int f65id = 1;
//    for(int i=0;i<answ.size();i++){
//        fprintf(ffff, "v %lf %lf %lf \n",answ[i][0].x(),
//                answ[i][0].y(),
//        answ[i][0].z());
//        fprintf(ffff, "v %lf %lf %lf \n",answ[i][1].x(),
//                answ[i][1].y(),
//                answ[i][1].z());
//        fprintf(ffff, "v %lf %lf %lf \n",answ[i][2].x(),
//                answ[i][2].y(),
//                answ[i][2].z());
//
//        fprintf(ffff,"f %d %d %d\n",f65id,f65id+1,f65id+2);
//        f65id+=3;
//    }
//
//
//
//
//
//
//
//
//
//    return 0;
//    fclose(file4);
    //fclose(file5);
    int f4id = 1;
    for(int i=0;i<mesh->VertexSize();i++){

        fprintf(file4, "v %lf %lf %lf\n", mesh->fast_iGameVertex[i].x(), mesh->fast_iGameVertex[i].y(),
                mesh->fast_iGameVertex[i].z());
        fprintf(file4, "v %lf %lf %lf\n", field_move_vertex[i].x(), field_move_vertex[i].y(),
                field_move_vertex[i].z());
        fprintf(file4, "l %d %d\n", f4id, f4id + 1);
        f4id+=2;

    }







    //fclose(file4);
    cout <<"build st "<< endl;
    std::vector <std::shared_ptr<std::thread> > build_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        build_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for(int i=0;i<mesh->FaceSize();i++) {
                if (i % thread_num != now_id)continue;
                if(i%500==0)
                    cout << "build "<< i << endl;
                faces_approximate_field[i] = ApproximateField(MeshKernel::iGameFaceHandle(i));

                //continue;

                MeshKernel::iGameVertex v0 = field_move_vertex[mesh->fast_iGameFace[i].vh(0)];
                MeshKernel::iGameVertex v1 = field_move_vertex[mesh->fast_iGameFace[i].vh(1)];
                MeshKernel::iGameVertex v2 = field_move_vertex[mesh->fast_iGameFace[i].vh(2)];

                MeshKernel::iGameVertex ov0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)];
                MeshKernel::iGameVertex ov1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)];
                MeshKernel::iGameVertex ov2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)];

                MeshKernel::iGameVertex normal = (v1 - v0) % (v2 - v0);
                MeshKernel::iGameVertex normal_o = (ov1 - ov0) % (ov2 - ov0);
                if(normal * normal_o <0){
                    field_move_face[i]=vector<MeshKernel::iGameVertex>{v0,v2,v1};
                    //field_move_face[i]=vector<MeshKernel::iGameVertex>{v0,v1,v2};
                    auto center = (v0 + v1 + v2)/3;
                }
                else{
                    field_move_face[i]=vector<MeshKernel::iGameVertex>{v0,v1,v2};
                }
                field_move_K2_triangle[i] = K2::Triangle_3(iGameVertex_to_Point_K2(field_move_face[i][0]),
                                                           iGameVertex_to_Point_K2(field_move_face[i][1]),
                                                           iGameVertex_to_Point_K2(field_move_face[i][2]));


            }
        },i);
    }


    for(int i=0;i<thread_num;i++)
        build_thread_pool[i]->join();





    std::vector <std::shared_ptr<std::thread> > one_ring_select_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        one_ring_select_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for(int i=0;i<mesh->FaceSize();i++) {
                if (i % thread_num != now_id)continue;
                if(i%500==0)
                    cout << "one_ring_select_thread_pool "<< i << endl;

                set<int>neighbor_field;
                for(int j=0;j<3;j++) {
                    for (auto neighbor_id: mesh->FastNeighborFhOfEdge_[mesh->fast_iGameFace[i].eh(j)]) {
                        if(neighbor_id == i)continue;
                        neighbor_field.insert(neighbor_id);
                    }
                }



                vector<K2::Triangle_3>neighbor_face;
                for(auto neighbor_id: neighbor_field){
                    for (int k = 0; k < faces_approximate_field[neighbor_id].bound_face_id.size(); k++) {
                        vector<MeshKernel::iGameVertex> tmp{
                                faces_approximate_field[neighbor_id].bound_face_vertex[faces_approximate_field[neighbor_id].bound_face_id[k][0]],
                                faces_approximate_field[neighbor_id].bound_face_vertex[faces_approximate_field[neighbor_id].bound_face_id[k][1]],
                                faces_approximate_field[neighbor_id].bound_face_vertex[faces_approximate_field[neighbor_id].bound_face_id[k][2]]};
                        K2::Triangle_3 tri_this(iGameVertex_to_Point_K2(tmp[0]),
                                                iGameVertex_to_Point_K2(tmp[1]),
                                                iGameVertex_to_Point_K2(tmp[2])
                        );

                        neighbor_face.push_back(tri_this);

                    }
                }



                for(int j=0;j<faces_approximate_field[i].bound_face_id.size();j++) {
                    if ((faces_approximate_field[i].bound_face_id[j][0] >= 3 ||
                         faces_approximate_field[i].bound_face_id[j][1] >= 3 ||
                         faces_approximate_field[i].bound_face_id[j][2] >= 3) &&
                        (faces_approximate_field[i].bound_face_id[j][0] < 3 ||
                         faces_approximate_field[i].bound_face_id[j][1] < 3 ||
                         faces_approximate_field[i].bound_face_id[j][2] < 3) ){
                        bool flag = false;
                        vector<MeshKernel::iGameVertex> tmp{
                                faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][0]],
                                faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][1]],
                                faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][2]]};
                        K2::Triangle_3 tri_this(iGameVertex_to_Point_K2(tmp[0]),
                                                iGameVertex_to_Point_K2(tmp[1]),
                                                iGameVertex_to_Point_K2(tmp[2]));
                        K2::Segment_3 e0(tri_this.vertex(0), tri_this.vertex(1));
                        K2::Segment_3 e1(tri_this.vertex(1), tri_this.vertex(2));
                        K2::Segment_3 e2(tri_this.vertex(2), tri_this.vertex(0));

                        vector<K2::Segment_3> vs_tmp;
                        for (auto other: neighbor_face) { //Point_3, or Segment_3, or Triangle_3, or std::vector < Point_3 >
                            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
                                    res_tt = intersection(tri_this, other);
                            if (res_tt) {
                                if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tt)) {
                                    vs_tmp.push_back(*s);
                                } else if (const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&*res_tt)) {
                                    vs_tmp.emplace_back(t->vertex(0), t->vertex(1));
                                    vs_tmp.emplace_back(t->vertex(1), t->vertex(2));
                                    vs_tmp.emplace_back(t->vertex(2), t->vertex(0));
                                } else if (std::vector<K2::Point_3> *vs = boost::get<std::vector<K2::Point_3 >>(&*res_tt)) {
                                    sort_by_polar_order(*vs, tri_this.supporting_plane().orthogonal_vector());
                                    for (int k = 0; k < vs->size(); k++) {
                                        vs_tmp.emplace_back(vs->operator[](k), vs->operator[]((k + 1) % vs->size()));
                                        // cerr << "run iiiiiiiiiiiiiiit" << endl;
                                    }
                                }
                            }
                        }
                        vector<K2::Segment_3> vs;
                        for (auto se: vs_tmp) {
                            if (!segment_in_line(se, e0) && !segment_in_line(se, e1) && !segment_in_line(se, e2)) {
                                vs.push_back(se);
                            }
                        }
                        vector<vector<K2::Point_3> > res = CGAL_CDT({iGameVertex_to_Point_K2(tmp[0]),
                                                                     iGameVertex_to_Point_K2(tmp[1]),
                                                                     iGameVertex_to_Point_K2(tmp[2])}, vs, tri_this);
                        //cout << res.size() <<endl;
                        for (auto each_tri: res) {
                            bool patch_flag = false;
                            K2::Point_3 center = CGAL::centroid(K2::Triangle_3(each_tri[0], each_tri[1], each_tri[2]));
                            for (auto j: neighbor_field) {
                                if (tri_this.supporting_plane().oriented_side(faces_approximate_field[i].center) !=
                                    tri_this.supporting_plane().oriented_side(faces_approximate_field[j].center) &&
                                    tri_this.supporting_plane().oriented_side(faces_approximate_field[i].center)+
                                    tri_this.supporting_plane().oriented_side(faces_approximate_field[j].center) ==0 &&
                                    faces_approximate_field[j].in_or_on_field(center)) {
                                    patch_flag = true;
                                    break;
                                }
                            }
                            if(!patch_flag){
                                flag = true;
                                break;
                            }
                        }
                        faces_approximate_field[i].bound_face_useful[j]=flag;
                        //cout <<i<<" "<<j <<" "<< (flag?"yes":"no") << endl;
                    }
                }

                //for()

                //for(auto neighbor_id : mesh->FastNeighborFhOfEdge_[])
                //faces_approximate_field[i];

//                //continue;
//
//                MeshKernel::iGameVertex v0 = field_move_vertex[mesh->fast_iGameFace[i].vh(0)];
//                MeshKernel::iGameVertex v1 = field_move_vertex[mesh->fast_iGameFace[i].vh(1)];
//                MeshKernel::iGameVertex v2 = field_move_vertex[mesh->fast_iGameFace[i].vh(2)];
//
//                MeshKernel::iGameVertex ov0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)];
//                MeshKernel::iGameVertex ov1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)];
//                MeshKernel::iGameVertex ov2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)];
//
//                MeshKernel::iGameVertex normal = (v1 - v0) % (v2 - v0);
//                MeshKernel::iGameVertex normal_o = (ov1 - ov0) % (ov2 - ov0);
//                if(normal * normal_o <0){
//                    field_move_face[i]=vector<MeshKernel::iGameVertex>{v0,v2,v1};
//                    //field_move_face[i]=vector<MeshKernel::iGameVertex>{v0,v1,v2};
//                    auto center = (v0 + v1 + v2)/3;
//                }
//                else{
//                    field_move_face[i]=vector<MeshKernel::iGameVertex>{v0,v1,v2};
//                }
//                field_move_K2_triangle[i] = K2::Triangle_3(iGameVertex_to_Point_K2(field_move_face[i][0]),
//                                                           iGameVertex_to_Point_K2(field_move_face[i][1]),
//                                                           iGameVertex_to_Point_K2(field_move_face[i][2]));


            }
        },i);
    }

    for(int i=0;i<thread_num;i++)
        one_ring_select_thread_pool[i]->join();


    time_path.push_back(std::chrono::system_clock::now());
//
//    int f15id = 0;
//    int xxx=0;
//    for(int i=0;i<mesh->FaceSize();i++) {
//        int ccc = 0;
//        set<int>secc;
//        for(int k=0;k<faces_approximate_field[i].bound_face_useful.size();k++) {
//            if ((faces_approximate_field[i].bound_face_id[k][0] >= 3 &&
//                 faces_approximate_field[i].bound_face_id[k][1] >= 3 &&
//                 faces_approximate_field[i].bound_face_id[k][2] >= 3) ||
//                (faces_approximate_field[i].bound_face_id[k][0] < 3 &&
//                 faces_approximate_field[i].bound_face_id[k][1] < 3 &&
//                 faces_approximate_field[i].bound_face_id[k][2] < 3))
//                continue;
//            if(faces_approximate_field[i].bound_face_useful[k] == true) {
//                ccc++;
//                secc.insert(k);
//            }
//        }
//        if(ccc!=1  || faces_approximate_field[i].bound_face_id.size() !=8 )continue;
//
//        cout << "ccc" << ccc << endl;
//
//
////        if(xxx <3 && faces_approximate_field[i].bound_face_id.size() !=8 )continue;
////        if(xxx >=3 &&faces_approximate_field[i].bound_face_id.size() !=6 )continue;
////        if(xxx == 6) break;
//        cout << faces_approximate_field[i].bound_face_id.size() << endl;
//        xxx++;
//        if(xxx <15)continue;
//        FILE *file200 = fopen( (input_filename + "_" + to_string(200) +".obj").c_str(), "w");
//        FILE *file201 = fopen( (input_filename + "_" + to_string(201) +".obj").c_str(), "w");
//        FILE *file202 = fopen( (input_filename + "_" + to_string(202) +".obj").c_str(), "w");
//        FILE *file203 = fopen( (input_filename + "_" + to_string(203) +".obj").c_str(), "w");
//        FILE *file204 = fopen( (input_filename + "_" + to_string(204) +".obj").c_str(), "w");
//        FILE *file205 = fopen( (input_filename + "_" + to_string(205) +".obj").c_str(), "w");
//
//        int tt200 = 1;
//        int tt202 = 1;
//        int tt203 = 1;
//        int tt205 = 1;
//        cout << *secc.begin() << " "<< faces_approximate_field[i].bound_face_id[*secc.begin()][0]
//        <<" "<<faces_approximate_field[i].bound_face_id[*secc.begin()][1]<<" "
//        <<faces_approximate_field[i].bound_face_id[*secc.begin()][2] <<" "<<endl;
//
//        K2::Triangle_3 face205;
//
//        for(int j=0;j<faces_approximate_field[i].bound_face_id.size();j++){
//            if(faces_approximate_field[i].bound_face_id[j][0] <3 and
//               faces_approximate_field[i].bound_face_id[j][1] <3 and
//               faces_approximate_field[i].bound_face_id[j][2] < 3) {
//                vector<MeshKernel::iGameVertex> tmp{
//                        faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][0]],
//                        faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][1]],
//                        faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][2]]};
//                fprintf(file200, "v %lf %lf %lf \n", CGAL::to_double(tmp[0].x()),
//                        CGAL::to_double(tmp[0].y()),
//                        CGAL::to_double(tmp[0].z()));
//                fprintf(file200, "v %lf %lf %lf \n", CGAL::to_double(tmp[1].x()),
//                        CGAL::to_double(tmp[1].y()),
//                        CGAL::to_double(tmp[1].z()));
//                fprintf(file200, "v %lf %lf %lf \n", CGAL::to_double(tmp[2].x()),
//                        CGAL::to_double(tmp[2].y()),
//                        CGAL::to_double(tmp[2].z()));
//                fprintf(file200, "f %d %d %d\n", 1, 1 + 1, 1 + 2);
//            }
//            else{
//                vector<MeshKernel::iGameVertex> tmp{
//                        faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][0]],
//                        faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][1]],
//                        faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][2]]};
//
//                fprintf(file201, "v %lf %lf %lf \n", CGAL::to_double(tmp[0].x()),
//                        CGAL::to_double(tmp[0].y()),
//                        CGAL::to_double(tmp[0].z()));
//                fprintf(file201, "v %lf %lf %lf \n", CGAL::to_double(tmp[1].x()),
//                        CGAL::to_double(tmp[1].y()),
//                        CGAL::to_double(tmp[1].z()));
//                fprintf(file201, "v %lf %lf %lf \n", CGAL::to_double(tmp[2].x()),
//                        CGAL::to_double(tmp[2].y()),
//                        CGAL::to_double(tmp[2].z()));
//                fprintf(file201, "f %d %d %d\n", tt200, tt200 + 1, tt200 + 2);
//                tt200 += 3;
//            }
//
//            if(secc.count(j)){
//                vector<MeshKernel::iGameVertex> tmp{
//                        faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][0]],
//                        faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][1]],
//                        faces_approximate_field[i].bound_face_vertex[faces_approximate_field[i].bound_face_id[j][2]]};
//                face205 = K2::Triangle_3 (iGameVertex_to_Point_K2(tmp[0]),
//                                    iGameVertex_to_Point_K2(tmp[1]),
//                                    iGameVertex_to_Point_K2(tmp[2]));
//
//
//                fprintf(file205, "v %lf %lf %lf \n", CGAL::to_double(tmp[0].x()),
//                        CGAL::to_double(tmp[0].y()),
//                        CGAL::to_double(tmp[0].z()));
//                fprintf(file205, "v %lf %lf %lf \n", CGAL::to_double(tmp[1].x()),
//                        CGAL::to_double(tmp[1].y()),
//                        CGAL::to_double(tmp[1].z()));
//                fprintf(file205, "v %lf %lf %lf \n", CGAL::to_double(tmp[2].x()),
//                        CGAL::to_double(tmp[2].y()),
//                        CGAL::to_double(tmp[2].z()));
//                fprintf(file205, "f %d %d %d\n", tt205, tt205 + 1, tt205 + 2);
//                tt205 += 3;
//            }
//        }
//
//
//
//        cout << "nei se" << mesh->NeighborFh(MeshKernel::iGameFaceHandle(i)).size() << endl;
//        set<int>se;
//        for(auto j : mesh->NeighborFh(MeshKernel::iGameFaceHandle(int(i)))){
//            se.insert(mesh->faces(j).vh(0));
//            se.insert(mesh->faces(j).vh(1));
//            se.insert(mesh->faces(j).vh(2));
//        }
//        int tt204 = 1;
//        cout <<"st do_quadratic_error_metric" << endl;
//        for(auto i: se){
//
//            fprintf(file204, "v %lf %lf %lf\n", mesh->fast_iGameVertex[i].x(), mesh->fast_iGameVertex[i].y(),
//                    mesh->fast_iGameVertex[i].z());
//            fprintf(file204, "v %lf %lf %lf\n", field_move_vertex[i].x(), field_move_vertex[i].y(),
//                    field_move_vertex[i].z());
//            fprintf(file204, "l %d %d\n", tt204, tt204 + 1);
//            tt204+=2;
//
//
//        }
//
//        bool ffffflag = false;
//        for(auto nei: mesh->NeighborFh(MeshKernel::iGameFaceHandle(i))){
//            if(nei != i){
//                for(int j=0;j<faces_approximate_field[nei].bound_face_id.size();j++){
//                    vector<MeshKernel::iGameVertex> tmp{
//                            faces_approximate_field[nei].bound_face_vertex[faces_approximate_field[nei].bound_face_id[j][0]],
//                            faces_approximate_field[nei].bound_face_vertex[faces_approximate_field[nei].bound_face_id[j][1]],
//                            faces_approximate_field[nei].bound_face_vertex[faces_approximate_field[nei].bound_face_id[j][2]]};
//
//                    K2::Triangle_3 this_face(iGameVertex_to_Point_K2(tmp[0]),
//                                             iGameVertex_to_Point_K2(tmp[1]),
//                                             iGameVertex_to_Point_K2(tmp[2]));
//                    CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
//                            res_tt = intersection(this_face, face205);
//                    if (res_tt) {
//                        if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tt)) {
//                            if(!segment_in_line(*s,K2::Segment_3(face205.vertex(0),face205.vertex(1))) &&
//                                    !segment_in_line(*s,K2::Segment_3(face205.vertex(2),face205.vertex(1))) &&
//                                    !segment_in_line(*s,K2::Segment_3(face205.vertex(0),face205.vertex(2))))
//                            ffffflag = true;
//                        }
//                    }
//
//                    if(faces_approximate_field[nei].bound_face_id[j][0] <3 and
//                       faces_approximate_field[nei].bound_face_id[j][1] <3 and
//                       faces_approximate_field[nei].bound_face_id[j][2] < 3) {
//
//                        fprintf(file203, "v %lf %lf %lf \n", CGAL::to_double(tmp[0].x()),
//                                CGAL::to_double(tmp[0].y()),
//                                CGAL::to_double(tmp[0].z()));
//                        fprintf(file203, "v %lf %lf %lf \n", CGAL::to_double(tmp[1].x()),
//                                CGAL::to_double(tmp[1].y()),
//                                CGAL::to_double(tmp[1].z()));
//                        fprintf(file203, "v %lf %lf %lf \n", CGAL::to_double(tmp[2].x()),
//                                CGAL::to_double(tmp[2].y()),
//                                CGAL::to_double(tmp[2].z()));
//                        fprintf(file203, "f %d %d %d\n", tt203, tt203 + 1, tt203 + 2);
//                        tt203 += 3;
//                    }
//                    else{
//                        fprintf(file202, "v %lf %lf %lf \n", CGAL::to_double(tmp[0].x()),
//                                CGAL::to_double(tmp[0].y()),
//                                CGAL::to_double(tmp[0].z()));
//                        fprintf(file202, "v %lf %lf %lf \n", CGAL::to_double(tmp[1].x()),
//                                CGAL::to_double(tmp[1].y()),
//                                CGAL::to_double(tmp[1].z()));
//                        fprintf(file202, "v %lf %lf %lf \n", CGAL::to_double(tmp[2].x()),
//                                CGAL::to_double(tmp[2].y()),
//                                CGAL::to_double(tmp[2].z()));
//                        fprintf(file202, "f %d %d %d\n", tt202, tt202 + 1, tt202 + 2);
//                        tt202 += 3;
//                    }
//                }
//            }
//        }
//        if(ffffflag)
//        break;
//    }
//    return 0;



    // exit(0);

    cout <<"build end "<< endl;

    cgal_polygon = make_shared<CGALPolygon>(mesh);

    std::list < K2::Triangle_3> triangles;

    mesh->initBBox();

    double max_move = 0;
    double min_move = 1e10;
    for (auto i: mesh->allfaces()) {
        max_move = max(max_move, abs(i.second.move_dist));
        min_move = min(min_move, abs(i.second.move_dist));
    }


    stx = mesh->BBoxMin.x() - max_move;
    sty = mesh->BBoxMin.y() - max_move;
    stz = mesh->BBoxMin.z() - max_move;




    printf("GL %lf\n", grid_len);

    unordered_map <grid, GridVertex, grid_hash, grid_equal> frame_grid_mp;

    vector <vector<int>> bfs_dir = {{0,  0,  1},
                                    {0,  1,  0},
                                    {1,  0,  0},
                                    {0,  0,  -1},
                                    {0,  -1, 0},
                                    {-1, 0,  0}};
    std::function < vector<grid>(grid) > get_neighbor = [&](grid g) {
        vector <grid> ret;
        for (auto i: bfs_dir) {
            int xx, yy, zz;
            xx = g.x + i[0];
            yy = g.y + i[1];
            zz = g.z + i[2];
            //if (xx >= 0 && yy >= 0 && zz >= 0)
            ret.push_back({xx, yy, zz});
        }
        return ret;
    };
    int fsize = mesh->FaceSize();

    MeshKernel::iGameVertex debug_v(2.004797,-10.943162, -6.274028);
    grid debug_g =  vertex_to_grid(debug_v);

    cout <<"v to g :" <<debug_g.x <<" "<< debug_g.y <<" "<<debug_g.z << endl;

    // return 0;
    cout <<"bfs start \n" << endl;
    std::mutex bfs_mutex;
    std::vector <std::shared_ptr<std::thread> > bfs_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        bfs_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for (int face_id = 0; face_id < fsize; face_id++) {
                if(face_id%thread_num !=  now_id)continue;
                if (face_id % 1000 == 0)
                    printf("%d/%d\n", face_id, fsize);
                auto fh = make_pair(MeshKernel::iGameFaceHandle(face_id),
                                    mesh->faces(MeshKernel::iGameFaceHandle(face_id)));
                MeshKernel::iGameVertex center = (mesh->fast_iGameVertex[fh.second.vh(0)] +
                                                  mesh->fast_iGameVertex[fh.second.vh(1)] +
                                                  mesh->fast_iGameVertex[fh.second.vh(2)]) / 3;
                grid now = vertex_to_grid(center);
                vector <MeshKernel::iGameFaceHandle> face_and_neighbor;
                face_and_neighbor.push_back(fh.first);
                for (auto i: mesh->NeighborFh(fh.first)) {
                    face_and_neighbor.push_back(i);
                }

                queue <grid> q;
                unordered_set <grid,grid_hash,grid_equal> is_visit;
                vector <grid> center_neighbor = get_neighbor(now);
                for (auto bfs_start_node: center_neighbor) {
                    is_visit.insert(bfs_start_node);
                    q.push(bfs_start_node);
                }
                while (!q.empty()) {
                    now = q.front();
                    q.pop();
                    //cout << now.x <<" "<< now.y <<""
                    std::unique_lock<std::mutex>lock1(bfs_mutex,std::defer_lock);
                    lock1.lock();
                    auto iter = frame_grid_mp.find(now);
                    if (iter == frame_grid_mp.end()) {
                        iter = frame_grid_mp.insert(make_pair(now, GridVertex())).first;
                    }
                    iter->second.face_list.push_back(fh.first);
                    lock1.unlock();
                    vector <grid> neighbor = get_neighbor(now);
                    for (auto j: neighbor) {
                        if (!is_visit.count(j)) {
                            double dist = cgal_vertex_triangle_dist(fh.second, (getGridVertex(j, 0)+getGridVertex(j, 7))/2, mesh);
                            double move_limit = *set<double>
                                                {(field_move_vertex[fh.second.vh(1)]-mesh->fast_iGameVertex[fh.second.vh(0)]).norm(),
                                                 (field_move_vertex[fh.second.vh(2)]-mesh->fast_iGameVertex[fh.second.vh(1)]).norm(),
                                                 (field_move_vertex[fh.second.vh(2)]-mesh->fast_iGameVertex[fh.second.vh(0)]).norm()
                                                }.rbegin();
                            if (dist <  (grid_len*1.74/2+move_limit) ) { //TODO : zheli youhua cheng pianyi juli de shiji jisuan
                                q.push(j);
                                is_visit.insert(j);
                            }
                        }
                    }
                }
            }
        },i);
    }

    for(int i=0;i<thread_num;i++)
        bfs_thread_pool[i]->join();

    cout <<"bfs end \n" << endl;
    auto frame_grid_mp_bk =  frame_grid_mp;
    for(auto i : frame_grid_mp_bk) {
        //    if(i.first.x == 33 && i.first.y == 51 && i.first.z == 7 )
        //cout <<" ?????????????33517????" << std::endl;
        for(auto j : container_grid_dir){
            int nx = i.first.x+j[0];
            int ny = i.first.y+j[1];
            int nz = i.first.z+j[2];

            // if(nx>=0 && ny>=0 && nz>=0){
            for(auto z : i.second.face_list)
                frame_grid_mp[{nx,ny,nz}].face_list.push_back(z);

        }
    }

    std::vector <std::shared_ptr<std::thread> > each_grid_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        each_grid_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            int each_grid_cnt =-1;
            for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
                each_grid_cnt++;
                if(each_grid_cnt % thread_num != now_id)continue;
                if(each_grid_cnt % (thread_num*1000) == now_id)
                    printf("each_grid_cnt %d/%d\n",each_grid_cnt,(int)frame_grid_mp.size());

                MeshKernel::iGameVertex grid_vertex = getGridVertex(each_grid->first, 0);

                sort(each_grid->second.face_list.begin(), each_grid->second.face_list.end());
                each_grid->second.face_list.resize(
                        unique(each_grid->second.face_list.begin(), each_grid->second.face_list.end()) -
                        each_grid->second.face_list.begin());
//                each_grid->second.grid_type = vertex_state(each_grid->second.face_list,
//                                                           getGridVertex(each_grid->first, 0));
            }
        }, i);
    }


    for(int i=0;i<thread_num;i++)
        each_grid_thread_pool[i]->join();

    cout << "each_grid_cnt succ " <<endl;

    std::function<bool(Plane_3 plane,MeshKernel::iGameVertex,double)> vertex_in_plane
            = [&](Plane_3 pla,MeshKernel::iGameVertex v,double eps){
                if(sqrt(squared_distance(pla,Point_3(v.x(),v.y(),v.z())))<eps){
                    return true;
                }
                return false;
            };// 20 10 10 5

    long long  sum_grid = 0;
    vector<vector<size_t> >face_type_012{{0,1,2}};
    map<int,vector<long long > > debug_time_use;

//28 6 5+
    FILE *file9 = fopen( (input_filename + "_9.obj").c_str(), "w");
    for (auto each_grid= frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
        //if(!(each_grid->first.x == 26 && each_grid->first.y == 28 && each_grid->first.z == 9  ))continue;
        //if(!(each_grid->first.x == 19 && each_grid->first.y == 19 && each_grid->first.z == 1 ))continue;
        // if(!(each_grid->first.x == 1 && each_grid->first.y == 0 && each_grid->first.z == 3 ))continue;
        auto small  = getGridVertex(each_grid->first,0);
        auto big  = getGridVertex(each_grid->first,7);
        static int f3_id = 1;
        for (int ii = 0; ii < 7; ii++) {
            for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
                int from = ii;
                int to = DirectedGridEdge[ii][jj];
                MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
                MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
                if(fv.x() < (mesh->BBoxMin.x() + mesh->BBoxMax.x())/2)continue;
                fprintf(file9, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
                fprintf(file9, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
                fprintf(file9, "l %d %d\n", f3_id, f3_id + 1);
                f3_id += 2;
            }
        }
    }
    std::vector<std::vector<std::size_t> > each_grid_face_list;
    for(auto  each_container_face : container_grid_face){
        each_grid_face_list.push_back({(size_t)each_container_face[0],(size_t)each_container_face[1],(size_t)each_container_face[2]});
        each_grid_face_list.push_back({(size_t)each_container_face[2],(size_t)each_container_face[3],(size_t)each_container_face[0]});
    }

    //  return 0;

    cout << "each_grid_cnt succ2 " <<endl;
    // 上述代码完成距离场建格子的过程 8 ;
    // atomic<int>sum_face_size(0);
    // atomic<int>maxx_face_size(0);
    //atomic<int> qq1(0);
    // atomic<int> qq2(0);
    std::mutex mu;
    std::vector<std::shared_ptr<std::thread> > each_frame_thread(thread_num);
    vector<K2::Point_3> final_gen_vertex;
    vector<vector<size_t> > final_gen_face;
    queue<unordered_map <grid, GridVertex, grid_hash, grid_equal>::iterator >q;
    for(unordered_map <grid, GridVertex, grid_hash, grid_equal>::iterator i = frame_grid_mp.begin();i!=frame_grid_mp.end();i++)
        q.push(i);


    //std::atomic<unsigned long long> sttime;
    //std::atomic<unsigned long long> endtime;
    for(int i=0;i<thread_num;i++){
        each_frame_thread[i] = make_shared<thread>([&](int id){
            int tt=0;
            while(1) {
                unique_lock<std::mutex> lock(mu);
                if (!q.empty()) {
                    unordered_map <grid, GridVertex, grid_hash, grid_equal>::iterator each_grid = q.front();
                    q.pop();
                    if (( frame_grid_mp.size() -q.size()) % (frame_grid_mp.size()/50) == 0) { //2520
                        cout <<id <<" "<< ( q.size()) << " // " <<" "<< frame_grid_mp.size() << endl;
                    }
                    lock.unlock();
                    {
                        set <MeshKernel::iGameFaceHandle > face_set;
                        vector<int> face_list;
                        for (int j = 0; j < 8; j++) {
                            grid g = getGridFrameVertex(each_grid->first, j);
                            unordered_map<grid, GridVertex, grid_hash, grid_equal>::iterator it = frame_grid_mp.find(g);
                            if (it == frame_grid_mp.end())continue;
                            for (auto k: it->second.face_list) {
                                if (!face_set.count(k)) {
                                    face_set.insert(k);
                                    face_list.push_back(k);
                                }
                            }
                        }

                        K2::Point_3 small = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 0));
                        K2::Point_3 big = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 7));

                        std::mutex mu2;
                        std::function<void(K2::Point_3,K2::Point_3,vector<int>,int)> dfs = [&](K2::Point_3 small,K2::Point_3 big,vector<int>face_list,int depth){
                            std::vector<K2::Point_3> ps;
                            for(int i=0;i<8;i++){
                                ps.push_back(getGridK2Vertex(small,big,i));
                            }
                            std::list<K2::Triangle_3 >  frame_faces_list;
                            for(auto i : each_grid_face_list){
                                frame_faces_list.emplace_back(ps[i[0]],ps[i[1]],ps[i[2]]);
                            }
//                    double xxx = CGAL::to_double(small.x());
//                    double yyy = CGAL::to_double(small.y());
//                    double zzz = CGAL::to_double(small.z());
//                    double xx = CGAL::to_double(big.x());
//                    double yy = CGAL::to_double(big.y());
//                    double zz = CGAL::to_double(big.z());

                            Tree frame_aabb_tree(frame_faces_list.begin(),frame_faces_list.end());
                            CGAL::Polyhedron_3<K2> frame_poly;
                            PMP::polygon_soup_to_polygon_mesh(ps, each_grid_face_list, frame_poly, CGAL::parameters::all_default());
                            std::function<bool(K2::Point_3)> vertex_in_frame = [&](K2::Point_3 v){
                                return small.x() <= v.x() && v.x() <= big.x() &&
                                       small.y() <= v.y() && v.y() <= big.y() &&
                                       small.z() <= v.z() && v.z() <= big.z();
                            };
                            vector<K2::Triangle_3 > maybe_used_face;
                            vector<int> maybe_used_face_belong_field;
                            vector<vector<K2::Segment_3> > maybe_used_face_seg_cutting;

                            std::vector<MeshKernel::iGameFaceHandle> debug_tet_list_belong_face;
                            vector<int> field_through_list;
                            int state = 0;
                            for (int i: face_list) {
                                bool useful = false;
                                for(int j=0;j<0;j++){
                                    if(faces_approximate_field[i].in_field(ps[j])){
                                        state|=(1<<j);
                                    }
                                }
                                if(faces_approximate_field[i].in_field(midpoint(K2::Segment_3(ps[0],ps[7])))){
                                    useful = true;
                                }
                                if(!useful)
                                {
                                    if(vertex_in_frame(faces_approximate_field[i].center)){
                                        useful = true;
                                    }
                                }
                                if(!useful)
                                    useful = CGAL::Polygon_mesh_processing::do_intersect(*faces_approximate_field[i].poly,frame_poly);
                                if(useful)
                                    field_through_list.push_back(i);
                            }
                            if(state == ((1<<8)-1)) {
                                cout <<"skip" << endl;
                                return;
                            }

                            //   cout <<"field_through_list.size():  " <<field_through_list.size() << endl;

                            function<bool(K2::Triangle_3)>  triangle_through_grid = [&](K2::Triangle_3 tri) {
                                if(vertex_in_frame(centroid(tri))){
                                    return true;
                                }
                                CGAL::Polyhedron_3<K2> this_face;
                                vector<K2::Point_3> vs{tri.vertex(0),tri.vertex(1),tri.vertex(2)};
                                PMP::polygon_soup_to_polygon_mesh(vs, face_type_012, this_face, CGAL::parameters::all_default());
                                return CGAL::Polygon_mesh_processing::do_intersect(frame_poly,this_face);
                            };

                            //std::unordered_map<std::size_t,int> face_belong_field_mp;
                            std::list<K2::Triangle_3> field_triangles;
                            std::unordered_map<std::size_t,vector<int>> face_belong_field_source_id;

                            std::unordered_map<std::size_t,int> face_belong_field_all_mp;
                            std::list<K2::Triangle_3> field_triangles_all;
                            std::list<K2::Triangle_3> field_triangles_part_aabbtree;
                            for (int i=0;i< field_through_list.size();i++) {
                                for(int j=0;j<faces_approximate_field[field_through_list[i]].bound_face_id.size();j++){
                                    vector<MeshKernel::iGameVertex> tmp{faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][0]],
                                                                        faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][1]],
                                                                        faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][2]]};
                                    K2::Triangle_3 tri_this(iGameVertex_to_Point_K2(tmp[0]),
                                                            iGameVertex_to_Point_K2(tmp[1]),
                                                            iGameVertex_to_Point_K2(tmp[2])
                                    );
                                    if(faces_approximate_field[field_through_list[i]].bound_face_id[j][0] >= 3 || faces_approximate_field[field_through_list[i]].bound_face_id[j][1] >= 3 || faces_approximate_field[field_through_list[i]].bound_face_id[j][2] >= 3) { // 去掉原表面
                                        if (faces_approximate_field[field_through_list[i]].bound_face_useful[j]  && triangle_through_grid(tri_this)) {
                                            field_triangles.push_back(tri_this);
                                            field_triangles_part_aabbtree.push_back(tri_this);
                                            //face_belong_field_mp[tri_this.id()] = i;
//                                if(faces_approximate_field[field_through_list[i]].bound_face_id[j][0] >= 3  && faces_approximate_field[field_through_list[i]].bound_face_id[j][1] >= 3 && faces_approximate_field[field_through_list[i]].bound_face_id[j][2] >= 3) {
//                                    static int f52id = 1;
//                                    fprintf(file5_2, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(0).x()),
//                                            CGAL::to_double(tri_this.vertex(0).y()),
//                                            CGAL::to_double(tri_this.vertex(0).z()));
//                                    fprintf(file5_2, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(1).x()),
//                                            CGAL::to_double(tri_this.vertex(1).y()),
//                                            CGAL::to_double(tri_this.vertex(1).z()));
//                                    fprintf(file5_2, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(2).x()),
//                                            CGAL::to_double(tri_this.vertex(2).y()),
//                                            CGAL::to_double(tri_this.vertex(2).z()));
//                                    fprintf(file5_2, "f %d %d %d\n", f52id, f52id + 1, f52id + 2);
//                                    f52id += 3;
//                                }
//                                else{
//                                    static int f53id = 1;
//                                    fprintf(file5_3, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(0).x()),
//                                            CGAL::to_double(tri_this.vertex(0).y()),
//                                            CGAL::to_double(tri_this.vertex(0).z()));
//                                    fprintf(file5_3, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(1).x()),
//                                            CGAL::to_double(tri_this.vertex(1).y()),
//                                            CGAL::to_double(tri_this.vertex(1).z()));
//                                    fprintf(file5_3, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(2).x()),
//                                            CGAL::to_double(tri_this.vertex(2).y()),
//                                            CGAL::to_double(tri_this.vertex(2).z()));
//                                    fprintf(file5_3, "f %d %d %d\n", f53id, f53id + 1, f53id + 2);
//                                    f53id += 3;
//                                }

                                            face_belong_field_source_id[tri_this.id()] = vector<int>{faces_approximate_field[field_through_list[i]].bound_face_id[j][0],
                                                                                                     faces_approximate_field[field_through_list[i]].bound_face_id[j][1],
                                                                                                     faces_approximate_field[field_through_list[i]].bound_face_id[j][2]};

                                        }
                                    }
                                    else
                                        field_triangles_part_aabbtree.push_back(tri_this);
                                    field_triangles_all.push_back(tri_this);
                                    face_belong_field_all_mp[tri_this.id()] = i;
                                }
                            }
//                    cout <<"field_triangles.size():  " <<field_triangles.size() << endl;
//
//                    int f51id = 1;
//                    for(auto this_face : field_triangles_all){
//                        fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                                CGAL::to_double(this_face.vertex(0).y()),
//                                CGAL::to_double(this_face.vertex(0).z()));
//                        fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                                CGAL::to_double(this_face.vertex(1).y()),
//                                CGAL::to_double(this_face.vertex(1).z()));
//                        fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                                CGAL::to_double(this_face.vertex(2).y()),
//                                CGAL::to_double(this_face.vertex(2).z()));
//                        fprintf(file5_1,"f %d %d %d\n",f51id,f51id+1,f51id+2);
//                        f51id+=3;
//                    }
//
                            int f52id = 1;
                            int f53id = 1;
                            int f54id = 1;
                            for (int i=0;i< field_through_list.size();i++) {
                                for (int j = 0; j < faces_approximate_field[field_through_list[i]].bound_face_id.size(); j++) {
                                    vector<MeshKernel::iGameVertex> tmp{faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][0]],
                                                                        faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][1]],
                                                                        faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][2]]};
                                    K2::Triangle_3 this_face(iGameVertex_to_Point_K2(tmp[0]),
                                                             iGameVertex_to_Point_K2(tmp[1]),
                                                             iGameVertex_to_Point_K2(tmp[2])
                                    );
//                            if (faces_approximate_field[field_through_list[i]].bound_face_useful[j]  && triangle_through_grid(this_face)) {
//                                cout << faces_approximate_field[field_through_list[i]].bound_face_id[j][0] <<" "<< faces_approximate_field[field_through_list[i]].bound_face_id[j][1] <<" "<<faces_approximate_field[field_through_list[i]].bound_face_id[j][2] << endl;
//                                                        fprintf(file5_2, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                                                                CGAL::to_double(this_face.vertex(0).y()),
//                                                                CGAL::to_double(this_face.vertex(0).z()));
//                                                        fprintf(file5_2, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                                                                CGAL::to_double(this_face.vertex(1).y()),
//                                                                CGAL::to_double(this_face.vertex(1).z()));
//                                                        fprintf(file5_2, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                                                                CGAL::to_double(this_face.vertex(2).y()),
//                                                                CGAL::to_double(this_face.vertex(2).z()));
//                                                        fprintf(file5_2,"f %d %d %d\n",f52id,f52id+1,f52id+2);
//                                                        f52id+=3;
//                            }
//                            else if(faces_approximate_field[field_through_list[i]].bound_face_useful[j]){
//                                cout << faces_approximate_field[field_through_list[i]].bound_face_id[j][0] <<" "<< faces_approximate_field[field_through_list[i]].bound_face_id[j][1] <<" "<<faces_approximate_field[field_through_list[i]].bound_face_id[j][2] << endl;
//                                fprintf(file5_3, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                                        CGAL::to_double(this_face.vertex(0).y()),
//                                        CGAL::to_double(this_face.vertex(0).z()));
//                                fprintf(file5_3, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                                        CGAL::to_double(this_face.vertex(1).y()),
//                                        CGAL::to_double(this_face.vertex(1).z()));
//                                fprintf(file5_3, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                                        CGAL::to_double(this_face.vertex(2).y()),
//                                        CGAL::to_double(this_face.vertex(2).z()));
//                                fprintf(file5_3,"f %d %d %d\n",f53id,f53id+1,f53id+2);
//                                f53id+=3;
//                            }
                                }
                            }
//                    for(auto this_face : field_triangles_all){
//                        fprintf(file5_4, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                                CGAL::to_double(this_face.vertex(0).y()),
//                                CGAL::to_double(this_face.vertex(0).z()));
//                        fprintf(file5_4, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                                CGAL::to_double(this_face.vertex(1).y()),
//                                CGAL::to_double(this_face.vertex(1).z()));
//                        fprintf(file5_4, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                                CGAL::to_double(this_face.vertex(2).y()),
//                                CGAL::to_double(this_face.vertex(2).z()));
//                        fprintf(file5_4,"f %d %d %d\n",f54id,f54id+1,f54id+2);
//                        f54id+=3;
//                    }



                            if(field_triangles.size() ==0)return ;

//                    int f55id = 1;
//                    for(auto this_face : field_triangles){
//                        fprintf(file5_5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                                CGAL::to_double(this_face.vertex(0).y()),
//                                CGAL::to_double(this_face.vertex(0).z()));
//                        fprintf(file5_5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                                CGAL::to_double(this_face.vertex(1).y()),
//                                CGAL::to_double(this_face.vertex(1).z()));
//                        fprintf(file5_5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                                CGAL::to_double(this_face.vertex(2).y()),
//                                CGAL::to_double(this_face.vertex(2).z()));
//                        fprintf(file5_5,"f %d %d %d\n",f55id,f55id+1,f55id+2);
//                        f55id+=3;
//                    }

                            //cout <<"field_triangles.size():  " <<field_triangles.size() << endl;


                            Tree aabb_tree(field_triangles_all.begin(), field_triangles_all.end());
                            Tree aabb_tree_part(field_triangles_part_aabbtree.begin(), field_triangles_part_aabbtree.end());





                            //Tree aabb_tree_cross_usage;

                            /*for(auto i : each_grid_face_list){
                                frame_faces_list.emplace_back(ps[i[0]],ps[i[1]],ps[i[2]]);
                            }*/



                            function<bool(K2::Triangle_3)> check_in_field = [&](K2::Triangle_3 tri){
                                bool flag = false;
                                K2::Point_3 this_center = CGAL::centroid(tri);
                                K2::Ray_3 ray(this_center,tri.supporting_plane().orthogonal_vector());
                                std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections;
                                aabb_tree.all_intersections(ray,std::back_inserter(intersections));
                                //vector<bool>cutting_field_id(field_through_list.size(),false);
                                vector<set<K2::Point_3 > > intersection_v(field_through_list.size());
                                set<int>is_special;
                                set<int> positive_side;
                                for(auto item : intersections) {
                                    if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
                                        //se.insert(*p);
                                        //cout <<"*********"<<endl;
                                        int this_field_belong_id = face_belong_field_all_mp[item.second->id()];

                                        if( *p == this_center){
                                            //flag = true;
                                            positive_side.insert(this_field_belong_id);
                                            //cout <<"f1" << endl;

                                        }
                                        if(intersection_v[this_field_belong_id].count(*p)){
                                            if(!is_special.count(this_field_belong_id)){
                                                if(faces_approximate_field[field_through_list[this_field_belong_id]].in_or_on_field(this_center)) {
                                                    flag = true;
                                                    //cout <<"f2" << endl;
                                                    break;
                                                }
                                            }
                                            else
                                                is_special.insert(this_field_belong_id);
                                        }
                                        intersection_v[this_field_belong_id].insert(*p);

                                    }
                                }



                                for(int j=0;j<field_through_list.size();j++){
                                    if(!is_special.count(j) && positive_side.count(j)==0){
                                        if(intersection_v[j].size()%2) {
                                            flag = true;
                                            //cout <<"f3" << endl;
                                            break;
                                        }
                                    }
                                }
                                bool flag_positive = false;
                                bool flag_negative = false;
                                for(int j=0;j<field_through_list.size();j++){
                                    if(positive_side.count(j) && intersection_v[j].size()%2==0){
                                        flag_positive = true;
                                    }
                                }
                                if(!flag && flag_positive){
                                    K2::Ray_3 ray_r(this_center,tri.supporting_plane().orthogonal_vector()*(-1));
                                    std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections_r;
                                    aabb_tree.all_intersections(ray,std::back_inserter(intersections_r));
                                    //vector<bool>cutting_field_id(field_through_list.size(),false);
                                    vector<set<K2::Point_3 > > intersection_v_r(field_through_list.size());
                                    set<int>is_special_r;
                                    set<int> negative_side;
                                    for(auto item : intersections_r) {
                                        if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
                                            int this_field_belong_id = face_belong_field_all_mp[item.second->id()];

                                            if( *p == this_center){
                                                negative_side.insert(this_field_belong_id);
                                            }
                                            if(intersection_v_r[this_field_belong_id].count(*p)){
                                                if(!is_special_r.count(this_field_belong_id)){
                                                    if(faces_approximate_field[field_through_list[this_field_belong_id]].in_or_on_field(this_center)) {
                                                        flag = true;
                                                        break;
                                                    }
                                                }
                                                else
                                                    is_special_r.insert(this_field_belong_id);
                                            }
                                            intersection_v_r[this_field_belong_id].insert(*p);

                                        }
                                    }
                                    for(int j=0;j<field_through_list.size();j++){
                                        if(negative_side.count(j) && intersection_v_r[j].size()%2==0){
                                            flag_negative = true;
                                        }
                                    }
                                }
                                if(flag_positive && flag_negative)
                                    flag = true;
                                return flag;
                            };


//**************************************oct dfs skip
//                    bool skip = false;
//                    K2::Point_3 skipcc;
//                    if(field_through_list.size() > 28 ) {
//                        skip = true;
//                        for(auto  each_container_face : container_grid_face){
//                            K2::Triangle_3 tri1(ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]]);
//                            K2::Triangle_3 tri2(ps[each_container_face[2]],ps[each_container_face[3]],ps[each_container_face[0]]);
//                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_1;
//                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_2;
//                            aabb_tree.all_intersections(tri1,std::back_inserter(intersections_1));
//                            aabb_tree.all_intersections(tri2,std::back_inserter(intersections_2));
//                            vector<K2::Segment_3>vs;
//
//                            for(auto i : intersections_1){
//                                if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
//                                    vs.push_back(*s);
//                                }
//                            }
//                            for(auto i : intersections_2){
//                                if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
//                                    vs.push_back(*s);
//                                }
//                            }
//                            auto tmp = CGAL_CDT({ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]],ps[each_container_face[3]]},vs,tri1);
//                            for(auto i:tmp){
//                                //cout <<"skip true" << endl;
//                                if(!check_in_field(K2::Triangle_3(i[0],i[1],i[2]))){
//                                    // cout <<"skip false" << endl;
//                                    skipcc = CGAL::centroid(K2::Triangle_3(i[0],i[1],i[2]));
//                                    skip = false;
//                                    break;
//                                }
//                            }
//
////                    for(K2::Triangle_3 this_face : field_triangles){
////                        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
////                                res_tt = intersection(tri1,this_face);
////                        if (res_tt) {
////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
////                                vs.push_back(*p);
////                            }
////                        }
////                        res_tt = intersection(tri2,this_face);
////                        if (res_tt) {
////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
////                                vs.push_back(*p);
////                            }
////                        }
////                    }
//
//                            if(!skip)
//                                break;
//                        }
//                    }
//                    if(skip) {
//                        cout <<"skip"<<endl;
//                        return;
//                    }


                            // cout << "depth" << depth<< endl;
//                            if(depth){
//                                bool skip = false;
//                                K2::Point_3 skipcc;
//                                skip = true;
//                                for(auto  each_container_face : container_grid_face){
//                                    K2::Triangle_3 tri1(ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]]);
//                                    K2::Triangle_3 tri2(ps[each_container_face[2]],ps[each_container_face[3]],ps[each_container_face[0]]);
//                                    std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_1;
//                                    std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_2;
//                                    aabb_tree.all_intersections(tri1,std::back_inserter(intersections_1));
//                                    aabb_tree.all_intersections(tri2,std::back_inserter(intersections_2));
//                                    vector<K2::Segment_3>vs;
//
//                                    for(auto i : intersections_1){
//                                        if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
//                                            vs.push_back(*s);
//                                        }
//                                    }
//                                    for(auto i : intersections_2){
//                                        if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
//                                            vs.push_back(*s);
//                                        }
//                                    }
//                                    auto tmp = CGAL_CDT({ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]],ps[each_container_face[3]]},vs,tri1);
//                                    for(auto i:tmp){
//                                        //cout <<"skip true" << endl;
//                                        if(!check_in_field(K2::Triangle_3(i[0],i[1],i[2]))){
//                                            // cout <<"skip false" << endl;
//                                            skipcc = CGAL::centroid(K2::Triangle_3(i[0],i[1],i[2]));
//                                            skip = false;
//                                            break;
//                                        }
//                                    }
//
////                    for(K2::Triangle_3 this_face : field_triangles){
////                        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
////                                res_tt = intersection(tri1,this_face);
////                        if (res_tt) {
////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
////                                vs.push_back(*p);
////                            }
////                        }
////                        res_tt = intersection(tri2,this_face);
////                        if (res_tt) {
////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
////                                vs.push_back(*p);
////                            }
////                        }
////                    }
//
//                                    if(!skip)
//                                        break;
//                                }
//
//                                if(skip) {
//                                    cout <<"skip"<<endl;
//                                    return;
//                                }
//                                cout <<"50!!" << endl;
//                            }
                            if(field_through_list.size() > 25+depth*2.5 && depth<3){ ///
                                //cout <<field_through_list.size() <<" "<<"dfs!!" << endl;
                                K2::Point_3 center = midpoint(K2::Segment_3(small,big));
                                K2::Vector_3 delta = center - small;
                                vector<K2::Point_3> new_small;
                                vector<K2::Point_3> new_big;
                                for(int i=0;i<8;i++){
                                    new_small.push_back(getGridK2Vertex(small,center,i));
                                }
                                for(int i=0;i<8;i++){
                                    new_big.push_back(new_small[i]+delta);
                                }
                                vector<std::shared_ptr<std::thread> >vt(9);
                                for(int i=0;i<8;i++){
                                    vt[i] =  make_shared<std::thread>([=](int now_id) {
                                        dfs(new_small[i] ,new_big[i],field_through_list,depth+1);
                                    }, i);
                                }
                                for(int i=0;i<8;i++){
                                    vt[i]->join();
                                }
                                return ;
                            }
                            //cout << "field_through_list.size()" << field_through_list.size()<< endl;

//                    if(field_through_list.size() > 20){
//                        cout << each_grid->first.x<<" "<< each_grid->first.y<<" "<< each_grid->first.z << endl;
//                        exit(0);
//
//                    }




//                    if(skip) {
//                        int f5id = 1;
//                        FILE * file5_1 = fopen( (input_filename + "_5_1.obj").c_str(), "w");
//                        for(auto this_face : field_triangles_all){
//                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                                    CGAL::to_double(this_face.vertex(0).y()),
//                                    CGAL::to_double(this_face.vertex(0).z()));
//                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                                    CGAL::to_double(this_face.vertex(1).y()),
//                                    CGAL::to_double(this_face.vertex(1).z()));
//                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                                    CGAL::to_double(this_face.vertex(2).y()),
//                                    CGAL::to_double(this_face.vertex(2).z()));
//                            fprintf(file5_1,"f %d %d %d\n",f5id,f5id+1,f5id+2);
//                            f5id+=3;
//                        }
//                        FILE * file10 = fopen( (input_filename + "_10.obj").c_str(), "w");
//                        int f3_id = 1;
//                        for (int ii = 0; ii < 7; ii++) {
//                            for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
//                                int from = ii;
//                                int to = DirectedGridEdge[ii][jj];
//                                MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
//                                MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
//                                fprintf(file10, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
//                                fprintf(file10, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
//                                fprintf(file10, "l %d %d\n", f3_id, f3_id + 1);
//                                f3_id += 2;
//                            }
//                        }
//                        fclose(file5_1);
//                        fclose(file10);
//                        cout << "skip" << endl;
//                        int xxx;
//                        cin>>xxx;
//                    }
//                    else{
//                        int f5id = 1;
//                        FILE * file5_1 = fopen( (input_filename + "_5_1.obj").c_str(), "w");
//                        for(auto this_face : field_triangles_all){
//                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                                    CGAL::to_double(this_face.vertex(0).y()),
//                                    CGAL::to_double(this_face.vertex(0).z()));
//                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                                    CGAL::to_double(this_face.vertex(1).y()),
//                                    CGAL::to_double(this_face.vertex(1).z()));
//                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                                    CGAL::to_double(this_face.vertex(2).y()),
//                                    CGAL::to_double(this_face.vertex(2).z()));
//                            fprintf(file5_1,"f %d %d %d\n",f5id,f5id+1,f5id+2);
//                            f5id+=3;
//                        }
//                        FILE * file10 = fopen( (input_filename + "_10.obj").c_str(), "w");
//                        int f3_id = 1;
//                        for (int ii = 0; ii < 7; ii++) {
//                            for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
//                                int from = ii;
//                                int to = DirectedGridEdge[ii][jj];
//                                MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
//                                MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
//                                fprintf(file10, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
//                                fprintf(file10, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
//                                fprintf(file10, "l %d %d\n", f3_id, f3_id + 1);
//                                f3_id += 2;
//                            }
//                        }
//                        fclose(file5_1);
//                        fclose(file10);
//                        cout << "no skip" << endl;
//                        int xxx;
//                        cin>>xxx;
//                    }






                            std::vector<std::vector<int> >maybe_used_face_source_id;



                            for(auto this_face : field_triangles){
                                unique_lock<std::mutex> lock(mu);
                                each_grid->second.generate_face_list_depth[depth].push_back(this_face);
                                lock.unlock();
                                std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
                                aabb_tree.all_intersections(this_face,std::back_inserter(intersections));
                                int this_face_belong_id = face_belong_field_all_mp[this_face.id()];
                                vector<bool>cutting_field_id(field_through_list.size(),false);
                                vector<K2::Segment_3> segment_cutting;
                                for(auto i : intersections){
                                    //   cout << "count exist???" << face_belong_field_all_mp.count(i.second->id()) << endl;
                                    int this_field_belong_id  = face_belong_field_all_mp[i.second->id()];
                                    if(this_field_belong_id != this_face_belong_id /*&& !cutting_field_id[this_field_belong_id]*/){
                                        //   cout<< "interse:"<<this_face.id()<<"  "<< this_face_belong_id<<" "<<this_field_belong_id << endl;
                                        if(const K2::Point_3* p = boost::get<K2::Point_3>(&(i.first))){

                                        }
                                        else if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
                                            bool same_edge = false;
                                            for (K2::Segment_3 edge: {K2::Segment_3(this_face.vertex(0),this_face.vertex(1)),
                                                                      K2::Segment_3(this_face.vertex(1),this_face.vertex(2)),
                                                                      K2::Segment_3(this_face.vertex(2),this_face.vertex(0))}) {
                                                if (segment_in_line(edge, *s))
                                                    same_edge = true;
                                            }
                                            if (!same_edge) {
                                                cutting_field_id[this_field_belong_id] = true;
                                                //segment_cutting.push_back(*s);
                                            }
                                        }
                                        else if(const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&(i.first))) {
                                            int cnt=0;
                                            for(int l = 0 ;l<3;l++){
                                                for(int m = 0 ;m<3;m++){
                                                    if(CGAL::squared_distance(this_face.vertex(l),t->vertex(m)) == CGAL::Epeck::FT(0)){
                                                        cnt++;
                                                        break;
                                                    }
                                                }
                                            }
                                            if(cnt!=3)
                                                cutting_field_id[this_field_belong_id] = true;
                                        }
                                        else{
                                            cutting_field_id[this_field_belong_id] = true;
                                        }
                                    }
                                    // cout <<"cutting_field_id[this_field_belong_id]: " <<cutting_field_id[this_field_belong_id] << endl;
                                }
                                bool useless = false;
                                //int ofid=-1;
                                for(int i = 0; i< cutting_field_id.size();i++){
                                    if(!cutting_field_id[i] && i!=this_face_belong_id ){
                                        if(faces_approximate_field[field_through_list[i]].in_field(CGAL::centroid(this_face))) {
                                            // cout << i <<" ?? "<< this_face_belong_id << endl;
                                            //ofid = field_through_list[i];
                                            useless = true;
                                            break;
                                        }
                                    }
                                    if(useless)
                                        break;
                                }
                                if(!useless) {
                                    maybe_used_face.push_back(this_face);
                                    maybe_used_face_belong_field.push_back(this_face_belong_id);

                                    std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_part;
                                    aabb_tree_part.all_intersections(this_face,std::back_inserter(intersections_part));
                                    for(auto ii : intersections_part) {
                                        //cout << "count exist???" << face_belong_field_all_mp.count(i.second->id()) << endl;
                                        int this_field_belong_id = face_belong_field_all_mp[ii.second->id()];
                                        if (this_field_belong_id != this_face_belong_id) {
                                            if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&(ii.first))) {
                                                segment_cutting.push_back(*s);
                                            }
                                        }
                                    }
                                    maybe_used_face_seg_cutting.push_back(segment_cutting);
                                    maybe_used_face_source_id.push_back(face_belong_field_source_id[this_face.id()]);
                                }
//                   else{ //TODO :DEBUG CODE
//                       static int f5id = 1;
//                       cout <<"f5id" <<f5id<<" "<<this_face.id()<<"  "<<this_face_belong_id <<" "<< ofid<<endl;
//                       FILE *file5_a = fopen( (input_filename + "_5.a"+std::to_string(f5id)+".obj").c_str(), "w");
//                       FILE *file5_b = fopen( (input_filename + "_5.b"+std::to_string(f5id)+".obj").c_str(), "w");
//                       fprintf(file5_a, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                               CGAL::to_double(this_face.vertex(0).y()),
//                               CGAL::to_double(this_face.vertex(0).z()));
//                       fprintf(file5_a, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                               CGAL::to_double(this_face.vertex(1).y()),
//                               CGAL::to_double(this_face.vertex(1).z()));
//                       fprintf(file5_a, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                               CGAL::to_double(this_face.vertex(2).y()),
//                               CGAL::to_double(this_face.vertex(2).z()));
//                       fprintf(file5_a,"f %d %d %d\n",1,1+1,1+2);
//                       int bid = 1;
//                       for (int k = 0; k < faces_approximate_field[ofid].bound_face_id.size(); k++) {
//                           vector<MeshKernel::iGameVertex> tmp{
//                                   faces_approximate_field[ofid].bound_face_vertex[faces_approximate_field[ofid].bound_face_id[k][0]],
//                                   faces_approximate_field[ofid].bound_face_vertex[faces_approximate_field[ofid].bound_face_id[k][1]],
//                                   faces_approximate_field[ofid].bound_face_vertex[faces_approximate_field[ofid].bound_face_id[k][2]]};
//                           K2::Triangle_3 tri_this(iGameVertex_to_Point_K2(tmp[0]),
//                                                   iGameVertex_to_Point_K2(tmp[1]),
//                                                   iGameVertex_to_Point_K2(tmp[2])
//                           );
//                           fprintf(file5_b, "v %lf %lf %lf \n",CGAL::to_double(tri_this.vertex(0).x()),
//                                   CGAL::to_double(tri_this.vertex(0).y()),
//                                   CGAL::to_double(tri_this.vertex(0).z()));
//                           fprintf(file5_b, "v %lf %lf %lf \n",CGAL::to_double(tri_this.vertex(1).x()),
//                                   CGAL::to_double(tri_this.vertex(1).y()),
//                                   CGAL::to_double(tri_this.vertex(1).z()));
//                           fprintf(file5_b, "v %lf %lf %lf \n",CGAL::to_double(tri_this.vertex(2).x()),
//                                   CGAL::to_double(tri_this.vertex(2).y()),
//                                   CGAL::to_double(tri_this.vertex(2).z()));
//                           fprintf(file5_b,"f %d %d %d\n",bid,bid+1,bid+2);
//                           bid +=3;
//                       }
//                       f5id++;
//                   }
                            }


//                    int f511id = 1;
//                    for(auto this_face : maybe_used_face){
//                        fprintf(file5_11, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                                CGAL::to_double(this_face.vertex(0).y()),
//                                CGAL::to_double(this_face.vertex(0).z()));
//                        fprintf(file5_11, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                                CGAL::to_double(this_face.vertex(1).y()),
//                                CGAL::to_double(this_face.vertex(1).z()));
//                        fprintf(file5_11, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                                CGAL::to_double(this_face.vertex(2).y()),
//                                CGAL::to_double(this_face.vertex(2).z()));
//                        fprintf(file5_11,"f %d %d %d\n",f511id,f511id+1,f511id+2);
//                        f511id+=3;
//                    }

//                int f5id = 1;
//                for(auto this_face : maybe_used_face){
//                    fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                            CGAL::to_double(this_face.vertex(0).y()),
//                            CGAL::to_double(this_face.vertex(0).z()));
//                    fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                            CGAL::to_double(this_face.vertex(1).y()),
//                            CGAL::to_double(this_face.vertex(1).z()));
//                    fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                            CGAL::to_double(this_face.vertex(2).y()),
//                            CGAL::to_double(this_face.vertex(2).z()));
//                    fprintf(file5,"f %d %d %d\n",f5id,f5id+1,f5id+2);
//                    f5id+=3;
//                }


                            //   continue;

                            // 第三关

                            std::function<vector<K2::Point_3>(K2::Triangle_3) > cutting_triangle_by_grid = [&](K2::Triangle_3 this_face) {
                                vector<K2::Point_3 > ret;
                                set<K2::Point_3> se;
                                std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
                                frame_aabb_tree.all_intersections(this_face,std::back_inserter(intersections));
                                for(int i=0;i<3;i++){
                                    if(vertex_in_frame(this_face.vertex(i))){
                                        se.insert(this_face.vertex(i));
                                    }
                                }

                                for(auto i : intersections) {
                                    if(const K2::Point_3* p = boost::get<K2::Point_3>(&(i.first))){
                                        se.insert(*p);
                                    }
                                    if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&(i.first))) {
                                        se.insert(s->vertex(0));
                                        se.insert(s->vertex(1));
                                    }
                                    else if(const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&(i.first))){
                                        se.insert(t->vertex(0));
                                        se.insert(t->vertex(1));
                                        se.insert(t->vertex(2));
                                    }
                                    else if(const std::vector<K2::Point_3> *v = boost::get<std::vector<K2::Point_3>>(&(i.first))){
                                        for(const auto& ii : *v)
                                            se.insert(ii);
                                    }
                                }
                                for(const auto& i : se)
                                    ret.push_back(i);

                                return ret;
                            };

                            vector<K2::Triangle_3 >possible_face_inner_part;

                            // cout <<"********1****" << endl;
                            vector<vector<K2::Point_3> > face_inner_grid_polygon(maybe_used_face.size());

                            // 下面这段是和边框切割的代码

                            //sttime += std::chrono::system_clock::now().time_since_epoch().count();

                            for(int maybe_used_face_id=0; maybe_used_face_id < maybe_used_face.size(); maybe_used_face_id++) {

                                vector<K2::Point_3 > triangle_vertex_list = cutting_triangle_by_grid(maybe_used_face[maybe_used_face_id]);

                                // cout <<"************" << endl;

                                K2::Vector_3 origin_direct = cross_product(
                                        (maybe_used_face[maybe_used_face_id][1] - maybe_used_face[maybe_used_face_id][0]),
                                        (maybe_used_face[maybe_used_face_id][2] - maybe_used_face[maybe_used_face_id][0]));

                                sort_by_polar_order(
                                        triangle_vertex_list, origin_direct);
                                //*****************************************************
                                // 修改triangle list
                                face_inner_grid_polygon[maybe_used_face_id] = triangle_vertex_list;

                            }


                            // 究极优化4 核心处搞搞搞！？！！！？！？！
                            //第四关
//                for(int i=0;i<maybe_used_face.size();i++) {
//
//                }
//
//
//                continue;

                            //cout << i <<" : "<< maybe_used_face.size() << endl;
                            vector<K2::Triangle_3 > generated_face_list;

//                vector<pair<K2::Point_2 ,int> > maybe_used_face_2d;
//                for(int i=0;i<maybe_used_face.size();i++) {
//                    maybe_used_face_2d.push_back(K2::)
//                }


//                vector<K2::Triangle_3 > maybe_used_face_sorted = maybe_used_face;
//                std::sort(maybe_used_face_sorted.begin(),maybe_used_face_sorted.end(),[&](K2::Triangle_3 a, K2::Triangle_3 b){
//                    return CGAL::squared_area(a.vertex(0),a.vertex(1),a.vertex(2)) <
//                            CGAL::squared_area(b.vertex(0),b.vertex(1),b.vertex(2));
//                });
//                vector<K2::Point_3>ray_x;
//
//                for(int i=0;i<maybe_used_face_sorted.size();i++){
//                    K2::Point_2 p0(maybe_used_face_sorted[i].vertex(0).y(),maybe_used_face_sorted[i].vertex(0).z());
//                    K2::Point_2 p1(maybe_used_face_sorted[i].vertex(1).y(),maybe_used_face_sorted[i].vertex(1).z());
//                    K2::Point_2 p2(maybe_used_face_sorted[i].vertex(2).y(),maybe_used_face_sorted[i].vertex(2).z());
//                    K2::Triangle_2 convert_2d(p0,p1,p2);
//
//                }

//                continue;

//
//                    int f65id = 1;
//                    for(auto this_face : field_triangles){
//                        fprintf(file5_6, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                                CGAL::to_double(this_face.vertex(0).y()),
//                                CGAL::to_double(this_face.vertex(0).z()));
//                        fprintf(file5_6, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                                CGAL::to_double(this_face.vertex(1).y()),
//                                CGAL::to_double(this_face.vertex(1).z()));
//                        fprintf(file5_6, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                                CGAL::to_double(this_face.vertex(2).y()),
//                                CGAL::to_double(this_face.vertex(2).z()));
//                        fprintf(file5_6,"f %d %d %d\n",f65id,f65id+1,f65id+2);
//                        f65id+=3;
//                    }





                            for(int i=0;i<maybe_used_face.size();i++) { // 处理单片面的切割

                                list<K2::Triangle_3 >now_tri_list;
                                K2::Segment_3 e0((maybe_used_face[i][0]),
                                                 (maybe_used_face[i][1]));
                                K2::Segment_3 e1((maybe_used_face[i][1]),
                                                 (maybe_used_face[i][2]));
                                K2::Segment_3 e2((maybe_used_face[i][2]),
                                                 (maybe_used_face[i][0]));
                                vector<pair<K2::Segment_3,K2::Triangle_3> >cutting_segment;
                                K2::Triangle_3  tri_i((maybe_used_face[i][0]),
                                                      (maybe_used_face[i][1]),
                                                      (maybe_used_face[i][2]));
                                //下面这个for 计算出所有的交线  直接用原来的面去切，不要用inner part 切 ，加快效率


                                // now_tri_list.push_back(field_move_K2_triangle[possible_face_list[i]]);

                                vector<K2::Segment_3>cs;
                                for(auto j : maybe_used_face_seg_cutting[i]){
                                    cs.push_back(j);
                                }

                                //下面这块用德劳内代替

                                //  cout<<face_inner_grid_part[i].size() <<" *** "<< face_inner_grid_polygon[i].size() << endl;

                                K2::Vector_3 origin_normal = CGAL::cross_product(tri_i.vertex(1) -tri_i.vertex(0),
                                                                                 tri_i.vertex(2) -tri_i.vertex(0));
                                // cout << "cs" << cs.size() << endl;
                                vector<vector<K2::Point_3> > cdt_res = CGAL_CDT(face_inner_grid_polygon[i],cs,tri_i);

                                //static int vid = 1;


//FIXME: 在这里启用cdt
                                //*********************

                                //  now_tri_list.clear();

                                for(int j=0;j<cdt_res.size();j++){
                                    K2::Vector_3 this_normal = CGAL::cross_product(cdt_res[j][1] - cdt_res[j][0],
                                                                                   cdt_res[j][2] - cdt_res[j][0]);
                                    if(this_normal * origin_normal < CGAL::Epeck::FT(0)){
                                        swap(cdt_res[j][1],cdt_res[j][2]);
                                    }
                                    now_tri_list.emplace_back(cdt_res[j][0],cdt_res[j][1],cdt_res[j][2]);
                                }

                                //*********************
                                int belong_field_id = maybe_used_face_belong_field[i];



                                //static int f57id = 1;


                                for(auto tri: now_tri_list) {


//                            fprintf(file5_7, "v %lf %lf %lf \n",CGAL::to_double(tri.vertex(0).x()),
//                                    CGAL::to_double(tri.vertex(0).y()),
//                                    CGAL::to_double(tri.vertex(0).z()));
//                            fprintf(file5_7, "v %lf %lf %lf \n",CGAL::to_double(tri.vertex(1).x()),
//                                    CGAL::to_double(tri.vertex(1).y()),
//                                    CGAL::to_double(tri.vertex(1).z()));
//                            fprintf(file5_7, "v %lf %lf %lf \n",CGAL::to_double(tri.vertex(2).x()),
//                                    CGAL::to_double(tri.vertex(2).y()),
//                                    CGAL::to_double(tri.vertex(2).z()));
//                            fprintf(file5_7,"f %d %d %d\n",f57id,f57id+1,f57id+2);
//                            f57id+=3;

                                    // todo : open this code this code is check weather is need generate

                                    //todo :   注意这里没有和外部面判断相交进行裁切 说不定有问题!!!!!!!!!；
/*****************************/
                                    // cout <<"start check "<< endl;

                                    // bool flag = check_in_approximate_field_list(maybe_used_face_field ,CGAL::centroid(tri));;
                                    //  bool flag2 = check_in_approximate_field_list(field_through_list ,CGAL::centroid(tri));
                                    bool flag = false;
                                    //Tree aabb_tree_final_round(field_triangles_final_round.begin(),field_triangles_final_round.end());

                                    //使用aabbtree 代替
                                    K2::Point_3 this_center = CGAL::centroid(tri);
                                    K2::Ray_3 ray(this_center,tri.supporting_plane().orthogonal_vector());
                                    std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections;
                                    aabb_tree.all_intersections(ray,std::back_inserter(intersections));
                                    //vector<bool>cutting_field_id(field_through_list.size(),false);
                                    vector<set<K2::Point_3 > > intersection_v(field_through_list.size());
                                    set<int>is_special;
                                    set<int> positive_side;
                                    for(auto item : intersections) {
                                        if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
                                            //se.insert(*p);
                                            //cout <<"*********"<<endl;
                                            int this_field_belong_id = face_belong_field_all_mp[item.second->id()];
                                            if (belong_field_id != this_field_belong_id){
                                                if( *p == this_center){
                                                    //flag = true;
                                                    positive_side.insert(this_field_belong_id);
                                                    //cout <<"f1" << endl;

                                                }
                                                if(intersection_v[this_field_belong_id].count(*p)){
                                                    if(!is_special.count(this_field_belong_id)){
                                                        if(faces_approximate_field[field_through_list[this_field_belong_id]].in_or_on_field(this_center)) {
                                                            flag = true;
                                                            //cout <<"f2" << endl;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                        is_special.insert(this_field_belong_id);
                                                }
                                                intersection_v[this_field_belong_id].insert(*p);
                                            }
                                        }
                                    }



                                    for(int j=0;j<field_through_list.size();j++){
                                        if(j==belong_field_id)continue;
                                        if(!is_special.count(j) && positive_side.count(j)==0){
                                            if(intersection_v[j].size()%2) {
                                                flag = true;
                                                //cout <<"f3" << endl;
                                                break;
                                            }
                                        }
                                    }
                                    bool flag_positive = false;
                                    bool flag_negative = false;
                                    for(int j=0;j<field_through_list.size();j++){
                                        if(positive_side.count(j) && intersection_v[j].size()%2==0 && tri.supporting_plane().has_on_positive_side(faces_approximate_field[field_through_list[j]].center)){
                                            flag_positive = true;
                                        }
                                    }
                                    if(!flag && flag_positive){
                                        K2::Ray_3 ray_r(this_center,tri.supporting_plane().orthogonal_vector()*(-1));
                                        std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections_r;
                                        aabb_tree.all_intersections(ray,std::back_inserter(intersections_r));
                                        //vector<bool>cutting_field_id(field_through_list.size(),false);
                                        vector<set<K2::Point_3 > > intersection_v_r(field_through_list.size());
                                        set<int>is_special_r;
                                        set<int> negative_side;
                                        for(auto item : intersections_r) {
                                            if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
                                                int this_field_belong_id = face_belong_field_all_mp[item.second->id()];
                                                if (belong_field_id != this_field_belong_id){
                                                    if( *p == this_center){
                                                        negative_side.insert(this_field_belong_id);
                                                    }
                                                    if(intersection_v_r[this_field_belong_id].count(*p)){
                                                        if(!is_special_r.count(this_field_belong_id)){
                                                            if(faces_approximate_field[field_through_list[this_field_belong_id]].in_or_on_field(this_center)) {
                                                                flag = true;
                                                                break;
                                                            }
                                                        }
                                                        else
                                                            is_special_r.insert(this_field_belong_id);
                                                    }
                                                    intersection_v_r[this_field_belong_id].insert(*p);
                                                }
                                            }
                                        }
                                        for(int j=0;j<field_through_list.size();j++){
                                            if(negative_side.count(j) && intersection_v_r[j].size()%2==0  && tri.supporting_plane().has_on_negative_side(faces_approximate_field[field_through_list[j]].center)){
                                                flag_negative = true;
                                            }
                                        }
                                    }
                                    if(flag_positive && flag_negative)
                                        flag = true;


                                    // 这个地方是不是可以去掉底相交呢？？？？？

                                    if(!flag) {
                                        int side_type = //get_side_type(Point_K2_to_iGameVertex(CGAL::centroid(tri)));
#ifdef MODEIN
                                                cgal_polygon->inMesh((CGAL::centroid(tri)));
#else
                                                cgal_polygon->outMesh((CGAL::centroid(tri)));
#endif
                                        if(side_type != 1) {
                                            flag = true;
                                            //cout <<"GGST" << endl;
                                        }
                                    }


                                    auto vvv0 = Point_K2_to_iGameVertex(tri.vertex(0));
                                    auto vvv1 = Point_K2_to_iGameVertex(tri.vertex(1));
                                    auto vvv2 = Point_K2_to_iGameVertex(tri.vertex(2));

                                    auto ddd = (vvv1 - vvv0) % (vvv2- vvv0);
                                    //if(ddd*direct <0)
                                    // swap(vvv1,vvv2);



//                       std::unique_lock<std::mutex>lock(mu);
//                       fprintf(file10,"v %lf %lf %lf\n",vvv0.x(),vvv0.y(),vvv0.z());
//                       fprintf(file10,"v %lf %lf %lf\n",vvv2.x(),vvv2.y(),vvv2.z());
//                       fprintf(file10,"v %lf %lf %lf\n",vvv1.x(),vvv1.y(),vvv1.z());
//                       fprintf(file10,"f %d %d %d\n",vid,vid+1,vid+2);
//                       vid+=3;

                                    if(flag) {
                                        continue;
                                    }

                                    //if(maybe_used_face_source_id[i][0] >= 3 && maybe_used_face_source_id[i][1] >= 3 && maybe_used_face_source_id[i][2] >= 3)
                                    generated_face_list.push_back(K2::Triangle_3(tri.vertex(0),tri.vertex(2),tri.vertex(1)));
//                        else{
//                            static int vid = 1;
//                            std::unique_lock<std::mutex>lock(mu);
//                            fprintf(file10,"v %lf %lf %lf\n",vvv0.x(),vvv0.y(),vvv0.z());
//                            fprintf(file10,"v %lf %lf %lf\n",vvv2.x(),vvv2.y(),vvv2.z());
//                            fprintf(file10,"v %lf %lf %lf\n",vvv1.x(),vvv1.y(),vvv1.z());
//                            fprintf(file10,"f %d %d %d\n",vid,vid+1,vid+2);
//                            vid+=3;
//                        }
/*****************************/
//                        auto vvv0 = Point_K2_to_iGameVertex(tri.vertex(0));
//                        auto vvv1 = Point_K2_to_iGameVertex(tri.vertex(1));
//                        auto vvv2 = Point_K2_to_iGameVertex(tri.vertex(2));
//
//                        auto ddd = (vvv1 - vvv0) % (vvv2- vvv0);
//                        //if(ddd*direct <0)
//                           // swap(vvv1,vvv2);
//
//
//
//                        std::unique_lock<std::mutex>lock(mu);
//                        fprintf(file10,"v %lf %lf %lf\n",vvv0.x(),vvv0.y(),vvv0.z());
//                        fprintf(file10,"v %lf %lf %lf\n",vvv2.x(),vvv2.y(),vvv2.z());
//                        fprintf(file10,"v %lf %lf %lf\n",vvv1.x(),vvv1.y(),vvv1.z());
//                        fprintf(file10,"f %d %d %d\n",vid,vid+1,vid+2);
//                        vid+=3;

                                }
                            }
                            //endtime += std::chrono::system_clock::now().time_since_epoch().count();

//                    std::chrono::time_point<std::chrono::system_clock> end_time = std::chrono::system_clock::now();
//
//
//                    std::chrono::milliseconds diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
//
//
//                    if(diff.count() >10000)
//                        cout <<"id: "<<id<<" field size: "<< field_through_list.size() <<"maybe use " << maybe_used_face.size()<< " gen size: "<< generated_face_list.size()<<" use time "<<" "<< diff.count()<<" skip?"<< skip<<endl;

//                if(generated_face_list.size() ==0 && !skip && field_through_list.size() >30){
//                    auto small  = getGridVertex(each_grid->first,0);
//                    auto big  = getGridVertex(each_grid->first,7);
//                    static int f3_id = 1;
//                    for (int ii = 0; ii < 7; ii++) {
//                        for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
//                            int from = ii;
//                            int to = DirectedGridEdge[ii][jj];
//                            MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
//                            MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
//                            fprintf(file9, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
//                            fprintf(file9, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
//                            fprintf(file9, "l %d %d\n", f3_id, f3_id + 1);
//                            f3_id += 2;
//                        }
//                    }
//                    int f5id = 1;
//                    for(auto this_face : field_triangles_final_round){
//                        fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
//                                CGAL::to_double(this_face.vertex(0).y()),
//                                CGAL::to_double(this_face.vertex(0).z()));
//                        fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
//                                CGAL::to_double(this_face.vertex(1).y()),
//                                CGAL::to_double(this_face.vertex(1).z()));
//                        fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
//                                CGAL::to_double(this_face.vertex(2).y()),
//                                CGAL::to_double(this_face.vertex(2).z()));
//                        fprintf(file5,"f %d %d %d\n",f5id,f5id+1,f5id+2);
//                        f5id+=3;
//                    }
//                    if(field_through_list.size() > 30){
//                        skip = true;
//                        for(auto  each_container_face : container_grid_face){
//                            K2::Triangle_3 tri1(ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]]);
//                            K2::Triangle_3 tri2(ps[each_container_face[2]],ps[each_container_face[3]],ps[each_container_face[0]]);
//                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_1;
//                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_2;
//                            aabb_tree.all_intersections(tri1,std::back_inserter(intersections_1));
//                            aabb_tree.all_intersections(tri2,std::back_inserter(intersections_2));
//                            vector<K2::Segment_3>vs;
//
//                            for(auto i : intersections_1){
//                                if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
//                                    vs.push_back(*s);
//                                }
//                            }
//                            for(auto i : intersections_2){
//                                if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
//                                    vs.push_back(*s);
//                                }
//                            }
//                            auto tmp = CGAL_CDT({ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]],ps[each_container_face[3]]},vs,tri1);
//                            for(auto i:tmp){
//                                //cout <<"skip true" << endl;
//                                if(!check_in_field(K2::Triangle_3(i[0],i[1],i[2]))){
//                                    // cout <<"skip false" << endl;
//                                    skipcc = CGAL::centroid(K2::Triangle_3(i[0],i[1],i[2]));
//                                    skip = false;
//                                    break;
//                                }
//                            }
//
////                    for(K2::Triangle_3 this_face : field_triangles){
////                        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
////                                res_tt = intersection(tri1,this_face);
////                        if (res_tt) {
////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
////                                vs.push_back(*p);
////                            }
////                        }
////                        res_tt = intersection(tri2,this_face);
////                        if (res_tt) {
////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
////                                vs.push_back(*p);
////                            }
////                        }
////                    }
//
//                            if(!skip)
//                                break;
//                        }
//                    }
//                    cout << CGAL::to_double(skipcc.x()) <<" "<< CGAL::to_double(skipcc.y())<<" "<< CGAL::to_double(skipcc.z()) << endl;
//                    exit(0);
//                }

//                vector<MeshKernel::iGameVertex >gen_vertex;
//                vector<vector<int> >gen_face;
                            unique_lock<std::mutex> lock2(mu2);
                            for(auto i : generated_face_list) {
                                each_grid->second.generate_face_list.push_back(i);

                            }
                            lock2.unlock();
                        };
                        dfs(small,big,face_list,0);
                    }


                }
                else {
                    lock.unlock();
                    break;
                }
            }



//            for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
//
//               // if(!(each_grid->first.x == 1 && each_grid->first.y == 0 && each_grid->first.z == 3 ))continue;
//
//                set <MeshKernel::iGameFaceHandle > face_set;
//                vector<int> face_list;
//                for (int j = 0; j < 8; j++) {
//                    grid g = getGridFrameVertex(each_grid->first, j);
//                    unordered_map<grid, GridVertex, grid_hash, grid_equal>::iterator it = frame_grid_mp.find(g);
//                    if (it == frame_grid_mp.end())continue;
//                    for (auto k: it->second.face_list) {
//                        if (!face_set.count(k)) {
//                            face_set.insert(k);
//                            face_list.push_back(k);
//                        }
//                    }
//                }
//
//                K2::Point_3 small = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 0));
//                K2::Point_3 big = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 7));
//                std::function<void(K2::Point_3,K2::Point_3,vector<int>,int)> dfs = [&](K2::Point_3 small,K2::Point_3 big,vector<int>face_list,int depth){
//                    std::vector<K2::Point_3> ps;
//                    for(int i=0;i<8;i++){
//                        ps.push_back(getGridK2Vertex(small,big,i));
//                    }
//                    std::list<K2::Triangle_3 >  frame_faces_list;
//                    for(auto i : each_grid_face_list){
//                        frame_faces_list.emplace_back(ps[i[0]],ps[i[1]],ps[i[2]]);
//                    }
////                    double xxx = CGAL::to_double(small.x());
////                    double yyy = CGAL::to_double(small.y());
////                    double zzz = CGAL::to_double(small.z());
////                    double xx = CGAL::to_double(big.x());
////                    double yy = CGAL::to_double(big.y());
////                    double zz = CGAL::to_double(big.z());
//
//                    Tree frame_aabb_tree(frame_faces_list.begin(),frame_faces_list.end());
//                    CGAL::Polyhedron_3<K2> frame_poly;
//                    PMP::polygon_soup_to_polygon_mesh(ps, each_grid_face_list, frame_poly, CGAL::parameters::all_default());
//                    std::function<bool(K2::Point_3)> vertex_in_frame = [&](K2::Point_3 v){
//                        return small.x() <= v.x() && v.x() <= big.x() &&
//                                small.y() <= v.y() && v.y() <= big.y() &&
//                                small.z() <= v.z() && v.z() <= big.z();
//                    };
//                    vector<K2::Triangle_3 > maybe_used_face;
//                    vector<int> maybe_used_face_belong_field;
//                    vector<vector<K2::Segment_3> > maybe_used_face_seg_cutting;
//
//                    std::vector<MeshKernel::iGameFaceHandle> debug_tet_list_belong_face;
//                    vector<int> field_through_list;
//                    int state = 0;
//                    for (int i: face_list) {
//                        bool useful = false;
//                        for(int j=0;j<0;j++){
//                            if(faces_approximate_field[i].in_field(ps[j])){
//                                state|=(1<<j);
//                            }
//                        }
//                        if(faces_approximate_field[i].in_field(midpoint(K2::Segment_3(ps[0],ps[7])))){
//                            useful = true;
//                        }
//                        if(!useful)
//                        {
//                            if(vertex_in_frame(faces_approximate_field[i].center)){
//                                useful = true;
//                            }
//                        }
//                        if(!useful)
//                            useful = CGAL::Polygon_mesh_processing::do_intersect(*faces_approximate_field[i].poly,frame_poly);
//                        if(useful)
//                            field_through_list.push_back(i);
//                    }
//                    if(state == ((1<<8)-1)) {
//                        cout <<"skip" << endl;
//                        return;
//                    }
//
//                 //   cout <<"field_through_list.size():  " <<field_through_list.size() << endl;
//
//                    function<bool(K2::Triangle_3)>  triangle_through_grid = [&](K2::Triangle_3 tri) {
//                        if(vertex_in_frame(centroid(tri))){
//                            return true;
//                        }
//                        CGAL::Polyhedron_3<K2> this_face;
//                        vector<K2::Point_3> vs{tri.vertex(0),tri.vertex(1),tri.vertex(2)};
//                        PMP::polygon_soup_to_polygon_mesh(vs, face_type_012, this_face, CGAL::parameters::all_default());
//                        return CGAL::Polygon_mesh_processing::do_intersect(frame_poly,this_face);
//                    };
//
//                    //std::unordered_map<std::size_t,int> face_belong_field_mp;
//                    std::list<K2::Triangle_3> field_triangles;
//                    std::unordered_map<std::size_t,vector<int>> face_belong_field_source_id;
//
//                    std::unordered_map<std::size_t,int> face_belong_field_all_mp;
//                    std::list<K2::Triangle_3> field_triangles_all;
//                    std::list<K2::Triangle_3> field_triangles_part_aabbtree;
//                    for (int i=0;i< field_through_list.size();i++) {
//                        for(int j=0;j<faces_approximate_field[field_through_list[i]].bound_face_id.size();j++){
//                            vector<MeshKernel::iGameVertex> tmp{faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][0]],
//                                                                faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][1]],
//                                                                faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][2]]};
//                            K2::Triangle_3 tri_this(iGameVertex_to_Point_K2(tmp[0]),
//                                                    iGameVertex_to_Point_K2(tmp[1]),
//                                                    iGameVertex_to_Point_K2(tmp[2])
//                            );
//                            if(faces_approximate_field[field_through_list[i]].bound_face_id[j][0] >= 3 || faces_approximate_field[field_through_list[i]].bound_face_id[j][1] >= 3 || faces_approximate_field[field_through_list[i]].bound_face_id[j][2] >= 3) { // 去掉原表面
//                                if (faces_approximate_field[field_through_list[i]].bound_face_useful[j]  && triangle_through_grid(tri_this)) {
//                                    field_triangles.push_back(tri_this);
//                                    field_triangles_part_aabbtree.push_back(tri_this);
//                                    //face_belong_field_mp[tri_this.id()] = i;
////                                if(faces_approximate_field[field_through_list[i]].bound_face_id[j][0] >= 3  && faces_approximate_field[field_through_list[i]].bound_face_id[j][1] >= 3 && faces_approximate_field[field_through_list[i]].bound_face_id[j][2] >= 3) {
////                                    static int f52id = 1;
////                                    fprintf(file5_2, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(0).x()),
////                                            CGAL::to_double(tri_this.vertex(0).y()),
////                                            CGAL::to_double(tri_this.vertex(0).z()));
////                                    fprintf(file5_2, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(1).x()),
////                                            CGAL::to_double(tri_this.vertex(1).y()),
////                                            CGAL::to_double(tri_this.vertex(1).z()));
////                                    fprintf(file5_2, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(2).x()),
////                                            CGAL::to_double(tri_this.vertex(2).y()),
////                                            CGAL::to_double(tri_this.vertex(2).z()));
////                                    fprintf(file5_2, "f %d %d %d\n", f52id, f52id + 1, f52id + 2);
////                                    f52id += 3;
////                                }
////                                else{
////                                    static int f53id = 1;
////                                    fprintf(file5_3, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(0).x()),
////                                            CGAL::to_double(tri_this.vertex(0).y()),
////                                            CGAL::to_double(tri_this.vertex(0).z()));
////                                    fprintf(file5_3, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(1).x()),
////                                            CGAL::to_double(tri_this.vertex(1).y()),
////                                            CGAL::to_double(tri_this.vertex(1).z()));
////                                    fprintf(file5_3, "v %lf %lf %lf \n", CGAL::to_double(tri_this.vertex(2).x()),
////                                            CGAL::to_double(tri_this.vertex(2).y()),
////                                            CGAL::to_double(tri_this.vertex(2).z()));
////                                    fprintf(file5_3, "f %d %d %d\n", f53id, f53id + 1, f53id + 2);
////                                    f53id += 3;
////                                }
//
//                                    face_belong_field_source_id[tri_this.id()] = vector<int>{faces_approximate_field[field_through_list[i]].bound_face_id[j][0],
//                                                                                             faces_approximate_field[field_through_list[i]].bound_face_id[j][1],
//                                                                                             faces_approximate_field[field_through_list[i]].bound_face_id[j][2]};
//
//                                }
//                            }
//                            else
//                                field_triangles_part_aabbtree.push_back(tri_this);
//                            field_triangles_all.push_back(tri_this);
//                            face_belong_field_all_mp[tri_this.id()] = i;
//                        }
//                    }
////                    cout <<"field_triangles.size():  " <<field_triangles.size() << endl;
////
////                    int f51id = 1;
////                    for(auto this_face : field_triangles_all){
////                        fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                                CGAL::to_double(this_face.vertex(0).y()),
////                                CGAL::to_double(this_face.vertex(0).z()));
////                        fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                                CGAL::to_double(this_face.vertex(1).y()),
////                                CGAL::to_double(this_face.vertex(1).z()));
////                        fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                                CGAL::to_double(this_face.vertex(2).y()),
////                                CGAL::to_double(this_face.vertex(2).z()));
////                        fprintf(file5_1,"f %d %d %d\n",f51id,f51id+1,f51id+2);
////                        f51id+=3;
////                    }
////
//                    int f52id = 1;
//                    int f53id = 1;
//                    int f54id = 1;
//                    for (int i=0;i< field_through_list.size();i++) {
//                        for (int j = 0; j < faces_approximate_field[field_through_list[i]].bound_face_id.size(); j++) {
//                            vector<MeshKernel::iGameVertex> tmp{faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][0]],
//                                                                faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][1]],
//                                                                faces_approximate_field[field_through_list[i]].bound_face_vertex[faces_approximate_field[field_through_list[i]].bound_face_id[j][2]]};
//                            K2::Triangle_3 this_face(iGameVertex_to_Point_K2(tmp[0]),
//                                                    iGameVertex_to_Point_K2(tmp[1]),
//                                                    iGameVertex_to_Point_K2(tmp[2])
//                            );
////                            if (faces_approximate_field[field_through_list[i]].bound_face_useful[j]  && triangle_through_grid(this_face)) {
////                                cout << faces_approximate_field[field_through_list[i]].bound_face_id[j][0] <<" "<< faces_approximate_field[field_through_list[i]].bound_face_id[j][1] <<" "<<faces_approximate_field[field_through_list[i]].bound_face_id[j][2] << endl;
////                                                        fprintf(file5_2, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                                                                CGAL::to_double(this_face.vertex(0).y()),
////                                                                CGAL::to_double(this_face.vertex(0).z()));
////                                                        fprintf(file5_2, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                                                                CGAL::to_double(this_face.vertex(1).y()),
////                                                                CGAL::to_double(this_face.vertex(1).z()));
////                                                        fprintf(file5_2, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                                                                CGAL::to_double(this_face.vertex(2).y()),
////                                                                CGAL::to_double(this_face.vertex(2).z()));
////                                                        fprintf(file5_2,"f %d %d %d\n",f52id,f52id+1,f52id+2);
////                                                        f52id+=3;
////                            }
////                            else if(faces_approximate_field[field_through_list[i]].bound_face_useful[j]){
////                                cout << faces_approximate_field[field_through_list[i]].bound_face_id[j][0] <<" "<< faces_approximate_field[field_through_list[i]].bound_face_id[j][1] <<" "<<faces_approximate_field[field_through_list[i]].bound_face_id[j][2] << endl;
////                                fprintf(file5_3, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                                        CGAL::to_double(this_face.vertex(0).y()),
////                                        CGAL::to_double(this_face.vertex(0).z()));
////                                fprintf(file5_3, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                                        CGAL::to_double(this_face.vertex(1).y()),
////                                        CGAL::to_double(this_face.vertex(1).z()));
////                                fprintf(file5_3, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                                        CGAL::to_double(this_face.vertex(2).y()),
////                                        CGAL::to_double(this_face.vertex(2).z()));
////                                fprintf(file5_3,"f %d %d %d\n",f53id,f53id+1,f53id+2);
////                                f53id+=3;
////                            }
//                        }
//                    }
////                    for(auto this_face : field_triangles_all){
////                        fprintf(file5_4, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                                CGAL::to_double(this_face.vertex(0).y()),
////                                CGAL::to_double(this_face.vertex(0).z()));
////                        fprintf(file5_4, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                                CGAL::to_double(this_face.vertex(1).y()),
////                                CGAL::to_double(this_face.vertex(1).z()));
////                        fprintf(file5_4, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                                CGAL::to_double(this_face.vertex(2).y()),
////                                CGAL::to_double(this_face.vertex(2).z()));
////                        fprintf(file5_4,"f %d %d %d\n",f54id,f54id+1,f54id+2);
////                        f54id+=3;
////                    }
//
//
//
//                    if(field_triangles.size() ==0)return ;
//
////                    int f55id = 1;
////                    for(auto this_face : field_triangles){
////                        fprintf(file5_5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                                CGAL::to_double(this_face.vertex(0).y()),
////                                CGAL::to_double(this_face.vertex(0).z()));
////                        fprintf(file5_5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                                CGAL::to_double(this_face.vertex(1).y()),
////                                CGAL::to_double(this_face.vertex(1).z()));
////                        fprintf(file5_5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                                CGAL::to_double(this_face.vertex(2).y()),
////                                CGAL::to_double(this_face.vertex(2).z()));
////                        fprintf(file5_5,"f %d %d %d\n",f55id,f55id+1,f55id+2);
////                        f55id+=3;
////                    }
//
//                    //cout <<"field_triangles.size():  " <<field_triangles.size() << endl;
//
//
//                    Tree aabb_tree(field_triangles_all.begin(), field_triangles_all.end());
//                    Tree aabb_tree_part(field_triangles_part_aabbtree.begin(), field_triangles_part_aabbtree.end());
//
//
//
//
//
//                    //Tree aabb_tree_cross_usage;
//
//                    /*for(auto i : each_grid_face_list){
//                        frame_faces_list.emplace_back(ps[i[0]],ps[i[1]],ps[i[2]]);
//                    }*/
//
//
//
//                    function<bool(K2::Triangle_3)> check_in_field = [&](K2::Triangle_3 tri){
//                        bool flag = false;
//                        K2::Point_3 this_center = CGAL::centroid(tri);
//                        K2::Ray_3 ray(this_center,tri.supporting_plane().orthogonal_vector());
//                        std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections;
//                        aabb_tree.all_intersections(ray,std::back_inserter(intersections));
//                        //vector<bool>cutting_field_id(field_through_list.size(),false);
//                        vector<set<K2::Point_3 > > intersection_v(field_through_list.size());
//                        set<int>is_special;
//                        set<int> positive_side;
//                        for(auto item : intersections) {
//                            if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
//                                //se.insert(*p);
//                                //cout <<"*********"<<endl;
//                                int this_field_belong_id = face_belong_field_all_mp[item.second->id()];
//
//                                if( *p == this_center){
//                                    //flag = true;
//                                    positive_side.insert(this_field_belong_id);
//                                    //cout <<"f1" << endl;
//
//                                }
//                                if(intersection_v[this_field_belong_id].count(*p)){
//                                    if(!is_special.count(this_field_belong_id)){
//                                        if(faces_approximate_field[field_through_list[this_field_belong_id]].in_or_on_field(this_center)) {
//                                            flag = true;
//                                            //cout <<"f2" << endl;
//                                            break;
//                                        }
//                                    }
//                                    else
//                                        is_special.insert(this_field_belong_id);
//                                }
//                                intersection_v[this_field_belong_id].insert(*p);
//
//                            }
//                        }
//
//
//
//                        for(int j=0;j<field_through_list.size();j++){
//                            if(!is_special.count(j) && positive_side.count(j)==0){
//                                if(intersection_v[j].size()%2) {
//                                    flag = true;
//                                    //cout <<"f3" << endl;
//                                    break;
//                                }
//                            }
//                        }
//                        bool flag_positive = false;
//                        bool flag_negative = false;
//                        for(int j=0;j<field_through_list.size();j++){
//                            if(positive_side.count(j) && intersection_v[j].size()%2==0){
//                                flag_positive = true;
//                            }
//                        }
//                        if(!flag && flag_positive){
//                            K2::Ray_3 ray_r(this_center,tri.supporting_plane().orthogonal_vector()*(-1));
//                            std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections_r;
//                            aabb_tree.all_intersections(ray,std::back_inserter(intersections_r));
//                            //vector<bool>cutting_field_id(field_through_list.size(),false);
//                            vector<set<K2::Point_3 > > intersection_v_r(field_through_list.size());
//                            set<int>is_special_r;
//                            set<int> negative_side;
//                            for(auto item : intersections_r) {
//                                if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
//                                    int this_field_belong_id = face_belong_field_all_mp[item.second->id()];
//
//                                    if( *p == this_center){
//                                        negative_side.insert(this_field_belong_id);
//                                    }
//                                    if(intersection_v_r[this_field_belong_id].count(*p)){
//                                        if(!is_special_r.count(this_field_belong_id)){
//                                            if(faces_approximate_field[field_through_list[this_field_belong_id]].in_or_on_field(this_center)) {
//                                                flag = true;
//                                                break;
//                                            }
//                                        }
//                                        else
//                                            is_special_r.insert(this_field_belong_id);
//                                    }
//                                    intersection_v_r[this_field_belong_id].insert(*p);
//
//                                }
//                            }
//                            for(int j=0;j<field_through_list.size();j++){
//                                if(negative_side.count(j) && intersection_v_r[j].size()%2==0){
//                                    flag_negative = true;
//                                }
//                            }
//                        }
//                        if(flag_positive && flag_negative)
//                            flag = true;
//                        return flag;
//                    };
//
//
////**************************************oct dfs skip
////                    bool skip = false;
////                    K2::Point_3 skipcc;
////                    if(field_through_list.size() > 28 ) {
////                        skip = true;
////                        for(auto  each_container_face : container_grid_face){
////                            K2::Triangle_3 tri1(ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]]);
////                            K2::Triangle_3 tri2(ps[each_container_face[2]],ps[each_container_face[3]],ps[each_container_face[0]]);
////                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_1;
////                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_2;
////                            aabb_tree.all_intersections(tri1,std::back_inserter(intersections_1));
////                            aabb_tree.all_intersections(tri2,std::back_inserter(intersections_2));
////                            vector<K2::Segment_3>vs;
////
////                            for(auto i : intersections_1){
////                                if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
////                                    vs.push_back(*s);
////                                }
////                            }
////                            for(auto i : intersections_2){
////                                if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
////                                    vs.push_back(*s);
////                                }
////                            }
////                            auto tmp = CGAL_CDT({ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]],ps[each_container_face[3]]},vs,tri1);
////                            for(auto i:tmp){
////                                //cout <<"skip true" << endl;
////                                if(!check_in_field(K2::Triangle_3(i[0],i[1],i[2]))){
////                                    // cout <<"skip false" << endl;
////                                    skipcc = CGAL::centroid(K2::Triangle_3(i[0],i[1],i[2]));
////                                    skip = false;
////                                    break;
////                                }
////                            }
////
//////                    for(K2::Triangle_3 this_face : field_triangles){
//////                        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
//////                                res_tt = intersection(tri1,this_face);
//////                        if (res_tt) {
//////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
//////                                vs.push_back(*p);
//////                            }
//////                        }
//////                        res_tt = intersection(tri2,this_face);
//////                        if (res_tt) {
//////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
//////                                vs.push_back(*p);
//////                            }
//////                        }
//////                    }
////
////                            if(!skip)
////                                break;
////                        }
////                    }
////                    if(skip) {
////                        cout <<"skip"<<endl;
////                        return;
////                    }
//
//
//                   // cout << "depth" << depth<< endl;
//                    if(field_through_list.size() > 25+depth*2.5 && depth<3){ ///
//                        cout <<"dfs!!" << endl;
//                        K2::Point_3 center = midpoint(K2::Segment_3(small,big));
//                        K2::Vector_3 delta = center - small;
//                        vector<K2::Point_3> new_small;
//                        vector<K2::Point_3> new_big;
//                        for(int i=0;i<8;i++){
//                            new_small.push_back(getGridK2Vertex(small,center,i));
//                        }
//                        for(int i=0;i<8;i++){
//                            new_big.push_back(new_small[i]+delta);
//                        }
//                        for(int i=0;i<8;i++){
//                            dfs(new_small[i] ,new_big[i],field_through_list,depth+1);
//                        }
//                        return ;
//                    }
//                    if(depth == 4){
//                        bool skip = false;
//                        K2::Point_3 skipcc;
//                        skip = true;
//                        for(auto  each_container_face : container_grid_face){
//                            K2::Triangle_3 tri1(ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]]);
//                            K2::Triangle_3 tri2(ps[each_container_face[2]],ps[each_container_face[3]],ps[each_container_face[0]]);
//                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_1;
//                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_2;
//                            aabb_tree.all_intersections(tri1,std::back_inserter(intersections_1));
//                            aabb_tree.all_intersections(tri2,std::back_inserter(intersections_2));
//                            vector<K2::Segment_3>vs;
//
//                            for(auto i : intersections_1){
//                                if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
//                                    vs.push_back(*s);
//                                }
//                            }
//                            for(auto i : intersections_2){
//                                if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
//                                    vs.push_back(*s);
//                                }
//                            }
//                            auto tmp = CGAL_CDT({ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]],ps[each_container_face[3]]},vs,tri1);
//                            for(auto i:tmp){
//                                //cout <<"skip true" << endl;
//                                if(!check_in_field(K2::Triangle_3(i[0],i[1],i[2]))){
//                                    // cout <<"skip false" << endl;
//                                    skipcc = CGAL::centroid(K2::Triangle_3(i[0],i[1],i[2]));
//                                    skip = false;
//                                    break;
//                                }
//                            }
//
////                    for(K2::Triangle_3 this_face : field_triangles){
////                        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
////                                res_tt = intersection(tri1,this_face);
////                        if (res_tt) {
////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
////                                vs.push_back(*p);
////                            }
////                        }
////                        res_tt = intersection(tri2,this_face);
////                        if (res_tt) {
////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
////                                vs.push_back(*p);
////                            }
////                        }
////                    }
//
//                            if(!skip)
//                                break;
//                        }
//
//                        if(skip) {
//                            cout <<"skip"<<endl;
//                            return;
//                        }
//                        cout <<"50!!" << endl;
//                    }
//                    //cout << "field_through_list.size()" << field_through_list.size()<< endl;
//
////                    if(field_through_list.size() > 20){
////                        cout << each_grid->first.x<<" "<< each_grid->first.y<<" "<< each_grid->first.z << endl;
////                        exit(0);
////
////                    }
//
//
//
//
////                    if(skip) {
////                        int f5id = 1;
////                        FILE * file5_1 = fopen( (input_filename + "_5_1.obj").c_str(), "w");
////                        for(auto this_face : field_triangles_all){
////                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                                    CGAL::to_double(this_face.vertex(0).y()),
////                                    CGAL::to_double(this_face.vertex(0).z()));
////                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                                    CGAL::to_double(this_face.vertex(1).y()),
////                                    CGAL::to_double(this_face.vertex(1).z()));
////                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                                    CGAL::to_double(this_face.vertex(2).y()),
////                                    CGAL::to_double(this_face.vertex(2).z()));
////                            fprintf(file5_1,"f %d %d %d\n",f5id,f5id+1,f5id+2);
////                            f5id+=3;
////                        }
////                        FILE * file10 = fopen( (input_filename + "_10.obj").c_str(), "w");
////                        int f3_id = 1;
////                        for (int ii = 0; ii < 7; ii++) {
////                            for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
////                                int from = ii;
////                                int to = DirectedGridEdge[ii][jj];
////                                MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
////                                MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
////                                fprintf(file10, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
////                                fprintf(file10, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
////                                fprintf(file10, "l %d %d\n", f3_id, f3_id + 1);
////                                f3_id += 2;
////                            }
////                        }
////                        fclose(file5_1);
////                        fclose(file10);
////                        cout << "skip" << endl;
////                        int xxx;
////                        cin>>xxx;
////                    }
////                    else{
////                        int f5id = 1;
////                        FILE * file5_1 = fopen( (input_filename + "_5_1.obj").c_str(), "w");
////                        for(auto this_face : field_triangles_all){
////                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                                    CGAL::to_double(this_face.vertex(0).y()),
////                                    CGAL::to_double(this_face.vertex(0).z()));
////                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                                    CGAL::to_double(this_face.vertex(1).y()),
////                                    CGAL::to_double(this_face.vertex(1).z()));
////                            fprintf(file5_1, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                                    CGAL::to_double(this_face.vertex(2).y()),
////                                    CGAL::to_double(this_face.vertex(2).z()));
////                            fprintf(file5_1,"f %d %d %d\n",f5id,f5id+1,f5id+2);
////                            f5id+=3;
////                        }
////                        FILE * file10 = fopen( (input_filename + "_10.obj").c_str(), "w");
////                        int f3_id = 1;
////                        for (int ii = 0; ii < 7; ii++) {
////                            for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
////                                int from = ii;
////                                int to = DirectedGridEdge[ii][jj];
////                                MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
////                                MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
////                                fprintf(file10, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
////                                fprintf(file10, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
////                                fprintf(file10, "l %d %d\n", f3_id, f3_id + 1);
////                                f3_id += 2;
////                            }
////                        }
////                        fclose(file5_1);
////                        fclose(file10);
////                        cout << "no skip" << endl;
////                        int xxx;
////                        cin>>xxx;
////                    }
//
//
//
//
//
//
//                    std::vector<std::vector<int> >maybe_used_face_source_id;
//
//
//
//                    for(auto this_face : field_triangles){
//                        std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
//                        aabb_tree.all_intersections(this_face,std::back_inserter(intersections));
//                        int this_face_belong_id = face_belong_field_all_mp[this_face.id()];
//                        vector<bool>cutting_field_id(field_through_list.size(),false);
//                        vector<K2::Segment_3> segment_cutting;
//                        for(auto i : intersections){
//                            //   cout << "count exist???" << face_belong_field_all_mp.count(i.second->id()) << endl;
//                            int this_field_belong_id  = face_belong_field_all_mp[i.second->id()];
//                            if(this_field_belong_id != this_face_belong_id /*&& !cutting_field_id[this_field_belong_id]*/){
//                                //   cout<< "interse:"<<this_face.id()<<"  "<< this_face_belong_id<<" "<<this_field_belong_id << endl;
//                                if(const K2::Point_3* p = boost::get<K2::Point_3>(&(i.first))){
//
//                                }
//                                else if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
//                                    bool same_edge = false;
//                                    for (K2::Segment_3 edge: {K2::Segment_3(this_face.vertex(0),this_face.vertex(1)),
//                                                              K2::Segment_3(this_face.vertex(1),this_face.vertex(2)),
//                                                              K2::Segment_3(this_face.vertex(2),this_face.vertex(0))}) {
//                                        if (segment_in_line(edge, *s))
//                                            same_edge = true;
//                                    }
//                                    if (!same_edge) {
//                                        cutting_field_id[this_field_belong_id] = true;
//                                        //segment_cutting.push_back(*s);
//                                    }
//                                }
//                                else if(const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&(i.first))) {
//                                    int cnt=0;
//                                    for(int l = 0 ;l<3;l++){
//                                        for(int m = 0 ;m<3;m++){
//                                            if(CGAL::squared_distance(this_face.vertex(l),t->vertex(m)) == CGAL::Epeck::FT(0)){
//                                                cnt++;
//                                                break;
//                                            }
//                                        }
//                                    }
//                                    if(cnt!=3)
//                                        cutting_field_id[this_field_belong_id] = true;
//                                }
//                                else{
//                                    cutting_field_id[this_field_belong_id] = true;
//                                }
//                            }
//                            // cout <<"cutting_field_id[this_field_belong_id]: " <<cutting_field_id[this_field_belong_id] << endl;
//                        }
//                        bool useless = false;
//                        //int ofid=-1;
//                        for(int i = 0; i< cutting_field_id.size();i++){
//                            if(!cutting_field_id[i] && i!=this_face_belong_id ){
//                                if(faces_approximate_field[field_through_list[i]].in_field(CGAL::centroid(this_face))) {
//                                    // cout << i <<" ?? "<< this_face_belong_id << endl;
//                                    //ofid = field_through_list[i];
//                                    useless = true;
//                                    break;
//                                }
//                            }
//                            if(useless)
//                                break;
//                        }
//                        if(!useless) {
//                            maybe_used_face.push_back(this_face);
//                            maybe_used_face_belong_field.push_back(this_face_belong_id);
//
//                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_part;
//                            aabb_tree_part.all_intersections(this_face,std::back_inserter(intersections_part));
//                            for(auto ii : intersections_part) {
//                                //cout << "count exist???" << face_belong_field_all_mp.count(i.second->id()) << endl;
//                                int this_field_belong_id = face_belong_field_all_mp[ii.second->id()];
//                                if (this_field_belong_id != this_face_belong_id) {
//                                    if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&(ii.first))) {
//                                        segment_cutting.push_back(*s);
//                                    }
//                                }
//                            }
//                            maybe_used_face_seg_cutting.push_back(segment_cutting);
//                            maybe_used_face_source_id.push_back(face_belong_field_source_id[this_face.id()]);
//                        }
////                   else{ //TODO :DEBUG CODE
////                       static int f5id = 1;
////                       cout <<"f5id" <<f5id<<" "<<this_face.id()<<"  "<<this_face_belong_id <<" "<< ofid<<endl;
////                       FILE *file5_a = fopen( (input_filename + "_5.a"+std::to_string(f5id)+".obj").c_str(), "w");
////                       FILE *file5_b = fopen( (input_filename + "_5.b"+std::to_string(f5id)+".obj").c_str(), "w");
////                       fprintf(file5_a, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                               CGAL::to_double(this_face.vertex(0).y()),
////                               CGAL::to_double(this_face.vertex(0).z()));
////                       fprintf(file5_a, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                               CGAL::to_double(this_face.vertex(1).y()),
////                               CGAL::to_double(this_face.vertex(1).z()));
////                       fprintf(file5_a, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                               CGAL::to_double(this_face.vertex(2).y()),
////                               CGAL::to_double(this_face.vertex(2).z()));
////                       fprintf(file5_a,"f %d %d %d\n",1,1+1,1+2);
////                       int bid = 1;
////                       for (int k = 0; k < faces_approximate_field[ofid].bound_face_id.size(); k++) {
////                           vector<MeshKernel::iGameVertex> tmp{
////                                   faces_approximate_field[ofid].bound_face_vertex[faces_approximate_field[ofid].bound_face_id[k][0]],
////                                   faces_approximate_field[ofid].bound_face_vertex[faces_approximate_field[ofid].bound_face_id[k][1]],
////                                   faces_approximate_field[ofid].bound_face_vertex[faces_approximate_field[ofid].bound_face_id[k][2]]};
////                           K2::Triangle_3 tri_this(iGameVertex_to_Point_K2(tmp[0]),
////                                                   iGameVertex_to_Point_K2(tmp[1]),
////                                                   iGameVertex_to_Point_K2(tmp[2])
////                           );
////                           fprintf(file5_b, "v %lf %lf %lf \n",CGAL::to_double(tri_this.vertex(0).x()),
////                                   CGAL::to_double(tri_this.vertex(0).y()),
////                                   CGAL::to_double(tri_this.vertex(0).z()));
////                           fprintf(file5_b, "v %lf %lf %lf \n",CGAL::to_double(tri_this.vertex(1).x()),
////                                   CGAL::to_double(tri_this.vertex(1).y()),
////                                   CGAL::to_double(tri_this.vertex(1).z()));
////                           fprintf(file5_b, "v %lf %lf %lf \n",CGAL::to_double(tri_this.vertex(2).x()),
////                                   CGAL::to_double(tri_this.vertex(2).y()),
////                                   CGAL::to_double(tri_this.vertex(2).z()));
////                           fprintf(file5_b,"f %d %d %d\n",bid,bid+1,bid+2);
////                           bid +=3;
////                       }
////                       f5id++;
////                   }
//                    }
//
//
////                    int f511id = 1;
////                    for(auto this_face : maybe_used_face){
////                        fprintf(file5_11, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                                CGAL::to_double(this_face.vertex(0).y()),
////                                CGAL::to_double(this_face.vertex(0).z()));
////                        fprintf(file5_11, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                                CGAL::to_double(this_face.vertex(1).y()),
////                                CGAL::to_double(this_face.vertex(1).z()));
////                        fprintf(file5_11, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                                CGAL::to_double(this_face.vertex(2).y()),
////                                CGAL::to_double(this_face.vertex(2).z()));
////                        fprintf(file5_11,"f %d %d %d\n",f511id,f511id+1,f511id+2);
////                        f511id+=3;
////                    }
//
////                int f5id = 1;
////                for(auto this_face : maybe_used_face){
////                    fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                            CGAL::to_double(this_face.vertex(0).y()),
////                            CGAL::to_double(this_face.vertex(0).z()));
////                    fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                            CGAL::to_double(this_face.vertex(1).y()),
////                            CGAL::to_double(this_face.vertex(1).z()));
////                    fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                            CGAL::to_double(this_face.vertex(2).y()),
////                            CGAL::to_double(this_face.vertex(2).z()));
////                    fprintf(file5,"f %d %d %d\n",f5id,f5id+1,f5id+2);
////                    f5id+=3;
////                }
//
//
//                    //   continue;
//
//                    // 第三关
//
//                    std::function<vector<K2::Point_3>(K2::Triangle_3) > cutting_triangle_by_grid = [&](K2::Triangle_3 this_face) {
//                        vector<K2::Point_3 > ret;
//                        set<K2::Point_3> se;
//                        std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
//                        frame_aabb_tree.all_intersections(this_face,std::back_inserter(intersections));
//                        for(int i=0;i<3;i++){
//                            if(vertex_in_frame(this_face.vertex(i))){
//                                se.insert(this_face.vertex(i));
//                            }
//                        }
//
//                        for(auto i : intersections) {
//                            if(const K2::Point_3* p = boost::get<K2::Point_3>(&(i.first))){
//                                se.insert(*p);
//                            }
//                            if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&(i.first))) {
//                                se.insert(s->vertex(0));
//                                se.insert(s->vertex(1));
//                            }
//                            else if(const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&(i.first))){
//                                se.insert(t->vertex(0));
//                                se.insert(t->vertex(1));
//                                se.insert(t->vertex(2));
//                            }
//                            else if(const std::vector<K2::Point_3> *v = boost::get<std::vector<K2::Point_3>>(&(i.first))){
//                                for(const auto& ii : *v)
//                                    se.insert(ii);
//                            }
//                        }
//                        for(const auto& i : se)
//                            ret.push_back(i);
//
//                        return ret;
//                    };
//
//                    vector<K2::Triangle_3 >possible_face_inner_part;
//
//                    // cout <<"********1****" << endl;
//                    vector<vector<K2::Point_3> > face_inner_grid_polygon(maybe_used_face.size());
//
//                    // 下面这段是和边框切割的代码
//                    for(int maybe_used_face_id=0; maybe_used_face_id < maybe_used_face.size(); maybe_used_face_id++) {
//
//                        vector<K2::Point_3 > triangle_vertex_list = cutting_triangle_by_grid(maybe_used_face[maybe_used_face_id]);
//
//                        // cout <<"************" << endl;
//
//                        K2::Vector_3 origin_direct = cross_product(
//                                (maybe_used_face[maybe_used_face_id][1] - maybe_used_face[maybe_used_face_id][0]),
//                                (maybe_used_face[maybe_used_face_id][2] - maybe_used_face[maybe_used_face_id][0]));
//
//                        sort_by_polar_order(
//                                triangle_vertex_list, origin_direct);
//                        //*****************************************************
//                        // 修改triangle list
//                        face_inner_grid_polygon[maybe_used_face_id] = triangle_vertex_list;
//
//                    }
//
//
//                    // 究极优化4 核心处搞搞搞！？！！！？！？！
//                    //第四关
////                for(int i=0;i<maybe_used_face.size();i++) {
////
////                }
////
////
////                continue;
//
//                    //cout << i <<" : "<< maybe_used_face.size() << endl;
//                    vector<K2::Triangle_3 > generated_face_list;
//
////                vector<pair<K2::Point_2 ,int> > maybe_used_face_2d;
////                for(int i=0;i<maybe_used_face.size();i++) {
////                    maybe_used_face_2d.push_back(K2::)
////                }
//
//
////                vector<K2::Triangle_3 > maybe_used_face_sorted = maybe_used_face;
////                std::sort(maybe_used_face_sorted.begin(),maybe_used_face_sorted.end(),[&](K2::Triangle_3 a, K2::Triangle_3 b){
////                    return CGAL::squared_area(a.vertex(0),a.vertex(1),a.vertex(2)) <
////                            CGAL::squared_area(b.vertex(0),b.vertex(1),b.vertex(2));
////                });
////                vector<K2::Point_3>ray_x;
////
////                for(int i=0;i<maybe_used_face_sorted.size();i++){
////                    K2::Point_2 p0(maybe_used_face_sorted[i].vertex(0).y(),maybe_used_face_sorted[i].vertex(0).z());
////                    K2::Point_2 p1(maybe_used_face_sorted[i].vertex(1).y(),maybe_used_face_sorted[i].vertex(1).z());
////                    K2::Point_2 p2(maybe_used_face_sorted[i].vertex(2).y(),maybe_used_face_sorted[i].vertex(2).z());
////                    K2::Triangle_2 convert_2d(p0,p1,p2);
////
////                }
//
////                continue;
//
////
////                    int f65id = 1;
////                    for(auto this_face : field_triangles){
////                        fprintf(file5_6, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                                CGAL::to_double(this_face.vertex(0).y()),
////                                CGAL::to_double(this_face.vertex(0).z()));
////                        fprintf(file5_6, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                                CGAL::to_double(this_face.vertex(1).y()),
////                                CGAL::to_double(this_face.vertex(1).z()));
////                        fprintf(file5_6, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                                CGAL::to_double(this_face.vertex(2).y()),
////                                CGAL::to_double(this_face.vertex(2).z()));
////                        fprintf(file5_6,"f %d %d %d\n",f65id,f65id+1,f65id+2);
////                        f65id+=3;
////                    }
//
//
//
//
//
//                    for(int i=0;i<maybe_used_face.size();i++) { // 处理单片面的切割
//
//                        list<K2::Triangle_3 >now_tri_list;
//                        K2::Segment_3 e0((maybe_used_face[i][0]),
//                                         (maybe_used_face[i][1]));
//                        K2::Segment_3 e1((maybe_used_face[i][1]),
//                                         (maybe_used_face[i][2]));
//                        K2::Segment_3 e2((maybe_used_face[i][2]),
//                                         (maybe_used_face[i][0]));
//                        vector<pair<K2::Segment_3,K2::Triangle_3> >cutting_segment;
//                        K2::Triangle_3  tri_i((maybe_used_face[i][0]),
//                                              (maybe_used_face[i][1]),
//                                              (maybe_used_face[i][2]));
//                        //下面这个for 计算出所有的交线  直接用原来的面去切，不要用inner part 切 ，加快效率
//
//
//                        // now_tri_list.push_back(field_move_K2_triangle[possible_face_list[i]]);
//
//                        vector<K2::Segment_3>cs;
//                        for(auto j : maybe_used_face_seg_cutting[i]){
//                            cs.push_back(j);
//                        }
//
//                        //下面这块用德劳内代替
//
//                        //  cout<<face_inner_grid_part[i].size() <<" *** "<< face_inner_grid_polygon[i].size() << endl;
//
//                        K2::Vector_3 origin_normal = CGAL::cross_product(tri_i.vertex(1) -tri_i.vertex(0),
//                                                                         tri_i.vertex(2) -tri_i.vertex(0));
//                        // cout << "cs" << cs.size() << endl;
//                        vector<vector<K2::Point_3> > cdt_res = CGAL_CDT(face_inner_grid_polygon[i],cs,tri_i);
//
//                        //static int vid = 1;
//
//
////FIXME: 在这里启用cdt
//                        //*********************
//
//                        //  now_tri_list.clear();
//
//                        for(int j=0;j<cdt_res.size();j++){
//                            K2::Vector_3 this_normal = CGAL::cross_product(cdt_res[j][1] - cdt_res[j][0],
//                                                                           cdt_res[j][2] - cdt_res[j][0]);
//                            if(this_normal * origin_normal < CGAL::Epeck::FT(0)){
//                                swap(cdt_res[j][1],cdt_res[j][2]);
//                            }
//                            now_tri_list.emplace_back(cdt_res[j][0],cdt_res[j][1],cdt_res[j][2]);
//                        }
//
//                        //*********************
//                        int belong_field_id = maybe_used_face_belong_field[i];
//
//
//
//                        //static int f57id = 1;
//
//
//                        for(auto tri: now_tri_list) {
//
//
////                            fprintf(file5_7, "v %lf %lf %lf \n",CGAL::to_double(tri.vertex(0).x()),
////                                    CGAL::to_double(tri.vertex(0).y()),
////                                    CGAL::to_double(tri.vertex(0).z()));
////                            fprintf(file5_7, "v %lf %lf %lf \n",CGAL::to_double(tri.vertex(1).x()),
////                                    CGAL::to_double(tri.vertex(1).y()),
////                                    CGAL::to_double(tri.vertex(1).z()));
////                            fprintf(file5_7, "v %lf %lf %lf \n",CGAL::to_double(tri.vertex(2).x()),
////                                    CGAL::to_double(tri.vertex(2).y()),
////                                    CGAL::to_double(tri.vertex(2).z()));
////                            fprintf(file5_7,"f %d %d %d\n",f57id,f57id+1,f57id+2);
////                            f57id+=3;
//
//                            // todo : open this code this code is check weather is need generate
//
//                            //todo :   注意这里没有和外部面判断相交进行裁切 说不定有问题!!!!!!!!!；
///*****************************/
//                            // cout <<"start check "<< endl;
//
//                            // bool flag = check_in_approximate_field_list(maybe_used_face_field ,CGAL::centroid(tri));;
//                            //  bool flag2 = check_in_approximate_field_list(field_through_list ,CGAL::centroid(tri));
//                            bool flag = false;
//                            //Tree aabb_tree_final_round(field_triangles_final_round.begin(),field_triangles_final_round.end());
//
//                            //使用aabbtree 代替
//                            K2::Point_3 this_center = CGAL::centroid(tri);
//                            K2::Ray_3 ray(this_center,tri.supporting_plane().orthogonal_vector());
//                            std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections;
//                            aabb_tree.all_intersections(ray,std::back_inserter(intersections));
//                            //vector<bool>cutting_field_id(field_through_list.size(),false);
//                            vector<set<K2::Point_3 > > intersection_v(field_through_list.size());
//                            set<int>is_special;
//                            set<int> positive_side;
//                            for(auto item : intersections) {
//                                if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
//                                    //se.insert(*p);
//                                    //cout <<"*********"<<endl;
//                                    int this_field_belong_id = face_belong_field_all_mp[item.second->id()];
//                                    if (belong_field_id != this_field_belong_id){
//                                        if( *p == this_center){
//                                            //flag = true;
//                                            positive_side.insert(this_field_belong_id);
//                                            //cout <<"f1" << endl;
//
//                                        }
//                                        if(intersection_v[this_field_belong_id].count(*p)){
//                                            if(!is_special.count(this_field_belong_id)){
//                                                if(faces_approximate_field[field_through_list[this_field_belong_id]].in_or_on_field(this_center)) {
//                                                    flag = true;
//                                                    //cout <<"f2" << endl;
//                                                    break;
//                                                }
//                                            }
//                                            else
//                                                is_special.insert(this_field_belong_id);
//                                        }
//                                        intersection_v[this_field_belong_id].insert(*p);
//                                    }
//                                }
//                            }
//
//
//
//                            for(int j=0;j<field_through_list.size();j++){
//                                if(j==belong_field_id)continue;
//                                if(!is_special.count(j) && positive_side.count(j)==0){
//                                    if(intersection_v[j].size()%2) {
//                                        flag = true;
//                                        //cout <<"f3" << endl;
//                                        break;
//                                    }
//                                }
//                            }
//                            bool flag_positive = false;
//                            bool flag_negative = false;
//                            for(int j=0;j<field_through_list.size();j++){
//                                if(positive_side.count(j) && intersection_v[j].size()%2==0 && tri.supporting_plane().has_on_positive_side(faces_approximate_field[field_through_list[j]].center)){
//                                    flag_positive = true;
//                                }
//                            }
//                            if(!flag && flag_positive){
//                                K2::Ray_3 ray_r(this_center,tri.supporting_plane().orthogonal_vector()*(-1));
//                                std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections_r;
//                                aabb_tree.all_intersections(ray,std::back_inserter(intersections_r));
//                                //vector<bool>cutting_field_id(field_through_list.size(),false);
//                                vector<set<K2::Point_3 > > intersection_v_r(field_through_list.size());
//                                set<int>is_special_r;
//                                set<int> negative_side;
//                                for(auto item : intersections_r) {
//                                    if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
//                                        int this_field_belong_id = face_belong_field_all_mp[item.second->id()];
//                                        if (belong_field_id != this_field_belong_id){
//                                            if( *p == this_center){
//                                                negative_side.insert(this_field_belong_id);
//                                            }
//                                            if(intersection_v_r[this_field_belong_id].count(*p)){
//                                                if(!is_special_r.count(this_field_belong_id)){
//                                                    if(faces_approximate_field[field_through_list[this_field_belong_id]].in_or_on_field(this_center)) {
//                                                        flag = true;
//                                                        break;
//                                                    }
//                                                }
//                                                else
//                                                    is_special_r.insert(this_field_belong_id);
//                                            }
//                                            intersection_v_r[this_field_belong_id].insert(*p);
//                                        }
//                                    }
//                                }
//                                for(int j=0;j<field_through_list.size();j++){
//                                    if(negative_side.count(j) && intersection_v_r[j].size()%2==0  && tri.supporting_plane().has_on_negative_side(faces_approximate_field[field_through_list[j]].center)){
//                                        flag_negative = true;
//                                    }
//                                }
//                            }
//                            if(flag_positive && flag_negative)
//                                flag = true;
//
//
//                            // 这个地方是不是可以去掉底相交呢？？？？？
//
//                            if(!flag) {
//                                int side_type = //get_side_type(Point_K2_to_iGameVertex(CGAL::centroid(tri)));
//                                        cgal_polygon->outMesh((CGAL::centroid(tri)));
//                                if(side_type != 1) {
//                                    flag = true;
//                                    //cout <<"GGST" << endl;
//                                }
//                            }
//                            auto vvv0 = Point_K2_to_iGameVertex(tri.vertex(0));
//                            auto vvv1 = Point_K2_to_iGameVertex(tri.vertex(1));
//                            auto vvv2 = Point_K2_to_iGameVertex(tri.vertex(2));
//
//                            auto ddd = (vvv1 - vvv0) % (vvv2- vvv0);
//                            //if(ddd*direct <0)
//                            // swap(vvv1,vvv2);
//
//
//
////                       std::unique_lock<std::mutex>lock(mu);
////                       fprintf(file10,"v %lf %lf %lf\n",vvv0.x(),vvv0.y(),vvv0.z());
////                       fprintf(file10,"v %lf %lf %lf\n",vvv2.x(),vvv2.y(),vvv2.z());
////                       fprintf(file10,"v %lf %lf %lf\n",vvv1.x(),vvv1.y(),vvv1.z());
////                       fprintf(file10,"f %d %d %d\n",vid,vid+1,vid+2);
////                       vid+=3;
//
//                            if(flag) {
//                                continue;
//                            }
//
//                            //if(maybe_used_face_source_id[i][0] >= 3 && maybe_used_face_source_id[i][1] >= 3 && maybe_used_face_source_id[i][2] >= 3)
//                            generated_face_list.push_back(K2::Triangle_3(tri.vertex(0),tri.vertex(2),tri.vertex(1)));
////                        else{
////                            static int vid = 1;
////                            std::unique_lock<std::mutex>lock(mu);
////                            fprintf(file10,"v %lf %lf %lf\n",vvv0.x(),vvv0.y(),vvv0.z());
////                            fprintf(file10,"v %lf %lf %lf\n",vvv2.x(),vvv2.y(),vvv2.z());
////                            fprintf(file10,"v %lf %lf %lf\n",vvv1.x(),vvv1.y(),vvv1.z());
////                            fprintf(file10,"f %d %d %d\n",vid,vid+1,vid+2);
////                            vid+=3;
////                        }
///*****************************/
////                        auto vvv0 = Point_K2_to_iGameVertex(tri.vertex(0));
////                        auto vvv1 = Point_K2_to_iGameVertex(tri.vertex(1));
////                        auto vvv2 = Point_K2_to_iGameVertex(tri.vertex(2));
////
////                        auto ddd = (vvv1 - vvv0) % (vvv2- vvv0);
////                        //if(ddd*direct <0)
////                           // swap(vvv1,vvv2);
////
////
////
////                        std::unique_lock<std::mutex>lock(mu);
////                        fprintf(file10,"v %lf %lf %lf\n",vvv0.x(),vvv0.y(),vvv0.z());
////                        fprintf(file10,"v %lf %lf %lf\n",vvv2.x(),vvv2.y(),vvv2.z());
////                        fprintf(file10,"v %lf %lf %lf\n",vvv1.x(),vvv1.y(),vvv1.z());
////                        fprintf(file10,"f %d %d %d\n",vid,vid+1,vid+2);
////                        vid+=3;
//
//                        }
//                    }
//
//
////                    std::chrono::time_point<std::chrono::system_clock> end_time = std::chrono::system_clock::now();
////
////
////                    std::chrono::milliseconds diff = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
////
////
////                    if(diff.count() >10000)
////                        cout <<"id: "<<id<<" field size: "<< field_through_list.size() <<"maybe use " << maybe_used_face.size()<< " gen size: "<< generated_face_list.size()<<" use time "<<" "<< diff.count()<<" skip?"<< skip<<endl;
//
////                if(generated_face_list.size() ==0 && !skip && field_through_list.size() >30){
////                    auto small  = getGridVertex(each_grid->first,0);
////                    auto big  = getGridVertex(each_grid->first,7);
////                    static int f3_id = 1;
////                    for (int ii = 0; ii < 7; ii++) {
////                        for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
////                            int from = ii;
////                            int to = DirectedGridEdge[ii][jj];
////                            MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
////                            MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
////                            fprintf(file9, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
////                            fprintf(file9, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
////                            fprintf(file9, "l %d %d\n", f3_id, f3_id + 1);
////                            f3_id += 2;
////                        }
////                    }
////                    int f5id = 1;
////                    for(auto this_face : field_triangles_final_round){
////                        fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(0).x()),
////                                CGAL::to_double(this_face.vertex(0).y()),
////                                CGAL::to_double(this_face.vertex(0).z()));
////                        fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(1).x()),
////                                CGAL::to_double(this_face.vertex(1).y()),
////                                CGAL::to_double(this_face.vertex(1).z()));
////                        fprintf(file5, "v %lf %lf %lf \n",CGAL::to_double(this_face.vertex(2).x()),
////                                CGAL::to_double(this_face.vertex(2).y()),
////                                CGAL::to_double(this_face.vertex(2).z()));
////                        fprintf(file5,"f %d %d %d\n",f5id,f5id+1,f5id+2);
////                        f5id+=3;
////                    }
////                    if(field_through_list.size() > 30){
////                        skip = true;
////                        for(auto  each_container_face : container_grid_face){
////                            K2::Triangle_3 tri1(ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]]);
////                            K2::Triangle_3 tri2(ps[each_container_face[2]],ps[each_container_face[3]],ps[each_container_face[0]]);
////                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_1;
////                            std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections_2;
////                            aabb_tree.all_intersections(tri1,std::back_inserter(intersections_1));
////                            aabb_tree.all_intersections(tri2,std::back_inserter(intersections_2));
////                            vector<K2::Segment_3>vs;
////
////                            for(auto i : intersections_1){
////                                if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
////                                    vs.push_back(*s);
////                                }
////                            }
////                            for(auto i : intersections_2){
////                                if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(i.first))){
////                                    vs.push_back(*s);
////                                }
////                            }
////                            auto tmp = CGAL_CDT({ps[each_container_face[0]],ps[each_container_face[1]],ps[each_container_face[2]],ps[each_container_face[3]]},vs,tri1);
////                            for(auto i:tmp){
////                                //cout <<"skip true" << endl;
////                                if(!check_in_field(K2::Triangle_3(i[0],i[1],i[2]))){
////                                    // cout <<"skip false" << endl;
////                                    skipcc = CGAL::centroid(K2::Triangle_3(i[0],i[1],i[2]));
////                                    skip = false;
////                                    break;
////                                }
////                            }
////
//////                    for(K2::Triangle_3 this_face : field_triangles){
//////                        CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Triangle_3)>::type
//////                                res_tt = intersection(tri1,this_face);
//////                        if (res_tt) {
//////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
//////                                vs.push_back(*p);
//////                            }
//////                        }
//////                        res_tt = intersection(tri2,this_face);
//////                        if (res_tt) {
//////                            if (const K2::Segment_3 *p = boost::get<K2::Segment_3>(&*res_tt)) {
//////                                vs.push_back(*p);
//////                            }
//////                        }
//////                    }
////
////                            if(!skip)
////                                break;
////                        }
////                    }
////                    cout << CGAL::to_double(skipcc.x()) <<" "<< CGAL::to_double(skipcc.y())<<" "<< CGAL::to_double(skipcc.z()) << endl;
////                    exit(0);
////                }
//
////                vector<MeshKernel::iGameVertex >gen_vertex;
////                vector<vector<int> >gen_face;
//                    for(auto i : generated_face_list) {
//                        each_grid->second.generate_face_list.push_back(i);
//                    }
//                };
//                dfs(small,big,face_list,0);
//            }

        },i);
    }
    time_path.push_back(std::chrono::system_clock::now());
    for(int i=0;i<thread_num;i++)
        each_frame_thread[i]->join();
    time_path.push_back(std::chrono::system_clock::now());


    FILE *fileans = fopen( (input_filename + "_times.off").c_str(), "w+");
    for(int i =1;i< time_path.size();i++){
        fprintf(fileans,"%d\t%lf\n",i, (time_path[i].time_since_epoch().count()-time_path[i-1].time_since_epoch().count())*1.0/1000000);
        cout << "time point "<<i<<" "<< (time_path[i].time_since_epoch().count()-time_path[i-1].time_since_epoch().count())*1.0/1000000 << endl;
    }
    fclose(fileans);
    //cout <<"sttime " << (endtime - sttime)*1.0/thread_num/1000000 << endl;


    for(auto i : debug_time_use){
        if(i.second.size()==0)continue;
        double sum = 0;
        for(auto j : i.second){
            sum+= j;
        }
        sum/=i.second.size();
        cout << i.first <<"\t" << sum << endl;
    }

    // cout << "qq1 qq2 "<<qq1 <<" "<< qq2 << endl;
    vector<K2::Triangle_3>generate_face_final;

    //FILE *file4 = fopen( (input_filename + "_4.obj").c_str(), "w");
    vector<FILE *> v_dep;
    vector<int>vfid;
    for(int i=0;i<5;i++){
        v_dep.emplace_back(fopen( (input_filename + "_depth"+to_string(i)+".obj").c_str(), "w"));
        vfid.push_back(1);
    }

    unordered_map<size_t ,int> vmp;
    double miny=(1<<30),maxy=-(1<<30);
    double minz=(1<<30),maxz=-(1<<30);
    for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++){
        for(int i = 0;i < each_grid->second.generate_face_list.size();i++){
            vmp[each_grid->second.generate_face_list[i].id()] = generate_face_final.size();
            generate_face_final.push_back(each_grid->second.generate_face_list[i]);
            final_gen_face.push_back({final_gen_vertex.size(),final_gen_vertex.size()+1,final_gen_vertex.size()+2});
            final_gen_vertex.push_back( each_grid->second.generate_face_list[i].vertex(0));
            final_gen_vertex.push_back( each_grid->second.generate_face_list[i].vertex(1));
            final_gen_vertex.push_back( each_grid->second.generate_face_list[i].vertex(2));
            for(int k=0;k<3;k++){
                miny=min(miny, CGAL::to_double(each_grid->second.generate_face_list[i].vertex(k).y()));
                maxy=max(maxy, CGAL::to_double(each_grid->second.generate_face_list[i].vertex(k).y()));
                minz=min(minz, CGAL::to_double(each_grid->second.generate_face_list[i].vertex(k).z()));
                maxz=max(maxz, CGAL::to_double(each_grid->second.generate_face_list[i].vertex(k).z()));
            }
        }
        for(int depth = 0;depth<4;depth++){
            for(int i = 0;i < each_grid->second.generate_face_list_depth[depth].size();i++) {
                if(*set<double>{ CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(0).x()),
                                 CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(1).x()),
                                 CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(2).x()),
                }.begin()< (mesh->BBoxMin.x() + mesh->BBoxMax.x() )/2)
                    continue;


                fprintf(v_dep[depth], "v %lf %lf %lf\n", CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(0).x()), CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(0).y()),
                        CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(0).z()));
                fprintf(v_dep[depth], "v %lf %lf %lf\n", CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(1).x()), CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(1).y()),
                        CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(1).z()));
                fprintf(v_dep[depth], "v %lf %lf %lf\n", CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(2).x()), CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(2).y()),
                        CGAL::to_double( each_grid->second.generate_face_list_depth[depth][i].vertex(2).z()));
                fprintf(v_dep[depth], "f %d %d %d\n", vfid[depth], vfid[depth] + 1, vfid[depth] + 2);
                vfid[depth]+=3;
            }
        }
    }
    cout <<"st build mesh" << endl;
    // MeshBuilder(generate_face_final).build();




    cout <<"f v : " << final_gen_face.size() <<" "<< final_gen_vertex.size() << endl;
    CGAL::Polyhedron_3<K2>  pmesh;
    PMP::repair_polygon_soup(final_gen_vertex, final_gen_face);

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(final_gen_vertex, final_gen_face, pmesh);
    PMP::duplicate_non_manifold_vertices(pmesh);
    PMP::stitch_borders(pmesh);
    PMP::merge_duplicated_vertices_in_boundary_cycles(pmesh);


    stringstream ss;
    ss<<pmesh;
    FILE *file13 = fopen( (input_filename + "_13.off").c_str(), "w+");
    fprintf(file13,"%s",ss.str().c_str());
    FILE *file6 = fopen( (input_filename + "_6.obj").c_str(), "w+");
    FILE *file7 = fopen( (input_filename + "_7.obj").c_str(), "w+");


    if(1){
        CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> * inside = new CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2>(pmesh);

        int xxid = 1;




        auto center = (mesh->BBoxMin + mesh->BBoxMax)/2;
//         K2::Point_3 vmin(center.x(),mesh->BBoxMin.y(),mesh->BBoxMin.z());
//         K2::Point_3 vmax(center.x(),mesh->BBoxMax.y(),mesh->BBoxMax.z());
//         K2::Point_3 vmid1(center.x(),mesh->BBoxMin.y(),mesh->BBoxMax.z());
//         K2::Point_3 vmid2(center.x(),mesh->BBoxMax.y(),mesh->BBoxMin.z());
        K2::Point_3 vmin(center.x(),miny,minz);
        K2::Point_3 vmax(center.x(),maxy,maxz);
        K2::Point_3 vmid1(center.x(),miny,maxz);
        K2::Point_3 vmid2(center.x(),maxy,minz);


        K2::Plane_3 plane(vmin,vmax,vmid1);
        vector<K2::Segment_3>vs;
        for(int i=0;i<mesh->FaceSize();i++) {

            MeshKernel::iGameVertex v0 = field_move_vertex[mesh->fast_iGameFace[i].vh(0)];
            MeshKernel::iGameVertex v1 = field_move_vertex[mesh->fast_iGameFace[i].vh(1)];
            MeshKernel::iGameVertex v2 = field_move_vertex[mesh->fast_iGameFace[i].vh(2)];

            MeshKernel::iGameVertex ov0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)];
            MeshKernel::iGameVertex ov1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)];
            MeshKernel::iGameVertex ov2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)];
            K2::Triangle_3 tri(iGameVertex_to_Point_K2(ov0),
                               iGameVertex_to_Point_K2(ov1),
                               iGameVertex_to_Point_K2(ov2)
            );
            K2::Vector_3 normal = tri.supporting_plane().orthogonal_vector();
            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Plane_3)>::type
                    res_sp = intersection(tri,plane);
            bool flag = false;
            vector<K2::Segment_3>now;
            if (res_sp) {
                if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_sp)) {
                    vs.push_back(*s);
                    now.push_back(*s);
                    flag = true;
                }
            }
            if(!flag){
                if(CGAL::centroid(tri).x()<= center.x()){
                    fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(tri.vertex(0).x()), CGAL::to_double(tri.vertex(0).y()),
                            CGAL::to_double(tri.vertex(0).z()));
                    fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(tri.vertex(1).x()), CGAL::to_double(tri.vertex(1).y()),
                            CGAL::to_double(tri.vertex(1).z()));
                    fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(tri.vertex(2).x()), CGAL::to_double(tri.vertex(2).y()),
                            CGAL::to_double(tri.vertex(2).z()));
                    fprintf(file6, "f %d %d %d\n", xxid, xxid + 1, xxid + 2);
                    xxid += 3;
                }
            }
            else{
                vector<vector<K2::Point_3> > res = CGAL_CDT({tri.vertex(0),tri.vertex(1),tri.vertex(2)},
                                                            now,
                                                            tri
                );
                for(int i=0;i<res.size();i++){
                    K2::Triangle_3 tri(res[i][0],res[i][1],res[i][2]);
                    if(CGAL::centroid(tri).x()<= center.x() ) {
                        K2::Vector_3 normal2 = CGAL::cross_product((res[i][1] - res[i][0]) ,  (res[i][2] - res[i][0]));
                        if(normal * normal2 <CGAL::Epeck::FT(0))
                            swap(res[i][1],res[i][2]);
                        fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(res[i][0].x()), CGAL::to_double(res[i][0].y()),
                                CGAL::to_double(res[i][0].z()));
                        fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(res[i][1].x()), CGAL::to_double(res[i][1].y()),
                                CGAL::to_double(res[i][1].z()));
                        fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(res[i][2].x()), CGAL::to_double(res[i][2].y()),
                                CGAL::to_double(res[i][2].z()));
                        fprintf(file6, "f %d %d %d\n", xxid, xxid + 1, xxid + 2);
                        xxid += 3;
                    }
                }
            }
        }

        for(auto tri :generate_face_final){
            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3 , K2::Plane_3)>::type
                    res_sp = intersection(tri,plane);
            K2::Vector_3 normal = tri.supporting_plane().orthogonal_vector();
            bool flag = false;
            vector<K2::Segment_3>now;
            if (res_sp) {
                if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_sp)) {
                    vs.push_back(*s);
                    now.push_back(*s);
                    flag = true;
                }
            }
            if(!flag){
                if(CGAL::centroid(tri).x()<= center.x()){
                    fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(tri.vertex(0).x()), CGAL::to_double(tri.vertex(0).y()),
                            CGAL::to_double(tri.vertex(0).z()));
                    fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(tri.vertex(1).x()), CGAL::to_double(tri.vertex(1).y()),
                            CGAL::to_double(tri.vertex(1).z()));
                    fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(tri.vertex(2).x()), CGAL::to_double(tri.vertex(2).y()),
                            CGAL::to_double(tri.vertex(2).z()));
                    fprintf(file6, "f %d %d %d\n", xxid, xxid + 1, xxid + 2);
                    xxid += 3;
                }
            }
            else{
                vector<vector<K2::Point_3> > res = CGAL_CDT({tri.vertex(0),tri.vertex(1),tri.vertex(2)},
                                                            now,
                                                            tri
                );
                for(int i=0;i<res.size();i++){
                    K2::Triangle_3 tri(res[i][0],res[i][1],res[i][2]);
                    if(CGAL::centroid(tri).x()<= center.x()) {
                        K2::Vector_3 normal2 = CGAL::cross_product((res[i][1] - res[i][0]) ,  (res[i][2] - res[i][0]));
                        if(normal * normal2 <CGAL::Epeck::FT(0))
                            swap(res[i][1],res[i][2]);
                        fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(res[i][0].x()), CGAL::to_double(res[i][0].y()),
                                CGAL::to_double(res[i][0].z()));
                        fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(res[i][1].x()), CGAL::to_double(res[i][1].y()),
                                CGAL::to_double(res[i][1].z()));
                        fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(res[i][2].x()), CGAL::to_double(res[i][2].y()),
                                CGAL::to_double(res[i][2].z()));
                        fprintf(file6, "f %d %d %d\n", xxid, xxid + 1, xxid + 2);
                        xxid += 3;
                    }
                }
            }
        }
        //vector<vector<K2::Point_3> > CGAL_CDT(vector<K2::Point_3> sorted_bound_vertex, vector<K2::Segment_3> cs,K2::Triangle_3 origin_face)
        vector<vector<K2::Point_3> > res = CGAL_CDT({vmin,vmid1,vmax,vmid2},
                                                    vs,
                                                    K2::Triangle_3(vmin,vmid1,vmax)
        );
        K2::Vector_3 normal = CGAL::cross_product((vmin-vmid1), (vmin-vmax));

        for(int i=0;i<res.size();i++){
            K2::Vector_3 normal2 = CGAL::cross_product((res[i][1] - res[i][0]) ,  (res[i][2] - res[i][0]));
            if(normal * normal2 >CGAL::Epeck::FT(0))
                swap(res[i][1],res[i][2]);
        }

        xxid = 1;
        for(int i=0;i<res.size();i++){
            K2::Triangle_3 tri(res[i][0],res[i][1],res[i][2]);
#ifdef MODEIN
            if( cgal_polygon->inMesh(centroid(tri)) && (*inside)(centroid(tri))!= CGAL::ON_BOUNDED_SIDE) {
#else
            if( cgal_polygon->outMesh(centroid(tri)) && (*inside)(centroid(tri))== CGAL::ON_BOUNDED_SIDE) {
#endif
                fprintf(file7, "v %lf %lf %lf\n", CGAL::to_double(res[i][0].x()), CGAL::to_double(res[i][0].y()),
                        CGAL::to_double(res[i][0].z()));
                fprintf(file7, "v %lf %lf %lf\n", CGAL::to_double(res[i][1].x()), CGAL::to_double(res[i][1].y()),
                        CGAL::to_double(res[i][1].z()));
                fprintf(file7, "v %lf %lf %lf\n", CGAL::to_double(res[i][2].x()), CGAL::to_double(res[i][2].y()),
                        CGAL::to_double(res[i][2].z()));
                fprintf(file7, "f %d %d %d\n", xxid, xxid + 1, xxid + 2);
                xxid += 3;
            }
        }





//        mesh->initBBox();
//        auto center = (mesh->BBoxMin + mesh->BBoxMax)/2;
//        list<K2::Triangle_3>l_inner;
//        for(auto i: generate_face_final)
//            l_inner.push_back(i);
//        Tree  tree_inner(l_inner.begin(),l_inner.end());
//        list<K2::Triangle_3>l_outer;
//        list<K2::Triangle_3>l_outer;

    }

    return 0;
}
// 1 2 3 4 5 6
// 5 1 2 3 7 8
//



