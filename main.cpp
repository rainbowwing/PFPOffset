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
#include "lib_impl.h"
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
//#include "BVH.h"


using namespace std;


shared_ptr <MeshKernel::SurfaceMesh> mesh;

shared_ptr<CGALPolygon>cgal_polygon;

Tree cgal_global_aabbtree;
double default_move = 0.1;
int thread_num = 16;
int file_id;
const double tolerance = 0.2;


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

vector <vector<int>> GridEdge = {{1, 2, 3}, //0
                                 {0, 4, 5}, //1
                                 {0, 4, 6}, //2
                                 {0, 5, 6}, //3
                                 {1, 2, 7}, //4
                                 {1, 3, 7}, //5
                                 {2, 3, 7}, //6
                                 {4, 5, 6}}; //7

vector <vector<int> > container_grid_face = {{0,3,5,1},
                                             {0,2,6,3},
                                             {0,1,4,2},
                                             {2,4,7,6},
                                             {3,6,7,5},
                                             {1,5,7,4}
                                             };
/*
f 1 4 6 2
f 1 3 7 4
f 1 2 5 3
f 3 5 8 7
f 4 7 8 6
f 2 6 8 5
 */
vector <grid> get_edge_neighbor(grid g1, grid g2) {

    grid now = min(g1,g2);
    int diff = -1;
    int nowx = now.x;
    int nowy = now.y;
    int nowz = now.z;
    if(g1.x != g2.x) {
        diff = 0;
    }
    else if(g1.y != g2.y) {
        diff = 1;
    }
    else {
        diff = 2;
    }

    if (diff == 0) {
        return {{nowx, nowy,     nowz},
                {nowx, nowy - 1, nowz},
                {nowx, nowy - 1, nowz - 1},
                {nowx, nowy,     nowz - 1}};
    } else if (diff == 1) {
        return {{nowx,     nowy, nowz},
                {nowx,     nowy, nowz - 1},
                {nowx - 1, nowy, nowz - 1},
                {nowx - 1, nowy, nowz}};
    } else {
        assert(diff == 2);
        return {{nowx,     nowy,     nowz},
                {nowx - 1, nowy,     nowz},
                {nowx - 1, nowy - 1, nowz},
                {nowx,     nowy - 1, nowz}};
    }
}




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

MeshKernel::iGameVertex getTinyGridVertex(const MeshKernel::iGameVertex& small, const MeshKernel::iGameVertex& big,int k) {

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





grid getGridFrameVertex(grid g, int k) {
    grid v = g;
    v.x += GridVertexDir[k][0];
    v.y += GridVertexDir[k][1];
    v.z += GridVertexDir[k][2];
    return v;
}

grid vertex_to_grid(MeshKernel::iGameVertex v) {
    int x = int((v.x() - stx) / grid_len + myeps);
    int y = int((v.y() - sty) / grid_len + myeps);
    int z = int((v.z() - stz) / grid_len + myeps);
    return grid{x, y, z};
}




bool vertex_in_grid(grid g, MeshKernel::iGameVertex v) {
    auto v1 = getGridVertex(g, 0);
    auto v2 = getGridVertex(g, 7);
    return v1.x() <= v.x() && v.x() <= v2.x() &&
           v1.y() <= v.y() && v.y() <= v2.y() &&
           v1.z() <= v.z() && v.z() <= v2.z();
}

struct GridVertex {
    int grid_type;
    vector <MeshKernel::iGameFaceHandle> face_list;

    GridVertex() {
        grid_type = -1;
    }
};


int vertex_state(std::vector<MeshKernel::iGameFaceHandle>face_list, MeshKernel::iGameVertex v,
                 std::shared_ptr<std::vector<MeshKernel::iGameFaceHandle> > field_face_list = nullptr){


    int  side_type=-2;
    vector<MeshKernel::iGameFace>vf;
    for(auto i : face_list){
        vf.push_back(mesh->fast_iGameFace.at(i));
    }
    cgal_aabbtree_query(vf, v, mesh, side_type, cgal_polygon);
    if(side_type == 0)
        return 0;
    bool in_field = false;
    //  cout <<"pos:"<< v.x()<<" "<< v.y()<<" "<< v.z()<<endl;
    for(MeshKernel::iGameFaceHandle fh : face_list) {
        double dist = cgal_vertex_triangle_dist(mesh->fast_iGameFace.at(fh), v, mesh);
        //   cout << "fh:"<< fh<<" dist:"<<dist <<endl;
        if(vertex_field_state(fh,v,mesh)){
            in_field = true;
            if(field_face_list){
                field_face_list->push_back(fh);
            }
        }

//        if(dist <= mesh->fast_iGameFace.at(fh).move_dist) {
//            in_field = true;
//            if(field_face_list){
//                field_face_list->push_back(fh);
//            }
//        }

    }

    if(in_field) {
        return 1;
    }
    return 2;
}


//struct ThickenCubeEdge{
//    grid from;
//    grid to;
//    int type;
//    MeshKernel::iGameVertex edge_mesh_intersect_vertex;
//    MeshKernel::iGameVertex field_far_vertex;
//};

//map<grid,vector<ThickenCubeEdge> > thicken_cube_edge_mp;
//
//double ternary_search(MeshKernel::iGameVertex v1,MeshKernel::iGameVertex v2,MeshKernel::iGameFaceHandle fh){
//    double l=0;
//    double r=1;
//    const double eps = 0.001;
//    Point_3 a(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)].x(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)].y(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)].z());
//    Point_3 b(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)].x(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)].y(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)].z());
//    Point_3 c(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)].x(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)].y(),
//              mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)].z());
//    Triangle tri (a,b,c);
//    while(r-l>eps)
//    {
//        double mid1=(2*l+r)/3;
//        double mid2=(l+2*r)/3;
//        MeshKernel::iGameVertex Vmid1 = v1*(1-mid1) + v2*mid1;
//        Point Pmid1(Vmid1.x(),Vmid1.y(),Vmid1.z());
//        MeshKernel::iGameVertex Vmid2 = v1*(1-mid2) + v2*mid2;
//        Point Pmid2(Vmid2.x(),Vmid2.y(),Vmid2.z());
//        auto dist1 =  CGAL::squared_distance(Pmid1,tri);
//        auto dist2 =  CGAL::squared_distance(Pmid2,tri);
//       if(dist1>dist2) l=mid1;
//       else r=mid2;
//    }
//    return l;
//}
//
//void solve_grid(grid g,std::vector<MeshKernel::iGameFaceHandle>face_handle_list){
//    sort(face_handle_list.begin(), face_handle_list.end());
//    face_handle_list.resize(unique(face_handle_list.begin(), face_handle_list.end()) - face_handle_list.begin());
//    vector<MeshKernel::iGameFace>face_list;
//    for(auto i : face_handle_list){
//        face_list.push_back(mesh->fast_iGameFace[i]);
//    }
//    std::vector<BVH::BVH_Face>bvh_faces;
//    for(MeshKernel::iGameFaceHandle i: face_handle_list){
//        BVH::BVH_Face tmp;
//        for(int id = 0 ; id < 3 ;id++) {
//            tmp.vertices.emplace_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(id)].x(),
//                                      mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(id)].y(),
//                                      mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(id)].z());
//        }
//        tmp.index = i;
//        bvh_faces.push_back(tmp);
//    }
//    BVH::BVH_Tree bvhTree;
//    bvhTree.buildBVH_Tree(bvh_faces);
//
//    for (int ev1 = 0; ev1 < GridEdge.size(); ev1++) {
//        for (auto ev2: GridEdge[ev1]) {
//            if (ev1 < ev2) {
//                for(MeshKernel::iGameFaceHandle i: face_handle_list) {
//                    grid gv1 = getGridFrameVertex(g, ev1);
//                    grid gv2 = getGridFrameVertex(g, ev2);
//                    MeshKernel::iGameVertex iv1 = getGridVertex(g, ev1);
//                    MeshKernel::iGameVertex iv2 = getGridVertex(g, ev2);
//                    double near_pos_div = ternary_search(iv1, iv2, i);
//                    MeshKernel::iGameVertex near_pos = iv1 * (1 - near_pos_div) + iv2 * near_pos_div;
//                    double dist = cgal_vertex_triangle_dist(mesh->fast_iGameFace[i],near_pos,mesh);
//                    if(dist > mesh->fast_iGameFace[i].move_dist)
//                        continue;
//
//
//
//                    MeshKernel::iGameVertex left_start;
//                    MeshKernel::iGameVertex right_start;
//                    bool need_left_start;
//                    bool need_right_start;
//                    int side = -1;
//                    cgal_aabbtree_query(face_list,near_pos,mesh,side,cgal_polygon);
//                    if(side ==0){
//
//                    }
//                    else{
//
//                    }
//                }
//            }
//        }
//    }
//
///*
//MeshKernel::iGameVertex cgal_aabbtree_query(std::vector<MeshKernel::iGameFace>& f_list,
//                                            MeshKernel::iGameVertex v,
//                                            std::shared_ptr<MeshKernel::SurfaceMesh>mesh,
//                                            int & side,std::shared_ptr<CGALPolygon>cgalinmesh);
// */
//
//
////    for(MeshKernel::iGameFaceHandle i: face_list){
////
////
////    }
//
//
//}
//
//vector<pair<MeshKernel::iGameVertex,MeshKernel::iGameVertex> >
//        solve_edge_field(MeshKernel::iGameVertex v1,MeshKernel::iGameVertex v2,std::vector<MeshKernel::iGameFaceHandle>face_list
//                ,BVH::BVH_Tree& local_bvh_tree){
//    vector<pair<MeshKernel::iGameVertex,MeshKernel::iGameVertex>  >segments_in_field;
//    for(MeshKernel::iGameFaceHandle fh : face_list){
//        // 先计算最近的点555
//
//
//
//
//        // 在计算平面与交点
//
//
//    }
//
//    return segments_in_field;
//};


vector<MeshKernel::iGameVertex> field_move_vertex;
vector<MeshKernel::iGameVertex> field_move_normal;

//vector<vector<int> >approximate_field_face_table = {{0,1,2},{3,4,5},{0,2,4},{4,3,0},{1,5,4},{4,2,1},{0,3,5},{5,1,0}};

inline K::Point_3 iGameVertex_to_Point(const MeshKernel::iGameVertex& v){
    return K::Point_3(v.x(),v.y(),v.z());
}

inline K2::Point_3 iGameVertex_to_Point_K2(const MeshKernel::iGameVertex& v){
    return K2::Point_3(v.x(),v.y(),v.z());
}

inline  MeshKernel::iGameVertex Point_K_to_iGameVertex(const K::Point_3& v){
    return MeshKernel::iGameVertex(CGAL::to_double(v.x()),CGAL::to_double(v.y()),CGAL::to_double(v.z()));
}

inline  MeshKernel::iGameVertex Point_K2_to_iGameVertex(const K2::Point_3& v){
    return MeshKernel::iGameVertex(CGAL::to_double(v.x()),CGAL::to_double(v.y()),CGAL::to_double(v.z()));
}

bool in_triangle_positive_side(MeshKernel::iGameFaceHandle fh,const MeshKernel::iGameVertex& v){
    MeshKernel::iGameVertex v0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)];
    MeshKernel::iGameVertex v1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)];
    MeshKernel::iGameVertex v2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)];
    MeshKernel::iGameVertex normal = ((v1 - v0) % (v2 - v0)).normalize();
    MeshKernel::iGameVertex center = (v0 + v1 + v2)/3;
    if(((v - center) * normal)>0){
        return true;
    }
    return false;
}

struct ApproximateField{
    vector<MeshKernel::iGameVertex> origin_vertices;
    vector<MeshKernel::iGameVertex> extend_vertices;
    vector<K2::Tetrahedron_3>tet_list;
    MeshKernel::iGameFaceHandle fh;
    ApproximateField(){}
    ApproximateField(MeshKernel::iGameFaceHandle fh){
        this->fh=fh;
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)]);
        //TODO 注意分解
        MeshKernel::iGameVertex normal = ((origin_vertices[1] - origin_vertices[0]) % (origin_vertices[2] - origin_vertices[0])).normalize();
        for(int i=0;i<3;i++){
            //extend_vertices.push_back(field_move_vertex[mesh->fast_iGameFace[fh].vh(i)]);
            if(!in_triangle_positive_side(this->fh,field_move_vertex[mesh->fast_iGameFace[fh].vh(i)])){
                extend_vertices.push_back(field_move_vertex[mesh->fast_iGameFace[fh].vh(i)]);
            }
            else{
                extend_vertices.push_back(origin_vertices[i] + normal * mesh->fast_iGameFace[fh].move_dist);
            }
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

    // 6个四面体法，防止相交处理麻烦 并且用tet 来判断内外这样每一个面的偏移就是6个tet，；
    }



    ApproximateField(MeshKernel::iGameFaceHandle fh,vector<MeshKernel::iGameVertex>moved_vertex){
        this->fh=fh;
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(0)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(1)]);
        origin_vertices.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[fh].vh(2)]);
        //TODO 注意分解
        MeshKernel::iGameVertex normal = ((origin_vertices[1] - origin_vertices[0]) % (origin_vertices[2] - origin_vertices[0])).normalize();
        for(int i=0;i<3;i++){
            extend_vertices.push_back(moved_vertex[i]);
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

        // 6个四面体法，防止相交处理麻烦 并且用tet 来判断内外这样每一个面的偏移就是6个tet，；
    }
};



vector<ApproximateField>faces_approximate_field;
vector<MeshKernel::iGameVertex> min_move_g;
vector<MeshKernel::iGameVertex> max_move_g;
MeshKernel::iGameVertex do_quadratic_error_metric(MeshKernel::iGameVertexHandle vh){
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
    MeshKernel::iGameVertex avg_move_vertex(0,0,0);
    for (auto f : mesh->FastNeighborFhOfVertex_[vh]){
        avg_move_dist += mesh->fast_iGameFace[f].move_dist;

    }
    avg_move_dist /= (1.0*mesh->FastNeighborFhOfVertex_[vh].size());

    for (auto f : mesh->FastNeighborFhOfVertex_[vh]) {
        MeshKernel::iGameVertex normal
        = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]
           - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]) %
          (mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]
           - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])).normalize();

        normal = normal * -1;
//        MeshKernel::iGameVertex new_v = v + normal * mesh->fast_iGameFace[f].move_dist;
//
//        MeshKernel::iGameVertex move_max_v = v + normal * mesh->fast_iGameFace[f].move_dist*1.35;
//        MeshKernel::iGameVertex move_min_v = v + normal * mesh->fast_iGameFace[f].move_dist*0.8;
        MeshKernel::iGameVertex new_v = v + normal * avg_move_dist;
        avg_move_vertex += normal * avg_move_dist;

        MeshKernel::iGameVertex move_max_v = v + normal * avg_move_dist * 1.35;
        MeshKernel::iGameVertex move_min_v = v + normal * avg_move_dist * 0.8;

        double d = -(normal.x() * new_v.x() + normal.y() * new_v.y() +  normal.z() * new_v.z());

        Eigen::Vector3d p(normal.x(),  normal.y(), normal.z());
        Eigen::Matrix3d A = p * p.transpose();
        Eigen::Vector3d D2(2*d*normal.x(),2*d*normal.y(),2*d*normal.z());
        Eigen::Vector3d LimitV(normal.x(),normal.y(),normal.z());

        double lower = normal.x()*move_min_v.x() + normal.y()*move_min_v.y() + normal.z()*move_min_v.z();
        double upper = normal.x()*move_max_v.x() + normal.y()*move_max_v.y() + normal.z()*move_max_v.z();



        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                hessian.coeffRef(i,j)+=A.coeff(i,j)*2;
        for(int i=0;i<3;i++)
            gradient.coeffRef(i) +=  D2.coeffRef(i);


        for(int i=0;i<3;i++)
            linearMatrix.coeffRef(cnt,i) = LimitV[i];
        lowerBound.coeffRef(cnt) = lower;
        upperBound.coeffRef(cnt) = upper;

        cnt++;
    }
    avg_move_vertex/= mesh->FastNeighborFhOfVertex_[vh].size();

    min_move_g[vh] =  v + avg_move_vertex;

    avg_move_vertex = v + avg_move_vertex;;
    hessian.coeffRef(0,0) += (2.0)/1000;
    hessian.coeffRef(1,1) += (2.0)/1000;
    hessian.coeffRef(2,2) += (2.0)/1000;

    gradient.coeffRef(0) -= (2.0/1000) * v.x();
    gradient.coeffRef(1) -= (2.0/1000) * v.y();
    gradient.coeffRef(2) -= (2.0/1000) * v.z();

//    MeshKernel::iGameVertex assist_vertex = avg_vertex/sum_move;

//    ErrorMat(0,0) += assist_vertex.x()/10;
//    ErrorMat(1,1) += assist_vertex.y()/10;
//    ErrorMat(2,2) += assist_vertex.z()/10;

//    Eigen::Matrix4d Q  =  ErrorMat;
//   // cout << Q << endl;
//
//    Eigen::Matrix4d A;
//    A << Q(0, 0), Q(0, 1), Q(0, 2), Q(0, 3),
//            Q(0, 1), Q(1, 1), Q(1, 2), Q(1, 3),
//            Q(0, 2), Q(1, 2), Q(2, 2), Q(2, 3),
//            0.f, 0.f, 0.f, 1.f;
//
//
//
//    if (abs(A.determinant())>myeps) {
//        Eigen::Vector4d b(0.f, 0.f, 0.f, 1.f);
//
//        //Eigen::Vector4d X = A.colPivHouseholderQr().solve(b);
//        Eigen::Vector4d X = A.inverse() * b;
//        //cout << "cost:" <<X.transpose() * Q * X << endl;
//        //cout << v.x() <<" "<< v.y()<<" "<< v.z() <<" "<<X[0]<<" "<< X[1]<<" "<< X[2] << endl;
//        return {X[0] , X[1] , X[2]};
//    }
//    else{
//        cout <<"err\n";
//        //return avg_vertex/sum_move;
//        return {1000,1000,1000};
//    }
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
    if(solver.solve()) {
        Eigen::VectorXd QPSolution = solver.getSolution();
        return {QPSolution.coeffRef(0), QPSolution.coeffRef(1), QPSolution.coeffRef(2)};
    }
    else{
        puts("use avg_move_vertex instead");
        return avg_move_vertex;
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




bool vertex_in_tiny_grid(const MeshKernel::iGameVertex& small,const MeshKernel::iGameVertex& big,
                         const MeshKernel::iGameVertex& v){
    return small.x() <= v.x()  && v.x() <= big.x() &&
    small.y() <= v.y()  && v.y() <= big.y() &&
    small.z() <= v.z()  && v.z() <= big.z();
}

bool face_through_grid(const MeshKernel::iGameVertex& small, const MeshKernel::iGameVertex& big,
                       const vector<MeshKernel::iGameVertex>& v){
    for(int i=0;i<v.size();i++){
        if(vertex_in_tiny_grid(small,big,v[i]))
            return true;
    }
    K::Triangle_3 tri(K::Point_3(v[0].x(),v[0].y(),v[0].z()),
                       K::Point_3(v[1].x(),v[1].y(),v[1].z()),
                       K::Point_3(v[2].x(),v[2].y(),v[2].z()));
    vector<K::Point_3>grid_vertex(8);
    for(int i=0;i<8;i++)
        grid_vertex[i] = iGameVertex_to_Point(getTinyGridVertex(small,big,i));

    for(auto  each_container_face : container_grid_face){
        K::Triangle_3 tri1(grid_vertex[each_container_face[0]],
                           grid_vertex[each_container_face[1]],
                           grid_vertex[each_container_face[2]]);
        K::Triangle_3 tri2(grid_vertex[each_container_face[1]],
                           grid_vertex[each_container_face[2]],
                           grid_vertex[each_container_face[3]]);
        if(sqrt(squared_distance(tri1,tri)) < myeps){
            return true;
        }
        if(sqrt(squared_distance(tri2,tri)) < myeps){
            return true;
        }
    }
    return false;

//    for(int i=0;i<6;i++){ // 用K1核就放宽成本，相差一点点的也算在里面，避免后期BFS误差，然后具体的用K2核心生成时修复;;;;;;
//        MeshKernel::iGameVertex v1 = getTinyGridVertex(small,big,i);
//        for(int j=0;j<DirectedGridEdge[i].size();j++){
//            MeshKernel::iGameVertex v2 = getTinyGridVertex(small,big,DirectedGridEdge[i][j]);
//            K::Segment_3 se(K::Point_3(v1.x(),v1.y(),v1.z()),
//                             K::Point_3(v2.x(),v2.y(),v2.z()));
//            if(sqrt(squared_distance(K::Point_3(v1.x(),v1.y(),v1.z()),tri)  ) < myeps ){
//                return true;
//            }
//            if(sqrt(squared_distance(K::Point_3(v2.x(),v2.y(),v2.z()),tri)  ) < myeps ){
//                return true;
//            }
//
//            CGAL::cpp11::result_of<K::Intersect_3(K::Segment_3, K::Triangle_3)>::type
//                    result = intersection(se, tri);
//            if(result)
//                return true;
//        }
//    }
    return false;
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


struct SharpPoint{
    K2::Point_3 p;
    vector<int> source_face_local_id;
};

int furthest_vertex_in_vector(int x , const vector<MeshKernel::iGameVertex>& v){
    int ans = x;
    for(int i=0;i<v.size();i++){
        if( (v[i] - v[x]).norm() > (v[ans] - v[x]).norm()){
            ans = i;
        }
    }

    return ans;
}

void delete_same_vertex(vector<MeshKernel::iGameVertex>& v){
    vector<MeshKernel::iGameVertex> tmp;
    for(auto i : v){
        bool flag = false;
        for(auto j : tmp){
            if( (i-j).norm() < myeps)
                flag = true;
        }
        if(!flag)
            tmp.push_back(i);
    }
    swap(tmp,v);
}


MeshKernel::iGameVertex sort_by_polar_order(vector<MeshKernel::iGameVertex>& v,MeshKernel::iGameVertex orthogonal_direction){
    MeshKernel::iGameVertex center(0,0,0);
    for(auto i : v){
        center = center + i;
    }
    center /= v.size();
    if(v.size() <=1)
        return center;

    function<int(double,double)> quadrant = [](double x,double y){
        if(x>0 && y > 0)return 1;
        else if(x<0 && y > 0)return 2;
        else if(x<0 && y < 0)return 3;
        return 4;
    };
    Plane_3 p;
    MeshKernel::iGameVertex x_axis = (v[0] - center).normalize();
    MeshKernel::iGameVertex y_axis = (orthogonal_direction % x_axis).normalize();
//    FILE *file11 = fopen(("../data/output" + to_string(file_id) + "_dianyun11.obj").c_str(), "w");
//    FILE *file12 = fopen(("../data/output" + to_string(file_id) + "_dianyun12.obj").c_str(), "w");
//    FILE *file13 = fopen(("../data/output" + to_string(file_id) + "_dianyun13.obj").c_str(), "w");
//    FILE *file14 = fopen(("../data/output" + to_string(file_id) + "_dianyun14.obj").c_str(), "w");
//    FILE *file15 = fopen(("../data/output" + to_string(file_id) + "_dianyun15.obj").c_str(), "w");
//    fprintf(file15,"v %lf %lf %lf\n",center.x(),center.y(),center.z());
//    fprintf(file15,"v %lf %lf %lf\n",(center + x_axis).x(),(center + x_axis).y(),(center + x_axis).z());
//    fprintf(file15,"v %lf %lf %lf\n",(center + y_axis).x(),(center + y_axis).y(),(center + y_axis).z());
//    fprintf(file15,"l 1 2\n");
//    fprintf(file15,"l 1 3\n");
//    for(int i=0;i<v.size();i++){
//        double x1 = (v[i] - center) * x_axis;
//        double y1 = (v[i] - center) * y_axis;
//        int q1 = quadrant(x1,y1);
//        if(q1 == 1){
//            fprintf(file11,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//        if(q1 == 2){
//            fprintf(file12,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//        if(q1 == 3){
//            fprintf(file13,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//        if(q1 == 4){
//            fprintf(file14,"v %lf %lf %lf\n",v[i].x(),v[i].y(),v[i].z());
//        }
//
//    }

    sort(v.begin(),v.end(),[&](MeshKernel::iGameVertex v1 , MeshKernel::iGameVertex v2){
        double x1 = (v1 - center) * x_axis;
        double y1 = (v1 - center) * y_axis;
        double x2 = (v2 - center) * x_axis;
        double y2 = (v2 - center) * y_axis;
        int q1 = quadrant(x1,y1);
        int q2 = quadrant(x2,y2);
        if(q1!=q2)return q1<q2;
        else
            return x1*y2 - x2*y1 > 0;

    });
//    for(auto i : v){
//        double x1 = (i - center) * x_axis;
//        double y1 = (i - center) * y_axis;
//        cout << quadrant(x1,y1) << "??\n";
//    }
    return center;
}

int main() {


    // freopen("../debugoutput.txt","w",stdout);
    default_move = 0.01;
    grid_len = 2.5;
    // 2.5 552531870234
    // 1.7 550170141470
    cout << grid_len <<endl;
    mix_factor = 0.5;
   // mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/debug4.obj2")); grid_len = 0.05;double default_move_dist = 0.1;

//    MeshKernel::iGameVertex v1(100,-1,-1);
//    MeshKernel::iGameVertex v2(1,-1,-1);
//    double t = ternary_search(v1,v2,MeshKernel::iGameFaceHandle(0)) ;
//    cout << t << endl;
//
//    return 0;


    // mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/Armadillo.obj")); grid_len = 2.5;


    mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/test_orgv2.obj2")); grid_len = 2.5; double default_move_dist = 0.8;

//    for(int i=0;i<mesh->FaceSize();i++){
//        mesh->faces(MeshKernel::iGameFaceHandle(i)).move_dist = 0.05;
//    }

    for(int i=0;i<mesh->FaceSize();i++){
        mesh->faces(MeshKernel::iGameFaceHandle(i)).move_dist = default_move_dist;
    }

    mesh->build_fast();

    //mix(30);


//    BVH::BVH_Tree bvhTree;
//    bvhTree.buildBVH_Tree(*mesh);
//    auto res = bvhTree.getIntersection(BVH::Ray(Vector3d(1,1,1),Vector3d(-1,-1,-1)));
//    cout << res.pos[0]<<" "<< res.pos[1]<<" "<< res.pos[2] << endl;


    //只动xy
    faces_approximate_field.resize(mesh->FaceSize());
    field_move_vertex.resize(mesh->VertexSize());
    min_move_g.resize(mesh->VertexSize());
    max_move_g.resize(mesh->VertexSize());

    field_move_normal.resize(mesh->FaceSize());


    for(int i=0;i<mesh->VertexSize();i++){
        field_move_vertex[i] = do_quadratic_error_metric(MeshKernel::iGameVertexHandle(i));
    }

    for(int i=0;i<mesh->FaceSize();i++){
        faces_approximate_field[i] = ApproximateField(MeshKernel::iGameFaceHandle(i));
        MeshKernel::iGam
        eVertex v0 = field_move_vertex[mesh->fast_iGameFace[i].vh(0)];
        MeshKernel::iGameVertex v1 = field_move_vertex[mesh->fast_iGameFace[i].vh(1)];
        MeshKernel::iGameVertex v2 = field_move_vertex[mesh->fast_iGameFace[i].vh(2)];

        MeshKernel::iGameVertex ov0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)];
        MeshKernel::iGameVertex ov1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)];
        MeshKernel::iGameVertex ov2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)];


        MeshKernel::iGameVertex normal = (v1 - v0) % (v2 - v0);
        MeshKernel::iGameVertex normal_o = (ov1 - ov0) % (ov2 - ov0);
        field_move_normal[i] = normal;
        if(normal * normal_o <0){
            field_move_normal[i] *= -1;
        }
    }


//

    file_id = 3018;
    FILE *file = fopen(("../data/output" + to_string(file_id) + ".obj").c_str(), "w");
    FILE *file0 = fopen(("../data/output" + to_string(file_id) + "_dianyun0.obj").c_str(), "w");
    FILE *file1 = fopen(("../data/output" + to_string(file_id) + "_dianyun1.obj").c_str(), "w");
    FILE *file2 = fopen(("../data/output" + to_string(file_id) + "_dianyun2.obj").c_str(), "w");
    FILE *file3 = fopen(("../data/output" + to_string(file_id) + "_dianyun3.obj").c_str(), "w");
    FILE *file4 = fopen(("../data/output" + to_string(file_id) + "_dianyun4.obj").c_str(), "w");
    FILE *file5 = fopen(("../data/output" + to_string(file_id) + "_dianyun5.obj").c_str(), "w");
    FILE *file6 = fopen(("../data/output" + to_string(file_id) + "_dianyun6.obj").c_str(), "w");
    FILE *file7 = fopen(("../data/output" + to_string(file_id) + "_dianyun7.obj").c_str(), "w");
    FILE *file8 = fopen(("../data/output" + to_string(file_id) + "_dianyun8.obj").c_str(), "w");
    FILE *file10 = fopen(("../data/output" + to_string(file_id) + "_dianyun10.obj").c_str(), "w");

    FILE *fileis1 = fopen(("../data/output" + to_string(file_id) + "_dianyunis1.obj").c_str(), "w");
    FILE *fileis2 = fopen(("../data/output" + to_string(file_id) + "_dianyunis2.obj").c_str(), "w");
    FILE *fileis3 = fopen(("../data/output" + to_string(file_id) + "_dianyunis3.obj").c_str(), "w");

    int cnt=1;

    for(int i=0;i<mesh->FaceSize();i++){
//        auto v0 = faces_approximate_field[i].extend_vertices[0];
//        auto v1 = faces_approximate_field[i].extend_vertices[1];
//        auto v2 = faces_approximate_field[i].extend_vertices[2];
        //field_move_vertex[i]
        auto v0 = field_move_vertex[mesh->fast_iGameFace[i].vh(0)];
        auto v1 = field_move_vertex[mesh->fast_iGameFace[i].vh(1)];
        auto v2 = field_move_vertex[mesh->fast_iGameFace[i].vh(2)];
        fprintf(file0, "v %lf %lf %lf\n", v0.x(), v0.y(),v0.z());
        fprintf(file0, "v %lf %lf %lf\n", v1.x(), v1.y(),v1.z());
        fprintf(file0, "v %lf %lf %lf\n", v2.x(), v2.y(),v2.z());
        //fprintf(file0, "f %d %d %d\n", cnt, cnt+1,cnt+2);
    }
    for(int i=0;i<mesh->FaceSize();i++){
        fprintf(file0, "f %d %d %d\n", cnt, cnt+1,cnt+2);
        cnt+=3;

    }
    cnt=0;
    for(int  i =0 ; i< mesh->VertexSize() ; i++){
    //    if(flag) {
//            cout << "origin " <<mesh->fast_iGameVertex[i].x()<<" "<< mesh->fast_iGameVertex[i].y()<<" "<<
//                    mesh->fast_iGameVertex[i].z()<<endl;
//            cout << "updated  " <<field_move_vertex[i].x()<<" "<< field_move_vertex[i].y()<<" "<<
//                    field_move_vertex[i].z()<<endl;
//            cout <<"***************"<< endl;

//            cout <<"***************"<< endl;
            fprintf(file1, "v %lf %lf %lf\n", mesh->fast_iGameVertex[i].x(), mesh->fast_iGameVertex[i].y(),
                   mesh->fast_iGameVertex[i].z());
            fprintf(file1, "v %lf %lf %lf\n", field_move_vertex[i].x(), field_move_vertex[i].y(),
                    field_move_vertex[i].z());

//            fprintf(file0, "v %lf %lf %lf\n", min_move_g[i].x(), min_move_g[i].y(),
//                    min_move_g[i].z());

            fprintf(file1,"l %d %d\n",cnt+1,cnt+2);
          //  fprintf(file0,"l %d %d\n",cnt+1,cnt+3);
            cnt+=2;

       // }
    }

    for(int  i =0 ; i< mesh->VertexSize() ; i++){

       if((mesh->fast_iGameVertex[i] - field_move_vertex[i]).norm()>100){
           fprintf(file1, "v %lf %lf %lf\n", mesh->fast_iGameVertex[i].x(), mesh->fast_iGameVertex[i].y(),mesh->fast_iGameVertex[i].z());
       }

    }


    cgal_polygon = make_shared<CGALPolygon>(mesh);

    std::list <Triangle> triangles;
    unordered_map <Triangle, MeshKernel::iGameFaceHandle, triangle_hash, triangle_equal> triangle_map;
    for (int i = 0; i < mesh->FaceSize(); i++) {
        Point_3 p0(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].x(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].y(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].z());

        Point_3 p1(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].x(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].y(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].z());

        Point_3 p2(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].x(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].y(),
                   mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].z());
        Triangle t(p0, p1, p2);
        triangle_map[t] = MeshKernel::iGameFaceHandle(i);
        triangles.push_back(t);
    }
    cgal_global_aabbtree = Tree(triangles.begin(), triangles.end());

    //Armadillo.obj
    // freopen("../log.txt","w",stdout);

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



    // mesh->init_plane_point();
    //  grid_len = (mesh->BBoxMax.x() - mesh->BBoxMin.x())/50;
    // grid_len = (mesh->BBoxMax.x() - mesh->BBoxMin.x())/50;
    // grid_len = min_move / 2.5;

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
            if (xx >= 0 && yy >= 0 && zz >= 0)
                ret.push_back({xx, yy, zz});
        }
        return ret;
    };
    int fsize = mesh->FaceSize();


//    double avg_area=0;
//    vector<bool>need_tiny_face;
//    need_tiny_face.resize(mesh->FaceSize());
//    for(auto i : mesh->allfaces()){
//        MeshKernel::iGameVertex v0 = mesh->fast_iGameVertex[i.second.vh(0)];
//        MeshKernel::iGameVertex v1 = mesh->fast_iGameVertex[i.second.vh(1)];
//        MeshKernel::iGameVertex v2 = mesh->fast_iGameVertex[i.second.vh(2)];
//        avg_area+=((v1 - v0)%(v2 - v0)).norm()/2;
//    }
//    avg_area/=mesh->FaceSize();
//    for(auto i : mesh->allfaces()){
//        MeshKernel::iGameVertex v0 = mesh->fast_iGameVertex[i.second.vh(0)];
//        MeshKernel::iGameVertex v1 = mesh->fast_iGameVertex[i.second.vh(1)];
//        MeshKernel::iGameVertex v2 = mesh->fast_iGameVertex[i.second.vh(2)];
//        if(((v1 - v0)%(v2 - v0)).norm()/2 < avg_area/10){
//            need_tiny_face[i.first] = true;
//        }
//        else
//            need_tiny_face[i.first] = false;
//    }



/*
    std::vector<std::shared_ptr<std::thread> > bfs_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++){
        bfs_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for(auto iter : build_mesh_vertex){
                if(iter.second%thread_num == now_id){
                    qef_result[iter.second] = qef_grid_pos(iter.first);
                    if(iter.second%100==0) // cut the output time usage
                        printf("succ %d/%d \n",iter.second, build_mesh_vertex.size());
                }
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        bfs_thread_pool[i]->join();
 */
    cout <<"bfs end \n" << endl;
    std::mutex bfs_mutex;
    std::vector <std::shared_ptr<std::thread> > bfs_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++){
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
                            double dist = cgal_vertex_triangle_dist(fh.second, getGridVertex(j, 0), mesh);
                            double move_limit = *set<double>
                                    {(field_move_vertex[fh.second.vh(0)]-mesh->fast_iGameVertex[fh.second.vh(0)]).norm(),
                                     (field_move_vertex[fh.second.vh(1)]-mesh->fast_iGameVertex[fh.second.vh(1)]).norm(),
                                     (field_move_vertex[fh.second.vh(2)]-mesh->fast_iGameVertex[fh.second.vh(2)]).norm()
                                     }.rbegin();
                            if (dist <  max(grid_len ,fh.second.move_dist)*1.1 ) { //TODO : zheli youhua cheng pianyi juli de shiji jisuan
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
    for(auto i : frame_grid_mp_bk){
        for(auto j : container_grid_dir){
            int nx = i.first.x+j[0];
            int ny = i.first.y+j[1];
            int nz = i.first.z+j[2];
            if(nx>=0 && ny>=0 && nz>=0){
                for(auto z : i.second.face_list)
                    frame_grid_mp[{nx,ny,nz}].face_list.push_back(z);
            }
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
                MeshKernel::iGameFaceHandle near_face_handle =
                        triangle_map[*cgal_global_aabbtree.closest_point_and_primitive(
                                Point_3(grid_vertex.x(), grid_vertex.y(), grid_vertex.z())).second];


                each_grid->second.face_list.push_back(near_face_handle);
                for (auto i: mesh->fast_iGameFace[near_face_handle].getSortedVertexHandle()) {
                    for (auto j: mesh->FastNeighborFhOfVertex_[i])
                        each_grid->second.face_list.push_back(j);
                }
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

    std::function<bool(Plane_3 plane,MeshKernel::iGameVertex,double)> vertex_in_plane
    = [&](Plane_3 pla,MeshKernel::iGameVertex v,double eps){
                if(sqrt(squared_distance(pla,Point_3(v.x(),v.y(),v.z())))<eps){
                    return true;
                }
                return false;
    };



    // 上述代码完成距离场建格子的过程 8 ;
    int tt=0;
    long long  sum_grid = 0;
    for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
        if(tt++%1000==0){
            cout << tt <<" // "<<  frame_grid_mp.size() << endl;
        }
      //  if(tt < 35000)continue;
        set<MeshKernel::iGameFaceHandle> face_set;
        vector<MeshKernel::iGameFaceHandle> face_list;
        vector<MeshKernel::iGameFaceHandle> possible_face_list;
        vector<MeshKernel::iGameFace>vf;

        function<bool(const K2::Point_3&)> need_generation_vertex = [&](const K2::Point_3& p){
            MeshKernel::iGameVertex igame_p2(CGAL::to_double(p.x()),
                                             CGAL::to_double(p.y()),
                                             CGAL::to_double(p.z()));

            bool in_tet_field = false;
            for(int field_face_id=0;field_face_id<face_list.size();field_face_id++){
                for(int each_tet_id =0; each_tet_id<6;each_tet_id++ ){
                    if(faces_approximate_field[face_list[field_face_id]].tet_list[each_tet_id].has_on_bounded_side(p)){
                        in_tet_field = true;
                        break;
                    }
                }
                if(in_tet_field) break;
            }
            int side_type = 0;
            //cout <<"******\n" << vf.size() << " "<<face_list.size() << endl;
            cgal_aabbtree_query(vf, igame_p2, mesh, side_type, cgal_polygon);

            if ( (!in_tet_field) && side_type ==1 )
                return true;
            else
                return false;
        };



        function<bool(const MeshKernel::iGameVertex& ,const vector<ApproximateField>&)> need_generation_face_in_tiny =
                [&](const MeshKernel::iGameVertex& v ,const vector<ApproximateField>& field){

            K2::Point_3 p = iGameVertex_to_Point_K2(v);

            bool in_tet_field = false;
            for(int field_face_id=0;field_face_id<field.size();field_face_id++){
                for(int each_tet_id =0; each_tet_id<6;each_tet_id++ ){
                    if(field[field_face_id].tet_list[each_tet_id].has_on_bounded_side(p)){
                        in_tet_field = true;
                        break;
                    }
                }
                if(in_tet_field) break;
            }
            int side_type = 0;
            //cout <<"******\n" << vf.size() << " "<<face_list.size() << endl;
            cgal_aabbtree_query(vf, v, mesh, side_type, cgal_polygon);

            if ( (!in_tet_field) && side_type ==1 )
                return true;
            else
                return false;
        };




        for (int j = 0; j < 8; j++) {
            grid g = getGridFrameVertex(each_grid->first, j);
            unordered_map<grid, GridVertex, grid_hash, grid_equal>::iterator it = frame_grid_mp.find(g);
            if (it == frame_grid_mp.end())continue;
            for (auto k: it->second.face_list) {
                if(!face_set.count(k)){
                    face_set.insert(k);
                    face_list.push_back(k);
                }
            }
        }
        MeshKernel::iGameVertex small = getGridVertex(each_grid->first,0);
        MeshKernel::iGameVertex big = getGridVertex(each_grid->first,7);

        for(MeshKernel::iGameFaceHandle i : face_list){
            vector<MeshKernel::iGameVertex> v;
            for(int j=0;j<3;j++)
                v.push_back(field_move_vertex[mesh->fast_iGameFace[i].vh(j)]);
            if(face_through_grid(small,big,v))
                possible_face_list.push_back(i);
        }

        for(auto i : face_list){
            vf.push_back(mesh->fast_iGameFace.at(i));
        }
        if(possible_face_list.size() ==0 )
            continue;
       // cout << possible_face_list.size() <<"/" << face_list.size() << endl;

        //face_through_grid


        DSU dsu(possible_face_list.size());
        vector<vector<int> >plane_cross;
        plane_cross.resize(possible_face_list.size());
        for(int i=0;i<possible_face_list.size();i++){
            for(int j=i+1;j<possible_face_list.size();j++){
                int root_i = dsu.find_root(i);
                MeshKernel::iGameVertex v0i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(0)];
                MeshKernel::iGameVertex v1i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(1)];
                MeshKernel::iGameVertex v2i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(2)];
                Plane_3 pla_i (Point_3(v0i.x(),v0i.y(),v0i.z()),
                             Point_3(v1i.x(),v1i.y(),v1i.z()),
                             Point_3(v2i.x(),v2i.y(),v2i.z()));
                MeshKernel::iGameVertex v0j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(0)];
                MeshKernel::iGameVertex v1j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(1)];
                MeshKernel::iGameVertex v2j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(2)];


                double edge_eps = sqrt(max({(v0i-v1i).norm(),(v2i-v1i).norm(),(v0i-v2i).norm()}))/15;
                if(vertex_in_plane(pla_i,v0j,edge_eps) &&
                        vertex_in_plane(pla_i,v1j,edge_eps) &&
                        vertex_in_plane(pla_i,v2j,edge_eps)
                        ) { //并且直接相连 ;;;;
                    K2::Triangle_3 tri_face1(K2::Point_3(v0i.x(),v0i.y(),v0i.z()),
                                             K2::Point_3(v1i.x(),v1i.y(),v1i.z()),
                                             K2::Point_3(v2i.x(),v2i.y(),v2i.z()));
                    K2::Triangle_3 tri_face2(K2::Point_3(v0j.x(),v0j.y(),v0j.z()),
                                             K2::Point_3(v1j.x(),v1j.y(),v1j.z()),
                                             K2::Point_3(v2j.x(),v2j.y(),v2j.z()));
                    bool link_flag = false;
                    CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
                            res = intersection(tri_face1, tri_face2);
                    if(res){
                        if (const K2::Segment_3* s = boost::get<K2::Segment_3>(&*res)) {
                            link_flag = true;
                        }else if (const K2:: Triangle_3* ti = boost::get< K2::Triangle_3>(&*res)) {
                            link_flag = true;
                        }else if(const std::vector <K2::Point_3>* vs = boost::get<std::vector <K2::Point_3>>(&*res)){
                            link_flag = true;
                        }
                    }
                    if(link_flag)
                        dsu.join(i,j);
                }
            }
        }

//        for(int i=0;i<face_list.size();i++){
//            cout <<"dsu.find_root(i):" << dsu.find_root(i) << endl;
//        }
        vector<SharpPoint> sharp_point_list;
        for(int i=0;i<possible_face_list.size();i++) {
            MeshKernel::iGameVertex v0i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(0)];
            MeshKernel::iGameVertex v1i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(1)];
            MeshKernel::iGameVertex v2i = field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(2)];
            K2::Triangle_3 tri_i (K2::Point_3(v0i.x(),v0i.y(),v0i.z()),
                                  K2::Point_3(v1i.x(),v1i.y(),v1i.z()),
                                  K2::Point_3(v2i.x(),v2i.y(),v2i.z()));

            for (int j = i + 1; j < possible_face_list.size(); j++) {
                if(dsu.find_root(i) != dsu.find_root(j)){

                    MeshKernel::iGameVertex v0j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(0)];
                    MeshKernel::iGameVertex v1j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(1)];
                    MeshKernel::iGameVertex v2j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(2)];
                    K2::Triangle_3 tri_j (K2::Point_3(v0j.x(),v0j.y(),v0j.z()),
                                          K2::Point_3(v1j.x(),v1j.y(),v1j.z()),
                                          K2::Point_3(v2j.x(),v2j.y(),v2j.z()));

                    CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
                            result = intersection(tri_i, tri_j);
                    if (result) {

                        if (const K2::Point_3 * p = boost::get<K2::Point_3>(&*result)) {
                            // std::cout << (*p) << std::endl;
                            // cout <<"poi: "<< *p << endl;
                            for (int k = j + 1; k < possible_face_list.size(); k++){
                                //cout << "PPPPPPPP " << endl;
                                if(dsu.find_root(i) == dsu.find_root(k) ||
                                    dsu.find_root(j) == dsu.find_root(k))continue;

                                MeshKernel::iGameVertex v0k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(0)];
                                MeshKernel::iGameVertex v1k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(1)];
                                MeshKernel::iGameVertex v2k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(2)];
                                K2::Triangle_3 tri_k (K2::Point_3(v0k.x(),v0k.y(),v0k.z()),
                                                    K2::Point_3(v1k.x(),v1k.y(),v1k.z()),
                                                    K2::Point_3(v2k.x(),v2k.y(),v2k.z()));
                                CGAL::cpp11::result_of<K2::Intersect_3(K2::Point_3, K2::Triangle_3)>::type
                                        result2 = intersection(*p, tri_k);
                                if (result2) {
                                    //  cout << "GGGGGGG2 " << endl;
                                    if (const K2::Point_3 * p2 = boost::get<K2::Point_3>(&*result2)) {
                                        if ( need_generation_vertex(*p2) ) {
                                            SharpPoint sp;
                                            sp.p = *p2;
                                            sp.source_face_local_id.push_back(i);
                                            sp.source_face_local_id.push_back(j);
                                            sp.source_face_local_id.push_back(k);
                                            sharp_point_list.push_back(sp);
                                        }
                                        //delete p2;
                                    }
                                }
                            }
                            //delete p;
                        }
                        if (const K2::Segment_3 * s = boost::get<K2::Segment_3 >(&*result)) {
                            plane_cross[dsu.find_root(i)].push_back(dsu.find_root(j));
                            plane_cross[dsu.find_root(j)].push_back(dsu.find_root(i));


                            // std::cout << (*p) << std::endl;
//                            cout <<"tri:" <<tri_i.has_on(s->point(0))<<" "<< tri_i.has_on(s->point(1))<<" "
//                            <<  tri_j.has_on(s->point(0))<<" "<< tri_j.has_on(s->point(1)) << endl;
//                             cout << "tri i "<< tri_i << endl;
//                            cout << "tri j "<< tri_j << endl;
                            for (int k = j + 1; k < possible_face_list.size(); k++){
                                if(dsu.find_root(i) == dsu.find_root(k) ||
                                   dsu.find_root(j) == dsu.find_root(k))continue;

                                //cout << "GGGGGGG " << endl;

                                MeshKernel::iGameVertex v0k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(0)];
                                MeshKernel::iGameVertex v1k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(1)];
                                MeshKernel::iGameVertex v2k = field_move_vertex[mesh->fast_iGameFace[possible_face_list[k]].vh(2)];
                                K2::Triangle_3 tri_k (K2::Point_3(v0k.x(),v0k.y(),v0k.z()),
                                                      K2::Point_3(v1k.x(),v1k.y(),v1k.z()),
                                                      K2::Point_3(v2k.x(),v2k.y(),v2k.z()));

                               // cout << "seg s " << *s << endl;
                               // cout << "tri_k  " << tri_k << endl;
                                CGAL::cpp11::result_of<K2::Intersect_3(K2::Segment_3, K2::Triangle_3)>::type
                                        result2 = intersection(*s, tri_k);
                                if (result2) {
                                  //  cout << "GGGGGGG2 " << endl;
                                    if (const K2::Point_3 * p2 = boost::get<K2::Point_3>(&*result2)) {

                                        if (need_generation_vertex(*p2) ) {
                                            SharpPoint sp;
                                            sp.p = *p2;
                                            sp.source_face_local_id.push_back(i);
                                            sp.source_face_local_id.push_back(j);
                                            sp.source_face_local_id.push_back(k);

//                                            if(sqrt(CGAL::to_double(squared_distance(tri_i,*p2)))>myeps ){
//                                                cout <<"errr!!???" << sqrt(CGAL::to_double(squared_distance(tri_i,*p2))) << endl;
//                                            }
//                                            if(sqrt(CGAL::to_double(squared_distance(tri_j,*p2)))>myeps ){
//                                                cout <<"errr!!???" << sqrt(CGAL::to_double(squared_distance(tri_j,*p2))) << endl;
//                                            }
//                                            if(sqrt(CGAL::to_double(squared_distance(tri_k,*p2)))>myeps ){
//                                                cout <<"errr!!???" << sqrt(CGAL::to_double(squared_distance(tri_k,*p2))) << endl;
//                                            }


                                            sharp_point_list.push_back(sp);
                                        }
                                       // delete p2;
                                    }
                                }
                            }
                           // delete s;
                        }
                    }
                }
            }
        }


        //cout << "sharp_point_list.size(): "<< sharp_point_list.size() << endl;

        vector<SharpPoint> sharp_point_list_merge;
        for(int i=0;i<sharp_point_list.size();i++){
            bool flag = false;
            for(int j=0;j<sharp_point_list_merge.size() && !flag; j++){
                if(sqrt(CGAL::to_double(squared_distance(sharp_point_list_merge[j].p,sharp_point_list[i].p)))<1e-4){
                    flag = true;
                    for(auto k : sharp_point_list[i].source_face_local_id)
                        sharp_point_list_merge[j].source_face_local_id.push_back(k);
                }
            }
            if(!flag){
                sharp_point_list_merge.push_back(sharp_point_list[i]);
            }
        }
        for(int i=0;i<sharp_point_list_merge.size();i++){
            sort(sharp_point_list_merge[i].source_face_local_id.begin(),sharp_point_list_merge[i].source_face_local_id.end());
            sharp_point_list_merge[i].source_face_local_id.resize(
                    unique(sharp_point_list_merge[i].source_face_local_id.begin(),
                           sharp_point_list_merge[i].source_face_local_id.end())-sharp_point_list_merge[i].source_face_local_id.begin());
        }

        swap(sharp_point_list_merge,sharp_point_list);

//        for(auto it : sharp_point_list){
//            for(auto cross_face_ph : it.source_face_local_id){
//                K2::Triangle_3 trii(K2::Point_3(field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(0)].x(),
//                                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(0)].y(),
//                                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(0)].z()),
//                                    K2::Point_3(field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(1)].x(),
//                                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(1)].y(),
//                                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(1)].z()),
//                                    K2::Point_3(field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(2)].x(),
//                                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(2)].y(),
//                                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(2)].z()));
//                cout << cross_face_ph<<" * "<<sqrt(CGAL::to_double(squared_distance(it.p,trii))) <<" "
//                     <<sqrt(CGAL::to_double(squared_distance(it.p,trii.supporting_plane()))) << endl;
//            }
//            cout <<"**************"<<endl;
//            exit(0);
//        }

        for(auto it : sharp_point_list){

            fprintf(file2,"v %lf %lf %lf\n",CGAL::to_double(it.p.x()),CGAL::to_double(it.p.y()),CGAL::to_double(it.p.z()));
        }


        // 注意合并处理 ;



        double maxlen = grid_len;
        for(int i=0;i<sharp_point_list.size();i++){
            for(int j=i+1;j<sharp_point_list.size();j++){
                double dist = sqrt(CGAL::to_double(squared_distance(sharp_point_list[i].p,sharp_point_list[j].p))/3.0);
                maxlen = min(maxlen , dist);
            }
            for (int j = 0; j < possible_face_list.size(); j++) {
                bool flag = true;
                for(int k=0;k<sharp_point_list[i].source_face_local_id.size();k++){
                    if( dsu.find_root(sharp_point_list[i].source_face_local_id[k]) == dsu.find_root(j))
                        flag = false;
                }
                if (flag) {
                    MeshKernel::iGameVertex v0j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(0)];
                    MeshKernel::iGameVertex v1j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(1)];
                    MeshKernel::iGameVertex v2j = field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(2)];
                    K2::Triangle_3 tri(K2::Point_3(v0j.x(),v0j.y(),v0j.z()),
                                   K2::Point_3(v1j.x(),v1j.y(),v1j.z()),
                                   K2::Point_3(v2j.x(),v2j.y(),v2j.z()));
                    double dist = sqrt(CGAL::to_double(squared_distance(sharp_point_list[i].p,tri))/3.0);
                    maxlen = min(maxlen , dist);
                }
            }
        }
        int num = 1 ;
        if(sharp_point_list.size() > 1 ){
            //cout << sharp_point_list.size() <<" "<< maxlen <<  endl;
            int div = int((grid_len / maxlen)+1);
            int d = 0;
            int tmp = div -1 ;
            while(tmp){
                tmp>>=1;
                d++;
            }
            num  = 1 << d;
           // cout << num << endl;
            sum_grid += 1LL* num*num*num;
        }
        else
            sum_grid += 1;

        if(num > 10){
            //cout << num << endl;
            //continue;
            num = 8;
        }
        num = 1;


        ;
        double tiny_grid_len = grid_len/num;
        // 边界先都算上




//        vector<MeshKernel::iGameVertex >sharp_point_list_vertex;
//        for(int i=0;i<sharp_point_list.size();i++){
//            double x = CGAL::to_double(sharp_point_list[i].p.x());
//            double y = CGAL::to_double(sharp_point_list[i].p.y());
//            double z = CGAL::to_double(sharp_point_list[i].p.z());
//            sharp_point_list_vertex.emplace_back(x,y,z);
//        }

        std::function<vector<pair<MeshKernel::iGameVertex,MeshKernel::iGameVertex > >(
                MeshKernel::iGameVertex , MeshKernel::iGameVertex , MeshKernel::iGameVertex ,
                MeshKernel::iGameVertex , MeshKernel::iGameVertex  ) > get_grid_intersect_triangle
                = [&](MeshKernel::iGameVertex v0, MeshKernel::iGameVertex v1, MeshKernel::iGameVertex v2,
                          MeshKernel::iGameVertex tinygrid_st, MeshKernel::iGameVertex  tinygrid_end){
            vector<pair<MeshKernel::iGameVertex,MeshKernel::iGameVertex > > ret;
            vector<MeshKernel::iGameVertex >vs;

            for(int j=0;j<8;j++){
                vs.push_back(getTinyGridVertex(tinygrid_st,tinygrid_end,j));
            }
            // 这里采用大三角形法： ;
            K2::Triangle_3 tri(iGameVertex_to_Point_K2(v0),
                               iGameVertex_to_Point_K2((v1-v0).normalize()*tiny_grid_len*2+ v0),
                               iGameVertex_to_Point_K2((v2-v0).normalize()*tiny_grid_len*2+ v0));
            for(auto i : container_grid_face){
                vector<MeshKernel::iGameVertex> cross_vertex;
                for(auto j : vector<vector<int> >{{0,1,2},{2,3,0}}) {
                    K2::Triangle_3 grid_tri(iGameVertex_to_Point_K2(vs[i[j[0]]]),
                                      iGameVertex_to_Point_K2(vs[i[j[1]]]),
                                      iGameVertex_to_Point_K2(vs[i[j[2]]]));

                    const auto result = intersection(tri, grid_tri);
                    if(result) {
                        if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*result)) {
                            MeshKernel::iGameVertex v1_f = Point_K2_to_iGameVertex(s->vertex(0));
                            MeshKernel::iGameVertex v2_f = Point_K2_to_iGameVertex(s->vertex(1));
                            MeshKernel::iGameVertex cross_product1 = (v1_f - v0) % (v2_f - v0);
                            MeshKernel::iGameVertex cross_product2 = (v1 - v0) % (v2 - v0) ;
                            if(cross_product1 * cross_product2 < 0)
                                swap(v1_f,v2_f);
                         //   ret.emplace_back(v1_f,v2_f);
                            cross_vertex.push_back(v1_f);
                            cross_vertex.push_back(v2_f);
                        }//todo : fix with point cut 这里还有bug !!!!!!!!!!!1;
                    }
                }
                //fixme： update
                if(cross_vertex.size()>=2){
                    int id1 =furthest_vertex_in_vector(0,cross_vertex);
                    int id2 =furthest_vertex_in_vector(id1,cross_vertex);
                    MeshKernel::iGameVertex v1_f = cross_vertex[id1];
                    MeshKernel::iGameVertex v2_f =  cross_vertex[id2];
                    MeshKernel::iGameVertex cross_product1 = (v1_f - v0) % (v2_f - v0);
                    MeshKernel::iGameVertex cross_product2 = (v1 - v0) % (v2 - v0) ;
                    if(cross_product1 * cross_product2 < 0)
                        swap(v1_f,v2_f);
                    ret.emplace_back(v1_f,v2_f);
                }

            }
                    //getTinyGridVertex
                    // container_grid_face


            return ret;
        };


        std::function<vector<pair<MeshKernel::iGameVertex,MeshKernel::iGameVertex > >(
                MeshKernel::iGameVertex , MeshKernel::iGameVertex , MeshKernel::iGameVertex ,
                MeshKernel::iGameVertex , MeshKernel::iGameVertex  ,MeshKernel::iGameVertex &) > get_grid_intersect_plane
                = [&](MeshKernel::iGameVertex v0, MeshKernel::iGameVertex v1, MeshKernel::iGameVertex v2,
                      MeshKernel::iGameVertex tinygrid_st, MeshKernel::iGameVertex  tinygrid_end,MeshKernel::iGameVertex & center){
                    vector<pair<MeshKernel::iGameVertex,MeshKernel::iGameVertex > > ret;
                    vector<MeshKernel::iGameVertex >vs;

                    for(int j=0;j<8;j++){
                        vs.push_back(getTinyGridVertex(tinygrid_st,tinygrid_end,j));
                    }
                    // 这里采用大三角形法： ;
                    K2::Plane_3 pla(iGameVertex_to_Point_K2(v0),
                                       iGameVertex_to_Point_K2(v1),
                                       iGameVertex_to_Point_K2(v2));
                    for(auto i : container_grid_face){
                        vector<MeshKernel::iGameVertex> cross_vertex;
                        for(auto j : vector<vector<int> >{{0,1,2},{2,3,0}}) {
                            K2::Triangle_3 grid_tri(iGameVertex_to_Point_K2(vs[i[j[0]]]),
                                                    iGameVertex_to_Point_K2(vs[i[j[1]]]),
                                                    iGameVertex_to_Point_K2(vs[i[j[2]]]));

                            const auto result = intersection(pla, grid_tri);
                            if(result) {
                                if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*result)) {
                                    MeshKernel::iGameVertex v1_f = Point_K2_to_iGameVertex(s->vertex(0));
                                    MeshKernel::iGameVertex v2_f = Point_K2_to_iGameVertex(s->vertex(1));
//                                    MeshKernel::iGameVertex cross_product1 = (v1_f - v0) % (v2_f - v0);
//                                    MeshKernel::iGameVertex cross_product2 = (v1 - v0) % (v2 - v0) ;
//                                    if(cross_product1 * cross_product2 < 0)
//                                        swap(v1_f,v2_f);
                                    //ret.emplace_back(v1_f,v2_f);
                                    cross_vertex.push_back(v1_f);
                                    cross_vertex.push_back(v2_f);
                                }//todo : fix with point cut 这里还有bug !!!!!!!!!!!1;
                            }
                        }
                        if(cross_vertex.size()>=2){
                            int id1 =furthest_vertex_in_vector(0,cross_vertex);
                            int id2 =furthest_vertex_in_vector(id1,cross_vertex);
                            MeshKernel::iGameVertex v1_f = cross_vertex[id1];
                            MeshKernel::iGameVertex v2_f =  cross_vertex[id2];
//                            MeshKernel::iGameVertex cross_product1 = (v1_f - v0) % (v2_f - v0);此处不需要;
//                            MeshKernel::iGameVertex cross_product2 = (v1 - v0) % (v2 - v0) ;
//                            if(cross_product1 * cross_product2 < 0)
//                                swap(v1_f,v2_f);
                            ret.emplace_back(v1_f,v2_f);
                        }
                    }
                    center = MeshKernel::iGameVertex(0,0,0);
                    for(auto i : ret){
                        center += i.first + i.second;
                    }
                    MeshKernel::iGameVertex cross_product1 = (v1 - v0) % (v2 - v0);
                    for(int i=0;i<ret.size();i++){
                          MeshKernel::iGameVertex cross_product2 = (ret[i].first - center) % (ret[i].second - center);
                    }
                    center /= ret.size()*2;

                    //getTinyGridVertex
                    // container_grid_face


                    return ret;
                };



        std::function<vector<vector<MeshKernel::iGameVertex> >(
                MeshKernel::iGameVertex , MeshKernel::iGameVertex , MeshKernel::iGameVertex ,
                        MeshKernel::iGameVertex , MeshKernel::iGameVertex) > get_grid_intersect_plane_full
                                                                  = [&](MeshKernel::iGameVertex v0, MeshKernel::iGameVertex v1, MeshKernel::iGameVertex v2,
                                                                        MeshKernel::iGameVertex tinygrid_st, MeshKernel::iGameVertex  tinygrid_end){
                    vector<vector<MeshKernel::iGameVertex> > ret;
                    vector<MeshKernel::iGameVertex >vs;

                    for(int j=0;j<8;j++){
                        vs.push_back(getTinyGridVertex(tinygrid_st,tinygrid_end,j));
                    }
                    // 这里采用大三角形法： ;
                    K2::Plane_3 pla(iGameVertex_to_Point_K2(v0),
                                    iGameVertex_to_Point_K2(v1),
                                    iGameVertex_to_Point_K2(v2));
                    for(auto i : container_grid_face){
                        ret.emplace_back();
                        vector<MeshKernel::iGameVertex>cross_vertex;
                        for(auto j : vector<vector<int> >{{0,1,2},{2,3,0}}) {
                            K2::Triangle_3 grid_tri(iGameVertex_to_Point_K2(vs[i[j[0]]]),
                                                    iGameVertex_to_Point_K2(vs[i[j[1]]]),
                                                    iGameVertex_to_Point_K2(vs[i[j[2]]]));

                            const auto result = intersection(pla, grid_tri);
                            if(result) {
                                if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*result)) {
                                    MeshKernel::iGameVertex v1_f = Point_K2_to_iGameVertex(s->vertex(0));
                                    MeshKernel::iGameVertex v2_f = Point_K2_to_iGameVertex(s->vertex(1));
//                                    MeshKernel::iGameVertex cross_product1 = (v1_f - v0) % (v2_f - v0);
//                                    MeshKernel::iGameVertex cross_product2 = (v1 - v0) % (v2 - v0) ;
//                                    if(cross_product1 * cross_product2 < 0)
//                                        swap(v1_f,v2_f);
                                    if((v1_f - v2_f).norm() > myeps){
                                        cross_vertex.push_back(v1_f);
                                        cross_vertex.push_back(v2_f);
                                    }
                                }//todo : fix with point cut 这里还有bug !!!!!!!!!!!1;
                            }
                        }
                        if(cross_vertex.size()>=2){
                            int id1 =furthest_vertex_in_vector(0,cross_vertex);
                            int id2 =furthest_vertex_in_vector(id1,cross_vertex);
                            MeshKernel::iGameVertex v1_f = cross_vertex[id1];
                            MeshKernel::iGameVertex v2_f =  cross_vertex[id2];
//                            MeshKernel::iGameVertex cross_product1 = (v1_f - v0) % (v2_f - v0);此处不需要;
//                            MeshKernel::iGameVertex cross_product2 = (v1 - v0) % (v2 - v0) ;
//                            if(cross_product1 * cross_product2 < 0)
//                                swap(v1_f,v2_f);
                            ret.rbegin()->push_back(v1_f);
                            ret.rbegin()->push_back(v2_f);
                        }
                        else if(cross_vertex.size()){
                            ret.rbegin()->push_back(*cross_vertex.begin());
                        }




                    }


                    //getTinyGridVertex
                    // container_grid_face


                    return ret;
                };

        num = 1;
        static int f5_id = 1;
        auto ccc = (small + big) /2 ;
        fprintf(file5,"v %lf %lf %lf\n",ccc.x(),ccc.y(),ccc.z());


        static int vid = 1;
        for(int x =0; x < num; x++){
            for(int y =0; y < num; y++){
                for(int z = 0; z < num; z++){

                    //int type3_id = -1;
                    MeshKernel::iGameVertex tinygrid_st = small
                            + MeshKernel::iGameVertex(x * tiny_grid_len, y * tiny_grid_len, z * tiny_grid_len);
                    MeshKernel::iGameVertex  tinygrid_end = small
                            + MeshKernel::iGameVertex( (x+1) * tiny_grid_len, (y+1) * tiny_grid_len,(z+1) * tiny_grid_len);

                    MeshKernel::iGameVertex tinygrid_center = (tinygrid_st + tinygrid_end)/2;

                    vector<SharpPoint>tiny_grid_sharp;
                   // map<MeshKernel::iGameFaceHandle,int> tiny_grid_face_state;
                    vector<int>local_face_state;
                    local_face_state.resize(possible_face_list.size());
                    fill(local_face_state.begin(),local_face_state.end(),1);

                    for(int i= 0 ;i < sharp_point_list.size();i++){
                        if(vertex_in_tiny_grid(tinygrid_st,tinygrid_end, Point_K2_to_iGameVertex(sharp_point_list[i].p))){
                            tiny_grid_sharp.push_back(sharp_point_list[i]);
                            for(auto j : sharp_point_list[i].source_face_local_id)
                                local_face_state[j] = 3;

//                            for(auto j : sharp_point_list[i].source_face_local_id){
//                                tiny_grid_face_state[possible_face_list[j]] = 3;
//                            }
                        }
                    }
                    //continue; ok
                    vector<int>built_face;
                    built_face.resize(possible_face_list.size());

                    if(tiny_grid_sharp.size() > 0) { // is type3 tiny grid
                        fprintf(fileis3,"v %lf %lf %lf\n",tinygrid_center.x(),tinygrid_center.y(),tinygrid_center.z());
                        MeshKernel::iGameVertex build_point(0,0,0);
                        vector<ApproximateField> approximates_field_lists;
                        for(int i = 0; i < tiny_grid_sharp.size();i++){
                            MeshKernel::iGameVertex this_sharp_vertex = Point_K2_to_iGameVertex(tiny_grid_sharp[i].p);
                            build_point = build_point + this_sharp_vertex;
                        }
                        build_point /= tiny_grid_sharp.size();

                        vector<pair<MeshKernel::iGameVertex, MeshKernel::iGameVertex > > triangles_point;
                        vector<vector<MeshKernel::iGameVertex> > origin_triangles_point;
                        function<void(MeshKernel::iGameVertex, MeshKernel::iGameVertex,MeshKernel::iGameVertex,MeshKernel::iGameFaceHandle)>
                                insert_into_triangles_point
                        = [&](MeshKernel::iGameVertex sharp, MeshKernel::iGameVertex v0, MeshKernel::iGameVertex v1,
                                MeshKernel::iGameFaceHandle origin_face_handle){
                                    auto frame_intersect_result = get_grid_intersect_triangle(
                                            sharp,
                                            v0,
                                            v1,
                                            tinygrid_st,
                                            tinygrid_end
                                    );

                            for(const auto& r : frame_intersect_result) {
                                        triangles_point.push_back(r);// 超新距离长应该加这里吧！！！！;还需要原始的平面三点;
                                        approximates_field_lists.push_back(ApproximateField(origin_face_handle,
                                                                                            {build_point,
                                                                                             r.first,
                                                                                             r.second}));
//                                        approximates_field_lists.emplace_back(origin_face_handle,
//                                                                              {build_point,
//                                                                               r.first,
//                                                                               r.second});

                                    }




                        };//改这个函数，实现功能.
                        // 此处先不考虑切割;



                        for(int i = 0; i < tiny_grid_sharp.size();i++) {
                            // 判断p是定点，边上的点，内部点;
                            MeshKernel::iGameVertex this_sharp_vertex = Point_K2_to_iGameVertex(tiny_grid_sharp[i].p);
                            for(auto local_source_face_id : tiny_grid_sharp[i].source_face_local_id){
                                MeshKernel::iGameFaceHandle fh = possible_face_list[local_source_face_id];
                                bool is_vertex = false;
                                for(int vid = 0;vid <3 && !is_vertex;vid++){
                                    if(this_sharp_vertex.dist(field_move_vertex[mesh->fast_iGameFace[fh].vh(vid)]) < myeps*10){
                                        insert_into_triangles_point(
                                                this_sharp_vertex,
                                                field_move_vertex[mesh->fast_iGameFace[fh].vh((vid+1)%3)],
                                                field_move_vertex[mesh->fast_iGameFace[fh].vh((vid+2)%3)],fh);
                                        is_vertex = true;
                                    }

                                }

                                if(is_vertex)continue;
                                bool is_line = false;
                                for(int vid = 0;vid <3 && !is_line;vid++) {
                                    MeshKernel::iGameVertex v1 = field_move_vertex[mesh->fast_iGameFace[fh].vh(vid)];
                                    MeshKernel::iGameVertex v2 = field_move_vertex[mesh->fast_iGameFace[fh].vh((vid+1)%3)];
                                    K::Line_3 l(iGameVertex_to_Point(v1),iGameVertex_to_Point(v2));
                                    if(sqrt(CGAL::squared_distance(l,iGameVertex_to_Point(this_sharp_vertex)) )< myeps*10){
                                        is_line = true;// 这里是变化面啊 ！！！！！;
                                        insert_into_triangles_point(this_sharp_vertex,
                                                                    field_move_vertex[mesh->fast_iGameFace[fh].vh((vid+1)%3)],
                                                                    field_move_vertex[mesh->fast_iGameFace[fh].vh((vid+2)%3)],fh);
                                        insert_into_triangles_point(this_sharp_vertex,
                                                                    field_move_vertex[mesh->fast_iGameFace[fh].vh((vid+2)%3)],
                                                                    field_move_vertex[mesh->fast_iGameFace[fh].vh(vid)],fh);
                                    }
                                }
                                if(is_line)continue;
                                insert_into_triangles_point(this_sharp_vertex,
                                                            field_move_vertex[mesh->fast_iGameFace[fh].vh(0)],
                                                            field_move_vertex[mesh->fast_iGameFace[fh].vh(1)],fh);
                                insert_into_triangles_point(this_sharp_vertex,
                                                            field_move_vertex[mesh->fast_iGameFace[fh].vh(1)],
                                                            field_move_vertex[mesh->fast_iGameFace[fh].vh(2)],fh);
                                insert_into_triangles_point(this_sharp_vertex,
                                                            field_move_vertex[mesh->fast_iGameFace[fh].vh(2)],
                                                            field_move_vertex[mesh->fast_iGameFace[fh].vh(0)],fh);


//                                static int vid = 1;
//                                auto v0 =  field_move_vertex[mesh->fast_iGameFace[fh].vh(0)];
//                                auto v1 =  field_move_vertex[mesh->fast_iGameFace[fh].vh(1)];
//                                auto v2 =  field_move_vertex[mesh->fast_iGameFace[fh].vh(2)];
//                                fprintf(file10,"v %lf %lf %lf \n",this_sharp_vertex.x(),this_sharp_vertex.y(),this_sharp_vertex.z());
//                                fprintf(file10,"v %lf %lf %lf \n",v0.x(),v0.y(),v0.z());
//                                fprintf(file10,"v %lf %lf %lf \n",v1.x(),v1.y(),v1.z());
//                                fprintf(file10,"f %d %d %d\n",vid,vid+1,vid+2);vid+=3;
//                                fprintf(file10,"v %lf %lf %lf \n",this_sharp_vertex.x(),this_sharp_vertex.y(),this_sharp_vertex.z());
//                                fprintf(file10,"v %lf %lf %lf \n",v1.x(),v1.y(),v1.z());
//                                fprintf(file10,"v %lf %lf %lf \n",v2.x(),v2.y(),v2.z());
//                                fprintf(file10,"f %d %d %d\n",vid,vid+1,vid+2);vid+=3;
//                                fprintf(file10,"v %lf %lf %lf \n",this_sharp_vertex.x(),this_sharp_vertex.y(),this_sharp_vertex.z());
//                                fprintf(file10,"v %lf %lf %lf \n",v2.x(),v2.y(),v2.z());
//                                fprintf(file10,"v %lf %lf %lf \n",v0.x(),v0.y(),v0.z());
//                                fprintf(file10,"f %d %d %d\n",vid,vid+1,vid+2);vid+=3;
                            }

                            // 这里还要2 轮修正 ！！！ ; 用大三角形状求交法 :!!!! 已经修正完了;;
                        }

                        vector<int>is_valid;
                        is_valid.resize(triangles_point.size());
                        for(int i=0;i<triangles_point.size();i++){
                            MeshKernel::iGameVertex i_v0 = build_point;
                            MeshKernel::iGameVertex i_v1 = triangles_point[i].first;
                            MeshKernel::iGameVertex i_v2 = triangles_point[i].second;
                            MeshKernel::iGameVertex i_center = (i_v0 + i_v1 + i_v2) /3;
                            for(int j=0;j<triangles_point.size();j++){
                                if(i==j)continue;
                                MeshKernel::iGameVertex j_v0 = build_point;
                                MeshKernel::iGameVertex j_v1 = triangles_point[j].first;
                                MeshKernel::iGameVertex j_v2 = triangles_point[j].second;
                                MeshKernel::iGameVertex j_center = (j_v0 + j_v1 + j_v2) /3;
                                MeshKernel::iGameVertex origin_dir = (j_v1 - j_v0) % (j_v2 - j_v0);
                                //todo 这里改成tet法

                            }
                        }


                        int cid = vid;
                        fprintf(file10,"v %lf %lf %lf\n",build_point.x(),build_point.y(),build_point.z());
                        vid++;
                        for(auto p : triangles_point){
                            MeshKernel::iGameVertex this_face_center = (p.first + p.second + build_point) /3;
                            bool flag = need_generation_face_in_tiny(this_face_center,approximates_field_lists);
                            if(!flag)continue;

                            fprintf(file10,"v %lf %lf %lf \n",p.first.x(),p.first.y(),p.first.z());
                            fprintf(file10,"v %lf %lf %lf \n",p.second.x(),p.second.y(),p.second.z());
                            fprintf(file10,"f %d %d %d \n",cid,vid,vid+1);
                            vid+=2;
                        }
                        continue;
                    }
                    static int f3_id = 1;
                    for(int ii=0;ii<7;ii++){
                        for(int jj=0;jj<DirectedGridEdge[ii].size();jj++){
                            int from = ii;
                            int to = DirectedGridEdge[ii][jj];
                            MeshKernel::iGameVertex fv = getTinyGridVertex(tinygrid_st,tinygrid_end,from);
                            MeshKernel::iGameVertex tv = getTinyGridVertex(tinygrid_st,tinygrid_end,to);
                            fprintf(file3,"v %lf %lf %lf\n",fv.x(),fv.y(),fv.z());
                            fprintf(file3,"v %lf %lf %lf\n",tv.x(),tv.y(),tv.z());
                            fprintf(file3,"l %d %d\n",f3_id,f3_id+1);
                            f3_id+=2;
                        }
                    }


                    set<int>se;
                    for(int i=0;i<possible_face_list.size();i++){
                        if(face_through_grid(tinygrid_st,tinygrid_end,{
                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(0)],
                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(1)],
                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(2)]
                        })) {
                            se.insert((int)dsu.find_root(i));
                        }
                    }
                    if(se.size()==1){
                        MeshKernel::iGameVertex center(0,0,0);
                        fprintf(fileis1,"v %lf %lf %lf\n",tinygrid_center.x(),tinygrid_center.y(),tinygrid_center.z());

                        int fid = *se.begin();
                        MeshKernel::iGameVertex origin_direct =
                                (field_move_vertex[mesh->fast_iGameFace[possible_face_list[fid]].vh(1)] -
                                        field_move_vertex[mesh->fast_iGameFace[possible_face_list[fid]].vh(0)]) %
                                        (field_move_vertex[mesh->fast_iGameFace[possible_face_list[fid]].vh(2)] -
                                        field_move_vertex[mesh->fast_iGameFace[possible_face_list[fid]].vh(0)]);
                        auto res = get_grid_intersect_plane(
                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[fid]].vh(0)],
                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[fid]].vh(1)],
                                field_move_vertex[mesh->fast_iGameFace[possible_face_list[fid]].vh(2)],
                                tinygrid_st,
                                tinygrid_end,
                                center
                        ); //此处是1类面的生成 !! ;;;;;;
                        //cout << "c info : "<<center.x() <<" "<< center.y()<<" "<< center.z() << endl;

                        int cid = vid;
                        fprintf(file10, "v %lf %lf %lf\n", center.x(), center.y(), center.z());
                        //fprintf(file8, "v %lf %lf %lf\n", center.x(), center.y(), center.z());
                        vid++;
                        for(auto p : res) {
                            MeshKernel::iGameVertex now_direct = (p.first - center) % (p.second - center);
                            if(now_direct * origin_direct < 0)
                                swap(p.first,p.second);
                            fprintf(file10, "v %lf %lf %lf \n", p.first.x(), p.first.y(), p.first.z());
                            fprintf(file10, "v %lf %lf %lf \n", p.second.x(), p.second.y(), p.second.z());
                            fprintf(file10, "f %d %d %d \n", cid, vid, vid + 1);
                            vid += 2;
                        }

                    }
                    else{
                        fprintf(fileis2,"v %lf %lf %lf\n",tinygrid_center.x(),tinygrid_center.y(),tinygrid_center.z());
                        std::vector<MeshKernel::iGameFaceHandle>tiny_grid_possible_face;
                        for(auto i :se){
                            tiny_grid_possible_face.push_back(possible_face_list[i]);
                        }

                        vector<MeshKernel::iGameVertex >vs;
                        for(int j=0;j<8;j++){
                            vs.push_back(getTinyGridVertex(tinygrid_st,tinygrid_end,j));
                        }

                        for(int local_id=0;local_id<tiny_grid_possible_face.size();local_id++) {
                            auto res = get_grid_intersect_plane_full(
                                    field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_id]].vh(0)],
                                    field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_id]].vh(1)],
                                    field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_id]].vh(2)],
                                    tinygrid_st,
                                    tinygrid_end
                            );


                            MeshKernel::iGameVertex origin_direct = (field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_id]].vh(1)] -
                                    field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_id]].vh(0)]) %
                                    (field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_id]].vh(2)] -
                                            field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_id]].vh(0)]);

                            K2::Plane_3 local_face_pla(iGameVertex_to_Point_K2(field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_id]].vh(0)]),
                                                       iGameVertex_to_Point_K2(field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_id]].vh(1)]),
                                                       iGameVertex_to_Point_K2(field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_id]].vh(2)])
                                                       );
                            for(int local_other_id=0;local_other_id<tiny_grid_possible_face.size();local_other_id++){
                                if(local_id == local_other_id)continue;
                                K2::Plane_3 local_other_face_pla(iGameVertex_to_Point_K2(field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_other_id]].vh(0)]),
                                                           iGameVertex_to_Point_K2(field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_other_id]].vh(1)]),
                                                           iGameVertex_to_Point_K2(field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_other_id]].vh(2)])
                                );
                                const auto result = intersection(local_face_pla, local_other_face_pla);
                                if(result) {
                                    int grid_face_id = -1;
                                    if (const K2::Line_3 *s = boost::get<K2::Line_3>(&*result)) {
                                        for(auto i : container_grid_face){
                                            vector<MeshKernel::iGameVertex> cross_vertex;
                                            grid_face_id++;
                                            for(auto j : vector<vector<int> >{{0,1,2},{2,3,0}}) {
                                                K2::Triangle_3 grid_tri(iGameVertex_to_Point_K2(vs[i[j[0]]]),
                                                                        iGameVertex_to_Point_K2(vs[i[j[1]]]),
                                                                        iGameVertex_to_Point_K2(vs[i[j[2]]]));

                                                const auto result2 = intersection(grid_tri, *s);
                                                if(result2){
                                                    if (const K2::Point_3 *p = boost::get<K2::Point_3>(&*result2)){
                                                        res[grid_face_id].push_back(Point_K2_to_iGameVertex(*p));
                                                    }
                                                    else if (const K2::Segment_3 *seg = boost::get<K2::Segment_3>(&*result2)){
                                                        res[grid_face_id].push_back(Point_K2_to_iGameVertex(seg->vertex(0)));
                                                        res[grid_face_id].push_back(Point_K2_to_iGameVertex(seg->vertex(1)));
                                                    }
                                                }
                                            }

                                        }

                                    }//todo : fix with point cut 这里还有bug !!!!!!!!!!!1;
                                }
                            }
                            vector<MeshKernel::iGameVertex >final_res;
                            for(auto i : res){
                                vector<MeshKernel::iGameVertex >build_vertex;
                                for(auto j : i){
                                    bool flag = true;
                                    for(int local_other_id=0;local_other_id<tiny_grid_possible_face.size();local_other_id++){
                                        if(local_id == local_other_id )continue;
                                        MeshKernel::iGameVertex v0 = field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_other_id]].vh(0)];
                                        MeshKernel::iGameVertex v1 = field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_other_id]].vh(1)];
                                        MeshKernel::iGameVertex v2 = field_move_vertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_other_id]].vh(2)];



                                        MeshKernel::iGameVertex ov0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_other_id]].vh(0)];
                                        MeshKernel::iGameVertex ov1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_other_id]].vh(1)];
                                        MeshKernel::iGameVertex ov2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[tiny_grid_possible_face[local_other_id]].vh(2)];





                                        MeshKernel::iGameVertex center = (v0 + v1 + v2)/3;
                                        MeshKernel::iGameVertex normal = (v1 - v0) % (v2 - v0);
                                        MeshKernel::iGameVertex normal_o = (ov1 - ov0) % (ov2 - ov0);
                                        if(normal * normal_o < myeps){
                                            normal *=-1;
                                            flag = false;
                                        }

                                        if (sqrt(CGAL::squared_distance(Plane_3(
                                                iGameVertex_to_Point(v0),
                                                iGameVertex_to_Point(v1),
                                                iGameVertex_to_Point(v2)),iGameVertex_to_Point(j))) > myeps){
                                            if( (j - center) * normal > 0) {
                                                flag = false;

                                            }
                                        }
                                    }
                                    if(flag )
                                        build_vertex.push_back(j);
                                }
                                if(build_vertex.size() > 1) {
                                    int id1 =furthest_vertex_in_vector(0,build_vertex);
                                    int id2 =furthest_vertex_in_vector(id1,build_vertex);
                                    final_res.emplace_back(build_vertex[id1]);
                                    final_res.emplace_back(build_vertex[id2]);
                                }
                            }


                            int cid = vid;
                            delete_same_vertex(final_res);

                            if(final_res.size() < 2 )continue;
                           // if(final_res.size() < 6 )continue;


                            MeshKernel::iGameVertex final_center = sort_by_polar_order(final_res,origin_direct);
                            fprintf(file10, "v %lf %lf %lf\n", final_center.x(), final_center.y(), final_center.z());
                            fprintf(file7, "v %lf %lf %lf\n", final_center.x(), final_center.y(), final_center.z());
                            vid++;
                            for(int i=0;i<final_res.size();i++) {

                                MeshKernel::iGameVertex v1 = final_res[i];
                                MeshKernel::iGameVertex v2 = final_res[(i+1)%final_res.size()];
                                MeshKernel::iGameVertex now_direct = (v1 - final_center) % (v2 - final_center);
                                if(now_direct * origin_direct < 0)
                                    swap(v1,v2);
                                fprintf(file10, "v %lf %lf %lf \n", v1.x(), v1.y(), v1.z());
                                fprintf(file10, "v %lf %lf %lf \n", v2.x(),v2.y(), v2.z());
                                fprintf(file6, "v %lf %lf %lf \n", v1.x(), v1.y(), v1.z());
                                fprintf(file6, "v %lf %lf %lf \n", v2.x(),v2.y(), v2.z());
                                fprintf(file10, "f %d %d %d \n", cid, vid, vid + 1);
                                vid += 2;
                            }

                            static int xxx = 2;
                            //printf("*************\n");
                            fprintf(file8, "v %lf %lf %lf\n", final_center.x(), final_center.y(), final_center.z());
                           // cout << final_res.size() << endl;
                            for(int i=0;i<final_res.size();i++){
                                fprintf(file8,"v %lf %lf %lf\n",final_res[i].x(),final_res[i].y(),final_res[i].z());
                                //printf("v %lf %lf %lf\n",final_res[i].x(),final_res[i].y(),final_res[i].z());
                            }
                            for(int i=0;i<final_res.size();i++){
                                fprintf(file8,"l %d %d\n",i-1 + xxx,i+xxx);
                               // printf("l %d %d\n",i-1 + xxx,i+xxx);
                            }
                            xxx+=final_res.size();
                           // return 0;


//                            MeshKernel::iGameVertex final_center = sort_by_polar_order(final_res);
//                            int cid = vid;
//                            fprintf(file10, "v %lf %lf %lf\n", final_center.x(), final_center.y(), final_center.z());
//                            vid++;
//                            for(int i=0;i<final_res.size();i++) {
//                                fprintf(file10, "v %lf %lf %lf \n", final_res[i].x(), final_res[i].y(), final_res[i].z());
//                            }
//                            for(int i=0;i<final_res.size();i++) {
//                                MeshKernel::iGameVertex v1 = final_res[i];
//                                MeshKernel::iGameVertex v2 = final_res[(i+1) % final_res.size()];
//                                MeshKernel::iGameVertex now_direct = (v1 - final_center) % (v2 - final_center);
//                                if(now_direct * origin_direct < 0)
//                                    swap(v1,v2);
//                                //fprintf(file10,"f %d %d %d",cid,vid+i,vid+((i+1)%final_res.size()));
//                            }
//                            vid+=final_res.size();


                            /*for(auto p : final_res) {
                                MeshKernel::iGameVertex now_direct = (p.first - final_center) % (p.second - final_center);
                                if(now_direct * origin_direct < 0)
                                    swap(p.first,p.second);
                                fprintf(file10, "v %lf %lf %lf \n", p.first.x(), p.first.y(), p.first.z());
                                fprintf(file10, "v %lf %lf %lf \n", p.second.x(), p.second.y(), p.second.z());
                                fprintf(file10, "f %d %d %d \n", cid, vid, vid + 1);
                                vid += 2;
                            }*/
                        }



                    }



//
//                    for(int i=0;i<possible_face_list.size();i++){
//                        if(local_face_state[i] != 3){
//                            for(int j=i+1;j<possible_face_list.size();j++){
//                                if(local_face_state[j] != 3){// xiangjiao
//                                    if(dsu.find_root(i) == dsu.find_root(j))continue;
//                                    vector<K2::Point_3> ivs;
//                                    vector<K2::Point_3> jvs;
//                                    for(int id = 0;id <3;id++){
//                                        ivs.push_back(iGameVertex_to_Point_K2(field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(id)]));
//                                        jvs.push_back(iGameVertex_to_Point_K2(field_move_vertex[mesh->fast_iGameFace[possible_face_list[j]].vh(id)]));
//                                    }
//                                    K2::Triangle_3 i_tri(ivs[0],ivs[1],ivs[2]);
//                                    K2::Triangle_3 j_tri(jvs[0],jvs[1],jvs[2]);
//                                    CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
//                                            ij_res = intersection(i_tri, j_tri);
//                                    bool link_flag = false;
//                                    if (ij_res) {
//                                        //  cout << "GGGGGGG2 " << endl;
//
//                                        if (const K2::Point_3 * p2 = boost::get<K2::Point_3>(&*ij_res)) {
//                                            if (const K2::Segment_3* s = boost::get<K2::Segment_3>(&*ij_res)) {
//                                                link_flag = true;
//                                            }else if (const K2:: Triangle_3* ti = boost::get< K2::Triangle_3>(&*ij_res)) {
//                                                link_flag = true;
//                                            }else if(const std::vector <K2::Point_3>* vs = boost::get<std::vector <K2::Point_3>>(&*ij_res)){
//                                                link_flag = true;
//                                            }
//                                            //delete p2;
//                                        }
//                                    }
//                                    if(link_flag) {
//                                        local_face_state[i] = 2;
//                                        local_face_state[j] = 2;
//                                    }
//                                }
//                            }
//                        }
//                    }
//
//
//                    for(int i=0;i<possible_face_list.size();i++){
//                        if(local_face_state[i] == 1){
//                            if(built_face[dsu.find_root(i)] == 0){
//                                built_face[dsu.find_root(i)] = 1;
//                                MeshKernel::iGameVertex center;
//                                auto res = get_grid_intersect_plane(
//                                        field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(0)],
//                                        field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(1)],
//                                        field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(2)],
//                                        tinygrid_st,
//                                        tinygrid_end,
//                                        center
//                                ); //此处是1类面的生成 !! ;;;;;;
//
//                                int cid = vid;
//                                fprintf(file10,"v %lf %lf %lf\n",center.x(),center.y(),center.z());
//                                vid++;
////                                for(auto p : res) {
////                                    fprintf(file10, "v %lf %lf %lf \n", p.first.x(), p.first.y(), p.first.z());
////                                    fprintf(file10, "v %lf %lf %lf \n", p.second.x(), p.second.y(), p.second.z());
////                                    fprintf(file10, "f %d %d %d \n", cid, vid, vid + 1);
////                                    vid += 2;
////                                }
//
//                            }
//                        }
//                    }
//

                    continue;

                    // 判断二类情况;




//                    for(int i= 0 ;i < sharp_point_list.size();i++){
//                        if(vertex_in_tiny_grid(tinygrid_st,tinygrid_end,sharp_point_list_vertex[i])){
//                            type3_id = i;
//                        }
//                    }
//                    vector<MeshKernel::iGameFaceHandle > this_grid_faces;
//                    for(auto i : possible_face_list) {
//                        vector<MeshKernel::iGameVertex > vlist;
//                        vlist.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)]);
//                        vlist.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)]);
//                        vlist.push_back(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)]);
//                        if(face_through_grid(tinygrid_st,tinygrid_end,vlist)) {
//                            this_grid_faces.push_back(i);
//                        }
//                    }

//TODO  这里没考虑同一bingchaji面片;


                }
            }
        }



//        std::unordered_map<grid,vector<int>,grid_hash, grid_equal >tiny_grid_map;
//        std::function<bool(grid)> is_in_bound = [&](grid g){
//            if(g.x <0 || g.y < 0 || g.z <0 || g.x >= num || g.y >= num || g.z >= num)
//                return false;
//            return true;
//        };
//        if(num != 1 )continue;
//        // FIXME : BFS STARTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT!!!!!!!!;
//        stringstream  ss;
//
//        for(int i=0;i<possible_face_list.size();i++){
//            std::vector<MeshKernel::iGameVertex>moved_faces_vertex;
//            moved_faces_vertex.push_back(field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(0)]);
//            moved_faces_vertex.push_back(field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(1)]);
//            moved_faces_vertex.push_back(field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(2)]);
//
//            MeshKernel::iGameVertex bfs2_start = (field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(0)]
//                    +  field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(1)]
//                    +  field_move_vertex[mesh->fast_iGameFace[possible_face_list[i]].vh(2)])/3;
//
//            queue<grid>q;
//            std::unordered_set<grid,grid_hash,grid_equal>is_visit;
//            grid start_grid(int((bfs2_start.x() - small.x())/tiny_grid_len+myeps),
//                            int((bfs2_start.y() - small.y())/tiny_grid_len+myeps),
//                            int((bfs2_start.z() - small.z())/tiny_grid_len+myeps)); //st 越界问题;;;;;;;;;
//            q.push(start_grid);
//            ss << "now turn " << i << endl;
//            ss <<"st? " <<start_grid.x <<" "<< start_grid.y <<" "<< start_grid.z <<" qsize"<<q.size() << endl;
//            ss << " stv ? : " << bfs2_start.x() <<" "<< bfs2_start.y() <<" "<< bfs2_start.z() << endl;
//            ss << " small ? : " << small.x() <<" "<< small.y() <<" "<< small.z() << endl;
//            ss << " big ? : " << big.x() <<" "<< big.y() <<" "<< big.z() << endl;
//            while(!q.empty()){
//                grid g = q.front();
//                q.pop();
//                ss << "inq "<<  g.x <<" "<< g.y <<" "<< g.z <<" "<<i<< endl;
//                if(is_visit.count(g))continue;
//                ss << "not vis "<<  g.x <<" "<< g.y <<" "<< g.z << " "<< i << endl;
//                is_visit.insert(g);
//                MeshKernel::iGameVertex st_igv(g.x * tiny_grid_len + small.x(),
//                                               g.y * tiny_grid_len + small.y(),
//                                               g.z * tiny_grid_len + small.z());
//                MeshKernel::iGameVertex end_igv((g.x+1) * tiny_grid_len + small.x(),
//                                                (g.y+1) * tiny_grid_len + small.y(),
//                                                (g.z+1) * tiny_grid_len + small.z());
//                ss <<tt<<" g" << g.x <<" "<< g.y <<" "<< g.z<< " to if : " <<i << endl;
//                ss << "st v "<< st_igv.x() <<" "<< st_igv.y() <<" "<< st_igv.z() << endl;
//                ss << "en v "<< end_igv.x() <<" "<< end_igv.y() <<" "<< end_igv.z() << endl;
//
//                ss << "0v "<< moved_faces_vertex[0].x() <<" "<< moved_faces_vertex[0].y() <<" "<< moved_faces_vertex[0].z() << endl;
//                ss << "1v "<< moved_faces_vertex[1].x() <<" "<< moved_faces_vertex[1].y() <<" "<< moved_faces_vertex[1].z() << endl;
//                ss << "2v "<< moved_faces_vertex[2].x() <<" "<< moved_faces_vertex[2].y() <<" "<< moved_faces_vertex[2].z() << endl;
//                ss << "ving0 : "<<vertex_in_tiny_grid(st_igv,end_igv,moved_faces_vertex[0]) << endl;
//                ss << "ving1 : "<<vertex_in_tiny_grid(st_igv,end_igv,moved_faces_vertex[1]) << endl;
//                ss << "ving2 : "<<vertex_in_tiny_grid(st_igv,end_igv,moved_faces_vertex[2]) << endl;
//                ss << "ving res :" << face_through_grid(st_igv,end_igv,moved_faces_vertex) << endl;
//                for(int kk = 0;kk<8;kk++)
//                {
//                    auto v = getTinyGridVertex(st_igv,end_igv,kk);
//                    ss << "v " << v.x()<<" "<< v.y()<<" "<< v.z() << endl;
//                }
//
//
//
//                if(face_through_grid(st_igv,end_igv,moved_faces_vertex)){
//                    ss <<tt <<" g"<< g.x <<" "<< g.y <<" "<< g.z<< " in if : " <<i << endl;
//                    if(is_in_bound(g))
//                        tiny_grid_map[g].push_back(i);
//                    for(auto d : bfs_dir){
//                        auto nei_grid = grid(g.x + d[0],g.y + d[1],g.z + d[2]);
//                        if(  !is_visit.count(nei_grid)){
//                            q.push(nei_grid);
//                        }
//                    }
//                }
//            }
//        }
//
//        for(auto i : sharp_point_list){
//            double x = CGAL::to_double(i.p.x());
//            double y = CGAL::to_double(i.p.y());
//            double z = CGAL::to_double(i.p.z());
//            int num_x = int((x-small.x())/tiny_grid_len + myeps);
//            int num_y = int((y-small.y())/tiny_grid_len + myeps);
//            int num_z = int((z-small.z())/tiny_grid_len + myeps);
//            if(tiny_grid_map[grid(num_x,num_y,num_z)].size() <3 && num == 1){
//                string tmp;
//                while(getline(ss,tmp))
//                    cout << tmp << endl;
//
//                cout << "???"<<num << " tt "<< tt <<endl;
//                cout << num_x <<" "<< num_y <<" "<< num_z << endl;
//                cout <<" ERRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRROR"<<endl;
//                auto it = tiny_grid_map.find(grid(num_x,num_y,num_z));
//                for(auto  bfs_face_ph : it->second){
//                    cout << "bfs_face_ph "<< bfs_face_ph << endl;
//                }
//                cout << CGAL::to_double(i.p.x()) <<" "<<  CGAL::to_double(i.p.y())<<" "<<  CGAL::to_double(i.p.z()) << endl;
//                for(auto  cross_face_ph : i.source_face_local_id){
//                    cout << "cross_face_ph "<< cross_face_ph << endl;
//                    cout <<"vh0 " <<field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(0)].x() <<" "
//                    <<field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(0)].y() <<" "
//                    <<field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(0)].z() << endl;
//                    cout <<"vh1 " <<field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(1)].x() <<" "
//                         <<field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(1)].y() <<" "
//                         <<field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(1)].z() << endl;
//                    cout <<"vh2 " <<field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(2)].x() <<" "
//                         <<field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(2)].y() <<" "
//                         <<field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(2)].z() << endl;
//                    K2::Triangle_3 trii(K2::Point_3(field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(0)].x(),
//                                          field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(0)].y(),
//                                          field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(0)].z()),
//                                        K2::Point_3(field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(1)].x(),
//                                          field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(1)].y(),
//                                          field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(1)].z()),
//                                        K2::Point_3(field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(2)].x(),
//                                          field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(2)].y(),
//                                          field_move_vertex[mesh->fast_iGameFace[possible_face_list[cross_face_ph]].vh(2)].z()));
//                    cout <<"squared_distance" <<sqrt(CGAL::to_double(squared_distance(i.p,trii))) << endl;
//                }
//                cout <<"***************************************************"<<endl;
//
//            }
//        }
//



        /*for(auto i : tiny_grid_map){

        }*/
        //sum_grid += tiny_grid_map.size();

        // 补充二轮 bfs 法;
        // 然后写完 1 2 3 类 ;
        //cout << possible_face_list.size() << endl;




        /*
         *  三类格子
         *  1356 行的地方先确定相交确定为两个
         *  两条然后相交线延长出去
         *  做格子外壳与三角面片求交
         *  然后连线
         *  正负方向通过在tet内外搞
         *  能写！
         *  二类格子
         *
         */







//        for(int i=0;i<sharp_point_list.size();i++){
//            for(int k1 = 0;k1<sharp_point_list[i].source_face_local_id.size();k1++)
//            {
//                for(int k2 = 0;k2<sharp_point_list[i].source_face_local_id.size();k2++){
//                    if(k1!=k2){
//
//                    }
//                }
//            }
//        }

//        for(int x=0;x<num;x++)
//            for(int y=0;y<num;y++)
//                for(int z=0;z<num;z++){
//                    MeshKernel::iGameVertex start_pos = small + MeshKernel::iGameVertex(x*tiny_grid_len,
//                                                                 y*tiny_grid_len,
//                                                                 z*tiny_grid_len);
//                    MeshKernel::iGameVertex end_pos = small + MeshKernel::iGameVertex(x * (tiny_grid_len+1),
//                                                                                     y * (tiny_grid_len+1),
//                                                                                     z * (tiny_grid_len+1));
//
//                    int tiny_grid_type = 0;
//                    for(int i = 0 ;i < sharp_point_list.size();i++){
//                        MeshKernel::iGameVertex sharp_vertex(CGAL::to_double(sharp_point_list[i].p.x()),
//                                                        CGAL::to_double(sharp_point_list[i].p.y()),
//                                                        CGAL::to_double(sharp_point_list[i].p.z()));
//                        if(vertex_in_tiny_grid(start_pos,end_pos,sharp_vertex)){
//                            tiny_grid_type = 3;
//                            for(int j=0;j<sharp_point_list[i].source_face_local_id.size();j++) {
//
//                            }
//                            break;
//                        }
//                    }
//                    if(tiny_grid_type == 3)
//                        continue;
//
//                    // check is type 3 !!!!!!
//
//
//
//
//                    /*
//                     *         for(MeshKernel::iGameFaceHandle i : face_list){
//            vector<MeshKernel::iGameVertex> v;
//            for(int j=0;j<3;j++)
//                v.push_back(field_move_vertex[mesh->fast_iGameFace[i].vh(j)]);
//            if(face_through_grid(small,big,v))
//                possible_face_list.push_back(i);
//        }
//                     */
//
//                }


//        if(num >60){
//           // cout <<"sharp_point_list.size()" <<sharp_point_list.size() << "\n";
//            for(auto it : sharp_point_list){
//                fprintf(file10,"v %lf %lf %lf\n",CGAL::to_double(it.p.x()),CGAL::to_double(it.p.y()),CGAL::to_double(it.p.z()));
//            }
//        }

        //vertex_in_grid


//        set<int>sse;
//        for(int i=0;i<face_list.size();i++){
//            sse.insert(dsu.find_root(i));
//        }

        //cout << "sse.size(): " << sse.size() << endl;
    }
    cout << sum_grid << endl;
    // 二类面判断交点在什么方向;














    return 0;


    vector<vector<int> >tet_face_order{{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
    // 用cgal



//    for(auto each_dir : dir){
//        ThickenCubeEdge edge_info;
//        edge_info.from = each_grid->first;
//        edge_info.to = grid(edge_info.from.x+each_dir[0],edge_info.from.y+each_dir[1],edge_info.from.z+each_dir[2]);
//        if(!frame_grid_mp.count(edge_info.to))continue;
//        edge_info.type = each_grid->second.grid_type;
//        vector<MeshKernel::iGameFaceHandle> face_list;
//        for(auto i : each_grid->second.face_list) {
//            if ()// zaili mian 不对 这个地方 123
//            {
//            }
//        }
//        if(edge_info.type == 0){
//
//        }
//        else if(edge_info.type ==1){
//
//        }
//
//    }



    //thicken_cube_edge_list.resize()

    //DSU
    //TODO 可改多线程 555




























    return 0;



    function<MeshKernel::iGameVertex(grid)>qef_grid_pos_all_consider = [&](grid g) {
        vector<double>b;
        vector<double>sx,sy,sz,nx,ny,nz;
        MeshKernel::iGameVertex center_of_weight(0,0,0);
        //cout<< "************************************************** "<< endl;
        //cout<< "grid g : "<<g.x<<" "<<g.y<<" "<< g.z << endl;
        // todo: 改这里，改成所有面能力的叠加 20220623 搞这里是时刻
        int center_cnt = 0;
        for (int ev1 = 0; ev1 < GridEdge.size(); ev1++) {
            for (auto ev2: GridEdge[ev1]) {
                if (ev1 < ev2) {
                    grid v1 = getGridFrameVertex(g,ev1);
                    grid v2 = getGridFrameVertex(g,ev2);
                    auto itv1 = frame_grid_mp.find(v1);
                    auto itv2 = frame_grid_mp.find(v2);
                    if(itv1 == frame_grid_mp.end() || itv2 == frame_grid_mp.end()){
                        continue;
                    }
                    if(itv1->second.grid_type > itv2->second.grid_type) {
                        swap(itv1, itv2);
                        swap(v1,v2);
                    }
                    // 这里确实有问题！ 4454545 先不管了 小丑 55555

                    if( (itv1->second.grid_type ==1||itv1->second.grid_type ==0) && itv2->second.grid_type ==2) { //判断1 2 的face list
                        //  cout << "itv1 : "<<itv1->first.x<<" "<<  itv1->first.y<<" "<<  itv1->first.z<<" "<<itv1->second.grid_type<< endl;
                        // cout << "itv2 : "<<itv2->first.x<<" "<<  itv2->first.y<<" "<<  itv2->first.z<<" "<< itv2->second.grid_type<<endl;

                        //cout <<"itv1fl size:"<< itv1->second.face_list.size() <<endl;


                        MeshKernel::iGameVertex v0_vertex = getGridVertex(v1, 0);
                        MeshKernel::iGameVertex v1_vertex = getGridVertex(v2, 0);
                        MeshKernel::iGameVertex pi = v0_vertex;
                        if(itv1->second.grid_type ==0)
                            pi = (v0_vertex + v1_vertex)/2;
                        MeshKernel::iGameVertex ni;
                        int tmp;
                        for (MeshKernel::iGameFaceHandle fht: itv1->second.face_list) {

//                            ni += (((mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(1)]
//                            - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)]) %
//                                    (mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(2)]
//                                     - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)]))*-1).normalize();

                            ni += (pi - cgal_aabbtree_query({mesh->fast_iGameFace.at(fht)}, pi, mesh, tmp,
                                                            nullptr)).normalize();
                        }
                        ni /= itv1->second.face_list.size();

                        double l = 0, r = 1;
                        double limit = 0.01;
                        //    cout << "st erfen"<< endl;
                        while (l + limit < r) {
                            double mid = (l + r) / 2;
                            MeshKernel::iGameVertex cut_point = (v0_vertex + (v1_vertex - v0_vertex) * mid);
                            bool flag = false;
                            std::shared_ptr<vector<MeshKernel::iGameFaceHandle> > field_face_list =
                                    std::make_shared<vector<MeshKernel::iGameFaceHandle> >();
                            int state = vertex_state(itv1->second.face_list,cut_point,field_face_list);
                            if (state == 2)
                                r = mid;
                            else
                                l = mid;

                            if (state == 1) {
                                pi = cut_point;
                                ni = (pi - cgal_aabbtree_query({mesh->fast_iGameFace.at(*field_face_list->begin())}, pi, mesh,tmp,
                                                               nullptr)).normalize();

//                                ni = (((mesh->fast_iGameVertex[mesh->fast_iGameFace.at(*field_face_list->begin()).vh(1)]
//                                         - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(*field_face_list->begin()).vh(0)]) %
//                                        (mesh->fast_iGameVertex[mesh->fast_iGameFace.at(*field_face_list->begin()).vh(2)]
//                                         - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(*field_face_list->begin()).vh(0)]))*-1).normalize();
                                flag = true;
                            }

//                            cout <<CGALPlaneDist(mesh->faces(fh),v0_vertex)<<" "<< CGALPlaneDist(mesh->faces(fh),v1_vertex)<<
//                            " "<<  CGALPlaneDist(mesh->faces(fh),pi)<<endl;
                            //   cout <<l<< "  en erfen"<< endl;
                            //   cout << "pi.:" << pi.x() << " " << pi.y() << " " << pi.z() << endl;
                            sx.push_back(pi.x());
                            sy.push_back(pi.y());
                            sz.push_back(pi.z());
                            //   cout << "ni.:" << ni.x() << " " << ni.y() << " " << ni.z() << endl;

                            nx.push_back(ni.x());
                            ny.push_back(ni.y());
                            nz.push_back(ni.z());
                            center_of_weight = center_of_weight + pi;

                            // cout <<"pi "<<pi.x() <<" "<< pi.y()<<" "<<pi.z() << endl;
//                            fprintf(file1, "v %lf %lf %lf\n",pi.x(),pi.y(),pi.z());

                            center_cnt++;
                            //    cout <<"sxsize:" <<sx.size() << endl;
                        }
                    }

                }
            }
        }
        // cout <<"center_cnt:" <<center_cnt << endl;

        vector<double>res(3);
        QEF_calculate({getGridVertex(g,0).x(),getGridVertex(g,0).y(),getGridVertex(g,0).z()},
                      {getGridVertex(g,7).x(),getGridVertex(g,7).y(),getGridVertex(g,7).z()},
                      sx.size(),sx,sy,sz,nx,ny,nz,res);
        //  cout << "res: "<<res[0]<<" "<<res[1]<<" "<<res[2]<<endl;
        // cout << "qefinner? "<< vertex_in_grid(g,MeshKernel::iGameVertex(res[0],res[1],res[2])) << endl;
        // cout << center_cnt << endl;
        //cout <<"jiuji bug :"<< vertex_in_grid(g,center_of_weight/center_cnt)<<endl;
        //  cout << (center_of_weight/center_cnt).x() <<" "<<  (center_of_weight/center_cnt).y()<<" "<<
        //                                                                                      (center_of_weight/center_cnt).z()<<endl;
        //  cout << getGirdVertex(g,0).x()<<" "<<getGirdVertex(g,0).y()<<" "<<getGirdVertex(g,0).z()<<endl;
        //   cout << getGirdVertex(g,7).x()<<" "<<getGirdVertex(g,7).y()<<" "<<getGirdVertex(g,7).z()<<endl;
        // return center_of_weight/center_cnt;


        //TODO : debug code need del

        if(vertex_in_grid(g,MeshKernel::iGameVertex(res[0],res[1],res[2]))){
            return MeshKernel::iGameVertex(res[0],res[1],res[2]);
        }
        else{
            // cout<< "out"<<endl;
            return center_of_weight/center_cnt;
        }
    };




    function<MeshKernel::iGameVertex(grid)>qef_grid_pos = [&](grid g) {

        vector<double>b;
        vector<double>sx,sy,sz,nx,ny,nz;
        MeshKernel::iGameVertex center_of_weight(0,0,0);
        //cout<< "************************************************** "<< endl;
        //cout<< "grid g : "<<g.x<<" "<<g.y<<" "<< g.z << endl;
        // todo: 改这里，改成所有面能力的叠加 20220623 搞这里是时刻
        int center_cnt = 0;
        for (int ev1 = 0; ev1 < GridEdge.size(); ev1++) {
            for (auto ev2: GridEdge[ev1]) {
                if (ev1 < ev2) {
                    grid v1 = getGridFrameVertex(g,ev1);
                    grid v2 = getGridFrameVertex(g,ev2);
                    auto itv1 = frame_grid_mp.find(v1);
                    auto itv2 = frame_grid_mp.find(v2);
                    if(itv1 == frame_grid_mp.end() || itv2 == frame_grid_mp.end()){
                        continue;
                    }
                    if(itv1->second.grid_type > itv2->second.grid_type) {
                        swap(itv1, itv2);
                        swap(v1,v2);
                    }


                    if( (itv1->second.grid_type ==1||itv1->second.grid_type ==0) && itv2->second.grid_type ==2) { //判断1 2 的face list
                        MeshKernel::iGameVertex v0_vertex = getGridVertex(v1, 0);
                        MeshKernel::iGameVertex v1_vertex = getGridVertex(v2, 0);
                        MeshKernel::iGameVertex edge_pi ;
                        MeshKernel::iGameVertex edge_ni ;
                        int edge_cnt = 0;
                        for (MeshKernel::iGameFaceHandle fht: itv1->second.face_list) {
                            int state_v0 = vertex_state({fht}, v0_vertex);
                            int state_v1 = vertex_state({fht}, v1_vertex);
                            if (state_v0 != state_v1 && state_v1 == 2) {
                                double l = 0, r = 1;
                                double limit = 0.01;
                                bool flag = false;
                                MeshKernel::iGameVertex pi,ni;
                                while (l + limit < r) {
                                    double mid = (l + r) / 2;
                                    MeshKernel::iGameVertex cut_point = (v0_vertex + (v1_vertex - v0_vertex) * mid);
                                    int state_mid = vertex_state({fht}, cut_point);
                                    if (state_mid == 2)
                                        r = mid;
                                    else
                                        l = mid;
                                    if (state_mid == 1){
                                        pi = (v0_vertex + (v1_vertex - v0_vertex) * l);
                                        ni = (((mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(1)]
                                                  - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)]) %
                                                 (mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(2)]
                                                  - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)])) *
                                                -1).normalize();
                                        flag = true;
                                    }
                                }
                                if(flag) {
                                    edge_pi+=pi;
                                    edge_ni+=ni;
                                    edge_cnt++;
                                }
                                else{
                                    edge_ni+= (((mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(1)]
                                            - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)]) %
                                           (mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(2)]
                                            - mesh->fast_iGameVertex[mesh->fast_iGameFace.at(fht).vh(0)])) *
                                          -1).normalize();
                                    edge_pi+=v1_vertex;
                                    edge_cnt++;
                                }
                            }
                        }
                        if(edge_cnt>0){
                            edge_pi/=edge_cnt;
                            edge_ni/=edge_cnt;
                            sx.push_back(edge_pi.x());
                            sy.push_back(edge_pi.y());
                            sz.push_back(edge_pi.z());
                            nx.push_back(edge_ni.x());
                            ny.push_back(edge_ni.y());
                            nz.push_back(edge_ni.z());
                            center_of_weight+=edge_pi;
                            center_cnt++;
                        }

                    }
                }//end  if (ev1 < ev2)
            }
        }
//        if(sx.size()==0 || center_cnt ==0){
//            return (MeshKernel::iGameVertex{getGridVertex(g,0).x(),getGridVertex(g,0).y(),getGridVertex(g,0).z()}+
//                    MeshKernel::iGameVertex{getGridVertex(g,7).x(),getGridVertex(g,7).y(),getGridVertex(g,7).z()})/2;
//
//        }
        // cout <<"center_cnt:" <<center_cnt << endl;

        vector<double>res(3);
        QEF_calculate({getGridVertex(g,0).x(),getGridVertex(g,0).y(),getGridVertex(g,0).z()},
                      {getGridVertex(g,7).x(),getGridVertex(g,7).y(),getGridVertex(g,7).z()},
                      sx.size(),sx,sy,sz,nx,ny,nz,res);



        if(vertex_in_grid(g,MeshKernel::iGameVertex(res[0],res[1],res[2]))){
            return MeshKernel::iGameVertex(res[0],res[1],res[2]);
        }
        else if(center_cnt>0){
            // cout<< "out"<<endl;
            return center_of_weight/center_cnt;
        }
        else{
            return qef_grid_pos_all_consider(g);
            /*cout <<"error"<<endl;
            return MeshKernel::iGameVertex{getGridVertex(g,0).x(),getGridVertex(g,0).y(),getGridVertex(g,0).z()}+
                   MeshKernel::iGameVertex{getGridVertex(g,7).x(),getGridVertex(g,7).y(),getGridVertex(g,7).z()};*/
        }
    };


    for(auto i : frame_grid_mp){
        if(i.second.grid_type==0){
            fprintf(file0, "v %lf %lf %lf\n", getGridVertex(i.first, 0).x(), getGridVertex(i.first, 0).y(),
                    getGridVertex(i.first, 0).z());
        }
    }

    for(auto i : frame_grid_mp){
        if(i.second.grid_type==1){
            fprintf(file1, "v %lf %lf %lf\n", getGridVertex(i.first, 0).x(), getGridVertex(i.first, 0).y(),
                    getGridVertex(i.first, 0).z());
        }
    }

    for(auto i : frame_grid_mp){
        if(i.second.grid_type==2){
            fprintf(file2, "v %lf %lf %lf\n", getGridVertex(i.first, 0).x(), getGridVertex(i.first, 0).y(),
                    getGridVertex(i.first, 0).z());
        }
    }

    unordered_map<grid,int ,grid_hash,grid_equal> build_mesh_vertex;
    vector<MeshKernel::iGameVertex>qef_result;
    std::vector<vector<int> > quad_face;

    string ss;
    int frame_grid_mp_cnt = 0;
    for(auto i : frame_grid_mp) {
        frame_grid_mp_cnt++;
        if(frame_grid_mp_cnt%1000==0)
            printf("%d/%d\n",frame_grid_mp_cnt,frame_grid_mp.size());
        auto v = getGridVertex(i.first, 0);
        if(i.second.grid_type==1 || i.second.grid_type==0) {
            auto tmp = get_neighbor(i.first);
            for(auto j : tmp){
                auto it = frame_grid_mp.find(j);
                if(it!= frame_grid_mp.end() && it->second.grid_type == 2){

                    auto v_gird = get_edge_neighbor(i.first,j);
                    for(auto vh : v_gird){
                        if(!build_mesh_vertex.count(vh)){
                            int id = (int)build_mesh_vertex.size()+1;
                            build_mesh_vertex.insert({vh ,id});
                        }
                    }
                    if (i.first < j) {
                        std::swap(v_gird[1], v_gird[3]);
                    }
                    quad_face.push_back({build_mesh_vertex[v_gird[0]],build_mesh_vertex[v_gird[1]] ,build_mesh_vertex[v_gird[2]],build_mesh_vertex[v_gird[3]]});
                }
            }
        }
    }
    qef_result.resize(build_mesh_vertex.size()+1);

    std::cout <<"start multi thread "<< endl;


    std::vector<std::shared_ptr<std::thread> > thread_pool(thread_num);
    for(int i=0;i<thread_num;i++){
        thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for(auto iter : build_mesh_vertex){
                if(iter.second%thread_num == now_id){
                    qef_result[iter.second] = qef_grid_pos(iter.first);
                    if(iter.second%500==0) // cut the output time usage
                        printf("succ %d/%d \n",iter.second, build_mesh_vertex.size());
                }
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        thread_pool[i]->join();

    for(int i=1 ;i<=build_mesh_vertex.size();i++){
        fprintf(file,"v %lf %lf %lf\n",qef_result[i].x(),qef_result[i].y(),qef_result[i].z());
    }



    for(auto i : quad_face){
        if((qef_result[i[0]] -
            qef_result[i[2]]).norm() <
           (qef_result[i[1]] -
            qef_result[i[3]]).norm()) {
            ss += "f";
            for (int k: {0, 1, 2})
                ss += (" " + to_string(i[k]));
            ss += "\n";
            ss += "f";
            for (int k: {2, 3, 0})
                ss += (" " + to_string(i[k]));
            ss += "\n";
        }
        else{
            ss += "f";
            for (int k: {1, 2, 3})
                ss += (" " + to_string(i[k]));
            ss += "\n";
            ss += "f";
            for (int k: {3, 0, 1})
                ss += (" " + to_string(i[k]));
            ss += "\n";
        }

    }


    fprintf(file,"%s",ss.c_str());



    fclose(file0);
    fclose(file1);
    fclose(file2);

    unordered_map<grid,int,grid_hash,grid_equal> build_debug_grid;
    for(auto i : frame_grid_mp){
        fprintf(file10, "v %lf %lf %lf\n", getGridVertex(i.first, 0).x(), getGridVertex(i.first, 0).y(),
                getGridVertex(i.first, 0).z());
        int id = build_debug_grid.size()+1;
        build_debug_grid[i.first] = id;
    }
    for(auto i : frame_grid_mp){
        auto nei = get_neighbor(i.first);
        for(auto j : nei){
            if(i.first < j && frame_grid_mp.count(j)) {
                if(build_debug_grid.count(j)) {
                    fprintf(file10, "l %d %d\n", build_debug_grid[i.first], build_debug_grid[j]);
                }
            }
        }
    }
    fclose(file10);

    return 0;
}
// 1 2 3 4 5 6
// 5 1 2 3 7 8
//