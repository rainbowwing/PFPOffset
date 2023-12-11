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
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <vector>
#include <iostream>
#include <CGAL/Search_traits_3.h>
#include <atomic>
#include "obj_input.h"
#include "grid.h"
#include "global_face.h"
#include "geometric_data.h"
#include "cdt.h"
#include "unique_hash.h"
#include "CoverageField.h"
#include "osqp.h"
#include "dp.h"
#include "sort_by_polar_order.h"
#include "remeshing.h"
#include "single_coverage_ray_detect.h"
#include "flag_parser.h"
#include "merge_initial.h"
#include "do_winding_number.h"
#include "disassemble_circle.h"
#include <omp.h>

int check_resolution = 3;
using namespace std;
int main(int argc, char* argv[]) {

    //test_wining_num();
//    exit(0);
    std::streambuf* original_cerr = std::cerr.rdbuf();

    // Redirect cerr to a file (or /dev/null)
    std::ofstream nullstream("/dev/null");
    std::cerr.rdbuf(nullstream.rdbuf());
    google::ParseCommandLineFlags(&argc, &argv, true);
    flag_parser();

    cout <<"CGAL_RELEASE_DATE:" << CGAL_RELEASE_DATE << endl;
    mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile(input_filename)); grid_len = 0.1;
    //update_model();
    start_x = 0;
    start_y = 0;
    start_z = 0;
    set_start();
    input_filename = input_filename.substr(0,input_filename.size()-4);
    input_filename+= string("_") + (running_mode==1 ? "offset_outward" : "offset_inward");
    FILE *file11 = fopen( (input_filename + "_grid.obj").c_str(), "w");
    FILE *file6 = fopen( (input_filename + "_tmp.obj").c_str(), "w");
    FILE *file7 = fopen( (input_filename + "_no_manifold.obj").c_str(), "w");
    FILE *file8 = fopen( (input_filename + "_result.obj").c_str(), "w");
    FILE *file13 = fopen( (input_filename + "_moveedge.obj").c_str(), "w");

    FILE *file14 = fopen( (input_filename + "_coveragefield.obj").c_str(), "w");
    FILE *file24 = fopen( (input_filename + "_coveragefield24.obj").c_str(), "w");
    FILE *file44 = fopen( (input_filename + "_check_resolution_delete.obj").c_str(), "w");
    FILE *file54 = fopen( (input_filename + "_check_resolution_reserver.obj").c_str(), "w");
    FILE *file34 = fopen( (input_filename + "_cutting_segment.obj").c_str(), "w");
    FILE *file676 = fopen( (input_filename + "_fffhole.obj").c_str(), "w");
    mesh->initBBox();
   // mesh->build_fast();
   // mesh->build_fast();
    auto start_clock = std::chrono::high_resolution_clock::now();
    cout <<"mesh->build_fast() succ" << endl;
    //double default_move_dist = 0.05;
    double x_len = (mesh->BBoxMax - mesh->BBoxMin).x();
    double y_len = (mesh->BBoxMax - mesh->BBoxMin).y();
    double z_len = (mesh->BBoxMax - mesh->BBoxMin).z();
    double min_bbox_len = min(min(x_len,y_len),z_len);
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

            sum += min(minx,maxx);

        }
        avg_edge_limit = sum/mesh->FaceSize();
        if(default_move <= 0) {
            default_move = min(min(x_len, y_len), z_len)*1e-3;
            cout <<"default_move_dist:" <<default_move << endl;
            //exit(0);
        }
        double min_len = min(min(x_len,y_len),z_len);
        //double min_len = sqrt(x_len*x_len+y_len*y_len+z_len*z_len)*1e-3;
        double rate = x_len/min_len*y_len/min_len*z_len/min_len;

        cout <<"??"<<(4*thread_num) <<" "<<rate << endl; // minlen rate 格了 4* thread_num 格子
        double need_div = max(cbrt((4*thread_num)/rate),1.0);
        grid_len = min_len / need_div;//(mesh->BBoxMax - mesh->BBoxMin).norm()/ thread_num * 2;
        cout << "grid_len "<< grid_len<<endl;

    }
//    FILE * file50 = fopen((input_filename.substr(0,input_filename.size()-4) + "_offset.obj2").c_str(),"w");
//
//    fprintf(file50,"#download from quad mesh, add offset distance of every by PFPOffset.\n");
//    for (int i = 0; i < mesh->VertexSize(); i++) {
//        //cout << mesh->fast_iGameVertex[i].x() <<" "<< mesh->fast_iGameVertex[i].y() <<" "<<  mesh->fast_iGameVertex[i].z()<<endl;
//        fprintf(file50,"v %.7lf %.7lf %.7lf\n",mesh->fast_iGameVertex[i].x(),
//                mesh->fast_iGameVertex[i].y(),
//                mesh->fast_iGameVertex[i].z());
//    }
    min_near_limit = 1e100;
    for (int i = 0; i < mesh->FaceSize(); i++) {
        auto v0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(0)];
        auto v1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(1)];
        auto v2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(2)];
        auto tmp = (v0 + v1 + v2)/3;
        if(mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist <= 0){// deal Illegal input
            mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist = default_move;
        }
        min_near_limit = min(min_near_limit,mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist);

        cout <<"facei "<<  mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist<<" "<<mesh->FastNeighborFhOfFace_[i].size()<< endl;
        //myeps = min(myeps,mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist);
//        fprintf(file50,"f %d %d %d %.7lf\n",mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(0)+1,
//                mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(1)+1,
//                mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(2)+1,
//                ((tmp.x() - mesh->BBoxMin.x())/(mesh->BBoxMax.x() - mesh->BBoxMin.x())/2+0.5)*default_move
//                );

        // mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist =  default_move_dist;
        //FILE *file14 = fopen( (input_filename + "_coveragefield.obj").c_str(), "w");
    }
    min_near_limit /= 10;
  //  exit(0);


    field_move_vertex.resize(mesh->VertexSize());
    field_move_vertices.resize(mesh->VertexSize());
    origin_mesh_vertices.resize(mesh->VertexSize());

    field_move_face.resize(mesh->FaceSize());
    field_move_K2_triangle.resize(mesh->FaceSize());

    cout <<"st do_quadratic_error_metric" << endl;

    std::vector <std::shared_ptr<std::thread> > build_neighbor(thread_num);
    for(int i=0;i<thread_num;i++) {
        //for (int i = 50; i < 51; i++) {
        build_neighbor[i] = make_shared<std::thread>([&](int now_id) {
            std::list<K2::Triangle_3> tri_list;
            std::unordered_map<unsigned long long ,int> mp;
            for (int i = 0; i < mesh->FaceSize(); i++) {
                auto v0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(0)];
                auto v1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(1)];
                auto v2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(2)];
                K2::Triangle_3 tri(iGameVertex_to_Point_K2(v0),
                                   iGameVertex_to_Point_K2(v1),
                                   iGameVertex_to_Point_K2(v2)
                );
                mp[tri.id()] = i;
                tri_list.push_back(tri);
            }
            Tree aabb(tri_list.begin(),tri_list.end());
            for (int i = 0; i < mesh->VertexSize(); i++) {
                if (i % thread_num != now_id)continue;
                if (i % 20 == 0)
                    cout << "build near: " << i << "/"<<mesh->VertexSize()<< endl;
                CGAL::Epeck::FT r(min_bbox_len/1000);
                std::list<Primitive> intersected_primitives;
                std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
                K2::Point_3 this_vertex = iGameVertex_to_Point_K2(mesh->fast_iGameVertex[i]);
//                cout << coverage_field_list[i].bbox_min <<" "<< coverage_field_list[i].bbox_max<<endl;
//                CGAL::Epeck::FT xx = coverage_field_list[i].bbox_min.x();
//                cout << xx - CGAL::Epeck::FT(0.05) << endl;
                K2::Iso_cuboid_3 bbox2(
                        this_vertex.x() - r,
                        this_vertex.y() - r,
                        this_vertex.z() - r,
                        this_vertex.x() + r,
                        this_vertex.y() + r,
                        this_vertex.z() + r
                );
                aabb.all_intersections(bbox2,std::back_inserter(intersections));
                cout << "intersectionssize:"<< intersections.size() << endl;
                for(auto item : intersections) {
//                    cout <<"v " <<item.second->vertex(0) << endl;
//                    cout <<"v " <<item.second->vertex(1) << endl;
//                    cout <<"v " <<item.second->vertex(2) << endl;
                    auto iter = mp.find(item.second->id());
                    if(iter == mp.end()) exit(0);// 异常错误
                    mesh->FastNeighborFhOfVertex_[i].insert(MeshKernel::iGameFaceHandle(iter->second));
                }
            }
        }, i);
    }

    for(int i=0;i<thread_num;i++)
        build_neighbor[i]->join();
    for (int i = 0; i < mesh->VertexSize(); i++) {
        cout << i <<" "<< mesh->FastNeighborFhOfVertex_[i].size() << endl;
    }

   // exit(0);
    for(int i =0;i< mesh->FastNeighborFhOfVertex_.size();i++){
        for(auto j : mesh->FastNeighborFhOfVertex_[i]){
            std::set<int>se;
            se.insert(mesh->fast_iGameFace[j].vh(0));
            se.insert(mesh->fast_iGameFace[j].vh(1));
            se.insert(mesh->fast_iGameFace[j].vh(2));
            for(auto k: mesh->FastNeighborFhOfVertex_[i])
            {
                if(j==k)continue;
                int cnt = se.count(mesh->fast_iGameFace[k].vh(0)) +
                        se.count(mesh->fast_iGameFace[k].vh(1)) +
                        se.count(mesh->fast_iGameFace[k].vh(2));
                if(cnt == 2){
                    mesh->FastNeighborFhOfFace_[j].insert(k);
                }
            }
        }
    }


    merge_limit.resize(mesh->VertexSize());



    std::vector <std::shared_ptr<std::thread> > dp_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        dp_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for (int i = 0; i < mesh->VertexSize(); i++) {
            //for (int i = 524; i < 525; i++) {
                if (i % thread_num != now_id)continue;
                if (i % 20 == 0)
                    cout << "dp: " << i << "/"<<mesh->VertexSize()<< endl;
                vector<MeshKernel::iGameFaceHandle> neighbor_list;
                double min_move_dist = 1e30;
                for(auto j : mesh->FastNeighborFhOfVertex_[i]){
                    neighbor_list.push_back(j);
                }
                if(neighbor_list.empty()) {
                    merge_limit[i] = 0 ;
                    continue;
                };
                // cout <<i<<"//"<<mesh->VertexSize()<<" dp_result.size(): " <<"start" << endl;
                vector<MeshKernel::iGameVertex> dp_result = solve_by_dp(MeshKernel::iGameVertexHandle(i),neighbor_list);
                if(dp_result.empty()){
                    merge_limit[i] = 0 ;
                    continue;
                }
                //cout <<i<<"//"<<mesh->VertexSize()<<" dp_result.size(): " <<dp_result.size() << endl;
                for(auto t: dp_result){
                    field_move_vertices[i].emplace_back(t.x(),t.y(),t.z());
                    min_move_dist = min(min_move_dist,mesh->fast_iGameVertex[i].dist(t));
                }

                merge_limit[i] = min_move_dist/10;
            }
        }, i);
    }
    for(int i=0;i<thread_num;i++)
        dp_thread_pool[i]->join();
    int cnte = 1;
    for (int i = 0; i < mesh->VertexSize(); i++) {
        for(int j = 0;j<field_move_vertices[i].size();j++){
            fprintf(file13,"v %lf %lf %lf\n",CGAL::to_double(field_move_vertices[i][j].x()),
                    CGAL::to_double(field_move_vertices[i][j].y()),
                    CGAL::to_double(field_move_vertices[i][j].z())
                    );
            fprintf(file13,"v %lf %lf %lf\n",mesh->fast_iGameVertex[i].x(),
                    mesh->fast_iGameVertex[i].y(),
                    mesh->fast_iGameVertex[i].z()
            );
            fprintf(file13,"l %d %d\n",cnte,cnte+1);
            cnte+=2;
        }
        //fprintf(file13,"v");
    }
    for(int i=0;i<mesh->VertexSize();i++){
        origin_mesh_vertices[i] = K2::Point_3 (mesh->fast_iGameVertex[i].x(),
                                               mesh->fast_iGameVertex[i].y(),
                                               mesh->fast_iGameVertex[i].z()
                                               );
    }
   // exit(0);

    merge_initial();


    for(int i=0;i<mesh->FaceSize();i++){
        coverage_field_list.push_back(CoverageField(MeshKernel::iGameFaceHandle(i)));

    }
    std::list<K2::Triangle_3>origin_face_list;
    for(int i=0;i<mesh->FaceSize();i++){
//                K2::Point_3 v0(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].x(),
//                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].y(),
//                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)].z());
//
//                K2::Point_3 v1(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].x(),
//                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].y(),
//                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)].z());
//
//                K2::Point_3 v2(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].x(),
//                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].y(),
//                               mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)].z());
        K2::Point_3 v0 = origin_mesh_vertices[mesh->fast_iGameFace[i].vh(0)];
        K2::Point_3 v1 = origin_mesh_vertices[mesh->fast_iGameFace[i].vh(1)];
        K2::Point_3 v2 = origin_mesh_vertices[mesh->fast_iGameFace[i].vh(2)];
        origin_face_list.emplace_back(v0,v1,v2);
    }

    Tree origin_face_tree(origin_face_list.begin(),origin_face_list.end());



    int pre_cnt = 0;
    for(int i=0;i<mesh->FaceSize();i++){
//        if(coverage_field_list[i].bound_face_vertex_inexact.size() == 7 ) {
        for (int j = 0; j < coverage_field_list[i].bound_face_vertex_exact.size(); j++) {
            fprintf(file14,"v %lf %lf %lf\n",CGAL::to_double(coverage_field_list[i].bound_face_vertex_exact[j].x()),
                    CGAL::to_double(coverage_field_list[i].bound_face_vertex_exact[j].y()),
                    CGAL::to_double(coverage_field_list[i].bound_face_vertex_exact[j].z()));
        }
        for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++){
            fprintf(file14,"f %d %d %d\n",coverage_field_list[i].bound_face_id[j][0]+1+pre_cnt,
                    coverage_field_list[i].bound_face_id[j][1]+1+pre_cnt,
                    coverage_field_list[i].bound_face_id[j][2]+1+pre_cnt);
            auto v0 = coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]];
            auto v1 = coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]];
            auto v2 = coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]];
//            if(is_same_triangle(Point_K2_to_Point_K(v0),
//                             Point_K2_to_Point_K(v1),
//                             Point_K2_to_Point_K(v2),
//                             K::Point_3(-464.815369,-287.047913,-743.737305),
//                             K::Point_3(464.136993,-287.023560,-744.651611),
//                             K::Point_3(464.045746,-17.413059,-745.819275),5.0
//                             )){
//                cout << v0 << endl;
//                cout << v1 << endl;
//                cout << v2 << endl;
//                cout <<"gaga "<< i<<" "<< j << endl;
//            }
        }
        pre_cnt += coverage_field_list[i].bound_face_vertex_exact.size();
//            break;
//        }
    }
    fclose(file14);
    //exit(0);
    vector<set<int> > coverage_intersection(mesh->FaceSize());
    //X: 全局删法1
    if(0){
        std::vector <std::shared_ptr<std::thread> > find_near(thread_num);
        for(int i=0;i<thread_num;i++)  {
            find_near[i] = make_shared<std::thread>([&](int now_id) {
                list<K2::Triangle_3> aabb_tree_face_list;
                map<unsigned long long ,int>aabb_tree_mp;
                for(int i=0;i<mesh->FaceSize();i++){
                    for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++){
                        K2::Triangle_3 tri(coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
                                           coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
                                           coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]);
                        aabb_tree_face_list.push_back(tri);
                        aabb_tree_mp[tri.id()] = i;
                    }
                }
                Tree aabb_tree(aabb_tree_face_list.begin(),aabb_tree_face_list.end());
                for (int i = 0; i < mesh->FaceSize(); i++) {
                    if (i % thread_num != now_id)continue;
                    cout <<"find near:"<<i<<":" <<mesh->FaceSize() << endl;
                    std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
                    aabb_tree.all_intersections(coverage_field_list[i].iso_cuboid_3,std::back_inserter(intersections));
                    for(auto item : intersections) {
                        auto iter = aabb_tree_mp.find(item.second->id());
                        if(iter->second != i)
                            coverage_intersection[i].insert(iter->second);
                    }
                }
            },i);
        }
        for(int i=0;i<thread_num;i++)
            find_near[i]->join();
        std::vector <std::shared_ptr<std::thread> > near_delete(thread_num);
        for(int i=0;i<thread_num;i++) {
            near_delete[i] = make_shared<std::thread>([&](int now_id) {
                for (int i = 0; i < mesh->FaceSize(); i++) {
                    if (i % thread_num != now_id)continue;
                    if(!coverage_field_list[i].useful)continue;
                    cout <<"near: "<< i <<":"<< mesh->FaceSize()<<endl;
                    coverage_field_list[i].self_face_delete_flag = true;
                    vector<K2::Triangle_3>neighbor_face;
                    for(auto neighbor_id: coverage_intersection[i]){
                        if(!coverage_field_list[neighbor_id].useful)continue;

                        for (int k = 0; k < coverage_field_list[neighbor_id].bound_face_id.size(); k++) {
                            K2::Triangle_3 tri_this(coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][0]],
                                                    coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][1]],
                                                    coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][2]]
                            );

                            neighbor_face.push_back(tri_this);
                        }
                    }

                    K2::Triangle_3 tri_this(coverage_field_list[i].bound_face_vertex_exact[0],
                                            coverage_field_list[i].bound_face_vertex_exact[1],
                                            coverage_field_list[i].bound_face_vertex_exact[2]);
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
                    vector<vector<K2::Point_3> > res = CGAL_CDT_NEW({tri_this.vertex(0),
                    tri_this.vertex(1),tri_this.vertex(2)}, vs, tri_this);
                    for(auto item : res){
                        K2::Point_3 center = centroid(K2::Triangle_3 (item[0],
                                                             item[1],
                                                             item[2]
                                                             ));
                        bool inner_flag = false;
                        for(auto neighbor_id: coverage_intersection[i]) {
                            if (!coverage_field_list[neighbor_id].useful)continue;
                            CGAL::Bounded_side bs = coverage_field_list[neighbor_id].bounded_side(center);
                            if(bs == CGAL::ON_BOUNDED_SIDE){
                                inner_flag = true;
                                break;
                            }
                            else if(bs == CGAL::ON_BOUNDARY){
                                if(tri_this.supporting_plane().oriented_side(coverage_field_list[i].center) !=
                                   tri_this.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) &&
                                   tri_this.supporting_plane().oriented_side(coverage_field_list[i].center)+
                                   tri_this.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) ==0 ){
                                    inner_flag = true;
                                    break;
                                }
                            }
                        }

                        if(!inner_flag){
                            coverage_field_list[i].self_face_delete_flag = false;
                            break;
                        }
                    }
                }
            },i);
        }

        for(int i=0;i<thread_num;i++)
            near_delete[i]->join();
        int cccc = 1;
        for(int i=0;i<mesh->FaceSize();i++){
            if( coverage_field_list[i].self_face_delete_flag) {
    //            int v0 = mesh->fast_iGameFace[i].vh(0);
    //            int v1 = mesh->fast_iGameFace[i].vh(1);
    //            int v2 = mesh->fast_iGameFace[i].vh(2);
    //
    //            fprintf(file24,"v %lf %lf %lf\n",mesh->fast_iGameVertex[v0].x(),
    //                    mesh->fast_iGameVertex[v0].y(),
    //                    mesh->fast_iGameVertex[v0].z()
    //                    );
    //            fprintf(file24,"v %lf %lf %lf\n",mesh->fast_iGameVertex[v1].x(),
    //                    mesh->fast_iGameVertex[v1].y(),
    //                    mesh->fast_iGameVertex[v1].z()
    //            );
    //            fprintf(file24,"v %lf %lf %lf\n",mesh->fast_iGameVertex[v2].x(),
    //                    mesh->fast_iGameVertex[v2].y(),
    //                    mesh->fast_iGameVertex[v2].z()
    //            );
    //            fprintf(file24,"f %d %d %d\n",cccc,cccc+1,cccc+2);
    //            cccc+=3;
    //            cout <<"new delete2" << endl;
                coverage_field_list[i].useful = false;
            }
        }
    }
    //X: 局部采样删除法
    if(1){
        std::vector <std::shared_ptr<std::thread> > find_near(thread_num);
        for(int i=0;i<thread_num;i++) {
            find_near[i] = make_shared<std::thread>([&](int now_id) {
                list<K2::Triangle_3> aabb_tree_face_list;
                map<unsigned long long ,int>aabb_tree_mp;
                for(int i=0;i<mesh->FaceSize();i++){
                    for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++){
                        K2::Triangle_3 tri(coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
                                           coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
                                           coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]);
                        aabb_tree_face_list.push_back(tri);
                        aabb_tree_mp[tri.id()] = i;
                    }
                }
                Tree aabb_tree(aabb_tree_face_list.begin(),aabb_tree_face_list.end());
                for (int i = 0; i < mesh->FaceSize(); i++) {
                    if (i % thread_num != now_id)continue;
                    cout <<"find near:"<<i<<":" <<mesh->FaceSize() << endl;
                    std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
                    aabb_tree.all_intersections(coverage_field_list[i].iso_cuboid_3,std::back_inserter(intersections));
                    for(auto item : intersections) {
                        auto iter = aabb_tree_mp.find(item.second->id());
                        if(iter->second != i)
                            coverage_intersection[i].insert(iter->second);
                    }
                }
            },i);
        }
        for(int i=0;i<thread_num;i++)
            find_near[i]->join();

        // 这里统一算一波绕数
        std::mutex winding_num_mutex;
        vector<K2::Point_3 > winding_vertex;
        std::vector <std::shared_ptr<std::thread> > arrange_winding(thread_num);
        for(int i=0;i<thread_num;i++) {
            arrange_winding[i] = make_shared<std::thread>([&](int now_id) {
                for (int i = 0; i < mesh->FaceSize(); i++) {
                    if (i % thread_num != now_id)continue;
                    if(!coverage_field_list[i].useful)continue;
                    cout <<"sample check "<< i <<"//"<<mesh->FaceSize()<<endl;
                    for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++) {
                        if ((coverage_field_list[i].bound_face_id[j][0] >= 3 ||
                             coverage_field_list[i].bound_face_id[j][1] >= 3 ||
                             coverage_field_list[i].bound_face_id[j][2] >= 3)) {
                            for(int k=0;k<coverage_field_list[i].bound_face_sampling_point[j].size();k++){
                                K2::Point_3 sample = coverage_field_list[i].bound_face_sampling_point[j][k];
                                if(CGAL::approximate_sqrt(origin_face_tree.squared_distance(sample)) <= CGAL::Epeck::FT(min_near_limit)){
                                    coverage_field_list[i].bound_face_sampling_point_state[j][k] = -1;
                                    //cout <<"-300:occur: " <<origin_face_tree.squared_distance(centroid(K2::Triangle_3(v0,v1,v2))) <<" "<< min_near_limit << " "<<default_move<<endl;
                                }
                                else {
                                    std::unique_lock<std::mutex> lock(winding_num_mutex);
                                    coverage_field_list[i].bound_face_sampling_point_state[j][k] = winding_vertex.size();
                                    winding_vertex.push_back(coverage_field_list[i].bound_face_sampling_point[j][k]);
                                    //这里不加mutex 没办法求绕数，那么直接不用omp了，还是继续手写吧
                                }
                            }
                        }
                        else{
                            coverage_field_list[i].bound_face_useful[j] = 0;
                        }
                    }

                }
            }, i);
        }
        for(int i=0;i<thread_num;i++)
            arrange_winding[i]->join();




        vector<int> ret  = winding_num(winding_vertex);



        std::vector <std::shared_ptr<std::thread> > near_delete(thread_num);
        for(int i=0;i<thread_num;i++) {
            near_delete[i] = make_shared<std::thread>([&](int now_id) {
                for (int i = 0; i < mesh->FaceSize(); i++) {
                    if (i % thread_num != now_id)continue;
                    if(!coverage_field_list[i].useful)continue;
                    cout <<"cal near: "<< i <<":"<< mesh->FaceSize()<<endl;
                   // coverage_field_list[i].self_face_delete_flag = true;
                    list<K2::Triangle_3>neighbor_face;
                    map<unsigned long long ,int> face_to_field;
                    for(auto neighbor_id: coverage_intersection[i]){
                        if(!coverage_field_list[neighbor_id].useful)continue;
                        if(neighbor_id == i)continue;
                        for (int k = 0; k < coverage_field_list[neighbor_id].bound_face_id.size(); k++) {
                            K2::Triangle_3 tri_this(coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][0]],
                                                    coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][1]],
                                                    coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][2]]
                            );

                            neighbor_face.push_back(tri_this);
                            face_to_field[tri_this.id()] = neighbor_id;
                        }
                    }
                    Tree aabb_local(neighbor_face.begin(),neighbor_face.end());

                    function<bool(K2::Point_3,K2::Triangle_3)>check_inner_vertex = [&](K2::Point_3 check_point,K2::Triangle_3 tri){
                        K2::Point_3 v0 = tri.vertex(0);
                        K2::Point_3 v1 = tri.vertex(1);
                        K2::Point_3 v2 = tri.vertex(2);
                        K2::Vector_3 d1 = (v1 - v0)/2;
                        K2::Vector_3 d2 = (v2 - v0)/2;
//                        cout <<"xxxxxxxxxxxx"<<endl;
//                        cout <<"v "<<v0<<endl;
//                        cout <<"v "<<v1<<endl;
//                        cout <<"v "<<v2<<endl;
//                        cout <<"f 1 2 3"<<endl;
                        int cnt = 0;


                            //cout << check_point<< endl;
                        for(auto neighbor_id: coverage_intersection[i]) {
                            if (!coverage_field_list[neighbor_id].useful)continue;
                            if(neighbor_id == i)continue;
                            auto side = coverage_field_list[neighbor_id].bounded_side(check_point);
                            if(side == CGAL::ON_BOUNDED_SIDE) {
                                return true;
                            }
                            if(side == CGAL::ON_BOUNDARY) {
                                if( tri.supporting_plane().oriented_side(coverage_field_list[i].center) !=
                                    tri.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) &&
                                    tri.supporting_plane().oriented_side(coverage_field_list[i].center)+
                                    tri.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) ==0){
                                    return true;
                                }
                            }
                        }
                        return false;


                            //下面写这段是新方法//todo 改造成爆掉一个点就炸了 ？ 或者有一个计数器
//                        map<int,set<K2::Point_3> >field_cnt;
//                        set<int>field_error;
//
//                           // cout << "v "<< check_point<<endl;
//                        K2::Vector_3 normal = -1 * tri.supporting_plane().orthogonal_vector();
//                        K2::Ray_3 ray(check_point,normal);
//                        std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
//
//                        aabb_local.all_intersections(ray,std::back_inserter(intersections));
//
//                        for(auto item : intersections) { //todo 这里改成批量插入
//                            int id = face_to_field.find(item.second->id())->second;
//                            if (const K2::Point_3 *p = boost::get<K2::Point_3>(&(item.first))) {
//                                if(*p == check_point)continue;
//                                if(field_cnt[id].count(*p)){
//                                    field_error.insert(id);
//                                }
//                                else
//                                    field_cnt[id].insert(*p);
//                            }
//                            else{
//                                field_error.insert(id);
//                            }
//                        }
//
//                        for(const auto& it : field_cnt){
//                            if(!field_error.count(it.first) && it.second.size() % 2 == 1) {
//                                return true;
//                            }
//                        }
//
//                        for(auto neighbor_id: field_error) {
//                            if (!coverage_field_list[neighbor_id].useful)continue;
//                            if(neighbor_id == i)continue;
//                            auto side = coverage_field_list[neighbor_id].bounded_side(check_point);
//                            if(side == CGAL::ON_BOUNDED_SIDE){
//                                return true;
//                            }
//                            if(side == CGAL::ON_BOUNDARY){
//                                if( tri.supporting_plane().oriented_side(coverage_field_list[i].center) !=
//                                    tri.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) &&
//                                    tri.supporting_plane().oriented_side(coverage_field_list[i].center)+
//                                    tri.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) ==0){
//                                    return true;
//                                }
//                            }
//                        }
//                        return false;
                    };

                    function<bool(vector<K2::Point_3>vv,K2::Triangle_3)>check_inner_vertex_all = [&](vector<K2::Point_3>vv ,K2::Triangle_3 tri){
                        K2::Point_3 v0 = tri.vertex(0);
                        K2::Point_3 v1 = tri.vertex(1);
                        K2::Point_3 v2 = tri.vertex(2);
                        K2::Vector_3 d1 = (v1 - v0)/2;
                        K2::Vector_3 d2 = (v2 - v0)/2;
//                        cout <<"xxxxxxxxxxxx"<<endl;
//                        cout <<"v "<<v0<<endl;
//                        cout <<"v "<<v1<<endl;
//                        cout <<"v "<<v2<<endl;
//                        cout <<"f 1 2 3"<<endl;


                        //cout << check_point<< endl;
                        for(auto neighbor_id: coverage_intersection[i]) {
                            int cnt = 0;
                            if (!coverage_field_list[neighbor_id].useful)continue;
                            if(neighbor_id == i)continue;
                            for(auto check_point : vv) {
                                auto side = coverage_field_list[neighbor_id].bounded_side(check_point);
                                if (side == CGAL::ON_BOUNDED_SIDE) {
                                    cnt++;
                                }
                                if (side == CGAL::ON_BOUNDARY) {
                                    if (tri.supporting_plane().oriented_side(coverage_field_list[i].center) !=
                                        tri.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) &&
                                        tri.supporting_plane().oriented_side(coverage_field_list[i].center) +
                                        tri.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) ==
                                        0) {
                                        cnt++;
                                    }
                                }
                            }
                            if(cnt == vv.size())return true;
                        }
                        return false;
                     };

                    for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++) {
                        K2::Triangle_3 tri_this(coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
                            coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
                            coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]);
                        if(check_inner_vertex_all(coverage_field_list[i].bound_face_sampling_point[j],tri_this)){
                            for(int k=0;k<coverage_field_list[i].bound_face_sampling_point[j].size();k++){
                                coverage_field_list[i].bound_face_sampling_point_state[j][k] = -1;
                            }
                        }
                        else{
                            for(int k=0;k<coverage_field_list[i].bound_face_sampling_point[j].size();k++){
                                if(coverage_field_list[i].bound_face_sampling_point_state[j][k]!=-1 && ret[coverage_field_list[i].bound_face_sampling_point_state[j][k]])
                                    coverage_field_list[i].bound_face_sampling_point_state[j][k] = 0;
                                else
                                    coverage_field_list[i].bound_face_sampling_point_state[j][k] = -1;
                            }
                        }


//                        for(int k=0;k<coverage_field_list[i].bound_face_sampling_point[j].size();k++){
//                            if(coverage_field_list[i].bound_face_sampling_point_state[j][k]!=-1){
//                                if(ret[coverage_field_list[i].bound_face_sampling_point_state[j][k]] ){
//                                    if(check_inner_vertex(coverage_field_list[i].bound_face_sampling_point[j][k],tri_this)){
//                                        coverage_field_list[i].bound_face_sampling_point_state[j][k] = -1;
//                                    }
//                                    else
//                                        coverage_field_list[i].bound_face_sampling_point_state[j][k] = 0;
//                                }
//                                else{
//                                    coverage_field_list[i].bound_face_sampling_point_state[j][k] = -1;
//                                }
//                            }
//                        }

                    }
                }
            },i);
        }

        for(int i=0;i<thread_num;i++)
            near_delete[i]->join();
        for(int i=0;i<mesh->FaceSize();i++) {
            for (int j = 0; j < coverage_field_list[i].bound_face_id.size(); j++) {
                if(coverage_field_list[i].bound_face_useful[j] && std::accumulate(coverage_field_list[i].bound_face_sampling_point_state[j].begin(),
                                   coverage_field_list[i].bound_face_sampling_point_state[j].end(),0
                                   ) == -1 * int(coverage_field_list[i].bound_face_sampling_point_state[j].size())){
                    coverage_field_list[i].bound_face_useful[j] = 2;
                }
            }
        }
        ofstream fsip("../occ2/find_inp.obj");
        ofstream fsop("../occ2/find_outp.obj");
        for(int i=0;i<mesh->FaceSize();i++) {
            for (int j = 0; j < coverage_field_list[i].bound_face_id.size(); j++) {
                for(int k=0;k<coverage_field_list[i].bound_face_sampling_point[j].size();k++){
                    cout << coverage_field_list[i].bound_face_sampling_point_state[j][k] << endl;
                    if(coverage_field_list[i].bound_face_sampling_point_state[j][k]==-1 || coverage_field_list[i].bound_face_useful[j]!=1 ) {
                        fsip <<"v "<< coverage_field_list[i].bound_face_sampling_point[j][k].x() << " "
                             << coverage_field_list[i].bound_face_sampling_point[j][k].y() << " "
                             << coverage_field_list[i].bound_face_sampling_point[j][k].z() << endl;
                    }
                    else{
                        fsop  <<"v "<< coverage_field_list[i].bound_face_sampling_point[j][k].x() << " "
                             << coverage_field_list[i].bound_face_sampling_point[j][k].y() << " "
                             << coverage_field_list[i].bound_face_sampling_point[j][k].z() << endl;
                    }
                }
            }
        }


//        int cnt44 = 1;
//        int cnt54 = 1;
//        for(int i=0;i<mesh->FaceSize();i++) {
//            for (int k = 0; k < coverage_field_list[i].bound_face_id.size(); k++) {
//                if(!coverage_field_list[i].bound_face_useful[k]) {
//
//                    K2::Triangle_3 tri_this(
//                            coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[k][0]],
//                            coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[k][1]],
//                            coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[k][2]]
//                    );
//                    fprintf(file44, "v %lf %lf %lf\n", CGAL::to_double(tri_this.vertex(0).x()),
//                            CGAL::to_double(tri_this.vertex(0).y()),
//                            CGAL::to_double(tri_this.vertex(0).z()));
//                    fprintf(file44, "v %lf %lf %lf\n", CGAL::to_double(tri_this.vertex(1).x()),
//                            CGAL::to_double(tri_this.vertex(1).y()),
//                            CGAL::to_double(tri_this.vertex(1).z()));
//                    fprintf(file44, "v %lf %lf %lf\n", CGAL::to_double(tri_this.vertex(2).x()),
//                            CGAL::to_double(tri_this.vertex(2).y()),
//                            CGAL::to_double(tri_this.vertex(2).z()));
//                    fprintf(file44, "f %d %d %d\n", cnt44, cnt44 + 1, cnt44 + 2);
//                    cnt44 += 3;
//                }
//                else{
//                    K2::Triangle_3 tri_this(
//                            coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[k][0]],
//                            coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[k][1]],
//                            coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[k][2]]
//                    );
//                    fprintf(file54, "v %lf %lf %lf\n", CGAL::to_double(tri_this.vertex(0).x()),
//                            CGAL::to_double(tri_this.vertex(0).y()),
//                            CGAL::to_double(tri_this.vertex(0).z()));
//                    fprintf(file54, "v %lf %lf %lf\n", CGAL::to_double(tri_this.vertex(1).x()),
//                            CGAL::to_double(tri_this.vertex(1).y()),
//                            CGAL::to_double(tri_this.vertex(1).z()));
//                    fprintf(file54, "v %lf %lf %lf\n", CGAL::to_double(tri_this.vertex(2).x()),
//                            CGAL::to_double(tri_this.vertex(2).y()),
//                            CGAL::to_double(tri_this.vertex(2).z()));
//                    fprintf(file54, "f %d %d %d\n", cnt54, cnt54 + 1, cnt54 + 2);
//                    cnt54 += 3;
//                }
//            }
//        }
        cout <<"find near show "<< endl;
    }
    ofstream fsd("../occ2/find_delete.obj");
    ofstream fsl("../occ2/find_last.obj");
    int cnt = 1;
    int cnt2 = 1;// 明天分析每一个步骤的结果
    for(int i=0;i<mesh->FaceSize();i++){
        for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++){
            if(coverage_field_list[i].bound_face_useful[j] != 1){
                int fid0 = coverage_field_list[i].bound_face_id[j][0];
                int fid1 = coverage_field_list[i].bound_face_id[j][1];
                int fid2 = coverage_field_list[i].bound_face_id[j][2];

                fsd <<"v " <<coverage_field_list[i].bound_face_vertex_exact[fid0] << endl;
                fsd <<"v "<< coverage_field_list[i].bound_face_vertex_exact[fid1] << endl;
                fsd <<"v "<< coverage_field_list[i].bound_face_vertex_exact[fid2] << endl;
                fsd <<"f "<< cnt <<" "<< cnt+1 <<" "<< cnt+2 << endl;
                cnt+=3;
            }
            else{
                int fid0 = coverage_field_list[i].bound_face_id[j][0];
                int fid1 = coverage_field_list[i].bound_face_id[j][1];
                int fid2 = coverage_field_list[i].bound_face_id[j][2];

                fsl <<"v " <<coverage_field_list[i].bound_face_vertex_exact[fid0] << endl;
                fsl <<"v "<< coverage_field_list[i].bound_face_vertex_exact[fid1] << endl;
                fsl <<"v "<< coverage_field_list[i].bound_face_vertex_exact[fid2] << endl;
                fsl <<"f "<< cnt2 <<" "<< cnt2+1 <<" "<< cnt2+2 << endl;
                cnt2+=3;
            }
        }
    }

//    cout <<"exit"<<endl;
//    exit(0);


//    std::vector <std::shared_ptr<std::thread> > one_ring_select_thread_pool(thread_num);
//    for(int i=0;i<thread_num;i++) {
//        one_ring_select_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
//            for(int i=0;i<mesh->FaceSize();i++) {
//                if (i % thread_num != now_id)continue;
//                if(i%500==0)
//                    cout << "one_ring_select_thread_pool "<< i << endl;
//                if(!coverage_field_list[i].useful)continue;
//
//
//                set<int>neighbor_field;
//                for (auto neighbor_id: mesh->FastNeighborFhOfFace_[i]) {
//                    if(neighbor_id == i)continue;
//                    if(coverage_field_list[neighbor_id].useful)
//                        neighbor_field.insert(neighbor_id);
//                }
//
//
//                vector<K2::Triangle_3>neighbor_face;
//                for(auto neighbor_id: neighbor_field){
//                    for (int k = 0; k < coverage_field_list[neighbor_id].bound_face_id.size(); k++) {
//                        K2::Triangle_3 tri_this(coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][0]],
//                                                coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][1]],
//                                                coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][2]]
//                        );
//
//                        neighbor_face.push_back(tri_this);
//                    }
//                }
//
//                for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++) {
//                    if ((coverage_field_list[i].bound_face_id[j][0] >= 3 ||
//                         coverage_field_list[i].bound_face_id[j][1] >= 3 ||
//                         coverage_field_list[i].bound_face_id[j][2] >= 3)){
//                        bool flag = false;
//
//                        K2::Triangle_3 tri_this(coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
//                                                coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
//                                                coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]);
//                        K2::Segment_3 e0(tri_this.vertex(0), tri_this.vertex(1));
//                        K2::Segment_3 e1(tri_this.vertex(1), tri_this.vertex(2));
//                        K2::Segment_3 e2(tri_this.vertex(2), tri_this.vertex(0));
//
//                        vector<K2::Segment_3> vs_tmp;
//                        for (auto other: neighbor_face) { //Point_3, or Segment_3, or Triangle_3, or std::vector < Point_3 >
//                            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
//                                    res_tt = intersection(tri_this, other);
//                            if (res_tt) {
//                                if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tt)) {
//                                    vs_tmp.push_back(*s);
//                                } else if (const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&*res_tt)) {
//                                    vs_tmp.emplace_back(t->vertex(0), t->vertex(1));
//                                    vs_tmp.emplace_back(t->vertex(1), t->vertex(2));
//                                    vs_tmp.emplace_back(t->vertex(2), t->vertex(0));
//                                } else if (std::vector<K2::Point_3> *vs = boost::get<std::vector<K2::Point_3 >>(&*res_tt)) {
//                                    sort_by_polar_order(*vs, tri_this.supporting_plane().orthogonal_vector());
//                                    for (int k = 0; k < vs->size(); k++) {
//                                        vs_tmp.emplace_back(vs->operator[](k), vs->operator[]((k + 1) % vs->size()));
//                                        // cerr << "run iiiiiiiiiiiiiiit" << endl;
//                                    }
//                                }
//                            }
//                        }
//                        vector<K2::Segment_3> vs;
//                        for (auto se: vs_tmp) {
//                            if (!segment_in_line(se, e0) && !segment_in_line(se, e1) && !segment_in_line(se, e2)) {
//                                vs.push_back(se);
//                            }
//                        }
//                        vector<vector<K2::Point_3> > res = CGAL_CDT_NEW({coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
//                                                                         coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
//                                                                         coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]}, vs, tri_this);
//                        //cout << res.size() <<endl;
//                        for (auto each_tri: res) {
//                            bool patch_flag = false;
//                            K2::Point_3 center = CGAL::centroid(K2::Triangle_3(each_tri[0], each_tri[1], each_tri[2]));
//                            for (auto j: neighbor_field) {
//                                if (tri_this.supporting_plane().oriented_side(coverage_field_list[i].center) !=
//                                    tri_this.supporting_plane().oriented_side(coverage_field_list[j].center) &&
//                                    tri_this.supporting_plane().oriented_side(coverage_field_list[i].center)+
//                                    tri_this.supporting_plane().oriented_side(coverage_field_list[j].center) ==0 &&
//                                    coverage_field_list[j].in_or_on_field(center)) {
//                                    patch_flag = true;
//                                    break;
//                                }
//                            }
//                            if(!patch_flag){
//                                flag = true;
//                                break;
//                            }
//                        }
//                        coverage_field_list[i].bound_face_useful[j]=(flag?1:0);
//                        //cout <<i<<" "<<j <<" "<< (flag?"yes":"no") << endl;
//                    }
//                    else{
//                        coverage_field_list[i].bound_face_useful[j]=0;
//                    }
//
//                }
//            }
//        },i);
//    }
//
//    for(int i=0;i<thread_num;i++)
//        one_ring_select_thread_pool[i]->join();

    for(auto  each_container_face : container_grid_face){
        each_grid_face_list.push_back({(size_t)each_container_face[0],(size_t)each_container_face[1],(size_t)each_container_face[2]});
        each_grid_face_list.push_back({(size_t)each_container_face[2],(size_t)each_container_face[3],(size_t)each_container_face[0]});
    }

    cout <<"build end "<< endl;

    //cgal_polygon = make_shared<CGALPolygon>(mesh);

    std::list < K2::Triangle_3> triangles;

    mesh->initBBox();

    double max_move = 0;
    double min_move = 1e10;
    for (auto i: mesh->fast_iGameFace) {
        max_move = max(max_move, abs(i.move_dist));
        min_move = min(min_move, abs(i.move_dist));
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
            if (xx >= 0 && yy >= 0 && zz >= 0)
                ret.push_back({xx, yy, zz});
        }
        return ret;
    };
    int fsize = mesh->FaceSize();

    cout <<"bfs start \n" << endl;
    for (int face_id = 0; face_id < fsize; face_id++) {
        if (face_id % 1000 == 0)
            printf("%d/%d\n", face_id, fsize);
        if(!coverage_field_list[face_id].useful)continue;
        MeshKernel::iGameVertex min_point(coverage_field_list[face_id].x_min,
                                          coverage_field_list[face_id].y_min,
                                          coverage_field_list[face_id].z_min
        );
        MeshKernel::iGameVertex max_point(coverage_field_list[face_id].x_max,
                                          coverage_field_list[face_id].y_max,
                                          coverage_field_list[face_id].z_max
        );
        grid g_min = vertex_to_grid(min_point);
        grid g_max = vertex_to_grid(max_point);
        for(int i=g_min.x;i<=g_max.x;i++){
            for(int j=g_min.y;j<=g_max.y;j++){
                for(int k=g_min.z;k<=g_max.z;k++){
                    grid now(i,j,k);
                    auto iter = frame_grid_mp.find(now);
                    if (iter == frame_grid_mp.end()) {
                        iter = frame_grid_mp.insert(make_pair(now, GridVertex())).first;
                    }
                    iter->second.field_list.push_back(MeshKernel::iGameFaceHandle(face_id));
                }
            }
        }
    }




    std::vector <std::shared_ptr<std::thread> > each_grid_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) { //todo 存在偶发性多线程异常
        each_grid_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            int each_grid_cnt =-1;
            for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) { // todo 逻辑修改
                each_grid_cnt++;
                if(each_grid_cnt % thread_num != now_id)continue;
                if(each_grid_cnt % (thread_num) == now_id)
                    printf("thread num :%d each_grid_cnt %d/%d\n",now_id,each_grid_cnt,(int)frame_grid_mp.size());

                MeshKernel::iGameVertex grid_vertex = getGridVertex(each_grid->first, 0);

                sort(each_grid->second.field_list.begin(), each_grid->second.field_list.end());
                each_grid->second.field_list.resize(
                        unique(each_grid->second.field_list.begin(), each_grid->second.field_list.end()) -
                        each_grid->second.field_list.begin());


                K2::Point_3 small = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 0));
                K2::Point_3 big = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 7));

                std::vector<K2::Point_3> ps;
                for(int i=0;i<8;i++){
                    ps.push_back(getGridK2Vertex(small,big,i));
                }
                std::list<K2::Triangle_3 >  frame_faces_list;

                for(auto i : each_grid_face_list){
                    frame_faces_list.emplace_back(ps[i[0]],ps[i[1]],ps[i[2]]);
                }
                Tree frame_aabb_tree(frame_faces_list.begin(),frame_faces_list.end());
                CGAL::Polyhedron_3<K2> frame_poly;
                PMP::polygon_soup_to_polygon_mesh(ps, each_grid_face_list, frame_poly, CGAL::parameters::all_default());

                for(int j = 0; j <  each_grid->second.field_list.size() ; j++){
                    int field_id = each_grid->second.field_list[j];
                    for(int k=0;k<coverage_field_list[field_id].bound_face_id.size();k++){

                        int v0_id = coverage_field_list[field_id].bound_face_id[k][0];
                        int v1_id = coverage_field_list[field_id].bound_face_id[k][1];
                        int v2_id = coverage_field_list[field_id].bound_face_id[k][2];

                        K2::Triangle_3 this_tri(coverage_field_list[field_id].bound_face_vertex_exact[v0_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v1_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v2_id]
                        );
                        //cout <<"useful: "<< coverage_field_list[field_id].bound_face_useful[k] << endl;
                        if(coverage_field_list[field_id].bound_face_useful[k] == 0)continue;
                        if(check_triangle_through_grid(small,big,frame_poly,this_tri)){
                            if(coverage_field_list[field_id].bound_face_useful[k] == 1){
                                each_grid->second.field_face_though_list[field_id].push_back(k);
                                each_grid->second.face_hash_id_map[this_tri.id()] = {field_id,k};
                                each_grid->second.build_aabb_tree_triangle_list.push_back(this_tri);
                            }
                        }

                    }
                }
                Tree this_grid_aabb_tree( each_grid->second.build_aabb_tree_triangle_list.begin(),
                                          each_grid->second.build_aabb_tree_triangle_list.end());
               // cout <<each_grid_cnt<<":::::::" <<each_grid->second.build_aabb_tree_triangle_list.size()<<endl;
                //X: 这里写查询 因为aabb树指针的bug，所以反着做，用每个格子去查询每个面

                for(auto item : each_grid->second.field_face_though_list){
                    int field_id = item.first;
                    for(int bound_face_id : item.second){
                        int v0_id = coverage_field_list[field_id].bound_face_id[bound_face_id][0];
                        int v1_id = coverage_field_list[field_id].bound_face_id[bound_face_id][1];
                        int v2_id = coverage_field_list[field_id].bound_face_id[bound_face_id][2];
                        K2::Triangle_3 this_tri(coverage_field_list[field_id].bound_face_vertex_exact[v0_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v1_id],
                                                coverage_field_list[field_id].bound_face_vertex_exact[v2_id]
                        );

                        std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;

                        this_grid_aabb_tree.all_intersections(this_tri,std::back_inserter(intersections));

                        unique_lock<mutex> lock(mutex);
                        for(auto item : intersections) { //todo 这里改成批量插入
                            if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(item.first))){

                                map<size_t,pair<int,int> >::iterator iter = each_grid->second.face_hash_id_map.find(item.second->id());
                                if(iter->second.first != field_id ) {
                                    coverage_field_list[field_id].bound_face_cutting_segment[bound_face_id].push_back(
                                            *s);
                                }
                            }
                            else if(const K2::Point_3 * p = boost::get<K2::Point_3>(&(item.first))){
                                map<size_t,pair<int,int> >::iterator iter = each_grid->second.face_hash_id_map.find(item.second->id());
                                if(iter->second.first != field_id ) {
                                    coverage_field_list[field_id].bound_face_cutting_point[bound_face_id].push_back(*p);
                                }
                            }
                            else if(const std::vector<K2::Point_3> * v = boost::get<std::vector<K2::Point_3> >(&(item.first)) ){
                                map<size_t,pair<int,int> >::iterator iter = each_grid->second.face_hash_id_map.find(item.second->id());
                                if(iter->second.first != field_id ){
                                    for(int j=0;j<v->size();j++){
                                        coverage_field_list[field_id].bound_face_cutting_point[bound_face_id].push_back(v->at(j));
                                    }
                                }
                            }
                        }
                        coverage_field_list[field_id].bound_face_cross_field_list[bound_face_id].push_back(each_grid->first);
                    }
                }
            }
            cout << "each_grid_cnt end thread num:"<< now_id<<endl;
        }, i);
    }

    for(int i=0;i<thread_num;i++)
        each_grid_thread_pool[i]->join();
    int ll = 0;
    for (int face_id = 0; face_id < fsize; face_id++) {
        if(!coverage_field_list[face_id].useful)continue;
        for(int j=0;j<coverage_field_list[face_id].bound_face_id.size();j++){
            if(coverage_field_list[face_id].bound_face_useful[j] == 1)
            for(int k=0;k<coverage_field_list[face_id].bound_face_cutting_segment[j].size();k++){
                fprintf(file34,"v %lf %lf %lf\n",CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(0).x()),
                        CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(0).y()),
                        CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(0).z())
                        );
                fprintf(file34,"v %lf %lf %lf\n",CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(1).x()),
                        CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(1).y()),
                        CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(1).z())
                );
                fprintf(file34,"l %d %d\n",ll+1,ll+2);
                ll+=2;
            }
            //fprintf(file34)
        }
    }
    fclose(file34);

    //exit(0);

    for (auto each_grid= frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
        //if(!(each_grid->first.x == 26 && each_grid->first.y == 28 && each_grid->first.z == 9  ))continue;
        //if(!(each_grid->first.x == 19 && each_grid->first.y == 19 && each_grid->first.z == 1 ))continue;
        //if(!(each_grid->first.x == 2 && each_grid->first.y == 3 && each_grid->first.z == 0 ))continue;
        auto small  = getGridVertex(each_grid->first,0);
        auto big  = getGridVertex(each_grid->first,7);
        static int f3_id = 1;
        for (int ii = 0; ii < 7; ii++) {
            for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
                int from = ii;
                int to = DirectedGridEdge[ii][jj];
                MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
                MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
                fprintf(file11, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
                fprintf(file11, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
                fprintf(file11, "l %d %d\n", f3_id, f3_id + 1);
                f3_id += 2;
            }
        }
    }


    std::vector <std::shared_ptr<std::thread> > field_vertex_numbering_thread_pool(thread_num);
    std::atomic<int>global_vertex_id_sum;
    for(int i=0;i<thread_num;i++) {
        field_vertex_numbering_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            //for (int field_id = 67; field_id < 68; field_id++) {
            for (int field_id = 0; field_id < fsize; field_id++) {
                if (field_id % thread_num != now_id)continue;
                if(field_id %(10*thread_num) ==0)
                    cout << "start do cdt "<< field_id <<"/"<<fsize << endl;
                coverage_field_list[field_id].field_id = field_id;
                coverage_field_list[field_id].do_cdt();
                coverage_field_list[field_id].renumber();
                global_vertex_id_sum+= coverage_field_list[field_id].renumber_bound_face_vertex.size();
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        field_vertex_numbering_thread_pool[i]->join();


//    for(int i=0;i<thread_num;i++) {
//        field_vertex_numbering_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
    //for (int field_id = 67; field_id < 68; field_id++) {
//    for (int field_id = 0; field_id < fsize; field_id++) {
////                if (field_id % thread_num != now_id)continue;
////                if(field_id %(10*thread_num) ==0)
//        cout << "start do cdt "<< field_id <<"/"<<fsize << endl;
//        coverage_field_list[field_id].field_id = field_id;
//        coverage_field_list[field_id].do_cdt();
//        coverage_field_list[field_id].renumber();
//        global_vertex_id_sum+= coverage_field_list[field_id].renumber_bound_face_vertex.size();
//    }
//        },i);
//    }
   // exit(0);




//    for (int field_id = 0; field_id < fsize; field_id++) {
//    //for (int field_id = 67; field_id < 68; field_id++) {
//    //for (int field_id = 4; field_id < 5; field_id++) {
//        for(int i=0;i<coverage_field_list[field_id].renumber_bound_face_vertex.size();i++){
//            auto pi = coverage_field_list[field_id].renumber_bound_face_vertex[i];
//            for(int j=i+1;j<coverage_field_list[field_id].renumber_bound_face_vertex.size();j++){
//                auto pj = coverage_field_list[field_id].renumber_bound_face_vertex[j];
//                if(pi == pj){
//                    cout <<"pi == pj gg"<<endl;
//                }
//            }
//        }
////        cout <<"2 is :"<< coverage_field_list[field_id].renumber_bound_face_vertex[2]<<endl;
////        cout <<"6 is :"<< coverage_field_list[field_id].renumber_bound_face_vertex[6]<<endl;
////        cout <<"dist:"<< CGAL::squared_distance(coverage_field_list[field_id].renumber_bound_face_vertex[2],
////                                                coverage_field_list[field_id].renumber_bound_face_vertex[6]
////                                                )<<endl;
////        cout <<"1 is :"<< coverage_field_list[field_id].renumber_bound_face_vertex[1]<<endl;
////        cout <<"12 is :"<< coverage_field_list[field_id].renumber_bound_face_vertex[12]<<endl;
////        cout <<"dist:"<< CGAL::squared_distance(coverage_field_list[field_id].renumber_bound_face_vertex[1],
////                                                coverage_field_list[field_id].renumber_bound_face_vertex[12]
////        )<<endl;
//        if(1)for(int i=0;i<coverage_field_list[field_id].renumber_bound_face_id.size();i++){
//            auto fi = K2::Triangle_3 (coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[i][0]],
//                                      coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[i][1]],
//                                      coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[i][2]]
//                                      );
//            for(int j=i+1;j<coverage_field_list[field_id].renumber_bound_face_id.size();j++){
//                auto fj = K2::Triangle_3 (coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[j][0]],
//                                          coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[j][1]],
//                                          coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[j][2]]
//                );
//                //continue;
//                CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
//                        res_tt = intersection(fi, fj);
//                if (res_tt) {
//                    if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tt)) {
//                        if(!segment_coincide_triangle(*s,fi) || !segment_coincide_triangle(*s,fj)){
//                            if(coverage_field_list[field_id].renumber_bound_face_useful[i]!=1 ||
//                                    coverage_field_list[field_id].renumber_bound_face_useful[j]!=1)
//                                continue;
//                            cout <<  "fieldid: "<< field_id<<" "<< i <<" "<< j << endl;
//                            cout <<  "usefuli"<<coverage_field_list[field_id].renumber_bound_face_useful[i]<<endl;
//                            cout <<  "usefulj"<<coverage_field_list[field_id].renumber_bound_face_useful[j]<<endl;
//                            cout <<"itttt"<< coverage_field_list[field_id].renumber_bound_face_id[i][0] <<" "
//                            << coverage_field_list[field_id].renumber_bound_face_id[i][1]<<" "
//                            <<coverage_field_list[field_id].renumber_bound_face_id[i][2]<<endl;
//                            cout <<"jtttt"<< coverage_field_list[field_id].renumber_bound_face_id[j][0] <<" "
//                                 << coverage_field_list[field_id].renumber_bound_face_id[j][1]<<" "
//                                 <<coverage_field_list[field_id].renumber_bound_face_id[j][2]<<endl;
//                            cout << is_same_triangle_cnt(fi,fj,CGAL::Epeck::FT(0))<<endl;
//                            cout << fi << endl;
//                            cout << fj << endl;
//                            cout << *s << endl;
//                            cout <<"se gg "<<endl;
//                            exit(0);
//                        }
//                    } else if (const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&*res_tt)) {
//                        cout <<"t gg"<<endl;
////                        vs_tmp.emplace_back(t->vertex(0), t->vertex(1));
////                        vs_tmp.emplace_back(t->vertex(1), t->vertex(2));
////                        vs_tmp.emplace_back(t->vertex(2), t->vertex(0));
//                    } else if (std::vector<K2::Point_3> *vs = boost::get<std::vector<K2::Point_3 >>(&*res_tt)) {
//                        cout <<"vs gg"<<endl;
////                        sort_by_polar_order(*vs, tri_this.supporting_plane().orthogonal_vector());
////                        for (int k = 0; k < vs->size(); k++) {
////                            vs_tmp.emplace_back(vs->operator[](k), vs->operator[]((k + 1) % vs->size()));
////                            // cerr << "run iiiiiiiiiiiiiiit" << endl;
////                        }
//                    }
//                }
//            }
//        }
//        //cout <<"fieldid" <<coverage_field_list[field_id].renumber_bound_face_vertex.size() <<" : "<<coverage_field_list[field_id].renumber_bound_face_id.size() << endl;
//    }
//    cout <<"check over1"<<endl;
//    vector<K2::Triangle_3 >trilist;
//    list<K2::Triangle_3>llist;
//    for (int field_id = 0; field_id < fsize; field_id++) {
//        for(int i=0;i<coverage_field_list[field_id].renumber_bound_face_id.size();i++){
//            auto fj = K2::Triangle_3 (coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[i][0]],
//                                      coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[i][1]],
//                                      coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[i][2]]);
//            if(coverage_field_list[field_id].renumber_bound_face_useful[i]==1) {
//                trilist.push_back(fj);
//                llist.push_back(fj);
//            }
//        }
//    }
//    Tree tr(llist.begin(),llist.end());
//    for(int i=0;i<trilist.size();i++){
//        cout << i <<"/"<<trilist.size()<<endl;
//        std::list<Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
//        tr.all_intersections(trilist[i], std::back_inserter(intersections));
//        for (auto item: intersections) {
//            if (const K2::Segment_3 *segment = boost::get<K2::Segment_3>(&(item.first))) {
//                if(!segment_coincide_triangle(*segment,trilist[i])){
//                    cout << *segment << endl;
//                    cout <<"se gg "<<endl;
//                    exit(0);
//                }
//            }
//        }
//    }

//    for(int i=0;i<trilist.size();i++){
//        for(int j=i+1;i<trilist.size();j++){
//            CGAL::cpp11::result_of<K2::Intersect_3(K2::Triangle_3, K2::Triangle_3)>::type
//                    res_tt = intersection(trilist[i], trilist[j]);
//            if (res_tt) {
//                if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&*res_tt)) {
//                    if(!segment_coincide_triangle(*s,trilist[i]) || !segment_coincide_triangle(*s,trilist[j])){
//                        cout << *s << endl;
//                        cout <<"se gg "<<endl;
//                        exit(0);
//                    }
//                } else if (const K2::Triangle_3 *t = boost::get<K2::Triangle_3>(&*res_tt)) {
//                    cout <<"t gg"<<endl;
////                        vs_tmp.emplace_back(t->vertex(0), t->vertex(1));
////                        vs_tmp.emplace_back(t->vertex(1), t->vertex(2));
////                        vs_tmp.emplace_back(t->vertex(2), t->vertex(0));
//                } else if (std::vector<K2::Point_3> *vs = boost::get<std::vector<K2::Point_3 >>(&*res_tt)) {
//                    cout <<"vs gg"<<endl;
////                        sort_by_polar_order(*vs, tri_this.supporting_plane().orthogonal_vector());
////                        for (int k = 0; k < vs->size(); k++) {
////                            vs_tmp.emplace_back(vs->operator[](k), vs->operator[]((k + 1) % vs->size()));
////                            // cerr << "run iiiiiiiiiiiiiiit" << endl;
////                        }
//                }
//            }
//        }
//    }
//    cout <<"check over2"<<endl;
//    exit(0);


//    ofstream fsnext("../occ2/fs_next.obj");
//    int next_cnt = 1;
//    for (int field_id = 0; field_id < fsize; field_id++) {
//        for (int i = 0; i < coverage_field_list[field_id].cdt_result.size(); i++) {
//            cout  << "coverage_field_list[field_id].cdt_result_useful[i]"<<coverage_field_list[field_id].cdt_result_useful[i]<<endl;
//            if(coverage_field_list[field_id].cdt_result_useful[i] == 1){
//                fsnext<<"v "<<coverage_field_list[field_id].cdt_result[i][0]<<endl;
//                fsnext<<"v "<<coverage_field_list[field_id].cdt_result[i][1]<<endl;
//                fsnext<<"v "<<coverage_field_list[field_id].cdt_result[i][2]<<endl;
//                fsnext<<"f "<< next_cnt <<" "<<next_cnt+1 <<" "<<next_cnt+2 << endl;
//                next_cnt+=3;
//            }
//        }
//    }
//    fsnext.close();

//    int next_cnt = 1;
//    for (int field_id = 0; field_id < fsize; field_id++) {
//        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_id.size(); i++) {
//            if(coverage_field_list[field_id].renumber_bound_face_useful[i] == 1){
//                fsnext<<"v "<<coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[i][0]]<<endl;
//                fsnext<<"v "<<coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[i][1]]<<endl;
//                fsnext<<"v "<<coverage_field_list[field_id].renumber_bound_face_vertex[coverage_field_list[field_id].renumber_bound_face_id[i][2]]<<endl;
//                fsnext<<"f "<< next_cnt <<" "<<next_cnt+1 <<" "<<next_cnt+2 << endl;
//                next_cnt+=3;
//            }
//
//        }
//    }
//    fsnext.close();
 //   exit(0);
    //global_vertex_list.resize(global_vertex_id_sum);
    global_vertex_list.clear();
    map<K2::Point_3,int> unique_vertices;

    for (int field_id = 0; field_id < fsize; field_id++) {
        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_vertex.size(); i++) {
            auto it = unique_vertices.find(coverage_field_list[field_id].renumber_bound_face_vertex[i]);
            if(it != unique_vertices.end()){
                coverage_field_list[field_id].renumber_bound_face_vertex_global_id[i] = it->second;
            }
            else{
                int c = global_vertex_list.size();
                global_vertex_list.push_back(coverage_field_list[field_id].renumber_bound_face_vertex[i]);
                unique_vertices[coverage_field_list[field_id].renumber_bound_face_vertex[i]] = c;
                coverage_field_list[field_id].renumber_bound_face_vertex_global_id[i] = c;
            }
        }
    }

    int global_face_cnt = 0;
    for (int field_id = 0; field_id < fsize; field_id++) {
        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_id.size(); i++) {
            int tv0 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][0]];
            int tv1 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][1]];
            int tv2 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][2]];
            //if(set<int>{tv0,tv1,tv2}.size() < 3)continue;
            global_face_cnt++;
        }
    }
    //exit(0);
    global_face_list.resize(global_face_cnt);
    global_face_cnt = 0;


    for (int field_id = 0; field_id < fsize; field_id++) {
        for (int i = 0; i < coverage_field_list[field_id].renumber_bound_face_id.size(); i++) {
            int tv0 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][0]];
            int tv1 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][1]];
            int tv2 =  coverage_field_list[field_id].renumber_bound_face_vertex_global_id[coverage_field_list[field_id].renumber_bound_face_id[i][2]];

            int useful = coverage_field_list[field_id].renumber_bound_face_useful[i];
            if(set<int>{tv0,tv1,tv2}.size() < 3)continue;
            for(int j =0 ; j < coverage_field_list[field_id].renumber_bound_face_cross_field_list[i].size();j++){
                grid cross_grid = coverage_field_list[field_id].renumber_bound_face_cross_field_list[i][j];
                frame_grid_mp[cross_grid].global_face_list.push_back(global_face_cnt);
            }
            coverage_field_list[field_id].renumber_bound_face_global_id[i] = global_face_cnt;
            global_face_list[global_face_cnt].field_id = field_id;
            global_face_list[global_face_cnt].idx0 = tv0;
            global_face_list[global_face_cnt].idx1 = tv1;
            global_face_list[global_face_cnt].idx2 = tv2;
            global_face_list[global_face_cnt].useful = useful;
            global_face_list[global_face_cnt].center = CGAL::centroid(K2::Triangle_3(global_vertex_list[tv0],global_vertex_list[tv1],global_vertex_list[tv2]));
            global_face_cnt++;
        }
    }

    //
    // exit(0);
    //exit(0);
    // 接下去就是求交

    cout << " 接下去时求交"<<endl;//deckel.obj
    queue<unordered_map <grid, GridVertex, grid_hash, grid_equal>::iterator > que;
    for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
        que.push(each_grid);
    }
    std::mutex que_mutex;
    std::vector <std::shared_ptr<std::thread> > face_generate_ray_detect_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        face_generate_ray_detect_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            while(true) {
                std::unique_lock<std::mutex> lock(que_mutex);
                if (que.empty()) {
                    lock.unlock();
                    return;
                }
                auto each_grid = que.front();
                que.pop();
                printf("face_generate_ray_detect_thread_pool %d: %d/%d\n", now_id,
                       (int) frame_grid_mp.size() - que.size(), (int) frame_grid_mp.size());

                lock.unlock();
                unordered_map<unsigned long long, pair<int, int> > triangle_mapping;
                std::list<K2::Triangle_3> l;
                K2::Point_3 small = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 0));
                K2::Point_3 big = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 7));
                std::vector<K2::Point_3> ps;
                for (int i = 0; i < 8; i++) {
                    ps.push_back(getGridK2Vertex(small, big, i));
                }

                CGAL::Polyhedron_3<K2> frame_poly;
                PMP::polygon_soup_to_polygon_mesh(ps, each_grid_face_list, frame_poly, CGAL::parameters::all_default());

                std::function<bool(K2::Point_3)> vertex_in_frame = [&](K2::Point_3 v) {
                    return small.x() <= v.x() && v.x() <= big.x() &&
                           small.y() <= v.y() && v.y() <= big.y() &&
                           small.z() <= v.z() && v.z() <= big.z();
                };

                for (int i = 0; i < each_grid->second.field_list.size(); i++) {
                    int field_id = each_grid->second.field_list[i];
                    bool useful = false;
                    if (coverage_field_list[field_id].in_field(midpoint(K2::Segment_3(ps[0], ps[7])))) {
                        useful = true;
                    }
                    if (!useful) {
                        if (vertex_in_frame(coverage_field_list[field_id].center)) {
                            useful = true;
                        }
                    }
                    if (!useful)
                        useful = CGAL::Polygon_mesh_processing::do_intersect(*coverage_field_list[field_id].poly,
                                                                             frame_poly);
                    if (useful) {
                        for (int j = 0; j < coverage_field_list[field_id].renumber_bound_face_global_id.size(); j++) {
                            int fid_global = coverage_field_list[field_id].renumber_bound_face_global_id[j];
                            K2::Point_3 v0 = global_vertex_list[global_face_list[fid_global].idx0];
                            K2::Point_3 v1 = global_vertex_list[global_face_list[fid_global].idx1];
                            K2::Point_3 v2 = global_vertex_list[global_face_list[fid_global].idx2];
                            K2::Triangle_3 tri(v0, v1, v2);
                            triangle_mapping[tri.id()] = {field_id, fid_global};
                            l.push_back(tri);
                        }
                    }
                }
                Tree aabb_tree(l.begin(), l.end());

                for (int i = 0; i < each_grid->second.global_face_list.size(); i++) {
                    int global_face_id = each_grid->second.global_face_list[i];
                    if (global_face_list[global_face_id].useful !=1 )continue;
                    K2::Point_3 v0 = global_vertex_list[global_face_list[global_face_id].idx0];
                    K2::Point_3 v1 = global_vertex_list[global_face_list[global_face_id].idx1];
                    K2::Point_3 v2 = global_vertex_list[global_face_list[global_face_id].idx2];

                    K2::Vector_3 ray_vec = K2::Triangle_3(v0, v1, v2).supporting_plane().orthogonal_vector();

                    K2::Ray_3 ray(global_face_list[global_face_id].center, -ray_vec);
                    std::list<Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections;
                    aabb_tree.all_intersections(ray, std::back_inserter(intersections));
                    for (auto item: intersections) {
                        pair<int, int> belong = triangle_mapping[item.second->id()];
                        int which_field = belong.first;
                        int which_id = belong.second;
                        if (which_field == global_face_list[global_face_id].field_id)continue;
                        if (const K2::Point_3 *p = boost::get<K2::Point_3>(&(item.first))) {
                            if (*p != global_face_list[global_face_id].center)
                                global_face_list[global_face_id].ray_detect_map[which_field].emplace_back(*p, which_id);
//                            else
//                                global_face_list[global_face_id].special_face_id.insert(which_field);
                        } else {
                            global_face_list[global_face_id].special_field_id.insert(which_field);
                        }
                    }
                }
            }

        },i);
    }
    for(int i=0;i<thread_num;i++)
        face_generate_ray_detect_thread_pool[i]->join();
    //exit(0);

    cout <<"start generate"<<endl;
    std::vector <std::shared_ptr<std::thread> > global_face_final_generate_thread_pool(thread_num);
    atomic<int> flag = 0;
    for(int i=0;i<thread_num;i++) {
        global_face_final_generate_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for (int i = 0; i < global_face_list.size(); i++) {
                if(i %(thread_num*1000) ==0)
                    cout << "global_face_list:" << i <<"/"<<global_face_list.size()<< endl;
                if (i % thread_num != now_id)continue;

                if(global_face_list[i].useful != 1 )continue;
                for(std::unordered_map<int,vector<pair<K2::Point_3,int> > >::iterator it = global_face_list[i].ray_detect_map.begin();
                    it != global_face_list[i].ray_detect_map.end(); it++
                        ){////field 交点 对应的面
                    unordered_set<int>se;
                    std::vector<K2::Point_3> intersection_v;
                    for(int j=0;j<it->second.size();j++){
                        if(se.count(it->second[j].second))continue;
                        se.insert(it->second[j].second);
                        intersection_v.push_back(it->second[j].first);
                    }
                    sort(intersection_v.begin(),intersection_v.end(),[&](const K2::Point_3 &a,const K2::Point_3 &b){
                        if(a.x() != b.x()){
                            return a.x() < b.x();
                        }
                        else if(a.y() != b.y()){
                            return a.y() < b.y();
                        }
                        return a.z() < b.z();
                    });
                    int old_size = intersection_v.size();
                    intersection_v.resize(unique(intersection_v.begin(),intersection_v.end())-intersection_v.begin());
                    if(old_size != intersection_v.size()){
                        global_face_list[i].special_field_id.insert(it->first);
                        continue;
                    }
                    if(intersection_v.size()%2 == 1) {
                        global_face_list[i].useful = -200;
                        break;
                    }
                    // 1 是自己  // 若全无交点，可想而知一定在外，若有交点，那不一定在外，则正方向被覆盖，然后需要探测负方向，和自己
                    // 这里加入自己射线正负方向的寻找
                    // 2 出现多次 加入
                }


                if(  global_face_list[i].useful == 1 ){
                    K2::Point_3 v0 = global_vertex_list[global_face_list[i].idx0];
                    K2::Point_3 v1 = global_vertex_list[global_face_list[i].idx1];
                    K2::Point_3 v2 = global_vertex_list[global_face_list[i].idx2];
                    if (running_mode == 2) {
                       // if (!cgal_polygon->inMesh(centroid(K2::Triangle_3(v0, v1, v2)))) {
                            //global_face_list[i].useful = -300;
                       // }
                    }
                    else{
                       // if (!cgal_polygon->outMesh(centroid(K2::Triangle_3(v0, v1, v2)))) {
                            //global_face_list[i].useful = -300;
                       // }
                    }

                    if(CGAL::approximate_sqrt(origin_face_tree.squared_distance(centroid(K2::Triangle_3(v0,v1,v2)))) <= CGAL::Epeck::FT(min_near_limit)){
                        //cout <<"-300:occur: " <<origin_face_tree.squared_distance(centroid(K2::Triangle_3(v0,v1,v2))) <<" "<< min_near_limit << " "<<default_move<<endl;
                        global_face_list[i].useful = -300;
                    }
                }
            }
        },i);
    }

    for(int i=0;i<thread_num;i++)
        global_face_final_generate_thread_pool[i]->join();

   // exit(0);


    std::vector <std::shared_ptr<std::thread> > special_case_detect(thread_num);
    for(int i=0;i<thread_num;i++) {
        special_case_detect[i] = make_shared<std::thread>([&](int now_id) {
            for (int i = 0; i < global_face_list.size(); i++) {
                if (i % thread_num != now_id)continue;
                bool inner_flag = true;
                if(global_face_list[i].useful == 1 && !global_face_list[i].special_field_id.empty() )[[unlikely]]{
                    cout <<"meeting rare case, do special method"<<std::endl;
                    K2::Triangle_3 tri(global_vertex_list[global_face_list[i].idx0],
                                       global_vertex_list[global_face_list[i].idx1],
                                       global_vertex_list[global_face_list[i].idx2]);
                    for(auto field_id : global_face_list[i].special_field_id){
                        if(in_single_coverage_field_ray_detect(field_id,tri)){
                            inner_flag = false;
                            break;
                        }
                    }
                    if(!inner_flag){
                        global_face_list[i].useful = -200;
                    }
                }
            }
        },i);
    }
    for(int i=0;i<thread_num;i++)
        special_case_detect[i]->join();

//    for(int i=0;i<global_vertex_list.size();i++){
//        fprintf(file6,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[i].x()),
//                CGAL::to_double(global_vertex_list[i].y()),
//                CGAL::to_double(global_vertex_list[i].z()));
//    }

    double sum_avg_edge = 0;

    vector<FinalFace> f_list;
    set<vector<int> >final_face_set;
    for (int i = 0; i < global_face_list.size(); i++) {
        //cout <<"global_face_list["<<i<<"].useful:" <<global_face_list[i].useful << endl;
        if(global_face_list[i].useful == 1) {
            vector<int> tmp{global_face_list[i].idx0,
                            global_face_list[i].idx1,
                            global_face_list[i].idx2};
            sort(tmp.begin(),tmp.end());
            if(final_face_set.count(tmp)){
                continue;
            }
            final_face_set.insert(tmp);

            FinalFace f(global_face_list[i].idx0,
                        global_face_list[i].idx1,
                        global_face_list[i].idx2,
                        centroid((K2::Triangle_3(global_vertex_list[global_face_list[i].idx0], global_vertex_list[global_face_list[i].idx1], global_vertex_list[global_face_list[i].idx2]))),
                        true);
            //if(!final_face_set.count(f))
            f_list.push_back(f);
        }
    }
    winding_num(f_list);
    map<int,int>is_vertex_useful;
    vector<int>final_vertex_useful;
    vector<K2::Point_3>final_vertex_useful_point;
    for (int i = 0; i < f_list.size(); i++) {
        if (f_list[i].flag) {
            for(auto id : {f_list[i].f0,f_list[i].f1,f_list[i].f2}){
                if(!is_vertex_useful.count(id)){
                    is_vertex_useful[id] = final_vertex_useful.size();
                    final_vertex_useful.push_back(id);
                }
            }
        }
    }
    for(int i=0;i<final_vertex_useful.size();i++){
        fprintf(file6,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[final_vertex_useful[i]].x()+start_x),
                CGAL::to_double(global_vertex_list[final_vertex_useful[i]].y()+start_y),
                CGAL::to_double(global_vertex_list[final_vertex_useful[i]].z()+start_z));
        final_vertex_useful_point.emplace_back(global_vertex_list[final_vertex_useful[i]].x()/*+CGAL::Epeck::FT(start_x)*/,
                                               global_vertex_list[final_vertex_useful[i]].y()/*+CGAL::Epeck::FT(start_y)*/,
                                               global_vertex_list[final_vertex_useful[i]].z()/*+CGAL::Epeck::FT(start_z)*/
                                               );
    }
    map<pair<int,int>,int >final_topo_check_mp;
    // 做一个冲突删面的逻辑

    vector<int>neighbor_build[is_vertex_useful.size()+10];
    map<pair<int,int>,vector<int> > check_manifold;
    function<pair<int,int>(int,int)> get_pair = [&](int a,int b){
        if(a>b)swap(a,b);
        return make_pair(a,b);
    };
    map<pair<int,int> ,FinalFace> hole_mp;


    // 这里新增检查逻辑，优先级： 相邻的面中存在更多的正确反向边。

    map<pair<int,int> ,vector<int> >check_conflict;
    map<int,vector<int> > conflict_with;
    vector<int> conflict_vec;
    for (int i = 0; i < f_list.size(); i++) {
        if (f_list[i].flag) {
            for(pair<int,int> edge : vector<pair<int,int> >{
                    std::make_pair(f_list[i].f0,f_list[i].f1),
                    std::make_pair(f_list[i].f1,f_list[i].f2),
                    std::make_pair(f_list[i].f2,f_list[i].f0)
            }
            ) {
                check_conflict[edge].push_back(i);
            }
        }
    }
    for(auto i : check_conflict){
        if(i.second.size() >= 2) {
            for(int j=0;j<i.second.size();j++){
                for(int k=0;k<i.second.size();k++){
                    if(j!=k){
                        conflict_with[i.second[j]].push_back(i.second[k]);
                        conflict_with[i.second[k]].push_back(i.second[j]);
                    }
                }
            }
        }
    }
    for(auto i : conflict_with){
        if(i.second.size() > 1){
            //conflict_vec.push_back(i.first);
            for(int j=0;j<i.second.size();j++){
                f_list[i.second[j]].flag = false;
            }
        }
    }
//    sort(conflict_vec.begin(),conflict_vec.end(),[&](int a,int b){
//        return  conflict_with[a].size() > conflict_with[b].size();
//    });
//
//    for(int i=0;i<conflict_vec.size();i++){
//        bool flag = false;
//        for(auto j: conflict_with[conflict_vec[i]]){
//            if(f_list[j].flag){
//                cout <<"deal conflict"<<endl;
//                flag = true;
//            }
//        }
//        if(flag){
//            f_list[conflict_vec[i]].flag = false;
//        }
//    }
//




    vector<vector<int> >final_face_list;



    for (int i = 0; i < f_list.size(); i++) {

        if(f_list[i].flag) {
            if(running_mode == 1) {
//                sum_avg_edge += sqrt(CGAL::to_double((CGAL::squared_distance(global_vertex_list[global_face_list[i].idx0] , global_vertex_list[global_face_list[i].idx1]))));
//                sum_avg_edge += sqrt(CGAL::to_double((CGAL::squared_distance(global_vertex_list[global_face_list[i].idx2] , global_vertex_list[global_face_list[i].idx1]))));
//                sum_avg_edge += sqrt(CGAL::to_double((CGAL::squared_distance(global_vertex_list[global_face_list[i].idx0] , global_vertex_list[global_face_list[i].idx2]))));

                fprintf(file6, "f %d %d %d\n", is_vertex_useful[f_list[i].f0] + 1, is_vertex_useful[f_list[i].f2] + 1,
                        is_vertex_useful[f_list[i].f1] + 1);
            }
            else {
                fprintf(file6, "f %d %d %d\n", is_vertex_useful[f_list[i].f0] + 1, is_vertex_useful[f_list[i].f1] + 1,
                        is_vertex_useful[f_list[i].f2] + 1);
//                check_manifold[get_pair(f_list[i].f0,f_list[i].f1)].push_back(f_list[i].f2);
//                check_manifold[get_pair(f_list[i].f0,f_list[i].f2)].push_back(f_list[i].f1);
//                check_manifold[get_pair(f_list[i].f1,f_list[i].f2)].push_back(f_list[i].f0);
//                check_manifold[get_pair(f_list[i].f1,f_list[i].f0)].push_back(f_list[i].f2);
//                check_manifold[get_pair(f_list[i].f2,f_list[i].f0)].push_back(f_list[i].f1);
//                check_manifold[get_pair(f_list[i].f2,f_list[i].f1)].push_back(f_list[i].f0);
                if(hole_mp.count(std::make_pair(f_list[i].f0,f_list[i].f1))){
                    auto old = hole_mp[std::make_pair(f_list[i].f0,f_list[i].f1)];
                    cout <<"hole_mp. error" << f_list[i].f0 <<" "<<f_list[i].f1 <<  endl;
                    cout <<"old:" <<old.f0 <<" "<< old.f1 <<" "<< old.f2 << endl;
                    cout <<"new:" << f_list[i].f0 <<" "<< f_list[i].f1 <<" "<< f_list[i].f2 << endl;
                    cout <<"old.f0 :"<< old.f0 <<" "<<CGAL::to_double(global_vertex_list[old.f0].x()) <<" "
                    <<CGAL::to_double(global_vertex_list[old.f0].y()) <<" "
                    <<CGAL::to_double(global_vertex_list[old.f0].z()) << endl;
                    cout <<"old.f1 :"<< old.f1 <<" "<<CGAL::to_double(global_vertex_list[old.f1].x()) <<" "
                    <<CGAL::to_double(global_vertex_list[old.f1].y()) <<" "
                    <<CGAL::to_double(global_vertex_list[old.f1].z()) << endl;
                    cout <<"old.f2 :"<< old.f2 <<" "<<CGAL::to_double(global_vertex_list[old.f2].x()) <<" "
                    <<CGAL::to_double(global_vertex_list[old.f2].y()) <<" "
                    <<CGAL::to_double(global_vertex_list[old.f2].z()) << endl;
                    cout << (K2::Triangle_3(global_vertex_list[old.f0],
                                            global_vertex_list[old.f1],
                                            global_vertex_list[old.f2])).is_degenerate()<<endl;

                    cout <<"new.f0 :"<< f_list[i].f0 <<" "<<CGAL::to_double(global_vertex_list[f_list[i].f0].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f0].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f0].z()) << endl;
                    cout <<"new.f1 :"<< f_list[i].f1 <<" "<<CGAL::to_double(global_vertex_list[f_list[i].f1].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f1].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f1].z()) << endl;
                    cout <<"new.f2 :"<< f_list[i].f2 <<" "<<CGAL::to_double(global_vertex_list[f_list[i].f2].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f2].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f2].z()) << endl;


                    cout << (K2::Triangle_3(global_vertex_list[f_list[i].f0],
                                            global_vertex_list[f_list[i].f1],
                                            global_vertex_list[f_list[i].f2])).is_degenerate()<<endl;
                }
                if(hole_mp.count(std::make_pair(f_list[i].f1,f_list[i].f2))){
                    auto old = hole_mp[std::make_pair(f_list[i].f1,f_list[i].f2)];
                    cout <<"hole_mp. error" << f_list[i].f1 <<" "<<f_list[i].f2 <<  endl;
                    cout <<"old:" <<old.f0 <<" "<< old.f1 <<" "<< old.f2 << endl;
                    cout <<"new:" << f_list[i].f0 <<" "<< f_list[i].f1 <<" "<< f_list[i].f2 << endl;
                    cout <<"old.f0 :"<< old.f0 <<" "<<CGAL::to_double(global_vertex_list[old.f0].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f0].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f0].z()) << endl;
                    cout <<"old.f1 :"<< old.f1 <<" "<<CGAL::to_double(global_vertex_list[old.f1].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f1].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f1].z()) << endl;
                    cout <<"old.f2 :"<< old.f2 <<" "<<CGAL::to_double(global_vertex_list[old.f2].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f2].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f2].z()) << endl;

                    cout << (K2::Triangle_3(global_vertex_list[old.f0],
                                            global_vertex_list[old.f1],
                                            global_vertex_list[old.f2])).is_degenerate()<<endl;

                    cout <<"new.f0 :"<< f_list[i].f0 <<" "<<CGAL::to_double(global_vertex_list[f_list[i].f0].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f0].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f0].z()) << endl;
                    cout <<"new.f1 :"<< f_list[i].f1 <<" "<<CGAL::to_double(global_vertex_list[f_list[i].f1].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f1].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f1].z()) << endl;
                    cout <<"new.f2 :"<< f_list[i].f2 <<" "<<CGAL::to_double(global_vertex_list[f_list[i].f2].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f2].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f2].z()) << endl;

                    cout << (K2::Triangle_3(global_vertex_list[f_list[i].f0],
                                            global_vertex_list[f_list[i].f1],
                                            global_vertex_list[f_list[i].f2])).is_degenerate()<<endl;

                }
                if(hole_mp.count(std::make_pair(f_list[i].f2,f_list[i].f0))){
                    auto old = hole_mp[std::make_pair(f_list[i].f2,f_list[i].f0)];
                    cout <<"hole_mp. error" << f_list[i].f2 <<" "<<f_list[i].f0 <<  endl;
                    cout <<"old:" <<old.f0 <<" "<< old.f1 <<" "<< old.f2 << endl;
                    cout <<"new:" << f_list[i].f0 <<" "<< f_list[i].f1 <<" "<< f_list[i].f2 << endl;
                    cout <<"old.f0 :"<< old.f0 <<" "<<CGAL::to_double(global_vertex_list[old.f0].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f0].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f0].z()) << endl;
                    cout <<"old.f1 :"<< old.f1 <<" "<<CGAL::to_double(global_vertex_list[old.f1].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f1].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f1].z()) << endl;
                    cout <<"old.f2 :"<< old.f2 <<" "<<CGAL::to_double(global_vertex_list[old.f2].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f2].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[old.f2].z()) << endl;
                    cout << (K2::Triangle_3(global_vertex_list[old.f0],
                                            global_vertex_list[old.f1],
                                            global_vertex_list[old.f2])).is_degenerate()<<endl;

                    cout <<"new.f0 :"<< f_list[i].f0 <<" "<<CGAL::to_double(global_vertex_list[f_list[i].f0].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f0].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f0].z()) << endl;
                    cout <<"new.f1 :"<< f_list[i].f1 <<" "<<CGAL::to_double(global_vertex_list[f_list[i].f1].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f1].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f1].z()) << endl;
                    cout <<"new.f2 :"<< f_list[i].f2 <<" "<<CGAL::to_double(global_vertex_list[f_list[i].f2].x()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f2].y()) <<" "
                         <<CGAL::to_double(global_vertex_list[f_list[i].f2].z()) << endl;
                    cout << (K2::Triangle_3(global_vertex_list[f_list[i].f0],
                                            global_vertex_list[f_list[i].f1],
                                            global_vertex_list[f_list[i].f2])).is_degenerate()<<endl;
                }
                hole_mp[std::make_pair(f_list[i].f0,f_list[i].f1)] = f_list[i];
                hole_mp[std::make_pair(f_list[i].f1,f_list[i].f2)] = f_list[i];
                hole_mp[std::make_pair(f_list[i].f2,f_list[i].f0)] = f_list[i];
                final_topo_check_mp[std::make_pair(is_vertex_useful[f_list[i].f0],is_vertex_useful[f_list[i].f1])] = final_face_list.size();
                final_topo_check_mp[std::make_pair(is_vertex_useful[f_list[i].f1],is_vertex_useful[f_list[i].f2])] = final_face_list.size();
                final_topo_check_mp[std::make_pair(is_vertex_useful[f_list[i].f2],is_vertex_useful[f_list[i].f0])] = final_face_list.size();
                final_face_list.push_back(vector<int>{is_vertex_useful[f_list[i].f0],
                                                      is_vertex_useful[f_list[i].f1],
                                                      is_vertex_useful[f_list[i].f2],
                                                      1});
            }
        }
    }
    map<int,vector<int> >hole_line;
    vector<vector<int> > each_hole;
    queue<int>hole_q;
    int debug_hole_cnt = 1;
    for(auto iter : hole_mp){
        int v0 = iter.first.first;
        int v1 = iter.first.second;
        if(!hole_mp.count(std::make_pair(v1,v0))){
//            cout << "v0" << v0 <<" "<< "v1"<<v1<<endl;
//            cout << "v0 :"<<v0<<" 具体 " <<global_vertex_list[v0]<< endl;
//            cout << "v1 :"<<v1<<" 具体 " <<global_vertex_list[v1]<< endl;
            hole_line[v1].push_back(v0);
//            cout <<"hole_line insert"<< v1 << endl;
//            cout <<"succ v1:"<<v1 <<"  v0:"<<v0 <<" f:"<<iter.second.f0 <<" "<<iter.second.f1 <<" "<< iter.second.f2 << endl;

        }
    }
    cout <<hole_line.size() << endl;
    for(auto i : hole_line){
        cout << i.first <<" "<< i.second.size()<<" ";
        for(auto j : i.second)
            cout << j <<" ";
        cout << endl;
    }

    int ss = is_vertex_useful.size();
    int cnt676 = 1;
    disassemble_circle(hole_line,each_hole);
    for(int i=0;i<each_hole.size();i++) {
        for(int j=0;j<each_hole[i].size();j++) {
            fprintf(file676,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[each_hole[i][j]].x()) + start_x,
                    CGAL::to_double(global_vertex_list[each_hole[i][j]].y()) + start_y,
                    CGAL::to_double(global_vertex_list[each_hole[i][j]].z()) + start_z
            );
            fprintf(file676,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[each_hole[i][(j+1)%each_hole[i].size()]].x()) + start_x,
                    CGAL::to_double(global_vertex_list[each_hole[i][(j+1)%each_hole[i].size()]].y()) + start_y,
                    CGAL::to_double(global_vertex_list[each_hole[i][(j+1)%each_hole[i].size()]].z()) + start_z
            );
            fprintf(file676,"l %d %d\n",cnt676,cnt676+1);
            cnt676+=2;
        }
    }
    // 一会抓出来反法向的面，检测所有边是否出现冲突，进行调整。

//    function<void(vector<int>)>


    fprintf(file6,"#*********\n");
    for(int i=0;i<each_hole.size();i++) {

        if(each_hole[i].size() == 3){
            fprintf(file6, "f %d %d %d\n", is_vertex_useful[each_hole[i][0]] + 1, is_vertex_useful[each_hole[i][1]] + 1,
                    is_vertex_useful[each_hole[i][2]] + 1);
//            final_topo_check_mp[std::make_pair(is_vertex_useful[each_hole[i][0]],is_vertex_useful[each_hole[i][1]])] = final_face_list.size();
//            final_topo_check_mp[std::make_pair(is_vertex_useful[each_hole[i][1]],is_vertex_useful[each_hole[i][2]])] = final_face_list.size();
//            final_topo_check_mp[std::make_pair(is_vertex_useful[each_hole[i][2]],is_vertex_useful[each_hole[i][0]])] = final_face_list.size();
            final_face_list.emplace_back(vector<int>{is_vertex_useful[each_hole[i][0]],
                                                     is_vertex_useful[each_hole[i][1]],
                                                     is_vertex_useful[each_hole[i][2]],
                                                     1});
            continue;
        }
        else{
            vector<K2::Point_3>position;
            for(int j=0;j<each_hole[i].size();j++) {
                position.push_back(global_vertex_list[each_hole[i][j]]);
                //origin_face_tree
            }
            vector<vector<int> > ret;
            fix_hole(each_hole[i],position,&origin_face_tree,ret);
            cout <<"ret.size():" <<ret.size() << endl;
            for(int j=0;j<ret.size();j++){
                if(ret[j].size() == 3){
                    final_face_list.emplace_back(vector<int>{is_vertex_useful[ret[j][0]],
                                                             is_vertex_useful[ret[j][1]],
                                                             is_vertex_useful[ret[j][2]],
                                                             1});
                    cout <<"pre_ret j"<< endl;
                }
                else{
                    cout <<"pre_ret_error"<<endl;
                }
            }
        }



//        K2::Vector_3 vec(0,0,0);
//        for(int j=0;j<each_hole[i].size();j++) {
//            vec += global_vertex_list[each_hole[i][j]] - K2::Point_3 (0,0,0);
//        }
//        vec /= each_hole[i].size();
//
//        fprintf(file6, "v %lf %lf %lf\n", CGAL::to_double(vec.x()) + start_x,
//                CGAL::to_double(vec.y()) + start_y,
//                CGAL::to_double(vec.z()) + start_z);
//        final_vertex_useful_point.emplace_back(vec.x()/*+CGAL::Epeck::FT(start_x)*/,
//                                               vec.y()/*+CGAL::Epeck::FT(start_y)*/,
//                                               vec.z()/*+CGAL::Epeck::FT(start_z)*/
//                                               );
//        for(int j=0;j<each_hole[i].size();j++) {
//            fprintf(file6, "f %d %d %d\n",ss + 1 ,is_vertex_useful[each_hole[i][j]] + 1,
//                    is_vertex_useful[each_hole[i][(j+1)%each_hole[i].size()]] + 1);
////            final_topo_check_mp[std::make_pair(ss, is_vertex_useful[each_hole[i][j]] ) ] = final_face_list.size();
////            final_topo_check_mp[std::make_pair(is_vertex_useful[each_hole[i][j]] , is_vertex_useful[each_hole[i][(j+1)%each_hole[i].size()]]) ] = final_face_list.size();
////            final_topo_check_mp[std::make_pair(is_vertex_useful[each_hole[i][(j+1)%each_hole[i].size()]], ss) ] = final_face_list.size();
//            final_face_list.emplace_back(vector<int>{ss,
//                                                     is_vertex_useful[each_hole[i][j]],
//                                                     is_vertex_useful[each_hole[i][(j+1)%each_hole[i].size()]],
//                                                     1});
//        }
//        ss++;
    }

    for(int times = 0; times < 3 ; times++) {
        for(int i=0;i<final_face_list.size();i++){
            if(final_face_list[i][3]) {
                final_topo_check_mp[std::make_pair(final_face_list[i][0], final_face_list[i][1] ) ] = i;
                final_topo_check_mp[std::make_pair(final_face_list[i][1], final_face_list[i][2] ) ] = i;
                final_topo_check_mp[std::make_pair(final_face_list[i][2], final_face_list[i][0] ) ] = i;
            }
        }

        bool final_flag = false;
        for(auto iter : final_topo_check_mp) {
            int v0 = iter.first.first;
            int v1 = iter.first.second;
            if(!final_topo_check_mp.count(std::make_pair(v1,v0))){
                cout <<"lack of "<< v1 <<" "<< v0 << endl;
                cout <<"deleta face" << iter.second<< endl;
                final_face_list[iter.second][3] = 0;
                final_flag = true;
            }
        }
        if(!final_flag) {
            cout <<"final ok" << endl;
            break;
        }
        final_topo_check_mp.clear();
        for(int i=0;i<final_face_list.size();i++){
            if(final_face_list[i][3]) {
                final_topo_check_mp[std::make_pair(final_face_list[i][0], final_face_list[i][1])] = i;
                final_topo_check_mp[std::make_pair(final_face_list[i][1], final_face_list[i][2])] = i;
                final_topo_check_mp[std::make_pair(final_face_list[i][2], final_face_list[i][0])] = i;
            }
        }
        hole_line.clear();
        each_hole.clear();
        for(auto iter : final_topo_check_mp) {
            int v0 = iter.first.first;
            int v1 = iter.first.second;
            if(!final_topo_check_mp.count(std::make_pair(v1,v0))){
                hole_line[v1].push_back(v0);
            }
        }
        disassemble_circle(hole_line,each_hole);
        for(int i=0;i<each_hole.size();i++) {
            if(each_hole[i].size() == 3){
                final_face_list.emplace_back(vector<int>{each_hole[i][0],
                                                         each_hole[i][1],
                                                         each_hole[i][2],
                                                         1});
                continue;
            }
            vector<K2::Point_3>position;
            for(int j=0;j<each_hole[i].size();j++) {
                position.push_back(final_vertex_useful_point[each_hole[i][j]]);
                //origin_face_tree
            }
            vector<vector<int> > ret;
            fix_hole(each_hole[i],position,&origin_face_tree,ret);
            for(int j=0;j<ret.size();j++){
                if(ret[j].size() == 3){
                    final_face_list.emplace_back(vector<int>{ret[j][0],
                                                             ret[j][1],
                                                             ret[j][2],
                                                         1});
                    cout <<"ret j"<< endl;
                }
                else{
                    cout <<"ret_error"<<endl;
                }
            }
        }





//        for(int i=0;i<each_hole.size();i++) {
//            if(each_hole[i].size() == 3){
//                final_face_list.emplace_back(vector<int>{each_hole[i][0],
//                                                         each_hole[i][1],
//                                                         each_hole[i][2],
//                                                         1});
//                continue;
//            }
//            cout <<"add p "<< endl;
//            K2::Vector_3 vec(0,0,0);
//            for(int j=0;j<each_hole[i].size();j++) {
//                vec += final_vertex_useful_point[each_hole[i][j]] - K2::Point_3 (0,0,0);
//            }
//            vec /= each_hole[i].size();
//
//            final_vertex_useful_point.emplace_back(vec.x(),
//                                                   vec.y(),
//                                                   vec.z()
//            );
//            for(int j=0;j<each_hole[i].size();j++) {
//                final_face_list.emplace_back(vector<int>{(int)final_vertex_useful_point.size()-1,
//                                                         each_hole[i][j],
//                                                         each_hole[i][(j+1)%each_hole[i].size()],
//                                                         1});
//            }
//        }

    }
    for(int i=0;i<final_vertex_useful_point.size();i++){
        fprintf(file8, "v %lf %lf %lf\n",
                CGAL::to_double(final_vertex_useful_point[i].x()+CGAL::Epeck::FT(start_x)),
                CGAL::to_double(final_vertex_useful_point[i].y()+CGAL::Epeck::FT(start_y)),
                CGAL::to_double(final_vertex_useful_point[i].z()+CGAL::Epeck::FT(start_z)));
    }
    for(int i=0;i<final_face_list.size();i++){
        if(final_face_list[i][3])
            fprintf(file8, "f %d %d %d\n",
                    final_face_list[i][0]+1,
                    final_face_list[i][1]+1,
                    final_face_list[i][2]+1
                    );
    }







//    for(int i=0;i<each_hole.size();i++){
//        cout <<"circle :"<< i <<" ";
//        for(int j=0;j<each_hole[i].size();j++){
//            cout << each_hole[i][j] <<" ";
//        }
//        cout << endl;
//    }


//    int c_id_cnt = is_vertex_useful.size();
//    while(0){
//        int now = hole_q.front();
//        hole_q.pop();
//        if(hole_line[now]!=-1){
//            int ss = hole_line[now];
//            vector<int>boundary;
//            //boundary.push_back(now);
//            boundary.push_back(ss);
//            while(ss != now){
//                ss = hole_line[ss];
//                cout <<"ss" << ss << endl;
//                boundary.push_back(ss);
//            }
//            cout <<"boundary:";
//            for(int i=0;i<boundary.size();i++){
//                cout <<boundary[i] <<" ";
//            }
//            cout <<endl;
//            if(boundary.size() == 3){
//                cout <<"add face in 3:"<< boundary.size()<< endl;
//                fprintf(file6,"f %d %d %d\n",is_vertex_useful[boundary[0]] + 1, is_vertex_useful[boundary[1]] + 1,
//                        is_vertex_useful[boundary[2]] + 1);
//            }
//            else {
//                cout <<"add face "<< boundary.size()<< endl;
//                K2::Vector_3 center_vec(0,0,0);
//                for (int i = 0; i < boundary.size(); i++) {
//                    center_vec += global_vertex_list[boundary[i]] - K2::Point_3(0, 0, 0);
//                }
//                center_vec /= boundary.size();
//                K2::Point_3 p = K2::Point_3(0, 0, 0) + center_vec;
//                int cid = c_id_cnt++;
//                fprintf(file6,"v %lf %lf %lf\n",CGAL::to_double(p.x()),
//                        CGAL::to_double(p.y()),
//                        CGAL::to_double(p.z()));
//
//                for(int i=0;i<boundary.size();i++){
//                    int fi = is_vertex_useful[boundary[i]];
//                    int se = is_vertex_useful[boundary[(i+1)%boundary.size()]];
//                    fprintf(file6,"f %d %d %d\n",fi,se,cid);
//                }
//
//            }
//            for(int i=0;i<boundary.size();i++){
//                hole_line[boundary[i]] = -1;
//            }
////            center_vec =
//
//
//        }
//    }







//    for(auto i : check_manifold){
//        if(i.second.size() != 4){
//            cout << i.first.first <<" "<< i.first.second <<" : ";
//            for(int j: i.second){
//                cout << j <<" ";
//            }
//            cout << endl;
//            cout << is_vertex_useful[i.first.first] <<" "<< is_vertex_useful[i.first.second] <<" : ";
//            for(int j: i.second){
//                cout << is_vertex_useful[j] <<" ";
//            }
//            cout <<"*******"<< endl;
//            cout << endl;
//
//            //exit(0);
//        }
//    }


    fclose(file6);
    //sum_avg_edge /=(global_face_list.size()*3*20);

    auto end_clock = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_clock - start_clock;

    string new_name = (input_filename + "_tmp.obj");
    string cmd;
    string tetwild_ex = "";
    if(tetwild_l>0){
        tetwild_ex += " -l "+to_string(tetwild_l) +" ";
    }
    if(tetwild_e>0){
        tetwild_ex += " -e "+to_string(tetwild_e) +" ";
    }
    //exit(0);
    if(result_mode == 1) {
        tetwild_ex+= "   --manifold-surface     ";
        cmd =
                ("../fTetWild/build/FloatTetwild_bin -i" + (input_filename + "_tmp.obj"))+ tetwild_ex+
                (" && mv " + new_name + "__sf.obj " + input_filename + "_final_result.obj");
    }
    else{
        for(int i=0;i<4;i++)new_name.pop_back();
        cmd = ("../TetWild/build/Tetwild " + (input_filename + "_tmp.obj") ) + tetwild_ex +
              (" && mv "+new_name+"__sf.obj " + input_filename+"_final_result.obj" );

    }
    exit(0);

    //system(("mv "+new_name+"__sf.obj " + input_filename+"_result.obj" ).c_str());
    //Remeshing().run((input_filename + "_result.obj").c_str());
    cout << cmd << endl;
    system(cmd.c_str());
    auto end2_clock = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff2 = end2_clock - end_clock;
    std::cout << input_filename<< "\t" << diff.count() <<"\t"<<diff2.count()<<"\t";
    auto mesh00 = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile(input_filename + "_final_result.obj"));
    cout <<"final result is saved in "<< input_filename + "_final_result.obj"<<endl;
    cout << mesh00->VertexSize() <<"\t"<< mesh00->FaceSize() << "\t"<< mesh->VertexSize()<<"\t"<<mesh->FaceSize()<<endl;
    stringstream time_use;
    FILE *file_ans = fopen(  (input_filename + "_time.txt").c_str(), "w");
    time_use  << input_filename<< "\t" << diff.count() <<"\t"<<diff2.count()<<"\t";
    time_use << mesh00->VertexSize() <<"\t"<< mesh00->FaceSize() << "\t"<< mesh->VertexSize()<<"\t"<<mesh->FaceSize()<<endl;
    fputs(time_use.str().c_str(),file_ans);
    fclose(file_ans);
    return 0;
}
// 1 2 3 4 5 6
// 5 1 2 3 7 8
//



