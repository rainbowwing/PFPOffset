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
#include "data_analyze.h"

int check_resolution = 3;
using namespace std;
//#define DEBUG

int main(int argc, char* argv[]) {
//    data_analyze();
//    exit(0);
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
#ifndef DEBUG
    set_start();
#endif
    input_filename = input_filename.substr(0,input_filename.size()-4);
    input_filename += string("_") + (running_mode==1 ? "offset_outward_"+ to_string(default_ratio) : "offset_inward_"+ to_string(default_ratio));
    FILE *file8 = fopen( (input_filename + "_result.obj").c_str(), "w");
    //FILE *file88 = fopen( (input_filename + "_result_have_degenerate.obj").c_str(), "w");
    FILE *file99 = fopen( (input_filename + "_resultdebug99.obj").c_str(), "w");
    FILE *file77 = fopen( (input_filename + "_result_have_without_build.obj").c_str(), "w");
#ifdef DEBUG
    FILE *file11 = fopen( (input_filename + "_grid.obj").c_str(), "w");
    FILE *file6 = fopen( (input_filename + "_tmp.obj").c_str(), "w");
    FILE *file7 = fopen( (input_filename + "_no_manifold.obj").c_str(), "w");

    FILE *file13 = fopen( (input_filename + "_moveedge.obj").c_str(), "w");
    FILE *file15 = fopen( (input_filename + "_extend_vertex.obj").c_str(), "w");

    FILE *file14 = fopen( (input_filename + "_coveragefield.obj").c_str(), "w");
    FILE *file24 = fopen( (input_filename + "_coveragefield24.obj").c_str(), "w");
//    FILE *file44 = fopen( (input_filename + "_check_resolution_delete.obj").c_str(), "w");
//    FILE *file54 = fopen( (input_filename + "_check_resolution_reserver.obj").c_str(), "w");
    FILE *file34 = fopen( (input_filename + "_cutting_segment.obj").c_str(), "w");
    FILE *file676 = fopen( (input_filename + "_fffhole.obj").c_str(), "w");
#endif
    mesh->initBBox();
   // mesh->build_fast();
   // mesh->build_fast();-f ../data_variable/914678.obj -i=2 -t=10 -E=0.0005 -r=0.003000
    auto start_clock = std::chrono::high_resolution_clock::now();
    cout <<"mesh->build_fast() succ" << endl;
    //double default_move_dist = 0.05;
    double x_len = (mesh->BBoxMax - mesh->BBoxMin).x();
    double y_len = (mesh->BBoxMax - mesh->BBoxMin).y();
    double z_len = (mesh->BBoxMax - mesh->BBoxMin).z();
    double min_bbox_len = min(min(x_len,y_len),z_len);

    {
        double sum = 0;
        for(int i=0;i<mesh->FaceSize();i++) {
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
            //default_move = min(min(x_len, y_len), z_len)*1e-3;
            default_move = sqrt(x_len*x_len+y_len*y_len+z_len*z_len)*default_ratio;
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
    //exit(0);
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
            double rate = (tmp.x() - mesh->BBoxMin.x())/(mesh->BBoxMax.x() - mesh->BBoxMin.x());


//            if(rate<0.25){
//                rate = 0;
//            }
//            else if(rate<0.5){
//                rate = 0.25;
//            }
//            else if(rate<0.75){
//                rate = 0.50;
//            }
//            else {
//                rate = 0.75;
//            }
          //  mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist = (rate*1.9+0.1)*default_move;

            cout << i <<" is: "<< mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist << endl;
        }
        min_near_limit = min(min_near_limit,mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist);

        cout <<"facei "<<  mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist<<" "<<mesh->FastNeighborFhOfFace_[i].size()<< endl;

       // mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist=0.133754;
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

//    std::vector <std::shared_ptr<std::thread> > build_neighbor(thread_num);
//    for(int i=0;i<thread_num;i++) {
//        //for (int i = 50; i < 51; i++) {
//        build_neighbor[i] = make_shared<std::thread>([&](int now_id) {
//            std::list<K2::Triangle_3> tri_list;
//            std::unordered_map<unsigned long long ,int> mp;
//            for (int i = 0; i < mesh->FaceSize(); i++) {
//                auto v0 = mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(0)];
//                auto v1 = mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(1)];
//                auto v2 = mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(2)];
//                K2::Triangle_3 tri(iGameVertex_to_Point_K2(v0),
//                                   iGameVertex_to_Point_K2(v1),
//                                   iGameVertex_to_Point_K2(v2)
//                );
//                mp[tri.id()] = i;
//                tri_list.push_back(tri);
//            }
//            Tree aabb(tri_list.begin(),tri_list.end());
//            for (int i = 0; i < mesh->VertexSize(); i++) {
//                if (i % thread_num != now_id)continue;
//                if (i % 20 == 0)
//                    cout << "build near: " << i << "/"<<mesh->VertexSize()<< endl;
//                CGAL::Epeck::FT r(min_bbox_len/1000);
//                std::list<Primitive> intersected_primitives;
//                std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
//                K2::Point_3 this_vertex = iGameVertex_to_Point_K2(mesh->fast_iGameVertex[i]);
////                cout << coverage_field_list[i].bbox_min <<" "<< coverage_field_list[i].bbox_max<<endl;
////                CGAL::Epeck::FT xx = coverage_field_list[i].bbox_min.x();
////                cout << xx - CGAL::Epeck::FT(0.05) << endl;
//                K2::Iso_cuboid_3 bbox2(
//                        this_vertex.x() - r,
//                        this_vertex.y() - r,
//                        this_vertex.z() - r,
//                        this_vertex.x() + r,
//                        this_vertex.y() + r,
//                        this_vertex.z() + r
//                );
//                aabb.all_intersections(bbox2,std::back_inserter(intersections));
//                cout << "intersectionssize:"<< intersections.size() << endl;
//                for(auto item : intersections) {
////                    cout <<"v " <<item.second->vertex(0) << endl;
////                    cout <<"v " <<item.second->vertex(1) << endl;
////                    cout <<"v " <<item.second->vertex(2) << endl;
//                    auto iter = mp.find(item.second->id());
//                    if(iter == mp.end()) exit(0);// 异常错误
//                    mesh->FastNeighborFhOfVertex_[i].insert(MeshKernel::iGameFaceHandle(iter->second));
//                }
//            }
//        }, i);
//    }
//
//    for(int i=0;i<thread_num;i++)
//        build_neighbor[i]->join();
    for (int i = 0; i < mesh->VertexSize(); i++) {
        cout << i <<" "<< mesh->FastNeighborFhOfVertex_[i].size() << endl;
    }

   // exit(0);
    for(int i =0;i< mesh->FastNeighborFhOfVertex_.size();i++){
        for(auto j : mesh->FastNeighborFhOfVertex_[i]) {
            std::set<int>se;
            se.insert(mesh->fast_iGameFace[j].vh(0));
            se.insert(mesh->fast_iGameFace[j].vh(1));
            se.insert(mesh->fast_iGameFace[j].vh(2));
            for(auto k: mesh->FastNeighborFhOfVertex_[i]) {
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
#ifdef DEBUG
    int cnte = 1;
    for (int i = 0; i < mesh->VertexSize(); i++) {
        for(int j = 0;j < field_move_vertices[i].size();j++){
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
#endif
    auto dp_clock = std::chrono::high_resolution_clock::now();
    for(int i=0;i<mesh->VertexSize();i++){
        origin_mesh_vertices[i] = K2::Point_3 (mesh->fast_iGameVertex[i].x(),
                                               mesh->fast_iGameVertex[i].y(),
                                               mesh->fast_iGameVertex[i].z()
                                               );

    }
   // exit(0);
    field_move_vertices_ex.resize(field_move_vertices.size());
    merge_initial();
//    for (int i = 0; i < mesh->VertexSize(); i++) {
//        bool flag = false;
//        if (field_move_vertices[i].size() >= 2) {
//            vector<MeshKernel::iGameFaceHandle> neighbor_list;
//            for(auto j : mesh->FastNeighborFhOfVertex_[i]){
//                neighbor_list.push_back(j);
//            }
//            if(neighbor_list.empty()) {
//                continue;
//            }
//            for(auto j : neighbor_list){
//                bool positive_side = false;
//                bool negative_side = false;
//                K2::Triangle_3 this_facet(origin_mesh_vertices[mesh->fast_iGameFace[j].vh(0)],
//                                          origin_mesh_vertices[mesh->fast_iGameFace[j].vh(1)],
//                                          origin_mesh_vertices[mesh->fast_iGameFace[j].vh(2)]);
//
//                for(auto k : field_move_vertices[i]) {
//                    positive_side |= this_facet.supporting_plane().has_on_positive_side(k);
//                    negative_side |= this_facet.supporting_plane().has_on_negative_side(k);
//                }
//                if(positive_side && negative_side) {
//                    K2::Point_3 center = centroid(this_facet);
//                    //flag = true;
//                    K2::Vector_3 vec = (origin_mesh_vertices[i] - center)
//                            / CGAL::approximate_sqrt((origin_mesh_vertices[i] - center).squared_length());
//
//                    K2::Point_3 extend_point = origin_mesh_vertices[i] + vec * mesh->fast_iGameFace[j].move_dist / 2;
//                    field_move_vertices_ex[i].push_back(extend_point);
//                }
//            }
//        }
//        for(auto k : field_move_vertices_ex[i]){
//            field_move_vertices[i].push_back(k);
//        }
//
//
//        if(flag) {
//            cout << i << "need extend "<< endl;
//            fprintf(file15, "v %lf %lf %lf\n",CGAL::to_double(origin_mesh_vertices[i].x())+start_x,
//                    CGAL::to_double(origin_mesh_vertices[i].y())+start_y,
//                    CGAL::to_double(origin_mesh_vertices[i].z())+start_z
//                    );
//        }
//        else {
//            cout << i << " not "<< endl;
//        }
//    }
   //exit(0);

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
        K2::Triangle_3 tri(v0,v1,v2);
        if(!tri.is_degenerate())
            origin_face_list.push_back(tri);
    }

    Tree origin_face_tree(origin_face_list.begin(),origin_face_list.end());

    for(int i=0;i<mesh->FaceSize();i++){
        //coverage_field_list.push_back(CoverageField(MeshKernel::iGameFaceHandle(i)));
        coverage_field_list.push_back(CoverageField(MeshKernel::iGameFaceHandle(i)));
    }

    auto polygon_clock = std::chrono::high_resolution_clock::now();


    int pre_cnt = 0;
    for(int i=0;i<mesh->FaceSize();i++){
//        if(coverage_field_list[i].bound_face_vertex_inexact.size() == 7 ) {
        for (int j = 0; j < coverage_field_list[i].bound_face_vertex_exact.size(); j++) {
#ifdef DEBUG
            fprintf(file14,"v %lf %lf %lf\n",CGAL::to_double(coverage_field_list[i].bound_face_vertex_exact[j].x())+start_x,
                    CGAL::to_double(coverage_field_list[i].bound_face_vertex_exact[j].y())+start_y,
                    CGAL::to_double(coverage_field_list[i].bound_face_vertex_exact[j].z())+start_z);
#endif
        }
        for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++){
#ifdef DEBUG
            fprintf(file14,"f %d %d %d\n",coverage_field_list[i].bound_face_id[j][0]+1+pre_cnt,
                    coverage_field_list[i].bound_face_id[j][1]+1+pre_cnt,
                    coverage_field_list[i].bound_face_id[j][2]+1+pre_cnt);
#endif
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
#ifdef DEBUG
    fclose(file14);
#endif
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
                            if(bs == CGAL::ON_BOUNDED_SIDE) {
                                inner_flag = true;
                                break;
                            }
                            else if(bs == CGAL::ON_BOUNDARY) {
                                if(tri_this.supporting_plane().oriented_side(coverage_field_list[i].center) !=
                                   tri_this.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) &&
                                   tri_this.supporting_plane().oriented_side(coverage_field_list[i].center) +
                                   tri_this.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) == 0 ) {
                                    inner_flag = true;
                                    break;
                                }
                            }
                        }

                        if(!inner_flag) {
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
                            for(int k=0;k<coverage_field_list[i].bound_face_sampling_point[j].size();k++) {
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

                    function<bool(K2::Point_3,K2::Triangle_3)>check_inner_vertex = [&](K2::Point_3 check_point,K2::Triangle_3 tri) {
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
//                        for(auto neighbor_id: coverage_intersection[i]) {
//                            if (!coverage_field_list[neighbor_id].useful)continue;
//                            if(neighbor_id == i)continue;
//                            auto side = coverage_field_list[neighbor_id].bounded_side(check_point);
//                            if(side == CGAL::ON_BOUNDED_SIDE) {
//                                return true;
//                            }
//                            if(side == CGAL::ON_BOUNDARY) {
//                                if( tri.supporting_plane().oriented_side(coverage_field_list[i].center) !=
//                                    tri.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) &&
//                                    tri.supporting_plane().oriented_side(coverage_field_list[i].center) +
//                                    tri.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) ==0){
//                                    return true;
//                                }
//                            }
//                        }
//                        return false;


                            //下面写这段是新方法//todo 改造成爆掉一个点就炸了 ？ 或者有一个计数器
                        map<int,set<K2::Point_3> >field_cnt;
                        set<int>field_error;

                           // cout << "v "<< check_point<<endl;
                        K2::Vector_3 normal = -1 * tri.supporting_plane().orthogonal_vector();
                        K2::Ray_3 ray(check_point,normal);
                        std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;

                        aabb_local.all_intersections(ray,std::back_inserter(intersections));

                        for(auto item : intersections) { //todo 这里改成批量插入
                            int id = face_to_field.find(item.second->id())->second;
                            if (const K2::Point_3 *p = boost::get<K2::Point_3>(&(item.first))) {
                                if(*p == check_point)continue;
                                if(field_cnt[id].count(*p)) {
                                    field_error.insert(id);
                                }
                                else
                                    field_cnt[id].insert(*p);
                            }
                            else {
                                field_error.insert(id);
                            }
                        }

                        for(const auto& it : field_cnt) {
                            if(!field_error.count(it.first) && it.second.size() % 2 == 1) {
                                return true;
                            }
                        }

                        for(auto neighbor_id: field_error) {
                            if (!coverage_field_list[neighbor_id].useful)continue;
                            if(neighbor_id == i)continue;
                            auto side = coverage_field_list[neighbor_id].bounded_side(check_point);
                            if(side == CGAL::ON_BOUNDED_SIDE){
                                return true;
                            }
                            if(side == CGAL::ON_BOUNDARY){
                                if( tri.supporting_plane().oriented_side(coverage_field_list[i].center) !=
                                    tri.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) &&
                                    tri.supporting_plane().oriented_side(coverage_field_list[i].center) +
                                    tri.supporting_plane().oriented_side(coverage_field_list[neighbor_id].center) ==0){
                                    return true;
                                }
                            }
                        }
                        return false;
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
                        if(check_inner_vertex_all(coverage_field_list[i].bound_face_sampling_point[j],tri_this)){ // 这里是决定到底全部删除，还是全部留下的逻辑所在之处
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

//                        for(int k=0;k<coverage_field_list[i].bound_face_sampling_point[j].size();k++) {
//                            if(coverage_field_list[i].bound_face_sampling_point_state[j][k] != -1){
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
//        FILE *file44 = fopen( (input_filename + "_check_resolution_delete.obj").c_str(), "w");
//        FILE *file54 = fopen( (input_filename + "_check_resolution_reserver.obj").c_str(), "w");
#ifdef DEBUG
        ofstream fsip((input_filename + "_check_resolution_delete.obj").c_str());
        ofstream fsop((input_filename + "_check_resolution_reserver.obj").c_str());
        for(int i=0;i<mesh->FaceSize();i++) {
            for (int j = 0; j < coverage_field_list[i].bound_face_id.size(); j++) {
                for(int k=0;k<coverage_field_list[i].bound_face_sampling_point[j].size();k++){
                    cout << coverage_field_list[i].bound_face_sampling_point_state[j][k] << endl;
                    if(coverage_field_list[i].bound_face_sampling_point_state[j][k]==-1 || coverage_field_list[i].bound_face_useful[j]!=1 ) {
                        fsip <<"v "<< coverage_field_list[i].bound_face_sampling_point[j][k].x()+start_x << " "
                             << coverage_field_list[i].bound_face_sampling_point[j][k].y()+start_y << " "
                             << coverage_field_list[i].bound_face_sampling_point[j][k].z()+start_z << endl;
                    }
                    else{
                        fsop  <<"v "<< coverage_field_list[i].bound_face_sampling_point[j][k].x()+start_x << " "
                             << coverage_field_list[i].bound_face_sampling_point[j][k].y()+start_y << " "
                             << coverage_field_list[i].bound_face_sampling_point[j][k].z()+start_z << endl;
                    }
                }
            }
        }
#endif

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
    auto speed_clock = std::chrono::high_resolution_clock::now();
    ofstream fsd(input_filename+"_find_delete.obj");
    ofstream fsl(input_filename+"_find_last.obj");
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
//    for (int face_id = 0; face_id < fsize; face_id++) {
//        MeshKernel::iGameVertex min_point(coverage_field_list[face_id].x_min,
//                                          coverage_field_list[face_id].y_min,
//                                          coverage_field_list[face_id].z_min
//        );
//        MeshKernel::iGameVertex max_point(coverage_field_list[face_id].x_max,
//                                          coverage_field_list[face_id].y_max,
//                                          coverage_field_list[face_id].z_max
//        );
//        cout << face_id <<": "<<mesh->fast_iGameFace[face_id].move_dist << endl;
//    }
    //exit(0);
    cout <<"bfs start \n" << endl;
    for (int face_id = 0; face_id < fsize; face_id++) {
        if (face_id % 1000 == 0)
            printf("%d/%d\n", face_id, fsize);
        if(!coverage_field_list[face_id].useful) continue;
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
    //exit(0);
    if(frame_grid_mp.size() > 5000){
        cout <<"len bug"<<": "<<frame_grid_mp.size()<< endl;
        exit(0);
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
                list<K2::Triangle_3> origin_face_list;
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
                            if(coverage_field_list[field_id].bound_face_useful[k] == 1) {
                                each_grid->second.field_face_though_list[field_id].push_back(k);
                                each_grid->second.face_hash_id_map[this_tri.id()] = {field_id,k};
                                each_grid->second.build_aabb_tree_triangle_list.push_back(this_tri);
                            }
                        }// 而是应该在这里直接插入，在这里直接插入可以形成较好的

                    }
                    K2::Triangle_3 origin_tri(origin_mesh_vertices[mesh->fast_iGameFace[field_id].vh(0)],
                    origin_mesh_vertices[mesh->fast_iGameFace[field_id].vh(1)],
                    origin_mesh_vertices[mesh->fast_iGameFace[field_id].vh(2)]);
                    origin_face_list.push_back(origin_tri);
                }
                origin_face_list.insert(origin_face_list.end(), each_grid->second.build_aabb_tree_triangle_list.begin(),
                                        each_grid->second.build_aabb_tree_triangle_list.end());
                Tree this_grid_aabb_tree( origin_face_list.begin(),
                                          origin_face_list.end()); // 这里插入新的面


//                Tree this_grid_aabb_tree( each_grid->second.build_aabb_tree_triangle_list.begin(),
//                                          each_grid->second.build_aabb_tree_triangle_list.end()); // 这里插入新的面


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
                                    coverage_field_list[field_id].bound_face_cutting_segment[bound_face_id].push_back(*s);
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
                                if(iter->second.first != field_id ) {
                                    for(int j=0;j<v->size();j++ ) {
                                        coverage_field_list[field_id].bound_face_cutting_segment[bound_face_id].push_back(K2::Segment_3(v->at(j),v->at((j+1)%v->size())));
                                    }
                                }
                            }
                            else if(const K2::Triangle_3 * tri = boost::get<K2::Triangle_3 >(&(item.first)) ){
                                map<size_t,pair<int,int> >::iterator iter = each_grid->second.face_hash_id_map.find(item.second->id());
                                if(iter->second.first != field_id ) {
                                    coverage_field_list[field_id].bound_face_cutting_segment[bound_face_id].push_back(K2::Segment_3(tri->vertex(0),tri->vertex(1)));
                                    coverage_field_list[field_id].bound_face_cutting_segment[bound_face_id].push_back(K2::Segment_3(tri->vertex(1),tri->vertex(2)));
                                    coverage_field_list[field_id].bound_face_cutting_segment[bound_face_id].push_back(K2::Segment_3(tri->vertex(2),tri->vertex(0)));
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
#ifdef DEBUG
    for (int face_id = 0; face_id < fsize; face_id++) {
        if(!coverage_field_list[face_id].useful)continue;
        for(int j=0;j<coverage_field_list[face_id].bound_face_id.size();j++){
            if(coverage_field_list[face_id].bound_face_useful[j] == 1)
            for(int k=0;k<coverage_field_list[face_id].bound_face_cutting_segment[j].size();k++){
                fprintf(file34,"v %lf %lf %lf\n",CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(0).x())+start_x,
                        CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(0).y()+start_y),
                        CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(0).z()+start_z)
                );
                fprintf(file34,"v %lf %lf %lf\n",CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(1).x()+start_x),
                        CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(1).y()+start_y),
                        CGAL::to_double(coverage_field_list[face_id].bound_face_cutting_segment[j][k].vertex(1).z()+start_z)
                );
                fprintf(file34,"l %d %d\n",ll+1,ll+2);
                ll+=2;
            }
            //fprintf(file34)
        }
    }
    fclose(file34);
#endif
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
#ifdef DEBUG
                fprintf(file11, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
                fprintf(file11, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
                fprintf(file11, "l %d %d\n", f3_id, f3_id + 1);
                f3_id += 2;
#endif
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
                global_vertex_id_sum += coverage_field_list[field_id].renumber_bound_face_vertex.size();
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
    //freopen("./out.txt","w",stdout);
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
                //cout << "thread" << now_id <<" st select"<<endl;
                for (int i = 0; i < each_grid->second.field_list.size(); i++) {
                    int field_id = each_grid->second.field_list[i];
//                    bool useful = false;
//                    if (coverage_field_list[field_id].in_field(midpoint(K2::Segment_3(ps[0], ps[7])))) {
//                        useful = true;
//                    }
//                    if (!useful) {
//                        if (vertex_in_frame(coverage_field_list[field_id].center)) {
//                            useful = true;
//                        }
//                    }
//                    if (!useful)
//                        useful = CGAL::Polygon_mesh_processing::do_intersect(*coverage_field_list[field_id].poly,
//                                                                             frame_poly);
                    if (1) {
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
                //cout << "thread" << now_id <<" st build aabb"<<endl;
                Tree aabb_tree(l.begin(), l.end());
                //cout << "thread" << now_id <<" st this thread"<<endl;
                for (int i = 0; i < each_grid->second.global_face_list.size(); i++) {
                    //cout << "thread" << now_id <<" st i"<<endl;
                    int global_face_id = each_grid->second.global_face_list[i];
                    //cout << "global_face_id" <<global_face_id << endl;
                    if (global_face_list[global_face_id].useful !=1 ) {
                        //cout << "thread" << now_id <<" continue i"<<endl;
                        continue;
                    }
//                    cout << "global_vertex_list check "<< global_vertex_list.size() <<" "
//                    <<global_face_list[global_face_id].idx0 <<" "<<global_face_list[global_face_id].idx1 <<" "
//                    <<global_face_list[global_face_id].idx2<<endl;
                    K2::Point_3 v0 = global_vertex_list[global_face_list[global_face_id].idx0];
                    //cout << "v0:"<<v0 << endl;
                    K2::Point_3 v1 = global_vertex_list[global_face_list[global_face_id].idx1];
                    //cout << "v1:"<<v1 << endl;
                    K2::Point_3 v2 = global_vertex_list[global_face_list[global_face_id].idx2];
                    //cout << "v2:"<<v2 << endl;
                    // 这里有没有退化，
                    /*
thread7 done 2.2
thread7 done 2.4
thread7 done 3
thread7 end i
thread7 st i

                     */
//                    cout << K2::Triangle_3(v0, v1, v2).is_degenerate() << endl;

                    K2::Vector_3 ray_vec = K2::Triangle_3(v0, v1, v2).supporting_plane().orthogonal_vector();
//                    cout << "ray_vec:"<<ray_vec << endl;
//                    cout << "centerbef :"<<global_face_list[global_face_id].center << endl;
                    K2::Ray_3 ray(global_face_list[global_face_id].center, -ray_vec);
//                    cout <<"center : " <<CGAL::to_double(global_face_list[global_face_id].center.x()) <<" "
//                    << CGAL::to_double(global_face_list[global_face_id].center.y())<<" "
//                    << CGAL::to_double(global_face_list[global_face_id].center.z()) << endl;
//                    cout <<"ray_vec : " <<CGAL::to_double(ray_vec.x()) <<" "
//                         << CGAL::to_double(ray_vec.y())<<" "
//                         << CGAL::to_double(ray_vec.z()) << endl;


                    std::list<Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections;
                    aabb_tree.all_intersections(ray, std::back_inserter(intersections));
                    //cout <<"aabbend"<<endl;
                    for (auto item: intersections) {
                        //cout << "thread" << now_id <<" done 0"<<endl;
                        pair<int, int> belong = triangle_mapping[item.second->id()];
                        int which_field = belong.first;
                        int which_id = belong.second;
                        //cout << "thread" << now_id <<" done 1"<<endl;
                        if (which_field == global_face_list[global_face_id].field_id)continue;
                        //cout << "thread" << now_id <<" done 2"<<endl;
                        if (const K2::Point_3 *p = boost::get<K2::Point_3>(&(item.first))) {
                            if (*p != global_face_list[global_face_id].center) {
                                //cout << "thread" << now_id <<" done 2.2"<<endl;
                                global_face_list[global_face_id].ray_detect_map[which_field].emplace_back(*p, which_id);
                                //cout << "thread" << now_id <<" done 2.4"<<endl;
                            }
//                            else
//                                global_face_list[global_face_id].special_face_id.insert(which_field);
                        } else {
                            // cout << "thread" << now_id <<" done 2.6"<<endl;
                            global_face_list[global_face_id].special_field_id.insert(which_field);
                        }
                        // cout << "thread" << now_id <<" done 3"<<endl;
                    }
                    // cout << "thread" << now_id <<" end i"<<endl;
                }
                // cout << "thread" << now_id <<" end this thread"<<endl;
            }

        },i);
    }
    for(int i=0;i<thread_num;i++)
        face_generate_ray_detect_thread_pool[i]->join();
    //exit(0); // 现在这里加exit 然后注释掉上半部分看看会不会出问题，不会就测测下面的是不是求交出问题了
    auto ray_clock = std::chrono::high_resolution_clock::now();
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

                    if(origin_face_tree.squared_distance(v0) <= CGAL::Epeck::FT(0) &&
                            origin_face_tree.squared_distance(v1) <= CGAL::Epeck::FT(0) &&
                            origin_face_tree.squared_distance(v2) <= CGAL::Epeck::FT(0)
                    ){
                        //cout <<"-300:occur: " <<origin_face_tree.squared_distance(centroid(K2::Triangle_3(v0,v1,v2))) <<" "<< min_near_limit << " "<<default_move<<endl;
                        global_face_list[i].useful = -300;
                        cout << "global_face_list[i].useful = -300" << endl;
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
    CGAL::Epeck::FT face_avg(0);
    vector<FinalFace> f_list;
    set<vector<int> >final_face_set;
    int xx77 = 1;
    for (int i = 0; i < global_face_list.size(); i++) {
//        if(is_same_triangle(global_vertex_list[global_face_list[i].idx0],
//                            global_vertex_list[global_face_list[i].idx1],
//                            global_vertex_list[global_face_list[i].idx2],
//                            K2::Point_3 (-2.322543,-23.230715,3.539219),
//                            K2::Point_3 (-2.053821,-23.227179,3.923453),
//                            K2::Point_3 (-2.321043,-20.235704,3.538151),
//                            0.3
//                            ) ){
//            cout << "find you " << global_face_list[i].useful<< endl;
//        }




        //cout <<"global_face_list["<<i<<"].useful:" <<global_face_list[i].useful << endl;
        if(global_face_list[i].useful == 1 ) {
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
            face_avg += (K2::Triangle_3(global_vertex_list[global_face_list[i].idx0], global_vertex_list[global_face_list[i].idx1], global_vertex_list[global_face_list[i].idx2])).squared_area();
            //if(!final_face_set.count(f))
            f_list.push_back(f);

            fprintf(file77,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[global_face_list[i].idx0].x()),
                    CGAL::to_double(global_vertex_list[global_face_list[i].idx0].y()),
                    CGAL::to_double(global_vertex_list[global_face_list[i].idx0].z())
            );
            fprintf(file77,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[global_face_list[i].idx1].x()),
                    CGAL::to_double(global_vertex_list[global_face_list[i].idx1].y()),
                    CGAL::to_double(global_vertex_list[global_face_list[i].idx1].z())
            );
            fprintf(file77,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[global_face_list[i].idx2].x()),
                    CGAL::to_double(global_vertex_list[global_face_list[i].idx2].y()),
                    CGAL::to_double(global_vertex_list[global_face_list[i].idx2].z())
            );
            fprintf(file77,"f %d %d %d\n",xx77,
                    xx77+1,
                    xx77+2
            );
            xx77+=3;

        }
    }
    if(f_list.size())
        face_avg /= (1.0*f_list.size());







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
    for(int i=0;i<final_vertex_useful.size();i++) {
#ifdef DEBUG
        fprintf(file6,"v %lf %lf %lf\n",CGAL::to_double(global_vertex_list[final_vertex_useful[i]].x()+start_x),
                CGAL::to_double(global_vertex_list[final_vertex_useful[i]].y()+start_y),
                CGAL::to_double(global_vertex_list[final_vertex_useful[i]].z()+start_z));
#endif
        final_vertex_useful_point.emplace_back(global_vertex_list[final_vertex_useful[i]].x()/*+CGAL::Epeck::FT(start_x)*/,
                                               global_vertex_list[final_vertex_useful[i]].y()/*+CGAL::Epeck::FT(start_y)*/,
                                               global_vertex_list[final_vertex_useful[i]].z()/*+CGAL::Epeck::FT(start_z)*/
        );

    }




    for(int i=0;i<final_vertex_useful_point.size();i++){
        fprintf(file99,"v %lf %lf %lf\n",CGAL::to_double(final_vertex_useful_point[i].x())+(start_x),
                CGAL::to_double(final_vertex_useful_point[i].y())+(start_y),
                CGAL::to_double(final_vertex_useful_point[i].z())+(start_z)
        );
//        fprintf(file88,"v %lf %lf %lf\n",CGAL::to_double(final_vertex_useful_point[i].x())+(start_x),
//                CGAL::to_double(final_vertex_useful_point[i].y())+(start_y),
//                CGAL::to_double(final_vertex_useful_point[i].z())+(start_z)
//        );
    }

//





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
            } ) {
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
                //f_list[i.second[j]].flag = false;
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

    for (int i = 0; i < f_list.size(); i++) {
        if(f_list[i].flag == 1 ) {
            if(running_mode == 1) {
                fprintf(file99, "f %d %d %d\n", is_vertex_useful[f_list[i].f0] + 1,
                        is_vertex_useful[f_list[i].f2] + 1,
                        is_vertex_useful[f_list[i].f1] + 1);
            }
            else{
                fprintf(file99, "f %d %d %d\n", is_vertex_useful[f_list[i].f0] + 1,
                        is_vertex_useful[f_list[i].f1] + 1,
                        is_vertex_useful[f_list[i].f2] + 1);
            }
        }
    }


    vector<vector<int> >final_face_list;
    list<K2::Triangle_3 > t_list;


    for (int i = 0; i < f_list.size(); i++) {
        if(f_list[i].flag) {
            if(running_mode == 1) {

//                sum_avg_edge += sqrt(CGAL::to_double((CGAL::squared_distance(global_vertex_list[global_face_list[i].idx0] , global_vertex_list[global_face_list[i].idx1]))));
//                sum_avg_edge += sqrt(CGAL::to_double((CGAL::squared_distance(global_vertex_list[global_face_list[i].idx2] , global_vertex_list[global_face_list[i].idx1]))));
//                sum_avg_edge += sqrt(CGAL::to_double((CGAL::squared_distance(global_vertex_list[global_face_list[i].idx0] , global_vertex_list[global_face_list[i].idx2]))));
#ifdef DEBUG
                fprintf(file6, "f %d %d %d\n", is_vertex_useful[f_list[i].f0] + 1, is_vertex_useful[f_list[i].f2] + 1,
                        is_vertex_useful[f_list[i].f1] + 1);
#endif
                t_list.push_back(K2::Triangle_3(global_vertex_list[f_list[i].f0],
                                global_vertex_list[f_list[i].f1],
                                global_vertex_list[f_list[i].f2]));
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
                                                      is_vertex_useful[f_list[i].f2],
                                                      is_vertex_useful[f_list[i].f1],
                                                      1});

            }
            else {
#ifdef DEBUG
                fprintf(file6, "f %d %d %d\n", is_vertex_useful[f_list[i].f0] + 1, is_vertex_useful[f_list[i].f1] + 1,
                        is_vertex_useful[f_list[i].f2] + 1);
#endif

                t_list.push_back(K2::Triangle_3(global_vertex_list[f_list[i].f0],
                                                global_vertex_list[f_list[i].f1],
                                                global_vertex_list[f_list[i].f2]));


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

    for (auto i: origin_face_list) {
        t_list.push_back(i);
    }
    Tree all_aabb(t_list.begin(),t_list.end());




    //vector<K2::Point_3>final_vertex_useful_point_d = final_vertex_useful_point;
//    sort(final_vertex_useful_point_d.begin(),final_vertex_useful_point_d.end(),[&](const K2::Point_3 &a,const K2::Point_3 &b){
//        if(a.x() != b.x())
//            return a.x() < b.x();
//        else if(a.y() != b.y())
//            return a.y() < b.y();
//        else
//            return a.z() < b.z();
//    });
//    for(int i=0;i<final_vertex_useful_point_d.size()-1;i++){
//        if(final_vertex_useful_point_d[i] == final_vertex_useful_point_d[i+1]){
//            cout << "coincide it" << endl;
//        }
//    }
//
//    Tree dtree(t_list.begin(),t_list.end());
//    for(auto i : t_list){
//        cnt++;
//
//        std::list<Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
//        dtree.all_intersections(i, std::back_inserter(intersections));
//        for (auto item: intersections) {
//            if (const K2::Point_3 * p = boost::get<K2::Point_3>(&(item.first))) {
//                continue;
//            }
//            else if(const K2::Segment_3 * s = boost::get<K2::Segment_3>(&(item.first))){
//                if(segment_coincide_triangle(*s,i)) {
//
//                }
//                else{
//                    cout <<"no1"<<endl;
//                }
//            }
//            else if(const K2::Triangle_3 * tri = boost::get<K2::Triangle_3>(&(item.first))){
//                if(i.id() != item.second->id()) {
//                    cout <<"no3:"<<cnt<<endl;
//                }
//                else{
//                   // cout <<"yes3:"<<cnt<<endl;
//                }
//                continue;
//            }
//            else if(const std::vector<K2::Point_3> * v = boost::get<std::vector<K2::Point_3> > (&(item.first))){
//                cout <<"no2"<<endl;
//            }
//        }
//    }
//    cout <<"dtree success"<< endl;

    // exit(0);
    map<int,vector<int> > hole_line;
    vector<vector<int> > each_hole;
    queue<int>hole_q;
    int debug_hole_cnt = 1;
    Tree dtree(t_list.begin(),t_list.end());

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
#ifdef DEBUG
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
#endif
        }
    }
    // 一会抓出来反法向的面，检测所有边是否出现冲突，进行调整。

//    function<void(vector<int>)>

#ifdef DEBUG
    fprintf(file6,"#*********\n");
#endif
    for(int i=0;i<each_hole.size();i++) {
        if(each_hole[i].size() == 3) {
#ifdef DEBUG
            fprintf(file6, "f %d %d %d\n", is_vertex_useful[each_hole[i][0]] + 1, is_vertex_useful[each_hole[i][1]] + 1,
                    is_vertex_useful[each_hole[i][2]] + 1);
#endif
//            final_topo_check_mp[std::make_pair(is_vertex_useful[each_hole[i][0]],is_vertex_useful[each_hole[i][1]])] = final_face_list.size();
//            final_topo_check_mp[std::make_pair(is_vertex_useful[each_hole[i][1]],is_vertex_useful[each_hole[i][2]])] = final_face_list.size();
//            final_topo_check_mp[std::make_pair(is_vertex_useful[each_hole[i][2]],is_vertex_useful[each_hole[i][0]])] = final_face_list.size();


            K2::Triangle_3 t(final_vertex_useful_point[is_vertex_useful[each_hole[i][0]]],
                             final_vertex_useful_point[is_vertex_useful[each_hole[i][1]]],
                             final_vertex_useful_point[is_vertex_useful[each_hole[i][2]]]
            );
            if(t.squared_area() <= face_avg/1000000)
           // if(t.squared_area() <= 1000000)
                final_face_list.emplace_back(vector<int>{is_vertex_useful[each_hole[i][0]],
                                                         is_vertex_useful[each_hole[i][1]],
                                                         is_vertex_useful[each_hole[i][2]],
                                                         1});



            if( !(t.is_degenerate()) ){
                cout <<"not zero"<<" "<< K2::Triangle_3 (final_vertex_useful_point[is_vertex_useful[each_hole[i][0]]],
                                                         final_vertex_useful_point[is_vertex_useful[each_hole[i][1]]],
                                                         final_vertex_useful_point[is_vertex_useful[each_hole[i][2]]]
                ).squared_area() <<" "<<face_avg/1000000<< endl;
            }
            continue;
        }
        else{
            vector<K2::Point_3>position;
            for(int j=0;j<each_hole[i].size();j++) {
                position.push_back(global_vertex_list[each_hole[i][j]]);
                //origin_face_tree
            }
            vector<vector<int> > ret;
            fix_hole(each_hole[i],position,&all_aabb,ret);
            cout <<"ret.size():" <<ret.size() << endl;
            for(int j=0;j<ret.size();j++){
                if(ret[j].size() == 3){
                    K2::Triangle_3 t(final_vertex_useful_point[is_vertex_useful[ret[j][0]]],
                                     final_vertex_useful_point[is_vertex_useful[ret[j][1]]],
                                     final_vertex_useful_point[is_vertex_useful[ret[j][2]]]
                    );
                    //if(t.squared_area() <= 1000000)
                    if(t.squared_area() <= face_avg/1000000)
                        final_face_list.emplace_back(vector<int>{is_vertex_useful[ret[j][0]],
                                                                 is_vertex_useful[ret[j][1]],
                                                                 is_vertex_useful[ret[j][2]],
                                                                 1});

                    if( !(t).is_degenerate()) {
                        cout <<" not zero "<<" "<< K2::Triangle_3 (final_vertex_useful_point[is_vertex_useful[ret[j][0]]],
                                                                 final_vertex_useful_point[is_vertex_useful[ret[j][1]]],
                                                                 final_vertex_useful_point[is_vertex_useful[ret[j][2]]]

                        ).squared_area() <<" "<<face_avg/1000000<< endl;
                    }
                    else
                        cout <<" zero " << endl;

                    cout <<"pre_ret j"<< endl;
                }
                else{
                    cout <<"pre_ret_error"<<endl;
                }
            }
        }
    }







//    for(int times = 0; times < 10 ; times++) {
//        for(int i=0;i<final_face_list.size();i++){
//            if(final_face_list[i][3]) {
//                final_topo_check_mp[std::make_pair(final_face_list[i][0], final_face_list[i][1] ) ] = i;
//                final_topo_check_mp[std::make_pair(final_face_list[i][1], final_face_list[i][2] ) ] = i;
//                final_topo_check_mp[std::make_pair(final_face_list[i][2], final_face_list[i][0] ) ] = i;
//            }
//        }
//
//        bool final_flag = false;
//        for(auto iter : final_topo_check_mp) {
//            int v0 = iter.first.first;
//            int v1 = iter.first.second;
//            if(!final_topo_check_mp.count(std::make_pair(v1,v0))){
//                cout <<"lack of "<< v1 <<" "<< v0 << endl;
//                cout <<"deleta face" << iter.second<< endl;
//                final_face_list[iter.second][3] = 0;
//                final_flag = true;
//            }
//        }
//        if(!final_flag) {
//            cout <<"final ok" << endl;
//            break;
//        }
//        final_topo_check_mp.clear();
//        for(int i=0;i<final_face_list.size();i++){
//            if(final_face_list[i][3]) {
//                final_topo_check_mp[std::make_pair(final_face_list[i][0], final_face_list[i][1])] = i;
//                final_topo_check_mp[std::make_pair(final_face_list[i][1], final_face_list[i][2])] = i;
//                final_topo_check_mp[std::make_pair(final_face_list[i][2], final_face_list[i][0])] = i;
//            }
//        }
//        hole_line.clear();
//        each_hole.clear();
//        for(auto iter : final_topo_check_mp) {
//            int v0 = iter.first.first;
//            int v1 = iter.first.second;
//            if(!final_topo_check_mp.count(std::make_pair(v1,v0))){
//                hole_line[v1].push_back(v0);
//            }
//        }
//        disassemble_circle(hole_line,each_hole);
//        for(int i=0;i<each_hole.size();i++) {
//            if(each_hole[i].size() == 3){
//                final_face_list.emplace_back(vector<int>{each_hole[i][0],
//                                                         each_hole[i][1],
//                                                         each_hole[i][2],
//                                                         1});
//                if( !(K2::Triangle_3 (final_vertex_useful_point[each_hole[i][0]],
//                                      final_vertex_useful_point[each_hole[i][1]],
//                                      final_vertex_useful_point[each_hole[i][2]]
//                ).is_degenerate()) ){
//                    cout <<"not zero x"<<" "<< K2::Triangle_3 (final_vertex_useful_point[each_hole[i][0]],
//                                                             final_vertex_useful_point[each_hole[i][1]],
//                                                             final_vertex_useful_point[each_hole[i][2]]
//                    ).squared_area() << endl;
//                }
//                else
//                    cout <<"is zero"<<endl;
//                continue;
//            }
//            vector<K2::Point_3>position;
//            for(int j=0;j<each_hole[i].size();j++) {
//                position.push_back(final_vertex_useful_point[each_hole[i][j]]);
//                //origin_face_tree
//            }
//            vector<vector<int> > ret;
//            fix_hole(each_hole[i],position,&origin_face_tree,ret);
//            for(int j=0;j<ret.size();j++){
//                if(ret[j].size() == 3){
//                    final_face_list.emplace_back(vector<int>{ret[j][0],
//                                                             ret[j][1],
//                                                             ret[j][2],
//                                                         1});
//                    K2::Triangle_3 t (final_vertex_useful_point[ret[j][0]],
//                    final_vertex_useful_point[ret[j][1]],
//                    final_vertex_useful_point[ret[j][2]]);
//                    if( !(t).is_degenerate() ){
//                        cout <<"not zero y"<<" "<< t.squared_area() << endl;
//                    }
//                    else
//                        cout <<"is zero"<<endl;
//
//                    cout <<"ret j"<< endl;
//                }
//                else{
//                    K2::Vector_3 vec(0,0,0);
//                    for(int k=0;k<ret[j].size();k++) {
//                        vec += final_vertex_useful_point[ret[j][k]] - K2::Point_3 (0,0,0);
//                    }
//                    vec /= each_hole[i].size();
//                    final_vertex_useful_point.emplace_back(vec.x(),
//                                                           vec.y(),
//                                                           vec.z()
//                    );
//                    for(int k=0;k<ret[j].size();k++) {
//                        final_face_list.emplace_back(vector<int>{(int)final_vertex_useful_point.size()-1,
//                                                                 ret[j][k],
//                                                                 ret[j][(k+1)%ret[j].size()],
//                                                                 1});
//                    }
//
//                    cout <<"ret_error dual"<<endl;
//                }
//            }
//        }
//
//
//
//
//
////        for(int i=0;i<each_hole.size();i++) {
////            if(each_hole[i].size() == 3){
////                final_face_list.emplace_back(vector<int>{each_hole[i][0],
////                                                         each_hole[i][1],
////                                                         each_hole[i][2],
////                                                         1});
////                continue;
////            }
////            cout <<"add p "<< endl;
////            K2::Vector_3 vec(0,0,0);
////            for(int j=0;j<each_hole[i].size();j++) {
////                vec += final_vertex_useful_point[each_hole[i][j]] - K2::Point_3 (0,0,0);
////            }
////            vec /= each_hole[i].size();
////
////            final_vertex_useful_point.emplace_back(vec.x(),
////                                                   vec.y(),
////                                                   vec.z()
////            );
////            for(int j=0;j<each_hole[i].size();j++) {
////                final_face_list.emplace_back(vector<int>{(int)final_vertex_useful_point.size()-1,
////                                                         each_hole[i][j],
////                                                         each_hole[i][(j+1)%each_hole[i].size()],
////                                                         1});
////            }
////        }
//
//    }




    for(int i=0;i<final_vertex_useful_point.size();i++){
        fprintf(file8, "v %lf %lf %lf\n",
                CGAL::to_double(final_vertex_useful_point[i].x()+CGAL::Epeck::FT(start_x)),
                CGAL::to_double(final_vertex_useful_point[i].y()+CGAL::Epeck::FT(start_y)),
                CGAL::to_double(final_vertex_useful_point[i].z()+CGAL::Epeck::FT(start_z)));
    }
    for (int i = 0; i <  final_face_list.size(); ++i) {
        if(running_mode == 1) {
            if (final_face_list[i][3]) {
//                fprintf(file88, "f %d %d %d\n",
//                        final_face_list[i][0] + 1,
//                        final_face_list[i][2] + 1,
//                        final_face_list[i][1] + 1);
                fprintf(file8, "f %d %d %d\n",
                        final_face_list[i][0] + 1,
                        final_face_list[i][1] + 1,
                        final_face_list[i][2] + 1);
            }
        }
        else{
            if (final_face_list[i][3]) {
//                fprintf(file8, "f %d %d %d\n",
//                        final_face_list[i][0] + 1,
//                        final_face_list[i][1] + 1,
//                        final_face_list[i][2] + 1);
                fprintf(file8, "f %d %d %d\n",
                        final_face_list[i][0] + 1,
                        final_face_list[i][1] + 1,
                        final_face_list[i][2] + 1);
            }
        }
    }


    // 构建拓扑

//    vector<vector<int> >topological_vertex_neighbor(final_vertex_useful_point.size());
//    for(int i=0;i<final_face_list.size();i++){
//        if(final_face_list[i][3]) {
//            topological_vertex_neighbor[final_face_list[i][0]].push_back(i);
//            topological_vertex_neighbor[final_face_list[i][1]].push_back(i);
//            topological_vertex_neighbor[final_face_list[i][2]].push_back(i);
//        }
//    }
//    DSU final_dsu;
//    map<pair<int,int>, vector<int> > record_cut;
//
//    for(int i=0;i<final_face_list.size();i++){
//        if(final_face_list[i][3]) {
//            if(set<int>{final_face_list[i][0],
//                        final_face_list[i][1],
//                        final_face_list[i][2],
//            }.size() != 3) {
//                final_face_list[i][3] = 0;
//                continue;
//            }
//            if(K2::Triangle_3 (final_vertex_useful_point[final_face_list[i][0]],
//                               final_vertex_useful_point[final_face_list[i][1]],
//                               final_vertex_useful_point[final_face_list[i][2]]
//            ).is_degenerate()){
//                if (final_vertex_useful_point[final_face_list[i][0]] == final_vertex_useful_point[final_face_list[i][1]]
//                    && final_vertex_useful_point[final_face_list[i][2]] == final_vertex_useful_point[final_face_list[i][1]]){
//                    final_face_list[i][3] = 0;
//                    final_dsu.join(final_face_list[i][0],final_face_list[i][1]);
//                    final_dsu.join(final_face_list[i][1],final_face_list[i][2]);
//                    final_dsu.join(final_face_list[i][2],final_face_list[i][0]);
//                    cout <<"out 3" << endl;
//                }
//                else if (final_vertex_useful_point[final_face_list[i][0]] == final_vertex_useful_point[final_face_list[i][1]]){
//                    final_face_list[i][3] = 0;
//                    final_dsu.join(final_face_list[i][0],final_face_list[i][1]);
//                    cout <<"out 2" << endl;
//                }
//                else if (final_vertex_useful_point[final_face_list[i][1]] == final_vertex_useful_point[final_face_list[i][2]]){
//                    final_face_list[i][3] = 0;
//                    final_dsu.join(final_face_list[i][1],final_face_list[i][2]);
//                    cout <<"out 2" << endl;
//                }
//                else if (final_vertex_useful_point[final_face_list[i][2]] == final_vertex_useful_point[final_face_list[i][0]]){
//                    final_face_list[i][3] = 0;
//                    final_dsu.join(final_face_list[i][2],final_face_list[i][0]);
//                    cout <<"out 2" << endl;
//                }
//                else if (K2::Segment_3(final_vertex_useful_point[final_face_list[i][0]],final_vertex_useful_point[final_face_list[i][1]])
//                        .has_on(final_vertex_useful_point[final_face_list[i][2]])){
//                    final_face_list[i][3] = 0;
//                    record_cut[{final_face_list[i][0],final_face_list[i][1]}].push_back(final_face_list[i][2]);
//                    record_cut[{final_face_list[i][1],final_face_list[i][0]}].push_back(final_face_list[i][2]);
//                    cout <<"out 2" << endl;
//                }
//                else if (K2::Segment_3(final_vertex_useful_point[final_face_list[i][1]],final_vertex_useful_point[final_face_list[i][2]])
//                        .has_on(final_vertex_useful_point[final_face_list[i][0]])){
//                    final_face_list[i][3] = 0;
//                    record_cut[{final_face_list[i][1],final_face_list[i][2]}].push_back(final_face_list[i][0]);
//                    record_cut[{final_face_list[i][2],final_face_list[i][1]}].push_back(final_face_list[i][0]);
//                    cout <<"out 2" << endl;
//                }
//                else if (K2::Segment_3(final_vertex_useful_point[final_face_list[i][0]],final_vertex_useful_point[final_face_list[i][2]])
//                        .has_on(final_vertex_useful_point[final_face_list[i][1]])){
//                    final_face_list[i][3] = 0;
//                    record_cut[{final_face_list[i][0],final_face_list[i][2]}].push_back(final_face_list[i][1]);
//                    record_cut[{final_face_list[i][2],final_face_list[i][0]}].push_back(final_face_list[i][1]);
//                    cout <<"out 2" << endl;
//                }
//                else{
//                    cout <<"error of degenerate"<< endl;
//                    exit(0);
//                }
//            }
//        }
//    }
//    int final_face_list_size = final_face_list.size();
//    for(int i=0;i<final_face_list_size;i++){
//        if(final_face_list[i][3]) {
//            int old_id0 = final_face_list[i][0];
//            int old_id1 = final_face_list[i][1];
//            int old_id2 = final_face_list[i][2];
//            if( record_cut.count({old_id0,old_id1}) ||
//                record_cut.count({old_id1,old_id2}) ||
//                record_cut.count({old_id2,old_id0}) ){
//                fprintf(file77, "#is deg\n");
//                auto center = centroid(K2::Triangle_3 (final_vertex_useful_point[final_face_list[i][0]],
//                                                       final_vertex_useful_point[final_face_list[i][1]],
//                                                       final_vertex_useful_point[final_face_list[i][2]]));
//                vector<int> new_tri;
//                new_tri.push_back(old_id0);
//                if(record_cut.count({old_id0,old_id1})){
//                    auto v = record_cut[{old_id0,old_id1}];
//                    sort(v.begin(),v.end(),[&](int a,int b) {
//                        return CGAL::squared_distance(final_vertex_useful_point[a],final_vertex_useful_point[old_id0])
//                               < CGAL::squared_distance(final_vertex_useful_point[b],final_vertex_useful_point[old_id0]);
//                    });
//                    for(int j=0;j<v.size();j++){
//                        new_tri.push_back(v[j]);
//                    }
//                }
//
//                new_tri.push_back(old_id1);
//                if(record_cut.count({old_id1,old_id2})){
//                    auto v = record_cut[{old_id1,old_id2}];
//                    sort(v.begin(),v.end(),[&](int a,int b) {
//                        return CGAL::squared_distance(final_vertex_useful_point[a],final_vertex_useful_point[old_id1])
//                               < CGAL::squared_distance(final_vertex_useful_point[b],final_vertex_useful_point[old_id1]);
//                    });
//                    for(int j=0;j<v.size();j++){
//                        new_tri.push_back(v[j]);
//                    }
//                }
//
//                new_tri.push_back(old_id2);
//                if(record_cut.count({old_id2,old_id0})){
//                    auto v = record_cut[{old_id2,old_id0}];
//                    sort(v.begin(),v.end(),[&](int a,int b) {
//                        return CGAL::squared_distance(final_vertex_useful_point[a],final_vertex_useful_point[old_id2])
//                               < CGAL::squared_distance(final_vertex_useful_point[b],final_vertex_useful_point[old_id2]);
//                    });
//                    for(int j=0;j<v.size();j++){
//                        new_tri.push_back(v[j]);
//                    }
//                }
//                new_tri.push_back(old_id0);
//                final_face_list[i][3] = 0;
//                int s = final_vertex_useful_point.size();
//                final_vertex_useful_point.push_back(center);
//                for(int j=0;j<(int)new_tri.size()-1;j++) {
//                    final_face_list.push_back({s,new_tri[j],new_tri[j+1],1});
//                }
//            }
//        }
//    }
//    for(int i=0;i<final_face_list_size;i++) {
//        if (final_face_list[i][3]) {
//            final_face_list[i][0] = final_dsu.find_root(final_face_list[i][0]);
//            final_face_list[i][1] = final_dsu.find_root(final_face_list[i][1]);
//            final_face_list[i][2] = final_dsu.find_root(final_face_list[i][2]);
//        }
//    }
//
//    for(int i=0;i<final_face_list.size();i++) {
//        if (final_face_list[i][3]) {
//            if (set < int > {final_face_list[i][0],
//                             final_face_list[i][1],
//                             final_face_list[i][2],
//            }.size() != 3) {
//                final_face_list[i][3] = 0;
//                continue;
//            }
//            if (K2::Triangle_3(final_vertex_useful_point[final_face_list[i][0]],
//                               final_vertex_useful_point[final_face_list[i][1]],
//                               final_vertex_useful_point[final_face_list[i][2]]
//            ).is_degenerate()) {
//                cout << "is_degenerate" << endl;
//            }
//            fprintf(file8, "f %d %d %d\n",
//                    final_face_list[i][0]+1,
//                    final_face_list[i][1]+1,
//                    final_face_list[i][2]+1
//            );
//        }
//    }









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

#ifdef DEBUG
    fclose(file6);
#endif
    //sum_avg_edge /=(global_face_list.size()*3*20);

    auto end_clock = std::chrono::high_resolution_clock::now();


    string new_name = (input_filename + "_tmp.obj");
    string cmd;
    string tetwild_ex = "";
    if(tetwild_l>0){
        tetwild_ex += " -l "+to_string(tetwild_l) +" ";
    }
    if(tetwild_e>0){
        tetwild_ex += " -e "+to_string(tetwild_e) +" ";
    }

    std::chrono::duration<double> diff0 = dp_clock - start_clock;
    std::chrono::duration<double> diff1 = polygon_clock - dp_clock;
    std::chrono::duration<double> diff2 = speed_clock - polygon_clock;
    std::chrono::duration<double> diff3 = ray_clock - speed_clock;
    std::chrono::duration<double> diff4 = end_clock - ray_clock;
    FILE *file_ans = fopen(  (input_filename + "_time.txt").c_str(), "w");
    stringstream time_use;
    time_use << diff0.count() <<" "<< diff1.count() <<" "<<diff2.count() <<" "<<diff3.count() <<" "<<diff4.count() << endl;
    fputs(time_use.str().c_str(),file_ans);

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

    return 0;
}
/*
# 查看swap内存大小
free -h
# 关闭所有交换空间
swapoff -a
# dd命令创建一个swap文件,/dev/zero不要改动，swapfile可以改名
dd if=/dev/zero of=/tmp/swapfile bs=1024 count=2048000
# 修改swapfile文件权限，就是600，不然会warning
chmod 600 /tmp/swapfile
# 格式化为swap分区
mkswap /tmp/swapfile
# 挂载并激活swap分区
swapon /tmp/swapfile
# 设置扩展的swap分区为自动挂载，这一步不做，swapoff和swapon的-a对这个文件swap分区不起效
echo /tmp/swapfile swap swap defaults 0 0 >> /etc/fstab

 */


