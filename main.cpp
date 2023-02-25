#define CGAL_HAS_THREADS
//#define MODEIN
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
#include "remeshing.h"
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>
#include "grid.h"
#include "io.h"
#include "MeshBuilder.h"
#include "LocalMesh.h"
#include "GridVertex.h"
#include "CGAL_CDT.h"
#include "ApproximateField.h"
#include "sort_by_polar_order.h"
#include "do_quadratic_error_metric.h"
#include <gflags/gflags.h>
using namespace std;

shared_ptr <MeshKernel::SurfaceMesh> mesh;

shared_ptr<CGALPolygon>cgal_polygon;
int thread_num = 16;
vector<MeshKernel::iGameVertex> field_move_vertex;
vector<vector<MeshKernel::iGameVertex> > field_move_face;
vector<K2::Triangle_3> field_move_K2_triangle;

vector<ApproximateField>faces_approximate_field;
double merge_eps = 1e-5;


DEFINE_int32(m, 3, "This arg is a integer means result mode which value can be chose in 1,2,3. \n mode 1  without step 6;mode 2 is with building mesh;mode 3 is with building mesh and remeshing.");
DEFINE_double(d, 0.5, "This arg is a double means the length running invariable thickening limited in 0.1~1.5. Example -d 0.5, Indicating that the offset distance is 0.75 times the average mesh edge length.");
DEFINE_double(l, 2.0, "This arg is a double which value indicates how many times the maximum offset distance is the ideal offset distance limited in 1.5~2.7. .You can set it is 2.0");
DEFINE_double(g, -1, "This arg is a a double means the length of edge length of each cell in grid length. If you can't calculate a length with better performance, it can be passed. Then it will use the default value.");
DEFINE_int32(t, 12, "thread num please set this value depend the cpu of you device.");
DEFINE_string(f, "", "file name, which can choose *.obj2 or *.obj.");
DEFINE_int32(i, 1,  "This arg is a integer means running mode which value can be chose in 1,2. \n 1 is offsetting to the outside of the mesh; 2 is offsetting to the inside of the mesh.");
DEFINE_double(e, 1e-5,  "This arg is a double means the eps. When the distance of two points is smaller than eps, we will these two point is coincide. ");
DEFINE_bool(s, false,  "This arg is a bool value which means the program will skip some cell which is most likely useless. It can improve performance, but may cause holes in the result. We suggest not use this function.");
int result_mode;
int skip_mode = 0;

int main(int argc, char* argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);

    result_mode = FLAGS_m;
    running_mode = FLAGS_i;
    double default_move_dist = FLAGS_d;
    grid_len = FLAGS_g;
    if(default_move_dist <0.1 || default_move_dist >1.5){
        cout <<"invariable thickening limit not in  0.1~1.5" << endl;
        exit(0);
    }
    max_l_times = FLAGS_l;
    if(max_l_times <1.5 || max_l_times >2.7){
        cout <<"i maximum offset distance not in  1.5~2.7" << endl;
        exit(0);
    }
    thread_num =  FLAGS_t;
    string input_filename(FLAGS_f);

    if(FLAGS_s) {
        skip_mode = 1;
        cout <<"skip mode on" << endl;
    }

    merge_eps = FLAGS_e;

    myeps = merge_eps / (1e-3);
    google::ShutDownCommandLineFlags();

    cout <<"CGAL_RELEASE_DATE:" << CGAL_RELEASE_DATE << endl;


    vector<chrono::time_point<chrono::system_clock> > time_path;// start = std::chrono::system_clock::now();
    mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile(input_filename));
    mesh->initBBox();
    mesh->build_fast();
    cout <<"mesh->build_fast() succ" << endl;
    // freopen("../debugoutput.txt","w",stdout);
    // 2.5 552531870234
    // 1.7 550170141470
    cout << grid_len <<endl;
    mix_factor = 0.5;
    double avg_input_mesh = 0;
    for(int i=0;i<mesh->FaceSize();i++){
        double minx = *set<double>{(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)]).norm(),
                                   (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)]).norm(),
                                   (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)]).norm()
        }.begin();
        double maxx = *set<double>{(mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)]).norm(),
                                   (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(1)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)]).norm(),
                                   (mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(2)] - mesh->fast_iGameVertex[mesh->fast_iGameFace[i].vh(0)]).norm()
        }.rbegin();
        avg_input_mesh += sqrt(minx*maxx);
    }
    avg_input_mesh/= mesh->FaceSize();
    if(grid_len <= 0){

        //grid_len = sum/mesh->FaceSize()*1.5;
        grid_len = avg_input_mesh*2.5;

        cout<<"GL:" <<grid_len << endl;
    }

  //3.59


    if(*input_filename.rbegin() != '2') {
        auto tmp = mesh->BBoxMax -  mesh->BBoxMin;
        double minx = mesh->BBoxMin.z();
        double maxx = mesh->BBoxMax.z();
        cout << "default_move_dist: "<< default_move_dist*avg_input_mesh << endl;
        for(int i=0;i<mesh->FaceSize();i++){
//            double this_z = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(0)]+
//                              mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(1)]+
//                              mesh->fast_iGameVertex[mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].vh(2)])/3).z();
//            double ratio = (maxx-this_z)/(maxx-minx);
            mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist = default_move_dist*avg_input_mesh;
//            mesh->fast_iGameFace[MeshKernel::iGameFaceHandle(i)].move_dist = max(default_move_dist*1.5*ratio*ratio,
//                                                                                 default_move_dist*0.5);
            //default_move_dist;
        }
    }




    //mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../data/test_orgv2.obj2")); grid_len = 12.5; double default_move_dist = 0.8;


    faces_approximate_field.resize(mesh->FaceSize());
    field_move_vertex.resize(mesh->VertexSize());
    min_move_g.resize(mesh->VertexSize());
    max_move_g.resize(mesh->VertexSize());

    field_move_face.resize(mesh->FaceSize());
    field_move_K2_triangle.resize(mesh->FaceSize());

    time_path.push_back(std::chrono::system_clock::now());
    for(int i=0;i<mesh->VertexSize();i++){
        bool is_succ = true;
        field_move_vertex[i] = do_quadratic_error_metric(mesh,MeshKernel::iGameVertexHandle(i),is_succ);
    }
    time_path.push_back(std::chrono::system_clock::now());

    cout <<"build st "<< endl;
    std::vector <std::shared_ptr<std::thread> > build_thread_pool(thread_num);
    for(int i=0;i<thread_num;i++) {
        build_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            for(int i=0;i<mesh->FaceSize();i++) {
                if (i % thread_num != now_id)continue;
                if(i%500==0)
                    cout << "build "<< i << endl;
                faces_approximate_field[i] = ApproximateField(MeshKernel::iGameFaceHandle(i),field_move_vertex,mesh);

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

            }
        },i);
    }

    for(int i=0;i<thread_num;i++)
        one_ring_select_thread_pool[i]->join();

    time_path.push_back(std::chrono::system_clock::now());


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

    MeshKernel::iGameVertex debug_v(-0.526370 ,2.852676 ,-0.017140);
    grid debug_g =  vertex_to_grid(debug_v);

    cout <<"v to g :" <<debug_g.x <<" "<< debug_g.y <<" "<<debug_g.z << endl;
    auto tmp = (getGridVertex(debug_g,0) + getGridVertex(debug_g,7))/2;
    cout << tmp.x()<<" "<< tmp.y()<<" "<< tmp.z() << endl;

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
//    FILE *file9 = fopen( (input_filename + "_9.obj").c_str(), "w");
//    for (auto each_grid= frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
//       //if(!(each_grid->first.x == 26 && each_grid->first.y == 28 && each_grid->first.z == 9  ))continue;
//        //if(!(each_grid->first.x == 19 && each_grid->first.y == 19 && each_grid->first.z == 1 ))continue;
//       // if(!(each_grid->first.x == 1 && each_grid->first.y == 0 && each_grid->first.z == 3 ))continue;
//        auto small  = getGridVertex(each_grid->first,0);
//        auto big  = getGridVertex(each_grid->first,7);
//        static int f3_id = 1;
//        for (int ii = 0; ii < 7; ii++) {
//            for (int jj = 0; jj < DirectedGridEdge[ii].size(); jj++) {
//                int from = ii;
//                int to = DirectedGridEdge[ii][jj];
//                MeshKernel::iGameVertex fv = getGridiGameVertex(small, big, from);
//                MeshKernel::iGameVertex tv = getGridiGameVertex(small, big, to);
//                fprintf(file9, "v %lf %lf %lf\n", fv.x(), fv.y(), fv.z());
//                fprintf(file9, "v %lf %lf %lf\n", tv.x(), tv.y(), tv.z());
//                fprintf(file9, "l %d %d\n", f3_id, f3_id + 1);
//                f3_id += 2;
//            }
//        }
//    }

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
                            for(int j=0;j<skip_mode*8;j++){ // close fast mode
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
                            }
                        }

                        if(field_triangles.size() ==0)return ;

                        Tree aabb_tree(field_triangles_all.begin(), field_triangles_all.end());
                        Tree aabb_tree_part(field_triangles_part_aabbtree.begin(), field_triangles_part_aabbtree.end());

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

                        std::vector<std::vector<int> >maybe_used_face_source_id;

                        for(auto this_face : field_triangles){
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
                        }



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

                        vector<K2::Triangle_3 > generated_face_list;

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

                            vector<K2::Segment_3>cs;
                            for(auto j : maybe_used_face_seg_cutting[i]){
                                cs.push_back(j);
                            }


                            K2::Vector_3 origin_normal = CGAL::cross_product(tri_i.vertex(1) -tri_i.vertex(0),
                                                                             tri_i.vertex(2) -tri_i.vertex(0));
                            // cout << "cs" << cs.size() << endl;
                            vector<vector<K2::Point_3> > cdt_res = CGAL_CDT(face_inner_grid_polygon[i],cs,tri_i);


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
                                    int side_type;
                                    if(running_mode == 2){
                                        side_type =cgal_polygon->inMesh((CGAL::centroid(tri)));
                                    }
                                    else{
                                        side_type = cgal_polygon->outMesh((CGAL::centroid(tri)));
                                    }

                                    if(side_type != 1) {
                                        flag = true;
                                        //cout <<"GGST" << endl;
                                    }
                                }


                                auto vvv0 = Point_K2_to_iGameVertex(tri.vertex(0));
                                auto vvv1 = Point_K2_to_iGameVertex(tri.vertex(1));
                                auto vvv2 = Point_K2_to_iGameVertex(tri.vertex(2));

                                auto ddd = (vvv1 - vvv0) % (vvv2- vvv0);

                                if(flag) {
                                    continue;
                                }

                                generated_face_list.push_back(K2::Triangle_3(tri.vertex(0),tri.vertex(2),tri.vertex(1)));


                            }
                        }

                        unique_lock<std::mutex> lock2(mu2);
                        for(auto i : generated_face_list) {
                            each_grid->second.generate_face_list.push_back(i);
                        }
                        lock2.unlock();
                    };
                    dfs(small,big,face_list,0);

                    //for(auto i :  each_grid->second.meshBuilder->generate_v)
                    //each_grid->second.res = MeshBuilder(each_grid->second.generate_face_list,small,big);

                    for(int i=0;i<each_grid->second.generate_face_list.size();i++){
                        for(int j=0;j<3;j++){
                            K2::Point_3 p = each_grid->second.generate_face_list[i].vertex(j);
                            if(abs(CGAL::to_double(p.x()-small.x())) < merge_eps*2){
                                each_grid->second.x_min_v.push_back(p);
                            }
                            if(abs(CGAL::to_double(p.x()-big.x())) < merge_eps*2){
                                each_grid->second.x_max_v.push_back(p);
                            }
                            if(abs(CGAL::to_double(p.y()-small.y())) < merge_eps*2){
                                each_grid->second.y_min_v.push_back(p);
                            }
                            if(abs(CGAL::to_double(p.y()-big.y())) < merge_eps*2){
                                each_grid->second.y_max_v.push_back(p);
                            }
                            if(abs(CGAL::to_double(p.z()-small.z())) < merge_eps*2){
                                each_grid->second.z_min_v.push_back(p);
                            }
                            if(abs(CGAL::to_double(p.z()-big.z())) < merge_eps*2){
                                each_grid->second.z_max_v.push_back(p);
                            }
                        }
                    }


                }
                else {
                    lock.unlock();
                    break;
                }
            }

        },i);
    }
    time_path.push_back(std::chrono::system_clock::now());
    for(int i=0;i<thread_num;i++)
        each_frame_thread[i]->join();
    time_path.push_back(std::chrono::system_clock::now());


    FILE *fileans = fopen( (input_filename + "_times.txt").c_str(), "w+");
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
//    double miny=(1<<30),maxy=-(1<<30);
//    double minz=(1<<30),maxz=-(1<<30);
    unordered_map<size_t ,int> vmp;
    for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++){
       // cout << each_grid->second.generate_face_list.size()*3 << endl;
        for(int i = 0;i < each_grid->second.generate_face_list.size();i++){// 局部id 全局id 按坐标轴排序 结果预存，每一个
            vmp[each_grid->second.generate_face_list[i].id()] = generate_face_final.size();
            generate_face_final.push_back(each_grid->second.generate_face_list[i]);
            final_gen_face.push_back({final_gen_vertex.size(),final_gen_vertex.size()+1,final_gen_vertex.size()+2});
            final_gen_vertex.push_back( each_grid->second.generate_face_list[i].vertex(0));
            final_gen_vertex.push_back( each_grid->second.generate_face_list[i].vertex(1));
            final_gen_vertex.push_back( each_grid->second.generate_face_list[i].vertex(2));
//            for(int k=0;k<3;k++){
//                miny=min(miny, CGAL::to_double(each_grid->second.generate_face_list[i].vertex(k).y()));
//                maxy=max(maxy, CGAL::to_double(each_grid->second.generate_face_list[i].vertex(k).y()));
//                minz=min(minz, CGAL::to_double(each_grid->second.generate_face_list[i].vertex(k).z()));
//                maxz=max(maxz, CGAL::to_double(each_grid->second.generate_face_list[i].vertex(k).z()));
//            }
        }
    }
    if(result_mode == 1){
        FILE *file12 = fopen( (input_filename + "_result.obj").c_str(), "w+");
        int vid = 1;
        for (auto each_grid = frame_grid_mp.begin(); each_grid != frame_grid_mp.end(); each_grid++) {
            for (int i = 0; i < each_grid->second.generate_face_list.size(); i++) {
                fprintf(file12,"v %lf %lf %lf\n",CGAL::to_double(each_grid->second.generate_face_list[i].vertex(0).x()),CGAL::to_double(each_grid->second.generate_face_list[i].vertex(0).y()),CGAL::to_double(each_grid->second.generate_face_list[i].vertex(0).z()));
                fprintf(file12,"v %lf %lf %lf\n",CGAL::to_double(each_grid->second.generate_face_list[i].vertex(1).x()),CGAL::to_double(each_grid->second.generate_face_list[i].vertex(1).y()),CGAL::to_double(each_grid->second.generate_face_list[i].vertex(1).z()));
                fprintf(file12,"v %lf %lf %lf\n",CGAL::to_double(each_grid->second.generate_face_list[i].vertex(2).x()),CGAL::to_double(each_grid->second.generate_face_list[i].vertex(2).y()),CGAL::to_double(each_grid->second.generate_face_list[i].vertex(2).z()));
                fprintf(file12,"f %d %d %d\n",vid,vid+1,vid+2);
                vid+=3;
            }
        }
        fclose(file12);
        exit(0);
    }


    unordered_map<grid, LocalMesh, grid_hash, grid_equal>local_mesh_mp;

    std::vector <std::shared_ptr<std::thread> > each_grid_merge_thread_pool(thread_num);
    std::vector <std::vector<pair<grid,LocalMesh> > > each_grid_merge_thread_pool_ans(thread_num);
    //thread_num = 1;
    for(int i=0;i<thread_num;i++) {
        each_grid_merge_thread_pool[i] = make_shared<std::thread>([&](int now_id) {
            int each_grid_cnt =-1;
            unordered_map <grid, GridVertex, grid_hash, grid_equal> frame_grid_mp_bk = frame_grid_mp;
            for (auto each_grid = frame_grid_mp_bk.begin(); each_grid != frame_grid_mp_bk.end(); each_grid++) {
                K2::Point_3 small = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 0));
                K2::Point_3 big = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 7));

                each_grid_cnt++;
                if(each_grid_cnt % thread_num != now_id)continue;
                if(each_grid_cnt % (thread_num*10) == now_id)
                    printf("each_grid_cnt2 %d/%d/%d\n",now_id,each_grid_cnt,(int)frame_grid_mp_bk.size());

                vector<K2::Point_3> neighbor_v;
                unordered_map <grid, GridVertex, grid_hash, grid_equal>::iterator it;
                grid xup = each_grid->first;
                xup.x += 1;
                it = frame_grid_mp_bk.find(xup);
                if(it!=frame_grid_mp_bk.end()){
                    for(auto j : it->second.x_min_v){
                        neighbor_v.push_back(j);
                    }
                }
                grid xdown = each_grid->first;
                xdown.x -= 1;
                it = frame_grid_mp_bk.find(xdown);
                if(it!=frame_grid_mp_bk.end()){
                    for(auto j : it->second.x_max_v){
                        neighbor_v.push_back(j);
                    }
                }

//********yyyyyyyyyyyyy
                grid yup = each_grid->first;
                yup.y += 1;
                it = frame_grid_mp_bk.find(yup);
                if(it!=frame_grid_mp_bk.end()){
                    for(auto j : it->second.y_min_v){
                        neighbor_v.push_back(j);
                    }
                }
                grid ydown = each_grid->first;
                ydown.y -= 1;
                it = frame_grid_mp_bk.find(ydown);
                if(it!=frame_grid_mp_bk.end()){
                    for(auto j : it->second.y_max_v){
                        neighbor_v.push_back(j);
                    }
                }
//********zzzzzzzzzzzzzz
                grid zup = each_grid->first;
                zup.z += 1;
                it = frame_grid_mp_bk.find(zup);
                if(it!=frame_grid_mp_bk.end()){
                    for(auto j : it->second.z_min_v){
                        neighbor_v.push_back(j);
                    }
                }
                grid zdown = each_grid->first;
                zdown.z -= 1;
                it = frame_grid_mp_bk.find(zdown);
                if(it!=frame_grid_mp_bk.end()){
                    for(auto j : it->second.z_max_v){
                        neighbor_v.push_back(j);
                    }
                }
                MeshBuilder builder(each_grid->second.generate_face_list,small,big,neighbor_v,merge_eps);
                LocalMesh localMesh;
                localMesh.final_v = builder.generate_v;
                localMesh.final_f = builder.generate_face;
                localMesh.built = false;
                each_grid_merge_thread_pool_ans[now_id].push_back({each_grid->first,localMesh});
            }
        }, i);
    }


    for(int i=0;i<thread_num;i++)
        each_grid_merge_thread_pool[i]->join();


    for(int i=0;i<thread_num;i++)
        for(auto j : each_grid_merge_thread_pool_ans[i]){
            local_mesh_mp[j.first] = j.second;
        }

    int global_v_id = 0;
    vector<K2::Point_3>global_v;
    DSU dsu;
    vector<vector<int>>global_face;
    int xxx=0;


    for (auto each_grid = local_mesh_mp.begin(); each_grid != local_mesh_mp.end(); each_grid++){
        cout << xxx++ <<"/?/"<< local_mesh_mp.size()<<endl;
        K2::Point_3 small = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 0));
        K2::Point_3 big = iGameVertex_to_Point_K2(getGridVertex(each_grid->first, 7));
        for(int j=0;j<each_grid->second.final_v.size();j++){
            dsu.append();
            each_grid->second.final_v_global_id.push_back(global_v_id);
            global_v.push_back(each_grid->second.final_v[j]);
            global_v_id++;
        }
        for(auto j : container_grid_dir){
            grid neighbor(each_grid->first.x+j[0],each_grid->first.y+j[1],each_grid->first.z+j[2]);
            auto nit = local_mesh_mp.find(neighbor);
            if(nit != local_mesh_mp.end()){
                if(nit->second.built){
                    for(int k=0;k<each_grid->second.final_v.size();k++){
                        K2::Point_3 p = each_grid->second.final_v[k];
//                        if(each_grid->second.final_v_global_id[k] != -1) { // 这里得改成dsu
//                            //cout << "each_grid->second.final_v_global_id[k]  "<<global_v_id <<" "<< each_grid->second.final_v_global_id[k] <<" ";
//                            int gid = each_grid->second.final_v_global_id[k];
//                            //cout <<  CGAL::to_double(global_v[gid].x()) <<" "<< CGAL::to_double(global_v[gid].y()) <<" "<<CGAL::to_double(global_v[gid].z()) << endl;
//                            continue;
//                        }

                        Fuzzy_circle fc0(p, merge_eps*2);
                        std::list<K2::Point_3> result0;
                        nit->second.kdtree->search(std::back_inserter(result0), fc0);
                        for(const auto l : result0){
                            if(sqrt(CGAL::to_double(CGAL::squared_distance(p,l))) < merge_eps){
//                                if(each_grid->second.final_v_global_id[k] != -1) {
//                                    each_grid->second.final_v_global_id[k] = nit->second.vmp[l.id()];
//                                }

                                    dsu.join(each_grid->second.final_v_global_id[k] , nit->second.vmp[l.id()]);
                                //}
                                break;
                            }
                        }

                    }
                }

            }
        }
        each_grid->second.built = true;

        for(int j=0;j<each_grid->second.final_f.size();j++){
            if(set<int>{ each_grid->second.final_v_global_id[each_grid->second.final_f[j][0]],
                         each_grid->second.final_v_global_id[each_grid->second.final_f[j][1]],
                         each_grid->second.final_v_global_id[each_grid->second.final_f[j][2]]
            }.size() == 3){
                global_face.push_back({ each_grid->second.final_v_global_id[each_grid->second.final_f[j][0]],
                                        each_grid->second.final_v_global_id[each_grid->second.final_f[j][1]],
                                        each_grid->second.final_v_global_id[each_grid->second.final_f[j][2]]
                                      });

            }
        }
        each_grid->second.kdtree = new KDTree();
        for(int j=0;j<each_grid->second.final_v.size();j++){
            int global_id = each_grid->second.final_v_global_id[j];
            K2::Point_3 p = global_v[global_id];
            if(abs(CGAL::to_double(p.x()-big.x())) < merge_eps*2){
                if(!each_grid->second.vmp.count(p.id())) {
                    each_grid->second.kdtree->insert(p);
                    each_grid->second.vmp[p.id()] = global_id;
                }
                //  fprintf(file16,"v %lf %lf %lf\n",CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
                //each_grid->second.x_max_final.push_back(global_id);
            }
            if(abs(CGAL::to_double(p.y()-big.y())) < merge_eps*2){
                if(!each_grid->second.vmp.count(p.id())) {
                    each_grid->second.kdtree->insert(p);
                    each_grid->second.vmp[p.id()] = global_id;
                }
                //  fprintf(file16,"v %lf %lf %lf\n",CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
                //each_grid->second.y_max_final.push_back(global_id);
            }
            if(abs(CGAL::to_double(p.z()-big.z())) < merge_eps*2){
                if(!each_grid->second.vmp.count(p.id())) {
                    each_grid->second.kdtree->insert(p);
                    each_grid->second.vmp[p.id()] = global_id;
                }
                //  fprintf(file16,"v %lf %lf %lf\n",CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
                //each_grid->second.z_max_final.push_back(global_id);
            }

            if(abs(CGAL::to_double(p.x()-small.x())) < merge_eps*2){
                if(!each_grid->second.vmp.count(p.id())) {
                    each_grid->second.kdtree->insert(p);
                    each_grid->second.vmp[p.id()] = global_id;
                }
                //  fprintf(file16,"v %lf %lf %lf\n",CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
                //each_grid->second.x_max_final.push_back(global_id);
            }
            if(abs(CGAL::to_double(p.y()-small.y())) < merge_eps*2){
                if(!each_grid->second.vmp.count(p.id())) {
                    each_grid->second.kdtree->insert(p);
                    each_grid->second.vmp[p.id()] = global_id;
                }
                // fprintf(file16,"v %lf %lf %lf\n",CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
                //each_grid->second.y_max_final.push_back(global_id);
            }
            if(abs(CGAL::to_double(p.z()-small.z())) < merge_eps*2){
                if(!each_grid->second.vmp.count(p.id())) {
                    each_grid->second.kdtree->insert(p);
                    each_grid->second.vmp[p.id()] = global_id;
                }
                //    fprintf(file16,"v %lf %lf %lf\n",CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
                //each_grid->second.z_max_final.push_back(global_id);
            }
        }
    }


    cout <<"st build mesh" << endl;

    FILE *file12 = fopen( (input_filename + "_result.obj").c_str(), "w+");
    vector<vector<int> >global_v_dsu(global_v.size());
    vector<K2::Point_3>global_true_v(global_v.size());
    for(int i=0;i<global_v.size();i++){
        global_v_dsu[dsu.find_root(i)].push_back(i);
    }
    for(int i=0;i<global_v.size();i++){
        K2::Point_3 zero(0,0,0);
        K2::Vector_3 now(0,0,0);
        for(auto j : global_v_dsu[i]){
            now += global_v[j] - zero;
        }
        now/= CGAL::Epeck::FT(global_v_dsu[i].size()*1.0);
        global_true_v[i] = zero + now;
    }


    vector<int>global_v_id_new(global_v.size());
    int id_new = 0;
    for(int i=0;i<global_v.size();i++){
        if(dsu.find_root(i) == i){
            fprintf(file12,"v %lf %lf %lf\n",CGAL::to_double(global_true_v[i].x()),
                    CGAL::to_double(global_true_v[i].y()),
                    CGAL::to_double(global_true_v[i].z())
            );
            global_v_id_new[i] = id_new++;
        }
    }
    for(int i=0;i<global_face.size();i++){
        int id0 = global_v_id_new[dsu.find_root(global_face[i][0])];
        int id1 = global_v_id_new[dsu.find_root(global_face[i][1])];
        int id2 = global_v_id_new[dsu.find_root(global_face[i][2])];
        if(set<int>{id0,id1,id2}.size() !=3)continue;

        fprintf(file12,"f %d %d %d\n",id0+ 1,
                id1 + 1,
                id2 + 1
        );

    }
    fclose(file12);
    if(result_mode == 3)
        Remeshing().run((input_filename + "_result.obj").c_str());


    return 0;
}
