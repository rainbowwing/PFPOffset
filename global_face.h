//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_GLOBAL_FACE_H
#define THICKEN2_GLOBAL_FACE_H
struct GlobalFace{
    int idx0,idx1,idx2;
    int field_id;
    int useful;
    K2::Point_3 center;
    set<int>special_field_id;
    // 特例： 交点是自己
    // 交点是二合一点
    //

    std::unordered_map<int,vector<pair<K2::Point_3,int> > > ray_detect_map;//field 交点 int
};



vector<K2::Point_3> get_sampling_point(K2::Triangle_3 tri){
    vector<K2::Point_3> ret;
//    for(int i=0;i<3;i++){
//        ret.push_back(tri.vertex(i));
//    }





    K2::Point_3 center = centroid(tri);
    ret.push_back(center);
    for(int i=0;i<3;i++){
        K2::Vector_3 vec = (center - tri.vertex(i))/CGAL::Epeck::FT(5000.0);
        ret.push_back(tri.vertex(i) + vec);
    }
    for(int i=0;i<3;i++){
        K2::Point_3 mid = midpoint((K2::Segment_3 (tri.vertex(i),tri.vertex((i+1)%3) )));
        K2::Vector_3 vec = (center - mid)/CGAL::Epeck::FT(5000.0);
        ret.push_back(mid + vec);
    }

    for(int i=0;i<3;i++){
        ret.push_back(centroid(K2::Triangle_3(tri.vertex(i),tri.vertex((i+1)%3),center)));
    }
    return ret;
}

vector<MeshKernel::iGameVertex> field_move_vertex;

vector<vector<K2::Point_3> > field_move_vertices;
vector<K2::Point_3>origin_mesh_vertices;

vector<vector<MeshKernel::iGameVertex> > field_move_face;
vector<K2::Triangle_3> field_move_K2_triangle;
vector<K2::Point_3> global_vertex_list;
vector<int> global_vertex_list_cnt;
vector<K2::Vector_3> global_vertex_list_avg;
vector<GlobalFace> global_face_list;
vector<double>merge_limit;
#endif //THICKEN2_GLOBAL_FACE_H
