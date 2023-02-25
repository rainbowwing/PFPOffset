//
// Created by rainbowwing on 2023/2/18.
//

#ifndef THICKEN2OUT_MESHBUILDER_H
#define THICKEN2OUT_MESHBUILDER_H

class MeshBuilder{
public:
    vector<K2::Triangle_3>face_list;
    vector<K2::Point_3> v;
    map<int,int> true_id;
    vector<int> sort_x;
    vector<int> sort_y;
    vector<int> sort_z;
    double variance_x;
    double variance_y;
    double variance_z;
    double avg_x;
    double avg_y;
    double avg_z;
    vector<K2::Point_3> generate_v;
    vector<vector<int> > generate_face;
    DSU dsu;
    unordered_map<size_t,int>vmp;
    unordered_map<size_t,int>ext_vmp;
    double merge_eps;
    int true_cnt = 0;
public:

    MeshBuilder(){};
    MeshBuilder(const vector<K2::Triangle_3>&generate_face_final, K2::Point_3 grid_min, K2::Point_3 grid_max,vector<K2::Point_3>ext_v,double merge_eps);

};

#endif //THICKEN2OUT_MESHBUILDER_H
