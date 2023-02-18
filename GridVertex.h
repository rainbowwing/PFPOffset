//
// Created by rainbowwing on 2023/2/18.
//

#ifndef THICKEN2OUT_GRIDVERTEX_H
#define THICKEN2OUT_GRIDVERTEX_H


struct GridVertex {
    int grid_type;
    vector <MeshKernel::iGameFaceHandle> face_list;
    vector <int> final_global_id;

    vector<K2::Triangle_3> generate_face_list;


    vector<K2::Point_3> x_max_v;
    vector<K2::Point_3> x_min_v;
    vector<K2::Point_3> y_max_v;
    vector<K2::Point_3> y_min_v;
    vector<K2::Point_3> z_max_v;
    vector<K2::Point_3> z_min_v;


    GridVertex() {

        grid_type = -1;
    }
};

#endif //THICKEN2OUT_GRIDVERTEX_H
