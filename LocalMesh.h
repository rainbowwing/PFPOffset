//
// Created by rainbowwing on 2023/2/18.
//

#ifndef THICKEN2OUT_LOCALMESH_H
#define THICKEN2OUT_LOCALMESH_H


struct LocalMesh {
    vector <K2::Point_3> final_v;
    vector <int> final_v_global_id;
    vector <vector<int> > final_f;
    vector<int> x_max_final;
    vector<int> x_min_final;
    vector<int> y_max_final;
    vector<int> y_min_final;
    vector<int> z_max_final;
    vector<int> z_min_final;

};

#endif //THICKEN2OUT_LOCALMESH_H
