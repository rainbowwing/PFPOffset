//
// Created by rainbowwing on 2023/2/18.
//

#ifndef THICKEN2OUT_LOCALMESH_H
#define THICKEN2OUT_LOCALMESH_H

#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>

typedef CGAL::Search_traits_3<K2> STraits;
typedef CGAL::Fuzzy_sphere<STraits> Fuzzy_circle;
typedef CGAL::Kd_tree<STraits> KDTree;

struct LocalMesh {
    vector<K2::Point_3> final_v;
    vector<int> final_v_global_id;
    vector<vector<int> > final_f;
    map<size_t, int> vmp;
    KDTree *kdtree;
    bool built;
};

#endif //THICKEN2OUT_LOCALMESH_H
