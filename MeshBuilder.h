//
// Created by rainbowwing on 2023/2/18.
//

#ifndef THICKEN2OUT_MESHBUILDER_H
#define THICKEN2OUT_MESHBUILDER_H
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>


typedef CGAL::Search_traits_3<K2>                      STraits;
typedef CGAL::Fuzzy_sphere<STraits>                    Fuzzy_circle;
typedef CGAL::Kd_tree<STraits>                         KDTree;

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
    K2::Point_3 grid_min;
    K2::Point_3 grid_max;
    vector<K2::Point_3> generate_v;
    vector<vector<int> > generate_face;
    DSU dsu;
    unordered_map<size_t,int>vmp;
    unordered_map<size_t,int>ext_vmp;



    double eps;
    int true_cnt = 0;
public:
//    vector<K2::Point_3> arrange(vector<K2::Point_3>ext_v){
//
//
//
////        vector<pair<K2::Point_3,int> >sorted_ex;
////        for(int i=0;i<generate_v.size();i++){
////            sorted_ex.push_back({generate_v[i],0});
////        }
////        for(int i=0;i<ext_v.size();i++){
////            sorted_ex.push_back({ext_v[i],1});
////        }
////        if(variance_x >= variance_y && variance_x >= variance_z) {
////            sort(sorted_ex.begin(),sorted_ex.end(),[&](pair<K2::Point_3,int> a,pair<K2::Point_3,int> b){
////                return a.first.x() < b.first.x();
////            });
////
////            for (int i = 0; i < sorted_ex.size(); i++) {
////                bool useful = true;
////                for (int j = i - 1; j >= 0; j--) {
////                    if (CGAL::to_double(sorted_ex[j].first.x() - sorted_ex[i].first.x()) *
////                        CGAL::to_double(sorted_ex[j].first.x() - sorted_ex[i].first.x()) > eps) {
////                        break;
////                    }
////                    if (CGAL::to_double(CGAL::squared_distance(sorted_ex[j].first, sorted_ex[i].first)) <= eps) {
////                        //OP
////                        useful = false;
////                    }
////                }
////                for (int j = i + 1; j < sorted_ex.size(); j++) {
////                    if (CGAL::to_double(sorted_ex[j].first.x() - sorted_ex[i].first.x()) *
////                        CGAL::to_double(sorted_ex[j].first.x() - sorted_ex[i].first.x()) > eps) {
////                        break;
////                    }
////                    if (CGAL::to_double(CGAL::squared_distance(sorted_ex[j].first, sorted_ex[i].first)) <= eps) {
////                        //OP
////                        if(sorted_ex[j].second == 0){
////                            useful = false;
////                        }
////                    }
////                }
////                if(sorted_ex[i].second == 1 && useful)
////                    sorted_ex[i].second = 2;
////            }
////        }
////        else if(variance_y >= variance_x && variance_y >= variance_z) {
////            sort(sorted_ex.begin(),sorted_ex.end(),[&](pair<K2::Point_3,int> a,pair<K2::Point_3,int> b){
////                return a.first.y() < b.first.y();
////            });
////
////            for (int i = 0; i < sorted_ex.size(); i++) {
////                bool useful = true;
////                for (int j = i - 1; j >= 0; j--) {
////                    if (CGAL::to_double(sorted_ex[j].first.y() - sorted_ex[i].first.y()) *
////                        CGAL::to_double(sorted_ex[j].first.y() - sorted_ex[i].first.y()) > eps) {
////                        break;
////                    }
////                    if (CGAL::to_double(CGAL::squared_distance(sorted_ex[j].first, sorted_ex[i].first)) <= eps) {
////                        //OP
////                        useful = false;
////                    }
////                }
////                for (int j = i + 1; j < sorted_ex.size(); j++) {
////                    if (CGAL::to_double(sorted_ex[j].first.y() - sorted_ex[i].first.y()) *
////                        CGAL::to_double(sorted_ex[j].first.y() - sorted_ex[i].first.y()) > eps) {
////                        break;
////                    }
////                    if (CGAL::to_double(CGAL::squared_distance(sorted_ex[j].first, sorted_ex[i].first)) <= eps) {
////                        //OP
////                        if(sorted_ex[j].second == 0){
////                            useful = false;
////                        }
////                    }
////                }
////                if(sorted_ex[i].second == 1 && useful)
////                    sorted_ex[i].second = 2;
////            }
////        }
////        else {
////            sort(sorted_ex.begin(),sorted_ex.end(),[&](pair<K2::Point_3,int> a,pair<K2::Point_3,int> b){
////                return a.first.z() < b.first.z();
////            });
////
////            for (int i = 0; i < sorted_ex.size(); i++) {
////                bool useful = true;
////                for (int j = i - 1; j >= 0; j--) {
////                    if (CGAL::to_double(sorted_ex[j].first.z() - sorted_ex[i].first.z()) *
////                        CGAL::to_double(sorted_ex[j].first.z() - sorted_ex[i].first.z()) > eps) {
////                        break;
////                    }
////                    if (CGAL::to_double(CGAL::squared_distance(sorted_ex[j].first, sorted_ex[i].first)) <= eps) {
////                        //OP
////                        useful = false;
////                    }
////                }
////                for (int j = i + 1; j < sorted_ex.size(); j++) {
////                    if (CGAL::to_double(sorted_ex[j].first.z() - sorted_ex[i].first.z()) *
////                        CGAL::to_double(sorted_ex[j].first.z() - sorted_ex[i].first.z()) > eps) {
////                        break;
////                    }
////                    if (CGAL::to_double(CGAL::squared_distance(sorted_ex[j].first, sorted_ex[i].first)) <= eps) {
////                        //OP
////                        if(sorted_ex[j].second == 0){
////                            useful = false;
////                        }
////                    }
////                }
////                if(sorted_ex[i].second == 1 && useful)
////                    sorted_ex[i].second = 2;
////            }
////        }
////        vector<K2::Point_3> ret;
////        for(auto i : sorted_ex){
////            if(i.second == 2){
////                ret.push_back(i.first);
////            }
////        }
//        return ret;
//
//
//
//    }
    MeshBuilder(){};
    MeshBuilder(const vector<K2::Triangle_3>&generate_face_final, K2::Point_3 grid_min, K2::Point_3 grid_max,vector<K2::Point_3>ext_v){
        double eps = merge_eps;
        this->grid_min=grid_min;
        this->grid_max=grid_max;
        this->face_list = generate_face_final;
        this->eps = eps;
        true_cnt = 0;
        for(auto i : this->face_list) {
            int sz = v.size();
            v.push_back(i.vertex(0));
            v.push_back(i.vertex(1));
            v.push_back(i.vertex(2));
        }
        dsu = DSU(v.size());
        for(int i=0;i<v.size();i++){
            sort_x.push_back(i);
            sort_y.push_back(i);
            sort_z.push_back(i);
        }
        variance_x=0;
        variance_y=0;
        variance_z=0;
        avg_x=0;
        avg_y=0;
        avg_z=0;

        for(int i=0;i<v.size();i++){
            vmp[v[i].id()] = i;
        }


        for(int i=0;i<v.size();i++){
            avg_x += CGAL::to_double(v[i].x());
            avg_y += CGAL::to_double(v[i].y());
            avg_z += CGAL::to_double(v[i].z());
        }
        avg_x/=v.size();
        avg_y/=v.size();
        avg_z/=v.size();
        for(int i=0;i<v.size();i++){
            variance_x += (CGAL::to_double(v[i].x()) - avg_x) * (CGAL::to_double(v[i].x()) - avg_x);
            variance_y += (CGAL::to_double(v[i].y()) - avg_y) * (CGAL::to_double(v[i].y()) - avg_y);
            variance_z += (CGAL::to_double(v[i].z()) - avg_z) * (CGAL::to_double(v[i].z()) - avg_z);
        }
//        sort(sort_x.begin(),sort_x.end(),[&](int a,int b){
//            return v[a].x() < v[b].x();
//        });
        KDTree kdtree;
        for(int i=0;i<v.size();i++){
            kdtree.insert(v[i]);
        }

        for(int i=0;i<v.size();i++){
            double r0 = merge_eps;
            Fuzzy_circle fc0(v[i], r0);
            std::list<K2::Point_3> result0;
            kdtree.search(std::back_inserter(result0), fc0);
            for(auto j : result0){
                int now_id = vmp[j.id()];
                if(now_id != i && CGAL::to_double(CGAL::squared_distance(j,v[i])) < merge_eps){
                    dsu.join(i,now_id);
                }
            }
        }
//        for(int i=0;i<v.size();i++){
//            for(int j=0;j<v.size();j++){
//                if(dsu.find_root(i) != dsu.find_root(j)){
//                    if(CGAL::to_double(CGAL::squared_distance(v[i],v[j])) < 1e-6)
//                        cout << "ij: "<< i<<" "<<j<< CGAL::to_double(CGAL::squared_distance(v[i],v[j])) << endl;
//                }
//            }
//        }
        //cout <<"prevs:" <<v.size() << endl; // 51 10




        // return ;
//        if(variance_x >= variance_y && variance_x >= variance_z) {
//            sort(sort_x.begin(),sort_x.end(),[&](int a,int b){
//                return v[a].x() < v[b].x();
//            });
//
//            for (int i = 0; i < sort_x.size(); i++) {
//                for (int j = i - 1; j >= 0; j--) {
//                    if (CGAL::to_double(v[sort_x[j]].x() - v[sort_x[i]].x()) *
//                        CGAL::to_double(v[sort_x[j]].x() - v[sort_x[i]].x()) > eps) {
//                        break;
//                    }
//                    if (CGAL::to_double(CGAL::squared_distance(v[sort_x[j]], v[sort_x[i]])) <= eps) {
//                        dsu.join(sort_x[j], sort_x[i]);
//                    }
//                }
//                for (int j = i + 1; j < sort_x.size(); j++) {
//                    if (CGAL::to_double(v[sort_x[j]].x() - v[sort_x[i]].x()) *
//                        CGAL::to_double(v[sort_x[j]].x() - v[sort_x[i]].x()) > eps) {
//                        break;
//                    }
//                    if (CGAL::to_double(CGAL::squared_distance(v[sort_x[j]], v[sort_x[i]])) <= eps) {
//                        dsu.join(sort_x[j], sort_x[i]);
//                    }
//                }
//            }
//        }
//        else if(variance_y >= variance_x && variance_y >= variance_z) {
//            sort(sort_y.begin(),sort_y.end(),[&](int a,int b){
//                return v[a].y() < v[b].y();
//            });
//
//            for (int i = 0; i < sort_y.size(); i++) {
//                for (int j = i - 1; j >= 0; j--) {
//                    if (CGAL::to_double(v[sort_y[j]].y() - v[sort_y[i]].y()) *
//                        CGAL::to_double(v[sort_y[j]].y() - v[sort_y[i]].y()) > eps) {
//                        break;
//                    }
//                    if (CGAL::to_double(CGAL::squared_distance(v[sort_y[j]], v[sort_y[i]])) <= eps) {
//                        dsu.join(sort_y[j], sort_y[i]);
//                    }
//                }
//                for (int j = i + 1; j < sort_x.size(); j++) {
//                    if (CGAL::to_double(v[sort_y[j]].y() - v[sort_y[i]].y()) *
//                        CGAL::to_double(v[sort_y[j]].y() - v[sort_y[i]].y()) > eps) {
//                        break;
//                    }
//                    if (CGAL::to_double(CGAL::squared_distance(v[sort_y[j]], v[sort_y[i]])) <= eps) {
//                        dsu.join(sort_y[j], sort_y[i]);
//                    }
//                }
//            }
//        }
//        else {
//            sort(sort_z.begin(),sort_z.end(),[&](int a,int b){
//                return v[a].z() < v[b].z();
//            });
//            for (int i = 0; i < sort_z.size(); i++) {
//                for (int j = i - 1; j >= 0; j--) {
//                    if (CGAL::to_double(v[sort_z[j]].z() - v[sort_z[i]].z()) *
//                        CGAL::to_double(v[sort_z[j]].z() - v[sort_z[i]].z()) > eps) {
//                        break;
//                    }
//                    if (CGAL::to_double(CGAL::squared_distance(v[sort_z[j]], v[sort_z[i]])) <= eps) {
//                        dsu.join(sort_z[j], sort_z[i]);
//                    }
//                }
//                for (int j = i + 1; j < sort_z.size(); j++) {
//                    if (CGAL::to_double(v[sort_z[j]].z() - v[sort_z[i]].z()) *
//                        CGAL::to_double(v[sort_z[j]].z() - v[sort_z[i]].z()) > eps) {
//                        break;
//                    }
//                    if (CGAL::to_double(CGAL::squared_distance(v[sort_z[j]], v[sort_z[i]])) <= eps) {
//                        dsu.join(sort_z[j], sort_z[i]);
//                    }
//                }
//            }
//        }




        // 加入格子情报  计算upmax

        kdtree.clear();
        for(int i=0;i<v.size();i++){
            int ri = dsu.find_root(i);
            if(i == ri){
                true_id[i] = true_cnt;
                true_cnt++;
                kdtree.insert(v[i]);
            }
        }


        int pre_vsize = v.size();
        vector<K2::Point_3> ext_v_bk = ext_v;
        ext_v.clear();
        for(int i=0;i<ext_v_bk.size();i++){
            double r0 = merge_eps;
            Fuzzy_circle fc0(ext_v_bk[i], r0);
            std::list<K2::Point_3> result0;
            kdtree.search(std::back_inserter(result0), fc0);
            bool flag = true;
            for(auto j : result0){
                if(CGAL::to_double(CGAL::squared_distance(j,ext_v_bk[i])) < merge_eps){
                    flag = false;
                    break;
                }
            }
            if(flag){
                ext_v.push_back(ext_v_bk[i]);
            }
        }


        map<size_t,int> ex_id;
        for(auto i : ext_v){
            ex_id[i.id()] = -1;
            kdtree.insert(i);
        }



        //这里跑双尺取，先内点自己合并，然后和外面尺取，最后剩下的点，开一个东西记录


        for(int i=0;i<pre_vsize;i+=3){
            int id0 = true_id[dsu.find_root(i)];
            int id1 = true_id[dsu.find_root(i+1)];
            int id2 = true_id[dsu.find_root(i+2)];// 直接中点法
            if(set<int>{id0,id1,id2}.size() != 3)
                continue;

            vector<int>cut0;
            vector<int>cut1;
            vector<int>cut2;

            K2::Point_3 mid0 = midpoint(K2::Segment_3(v[dsu.find_root(i)],v[dsu.find_root(i+1)]));
            double r0 = sqrt(CGAL::to_double(K2::Segment_3(v[dsu.find_root(i)],v[dsu.find_root(i+1)]).squared_length())/2);
            Fuzzy_circle fc0(mid0, r0);
            std::list<K2::Point_3> result0;
            kdtree.search(std::back_inserter(result0), fc0);

            K2::Point_3 mid1 = midpoint(K2::Segment_3(v[dsu.find_root(i+1)],v[dsu.find_root(i+2)]));
            double r1 = sqrt(CGAL::to_double(K2::Segment_3(v[dsu.find_root(i+1)],v[dsu.find_root(i+2)]).squared_length())/2);
            Fuzzy_circle fc1(mid1, r1);
            std::list<K2::Point_3> result1;
            kdtree.search(std::back_inserter(result1), fc1);

            K2::Point_3 mid2 = midpoint(K2::Segment_3(v[dsu.find_root(i+2)],v[dsu.find_root(i)]));
            double r2= sqrt(CGAL::to_double(K2::Segment_3(v[dsu.find_root(i+2)],v[dsu.find_root(i)]).squared_length())/2);
            Fuzzy_circle fc2(mid2, r2);
            std::list<K2::Point_3> result2;
            kdtree.search(std::back_inserter(result2), fc2);

            for(auto j: result0){
                if(CGAL::to_double(CGAL::squared_distance(j, K2::Segment_3(v[dsu.find_root(i)],v[dsu.find_root(i+1)]) )) <= eps){
                    if(ex_id.count(j.id())){
                        if(ex_id[j.id()] == -1){
                            int new_id = v.size();
                            ex_id[j.id()] = new_id;
                            v.push_back(j);
                            dsu.append();
                            true_id[new_id] = true_cnt;
                            true_cnt++;
                            cut0.push_back(new_id);
                        }
                        else
                            cut0.push_back(ex_id[j.id()]);
                        continue;
                    }
                    if(dsu.find_root(vmp[j.id()]) == dsu.find_root(i+2)){
                        goto end;
                    }
                    else if(dsu.find_root(vmp[j.id()]) == dsu.find_root(i) || dsu.find_root(vmp[j.id()]) == dsu.find_root(i+1))
                        continue;
                    else{
                        cut0.push_back(dsu.find_root(vmp[j.id()]));
                    }
                }
            }

            for(auto j: result1){
                if(CGAL::to_double(CGAL::squared_distance(j, K2::Segment_3(v[dsu.find_root(i+1)],v[dsu.find_root(i+2)]) )) <= eps){
                    if(ex_id.count(j.id())){
                        if(ex_id[j.id()] == -1){
                            int new_id = v.size();
                            ex_id[j.id()] = new_id;
                            v.push_back(j);
                            dsu.append();
                            true_id[new_id] = true_cnt;
                            true_cnt++;
                            cut1.push_back(new_id);
                        }
                        else
                            cut1.push_back(ex_id[j.id()]);
                        continue;
                    }

                    if(dsu.find_root(vmp[j.id()]) == dsu.find_root(i)){
                        goto end;
                    }
                    else if(dsu.find_root(vmp[j.id()]) == dsu.find_root(i+1) || dsu.find_root(vmp[j.id()]) == dsu.find_root(i+2))
                        continue;
                    else{
                        cut1.push_back(dsu.find_root(vmp[j.id()]));
                    }
                }
            }

            for(auto j: result2){
                if(CGAL::to_double(CGAL::squared_distance(j, K2::Segment_3(v[dsu.find_root(i+2)],v[dsu.find_root(i)]) )) <= eps){
                    if(ex_id.count(j.id())){
                        if(ex_id[j.id()] == -1){
                            int new_id = v.size();
                            ex_id[j.id()] = new_id;
                            v.push_back(j);
                            dsu.append();
                            true_id[new_id] = true_cnt;
                            true_cnt++;
                            cut2.push_back(new_id);
                        }
                        else
                            cut2.push_back(ex_id[j.id()]);
                        continue;
                    }

                    if(dsu.find_root(vmp[j.id()]) == dsu.find_root(i+1)){
                        goto end;
                    }
                    else if(dsu.find_root(vmp[j.id()]) == dsu.find_root(i+2) || dsu.find_root(vmp[j.id()]) == dsu.find_root(i))
                        continue;
                    else{
                        cut2.push_back(dsu.find_root(vmp[j.id()]));
                    }
                }
            }

            // todo: 究极debug
//            if(cut0.size() + cut1.size() + cut2.size() != 0) {
//                int new_id = v.size();
//                v.push_back(CGAL::centroid(K2::Triangle_3(v[dsu.find_root(i)],v[dsu.find_root(i+1)],v[dsu.find_root(i+2)])));
//                dsu.append();
//
//                true_id[new_id] = true_cnt;
//                true_cnt++;
//
//                sort(cut0.begin(),cut0.end());
//                cut0.resize(unique(cut0.begin(),cut0.end())-cut0.begin());
//                sort(cut0.begin(),cut0.end(),[&](int a,int b){
//                    return CGAL::squared_distance(v[dsu.find_root(i)],v[a]) < CGAL::squared_distance(v[dsu.find_root(i)],v[b]);
//                });
//
//                sort(cut1.begin(),cut1.end());
//                cut1.resize(unique(cut1.begin(),cut1.end())-cut1.begin());
//                sort(cut1.begin(),cut1.end(),[&](int a,int b){
//                    return CGAL::squared_distance(v[dsu.find_root(i+1)],v[a]) < CGAL::squared_distance(v[dsu.find_root(i+1)],v[b]);
//                });
//
//                sort(cut2.begin(),cut2.end());
//                cut2.resize(unique(cut2.begin(),cut2.end())-cut2.begin());
//                sort(cut2.begin(),cut2.end(),[&](int a,int b){
//                    return CGAL::squared_distance(v[dsu.find_root(i+2)],v[a]) < CGAL::squared_distance(v[dsu.find_root(i+2)],v[b]);
//                });
//                cut0.insert(cut0.begin(),dsu.find_root(i));
//                cut0.push_back(dsu.find_root(i+1));
//
//                cut1.insert(cut1.begin(),dsu.find_root(i+1));
//                cut1.push_back(dsu.find_root(i+2));
//
//                cut2.insert(cut2.begin(),dsu.find_root(i+2));
//                cut2.push_back(dsu.find_root(i));
//
//                for(int j=0;j<(int)cut0.size()-1;j++){
//                    generate_face.push_back({true_id[new_id],true_id[cut0[j]],true_id[cut0[j+1]]});
//                }
//                for(int j=0;j<(int)cut1.size()-1;j++){
//                    generate_face.push_back({true_id[new_id],true_id[cut1[j]],true_id[cut1[j+1]]});
//                }
//                for(int j=0;j<(int)cut2.size()-1;j++){
//                    generate_face.push_back({true_id[new_id],true_id[cut2[j]],true_id[cut2[j+1]]});
//                }
//            }
//            else
                generate_face.push_back({id0,id1,id2});
            end:;
        }
        int cnt=0;
        for(int i=0;i<v.size();i++){
            if(dsu.find_root(i) == i){
                cnt++;
            }
        }
//        cout << cnt<<"??" << true_cnt << endl;
//        cout <<"pre_vsize:: " <<pre_vsize << endl;
//        if(cnt != true_cnt){
//            for(int i=0;i<v.size();i++){
//                cout << i <<" "<< dsu.find_root(i) <<" "<< true_id[dsu.find_root(i)] << endl;
//            }
//
//
//        }

//        assert(cnt == true_cnt);


        generate_v.resize(cnt);
        for(int i=0;i<v.size();i++){
            if(dsu.find_root(i) == i){
                generate_v[true_id[i]] = v[i];
            }
        }

//        for(int i=0;i<v.size();i++){
//            cout << "eachv"<<i<<" "<< true_id[i] <<" "<< dsu.find_root(i) <<" " <<CGAL::to_double(v[i].x()) <<" "<<CGAL::to_double(v[i].y())<<" "<< CGAL::to_double(v[i].z())<<endl;;
//        }
//
//
//        for(int i=0;i<generate_v.size();i++){
//            cout <<"generate_v：" <<i <<" "<< CGAL::to_double(generate_v[i].x()) <<" "<<CGAL::to_double(generate_v[i].y())<<" "<< CGAL::to_double(generate_v[i].z())<<endl;
//        }
//
//
//        for(auto item : generate_face){
//            cout <<"??" <<CGAL::to_double(CGAL::squared_distance(generate_v[item[0]],generate_v[item[1]])) << " "<<item[0]<<" "<<item[1]<<" "<< (CGAL::squared_distance(generate_v[item[0]],generate_v[item[1]])<merge_eps) << endl;
//            cout <<"??" <<CGAL::to_double(CGAL::squared_distance(generate_v[item[1]],generate_v[item[2]])) << " "<<item[1]<<" "<<item[2]<<" "<< (CGAL::squared_distance(generate_v[item[1]],generate_v[item[2]])<merge_eps)<<  endl;
//            cout <<"??" <<CGAL::to_double(CGAL::squared_distance(generate_v[item[2]],generate_v[item[0]])) <<" " <<item[2]<<" "<<item[0]<<" "<< (CGAL::squared_distance(generate_v[item[2]],generate_v[item[0]])<merge_eps)<<  endl;
//        }
//



    } // 加入格子情报

};

#endif //THICKEN2OUT_MESHBUILDER_H
