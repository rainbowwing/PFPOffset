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

    int true_cnt = 0;
public:

    MeshBuilder(){};
    MeshBuilder(const vector<K2::Triangle_3>&generate_face_final, K2::Point_3 grid_min, K2::Point_3 grid_max,vector<K2::Point_3>ext_v,int a,int b,int c){

        this->grid_min=grid_min;
        this->grid_max=grid_max;
        this->face_list = generate_face_final;
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
                if(now_id != i && sqrt(CGAL::to_double(CGAL::squared_distance(j,v[i]))) < merge_eps){
                    dsu.join(i,now_id);
                }
            }
        }


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
                if(sqrt(CGAL::to_double(CGAL::squared_distance(j,ext_v_bk[i]))) < merge_eps){
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
                if(sqrt(CGAL::to_double(CGAL::squared_distance(j, K2::Segment_3(v[dsu.find_root(i)],v[dsu.find_root(i+1)]) ))) <= merge_eps){
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
                if(sqrt(CGAL::to_double(CGAL::squared_distance(j, K2::Segment_3(v[dsu.find_root(i+1)],v[dsu.find_root(i+2)]) ))) <= merge_eps){
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
                if(sqrt(CGAL::to_double(CGAL::squared_distance(j, K2::Segment_3(v[dsu.find_root(i+2)],v[dsu.find_root(i)]) ))) <= merge_eps){
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

            if(cut0.size() + cut1.size() + cut2.size() != 0) {
                int new_id = v.size();
                v.push_back(CGAL::centroid(K2::Triangle_3(v[dsu.find_root(i)],v[dsu.find_root(i+1)],v[dsu.find_root(i+2)])));
                dsu.append();

                true_id[new_id] = true_cnt;
                true_cnt++;

                sort(cut0.begin(),cut0.end());
                cut0.resize(unique(cut0.begin(),cut0.end())-cut0.begin());
                sort(cut0.begin(),cut0.end(),[&](int a,int b){
                    return CGAL::squared_distance(v[dsu.find_root(i)],v[a]) < CGAL::squared_distance(v[dsu.find_root(i)],v[b]);
                });

                sort(cut1.begin(),cut1.end());
                cut1.resize(unique(cut1.begin(),cut1.end())-cut1.begin());
                sort(cut1.begin(),cut1.end(),[&](int a,int b){
                    return CGAL::squared_distance(v[dsu.find_root(i+1)],v[a]) < CGAL::squared_distance(v[dsu.find_root(i+1)],v[b]);
                });

                sort(cut2.begin(),cut2.end());
                cut2.resize(unique(cut2.begin(),cut2.end())-cut2.begin());
                sort(cut2.begin(),cut2.end(),[&](int a,int b){
                    return CGAL::squared_distance(v[dsu.find_root(i+2)],v[a]) < CGAL::squared_distance(v[dsu.find_root(i+2)],v[b]);
                });
                cut0.insert(cut0.begin(),dsu.find_root(i));
                cut0.push_back(dsu.find_root(i+1));

                cut1.insert(cut1.begin(),dsu.find_root(i+1));
                cut1.push_back(dsu.find_root(i+2));

                cut2.insert(cut2.begin(),dsu.find_root(i+2));
                cut2.push_back(dsu.find_root(i));

                for(int j=0;j<(int)cut0.size()-1;j++){
                    generate_face.push_back({true_id[new_id],true_id[cut0[j]],true_id[cut0[j+1]]});
                }
                for(int j=0;j<(int)cut1.size()-1;j++){
                    generate_face.push_back({true_id[new_id],true_id[cut1[j]],true_id[cut1[j+1]]});
                }
                for(int j=0;j<(int)cut2.size()-1;j++){
                    generate_face.push_back({true_id[new_id],true_id[cut2[j]],true_id[cut2[j+1]]});
                }
            }
            else
                generate_face.push_back({id0,id1,id2});
            end:;
        }
        int cnt=0;
        for(int i=0;i<v.size();i++){
            if(dsu.find_root(i) == i){
                cnt++;
            }
        }
        generate_v.resize(cnt);
        for(int i=0;i<v.size();i++){
            if(dsu.find_root(i) == i){
                generate_v[true_id[i]] = v[i];
            }

        }

    }

};

#endif //THICKEN2OUT_MESHBUILDER_H
