//
// Created by rainbowwing on 2023/10/3.
//

#ifndef THICKEN2_ONE_RING_SELECT_H
#define THICKEN2_ONE_RING_SELECT_H
void one_ring_select(){
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
                    for (int k = 0; k < coverage_field_list[neighbor_id].bound_face_id.size(); k++) {
                        K2::Triangle_3 tri_this(coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][0]],
                                                coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][1]],
                                                coverage_field_list[neighbor_id].bound_face_vertex_exact[coverage_field_list[neighbor_id].bound_face_id[k][2]]
                        );

                        neighbor_face.push_back(tri_this);
                    }
                }

                for(int j=0;j<coverage_field_list[i].bound_face_id.size();j++) {
                    if ((coverage_field_list[i].bound_face_id[j][0] >= 3 ||
                         coverage_field_list[i].bound_face_id[j][1] >= 3 ||
                         coverage_field_list[i].bound_face_id[j][2] >= 3)){
                        bool flag = false;

                        K2::Triangle_3 tri_this(coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
                                                coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
                                                coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]);
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
                        vector<vector<K2::Point_3> > res = CGAL_CDT_NEW({coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][0]],
                                                                         coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][1]],
                                                                         coverage_field_list[i].bound_face_vertex_exact[coverage_field_list[i].bound_face_id[j][2]]}, vs, tri_this);
                        //cout << res.size() <<endl;
                        for (auto each_tri: res) {
                            bool patch_flag = false;
                            K2::Point_3 center = CGAL::centroid(K2::Triangle_3(each_tri[0], each_tri[1], each_tri[2]));
                            for (auto j: neighbor_field) {
                                if (tri_this.supporting_plane().oriented_side(coverage_field_list[i].center) !=
                                    tri_this.supporting_plane().oriented_side(coverage_field_list[j].center) &&
                                    tri_this.supporting_plane().oriented_side(coverage_field_list[i].center)+
                                    tri_this.supporting_plane().oriented_side(coverage_field_list[j].center) ==0 &&
                                    coverage_field_list[j].in_or_on_field(center)) {
                                    patch_flag = true;
                                    break;
                                }
                            }
                            if(!patch_flag){
                                flag = true;
                                break;
                            }
                        }
                        coverage_field_list[i].bound_face_useful[j]=flag;
                        //cout <<i<<" "<<j <<" "<< (flag?"yes":"no") << endl;
                    }
                    else{
                        coverage_field_list[i].bound_face_useful[j]=false;
                    }

                }
            }
        },i);
    }

    for(int i=0;i<thread_num;i++)
        one_ring_select_thread_pool[i]->join();
}
#endif //THICKEN2_ONE_RING_SELECT_H
