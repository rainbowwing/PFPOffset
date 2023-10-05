//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_COVERAGEFIELD_H
#define THICKEN2_COVERAGEFIELD_H
struct CoverageField {
    vector<K2::Point_3 > bound_face_vertex_exact;
    vector<K::Point_3 > bound_face_vertex_inexact;
    vector<RoughVertex> bound_face_vertex_rough;
    vector<vector<int> > bound_face_id;

    vector<vector<grid> >  bound_face_cross_field_list;
    vector<vector<K2::Segment_3> > bound_face_cutting_segment;
    vector<vector<K2::Point_3> > bound_face_cutting_point;
    vector<bool> bound_face_useful;

    K2::Point_3 center;
    std::vector<std::vector<K2::Point_3> > cdt_result;
    vector<int>cdt_result_cross_field_list_id;
    vector<bool>cdt_result_cross_field_list_useful;

    vector<K2::Point_3 > renumber_bound_face_vertex;
    vector<vector<int> > renumber_bound_face_id;
    vector<vector<grid> > renumber_bound_face_cross_field_list;
    int field_id;
    void do_cdt(){
        //            vector<K2::Point_3> sorted_bound_vertex;
//            vector<K2::Segment_3> cs;
//            set<pair<int,int> >cs_set;
        for(int i=0;i<bound_face_id.size();i++){
            vector<K2::Point_3> sorted_bound_vertex{bound_face_vertex_exact[bound_face_id[i][0]],
                                                    bound_face_vertex_exact[bound_face_id[i][1]],
                                                    bound_face_vertex_exact[bound_face_id[i][2]]
            };
            if(!bound_face_useful[i]) {
                cdt_result.push_back(sorted_bound_vertex);
                cdt_result_cross_field_list_id.push_back(i);
                cdt_result_cross_field_list_useful.push_back(false);
                continue;
            }

            vector<K2::Segment_3> cs;
            for(auto j : bound_face_cutting_segment[i]){
                int cnt = 0;
                cnt += (j.vertex(0) == bound_face_vertex_exact[bound_face_id[i][0]] );
                cnt += (j.vertex(0) == bound_face_vertex_exact[bound_face_id[i][1]] );
                cnt += (j.vertex(0) == bound_face_vertex_exact[bound_face_id[i][2]] );
                cnt += (j.vertex(1) == bound_face_vertex_exact[bound_face_id[i][0]] );
                cnt += (j.vertex(1) == bound_face_vertex_exact[bound_face_id[i][1]] );
                cnt += (j.vertex(1) == bound_face_vertex_exact[bound_face_id[i][2]] );
                if(cnt <2 ){
                    cs.push_back(j);
                }
            }
            sort(cs.begin(),cs.end(),[&](const K2::Segment_3& a,const K2::Segment_3& b){
                if(a.vertex(0).x() != b.vertex(0).x()){
                    return a.vertex(0).x() < b.vertex(0).x();
                }
                else if(a.vertex(0).y() != b.vertex(0).y()){
                    return a.vertex(0).y() < b.vertex(0).y();
                }
                else if(a.vertex(0).z() != b.vertex(0).z()){
                    return a.vertex(0).z() < b.vertex(0).z();
                }
                else if(a.vertex(1).x() != b.vertex(1).x()){
                    return a.vertex(1).x() < b.vertex(1).x();
                }
                else if(a.vertex(1).y() != b.vertex(1).y()){
                    return a.vertex(1).y() < b.vertex(1).y();
                } else
                    return a.vertex(1).z() < b.vertex(1).z();
            });


            cs.resize(std::unique(cs.begin(),cs.end())-cs.begin());

            K2::Triangle_3 tri(bound_face_vertex_exact[bound_face_id[i][0]],
                               bound_face_vertex_exact[bound_face_id[i][1]],
                               bound_face_vertex_exact[bound_face_id[i][2]]);

            for(auto j : bound_face_cutting_point[i]){
                if(j != bound_face_vertex_exact[bound_face_id[i][0]] &&
                   j != bound_face_vertex_exact[bound_face_id[i][1]] &&
                   j != bound_face_vertex_exact[bound_face_id[i][2]]
                        )
                    sorted_bound_vertex.push_back(j);
            }
            sort(sorted_bound_vertex.begin(),sorted_bound_vertex.end(),[&](const K2::Point_3& a, const K2::Point_3& b){
                if(a.x() != b.x()){
                    return a.x() < b.x();
                }
                else if(a.y() != b.y()){
                    return a.y() < b.y();
                }
                return a.z() < b.z();
            });
            sorted_bound_vertex.resize(std::unique(sorted_bound_vertex.begin(),sorted_bound_vertex.end())-sorted_bound_vertex.begin());

            vector<vector<K2::Point_3> > res = CGAL_CDT_NEW(sorted_bound_vertex,cs,tri);

            for(int j=0;j<res.size();j++){
                cdt_result.push_back(res[j]);
                cdt_result_cross_field_list_id.push_back(i);
                cdt_result_cross_field_list_useful.push_back(true);
            }

        }

    }

    Kd_tree * tree;
    vector<K::Vector_3>encode_vec;
    vector<int>encode_num;
    vector<int>renumber_bound_face_vertex_global_id;
    vector<int>renumber_bound_face_global_id;
    vector<bool>renumber_bound_face_useful;

    std::unordered_map<unsigned long long,int> encode_map;
    int get_kdtree_id(K::Point_3 p){
        std::vector<Point> result;
        Fuzzy_circle fs(p, myeps);
        tree->search(std::back_inserter(result), fs);
        //  cout <<"result.size() : " << result.size() << endl;
        return encode_map[unique_hash_value(*result.begin())];
    };
    MeshKernel::iGameVertex bbox_min,bbox_max;
    void renumber(){
        std::vector<K::Point_3> kd_tree_points;
        DSU dsu;
        for(int i=0;i<cdt_result.size();i++){
            for(int j=0;j<3;j++){
                kd_tree_points.emplace_back(CGAL::to_double(cdt_result[i][j].x()),
                                            CGAL::to_double(cdt_result[i][j].y()),
                                            CGAL::to_double(cdt_result[i][j].z())
                );
            }
        }

        tree = new Kd_tree(kd_tree_points.begin(),kd_tree_points.end());
        for(int i=0;i<kd_tree_points.size();i++){
            std::vector<Point> result;
            Fuzzy_circle fs(kd_tree_points[i], myeps);
            tree->search(std::back_inserter(result), fs);
            for (const Point& p : result) {
                dsu.join(unique_hash_value(kd_tree_points[i]),unique_hash_value(p));
            }
        }
        int encode_cnt = 0;
        encode_map = dsu.encode(encode_cnt);
        renumber_bound_face_vertex.resize(encode_cnt);
        encode_vec.resize(encode_cnt);
        encode_num.resize(encode_cnt);
        renumber_bound_face_vertex_global_id.resize(encode_cnt);
        for(int i=0;i<encode_cnt;i++){
            encode_vec[i] = K::Vector_3(0,0,0);
            encode_num[i] = 0;
        }

        for(int i=0;i<kd_tree_points.size();i++){
            //  cout <<"getbelong:" <<i<<" "<< encode_map[CGAL::hash_value(kd_tree_points[i])] <<endl;
            encode_vec[encode_map[unique_hash_value(kd_tree_points[i])]] += kd_tree_points[i] - K::Point_3 (0,0,0);
            encode_num[encode_map[unique_hash_value(kd_tree_points[i])]]++;
        }


        for(int i=0;i<encode_cnt;i++){
            K::Point_3 avg = K::Point_3 (0,0,0) + encode_vec[i] / encode_num[i];
            renumber_bound_face_vertex[i] = K2::Point_3(avg.x(),avg.y(),avg.z());
        }
        //重构搜索要用原来的点
        for(int i=0;i<cdt_result.size();i++){

            K::Point_3 v0(CGAL::to_double(cdt_result[i][0].x()),CGAL::to_double(cdt_result[i][0].y()),CGAL::to_double(cdt_result[i][0].z()));
            K::Point_3 v1(CGAL::to_double(cdt_result[i][1].x()),CGAL::to_double(cdt_result[i][1].y()),CGAL::to_double(cdt_result[i][1].z()));
            K::Point_3 v2(CGAL::to_double(cdt_result[i][2].x()),CGAL::to_double(cdt_result[i][2].y()),CGAL::to_double(cdt_result[i][2].z()));
            int id0 = get_kdtree_id(v0);
            int id1 = get_kdtree_id(v1);
            int id2 = get_kdtree_id(v2);
            if(set<int>{id0,id1,id2}.size() != 3)continue;
            renumber_bound_face_id.push_back({id0,id1,id2});
            renumber_bound_face_cross_field_list.push_back(bound_face_cross_field_list[cdt_result_cross_field_list_id[i]]);
            renumber_bound_face_useful.push_back(cdt_result_cross_field_list_useful[i]);

        }
        std::list<K2::Triangle_3>tri_list;
        for(int i=0;i<renumber_bound_face_id.size();i++){
            tri_list.emplace_back(renumber_bound_face_vertex[renumber_bound_face_id[i][0]],
                                  renumber_bound_face_vertex[renumber_bound_face_id[i][1]],
                                  renumber_bound_face_vertex[renumber_bound_face_id[i][2]]);
        }
        Tree aabb_tree(tri_list.begin(),tri_list.end());
        auto iter = tri_list.begin();
        for(int i=0;i<renumber_bound_face_id.size();i++,iter++){
            K2::Point_3 tri_center = CGAL::centroid(*iter);
            K2::Ray_3 ray(tri_center,iter->supporting_plane().orthogonal_vector());
            std::list< Tree::Intersection_and_primitive_id<K2::Ray_3>::Type> intersections;
            aabb_tree.all_intersections(ray,std::back_inserter(intersections));
            vector<K2::Point_3> intersection_v;
            for(auto item : intersections) {
                if(const K2::Point_3* p = boost::get<K2::Point_3>(&(item.first))){
                    if(*p != tri_center){
                        intersection_v.push_back(*p);
                    }
                }
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
            intersection_v.resize(unique(intersection_v.begin(),intersection_v.end())-intersection_v.begin());
            if(intersection_v.size() % 2 == 0) { // 调整方向向外
                swap(renumber_bound_face_id[i][1],renumber_bound_face_id[i][2]);
            }
        }
        renumber_bound_face_global_id.resize(renumber_bound_face_id.size());
    }
    vector<RoughVertex> extend_vertex;
    unordered_set<RoughVertex,rough_vertex_hash,rough_vertex_equal>rough_vertex_set;
    bool done;
    MeshKernel::iGameFaceHandle fh;
    CoverageField(MeshKernel::iGameFaceHandle fh){
        this->fh = fh;
        rebuild();
    }

    void rebuild() {
        done = true;
        bound_face_vertex_rough.clear();
        bound_face_vertex_exact.clear();
        bound_face_vertex_inexact.clear();
        rough_vertex_set.clear();
        for(int i=0;i<3;i++){
           bound_face_vertex_rough.push_back(mesh_vertex_rough[mesh->fast_iGameFace[fh].vh(i)]);
        }
        for(int i=0;i<3;i++){
            for(auto v :field_move_vertices_rough[mesh->fast_iGameFace[fh].vh(i)]){
                bound_face_vertex_rough.push_back(v);
            }
        }


        unordered_map<RoughVertex,int,rough_vertex_hash,rough_vertex_equal> mp;
        for(int i=0;i<bound_face_vertex_rough.size();i++){

            bound_face_vertex_exact.emplace_back(bound_face_vertex_rough[i].exact());
            bound_face_vertex_inexact.emplace_back(CGAL::to_double(bound_face_vertex_exact[i].x()),
                                                   CGAL::to_double(bound_face_vertex_exact[i].y()),
                                                   CGAL::to_double(bound_face_vertex_exact[i].z())
                                                   );
            mp[bound_face_vertex_rough[i]] = i;
        }
        bbox_min = bbox_max = Point_K_to_iGameVertex(bound_face_vertex_inexact[0]);


        int base_vertex_size = bound_face_vertex_rough.size();
        for(int i=0;i<extend_vertex.size();i++){
            bound_face_vertex_rough.push_back(extend_vertex[i]);
            bound_face_vertex_exact.push_back(bound_face_vertex_rough.rbegin()->exact());
            bound_face_vertex_inexact.emplace_back(CGAL::to_double(bound_face_vertex_exact.rbegin()->x()),
                                      CGAL::to_double(bound_face_vertex_exact.rbegin()->y()),
                                      CGAL::to_double(bound_face_vertex_exact.rbegin()->z()));

            mp[extend_vertex[i]] = (int)bound_face_vertex_rough.size()-1;

        }
        for(int i=1;i<extend_vertex.size();i++){
            bbox_min.x() = min(bbox_min.x(),bound_face_vertex_inexact[i].x());
            bbox_min.y() = min(bbox_min.y(),bound_face_vertex_inexact[i].y());
            bbox_min.z() = min(bbox_min.z(),bound_face_vertex_inexact[i].z());
            bbox_max.x() = max(bbox_max.x(),bound_face_vertex_inexact[i].x());
            bbox_max.y() = max(bbox_max.y(),bound_face_vertex_inexact[i].y());
            bbox_max.z() = max(bbox_max.z(),bound_face_vertex_inexact[i].z());
        }

        Delaunay3D_K2 dt;
        dt.insert(bound_face_vertex_exact.begin(), bound_face_vertex_exact.end());
        std::vector<K2::Triangle_3> surface_triangles; // 这里存在一个精度问题
        for (auto fit = dt.finite_cells_begin(); fit != dt.finite_cells_end(); ++fit) {
            for (int i = 0; i < 4; ++i) {
                if (dt.is_infinite(fit->neighbor(i))) {
                    surface_triangles.push_back(dt.triangle(fit, i));
                }
            }
        }
        bound_face_id.clear();
        bound_face_useful.clear();
        extend_vertex.clear();

        K2::Vector_3 center_vec = {0,0,0};

        for (const auto& triangle : surface_triangles) {
            if(mp.count(transfer(triangle.vertex(0)))+ mp.count(transfer(triangle.vertex(1)))+mp.count(transfer(triangle.vertex(2)))!=3)
                cout <<"mdzz3"<<endl;
            int v0_id = mp[transfer(triangle.vertex(0))];
            int v1_id = mp[transfer(triangle.vertex(1))];
            int v2_id = mp[transfer(triangle.vertex(2))];
            if(triangle.is_degenerate())cout <<"mdzz"<<endl;
            rough_vertex_set.insert(bound_face_vertex_rough[v0_id]);
            rough_vertex_set.insert(bound_face_vertex_rough[v1_id]);
            rough_vertex_set.insert(bound_face_vertex_rough[v2_id]);
            bound_face_id.push_back({v0_id, v1_id, v2_id});
            if(K2::Triangle_3(bound_face_vertex_exact[v0_id],
                              bound_face_vertex_exact[v1_id],
                              bound_face_vertex_exact[v2_id]).is_degenerate()){
                cout <<"mdzz2"<<endl;
            }
            center_vec += (centroid(K2::Triangle_3(bound_face_vertex_exact[v0_id],
                                                   bound_face_vertex_exact[v1_id],
                                                   bound_face_vertex_exact[v2_id])) - K2::Point_3(0,0,0)) ;
        }


        for(int i=base_vertex_size;i<bound_face_vertex_rough.size();i++){
            if(rough_vertex_set.count(bound_face_vertex_rough[i])){
                extend_vertex.push_back(bound_face_vertex_rough[i]);
            }
        }
        if(!surface_triangles.size()){
            cout <<"okkk?"<<endl;
            for(auto i : bound_face_vertex_inexact){
                cout << "v" <<" "<< i <<endl;
                /*
                 * v 49.969 314.912 57.9381
v 49.969 322.912 52.9381
v 49.969 318.912 55.4381
v 53.969 314.912 57.9381
v 53.969 322.912 52.9381
v 53.969 318.912 55.4381
                 */
            }
            cout << bound_face_vertex_inexact.size() << endl;
            cout <<"okkk"<<endl;
            //exit(0); // 找这里的bug
        }
        center =  K2::Point_3(0,0,0) + (center_vec / surface_triangles.size());

        std::vector<std::vector<std::size_t> > faces_list;

        for (auto i: bound_face_id) {
            faces_list.push_back({std::size_t(i[0]), std::size_t(i[2]), std::size_t(i[1])});
            bound_face_useful.push_back(true);
        }
        bound_face_cross_field_list.resize(bound_face_id.size());
        bound_face_cutting_segment.resize(bound_face_id.size());
        bound_face_cutting_point.resize(bound_face_id.size());
        poly = new CGAL::Polyhedron_3<K2>();

        PMP::polygon_soup_to_polygon_mesh(bound_face_vertex_exact, faces_list, *poly, CGAL::parameters::all_default());
        inside_ptr = new CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2>(*poly);

    }



public:
    bool in_field(K2::Point_3 v) {
        //CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> inside(*poly);
        if ((*inside_ptr)(v) == CGAL::ON_BOUNDED_SIDE)
            return true;
        return false;
    }

    bool in_or_on_field(K2::Point_3 v) {
        //CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> inside(*poly);
        auto side = (*inside_ptr)(v);
        if (side== CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY)
            return true;
        return false;
    }

    int side(K2::Point_3 v) {
        //CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> inside(*poly);
        auto side = (*inside_ptr)(v);
        return side;
    }
    CGAL::Polyhedron_3<K2> * poly;
    vector<int>nearly_field;
    CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<K2>, K2> * inside_ptr;


};



vector<CoverageField> coverage_field_list;
#endif //THICKEN2_COVERAGEFIELD_H
