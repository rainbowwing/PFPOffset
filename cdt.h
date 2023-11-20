//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_CDT_H
#define THICKEN2_CDT_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#include <vector>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Fuzzy_iso_box.h>
typedef CGAL::Arr_segment_traits_2<K2> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

typedef CGAL::Constrained_Delaunay_triangulation_2<K2> CDT;
typedef CGAL::Polygon_2<K2> Polygon_2;


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3 Point;
typedef CGAL::Search_traits_3<K> STraits;
typedef CGAL::Fuzzy_sphere<STraits> Fuzzy_circle;
typedef CGAL::Kd_tree<STraits> Kd_tree;


typedef CGAL::Search_traits_3<K2> STraitsK2;
typedef CGAL::Fuzzy_sphere<STraitsK2> Fuzzy_circle_K2;
typedef CGAL::Kd_tree<STraitsK2> Kd_tree_K2;
typedef CGAL::Fuzzy_iso_box<STraitsK2> FuzzyBoxK2;

typedef CGAL::Delaunay_triangulation_3<K> Delaunay3D;
typedef CGAL::Delaunay_triangulation_3<K2> Delaunay3DK2;

typedef CGAL::Constrained_Delaunay_triangulation_2<K2> CDT;



vector<vector<K2::Point_3> > CGAL_CDT_NEW(vector<K2::Point_3> sorted_bound_vertex, vector<K2::Segment_3> cs,K2::Triangle_3 origin_face) {

    //CGAL::make_conforming_Delaunay_2();
    // cout << sorted_bound_vertex.size() <<" "<< cs.size() << endl;
    K2::Point_3 base_point_3d = origin_face.vertex(0);
    CDT cdt;// 为什么不使用射线法投影回去？

    K2::Vector_3 origin_face_v0 = origin_face.vertex(1) - origin_face.vertex(0);
    K2::Vector_3 origin_face_v1 = origin_face.vertex(2) - origin_face.vertex(0);
    K2::Plane_3 pla = origin_face.supporting_plane();
    K2::Vector_3 X_axis = pla.base1();
    K2::Vector_3 Y_axis = pla.base2();

    std::vector<K2::Point_2> polygon_vertices;
    Arrangement_2 arrangement;
    for(int i=0;i<sorted_bound_vertex.size();i++){
        CGAL::Epeck::FT x = (sorted_bound_vertex[i]-base_point_3d) * X_axis;
        CGAL::Epeck::FT y = (sorted_bound_vertex[i]-base_point_3d) * Y_axis;
        //CGAL::insert(arrangement,K2::Point_2(x,y));
        cdt.insert(K2::Point_2(x,y));
        polygon_vertices.emplace_back(x,y);
    }

    Polygon_2 polygon(polygon_vertices.begin(), polygon_vertices.end());

    std::vector<K2::Point_2> bounded_polygon_vertices;
    for(int i=0;i<3;i++){
        CGAL::Epeck::FT x = (origin_face.vertex(i)-base_point_3d) * X_axis;
        CGAL::Epeck::FT y = (origin_face.vertex(i)-base_point_3d) * Y_axis;
        bounded_polygon_vertices.emplace_back(x,y);
    }
    K2::Triangle_2 origin_tri_2d(bounded_polygon_vertices[0],
                    bounded_polygon_vertices[1],
                    bounded_polygon_vertices[2]
                    );
    Polygon_2 bounded_polygon(bounded_polygon_vertices.begin(),bounded_polygon_vertices.end());


    for(int i=0;i<3;i++){
        CGAL::Epeck::FT x0 = (origin_face.vertex(i)-base_point_3d) * X_axis;
        CGAL::Epeck::FT y0 = (origin_face.vertex(i)-base_point_3d) * Y_axis;
        CGAL::Epeck::FT x1 = (origin_face.vertex((i+1)%3)-base_point_3d) * X_axis;
        CGAL::Epeck::FT y1 = (origin_face.vertex((i+1)%3)-base_point_3d) * Y_axis;
        K2::Segment_2 seg(K2::Point_2(x0,y0) ,K2::Point_2(x1,y1));
        CGAL::insert(arrangement,seg);
        //cdt.insert_constraint(K2::Point_2(x0,y0) ,K2::Point_2(x1,y1));
    }
    for(auto s: cs){
        CGAL::Epeck::FT x0 = (s.vertex(0)-base_point_3d) * X_axis;
        CGAL::Epeck::FT y0 = (s.vertex(0)-base_point_3d) * Y_axis;
        CGAL::Epeck::FT x1 = (s.vertex(1)-base_point_3d) * X_axis;
        CGAL::Epeck::FT y1 = (s.vertex(1)-base_point_3d) * Y_axis;
        K2::Segment_2 seg(K2::Point_2(x0,y0) ,K2::Point_2(x1,y1));
        CGAL::insert(arrangement,seg);
        //cdt.insert_constraint(K2::Point_2(x0,y0) ,K2::Point_2(x1,y1));
    }
    for (auto eit = arrangement.edges_begin(); eit != arrangement.edges_end(); ++eit) {
        K2::Segment_2 seg = eit->curve();
        K2::Point_2 source = seg.source();
        K2::Point_2 target = seg.target();
        cdt.insert_constraint(source,target);
//        std::cout << "Segment: ((" << source.x() << ", " << source.y() << "), ("
//                  << target.x() << ", " << target.y() << "))\n";
    }


    K2::Vector_3 vec = cross_product(origin_face.vertex(1)-origin_face.vertex(0),
                                     origin_face.vertex(2)-origin_face.vertex(0)
    );

    vector<vector<K2::Point_3> > ret;

    for(int i=0;i<3;i++){
        K2::Point_3 newp = base_point_3d + origin_tri_2d.vertex(i).x()*X_axis +  origin_tri_2d.vertex(i).y()*Y_axis;
        if(CGAL::squared_distance(newp,origin_face)!=CGAL::Epeck::FT(0) ){
            cout <<i<<"newp_dist_err"<< CGAL::squared_distance(newp,origin_face)<<endl;
        }
    }

    exit(0);

    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        K2::Triangle_2 t = cdt.triangle(fit);
        if(bounded_polygon.has_on_bounded_side(CGAL::centroid(t))) {
            vector<K2::Point_3>tmp(3);
            for(int j=0;j<3;j++){
                tmp[j] = base_point_3d + t[j].x()*X_axis + t[j].y()*Y_axis;
            }
            K2::Vector_3 vec1 = cross_product(tmp[1] - tmp[0],tmp[2] - tmp[0]);
            if(CGAL::squared_distance(tmp[0],origin_face) !=CGAL::Epeck::FT(0)){
                cout  <<"REMAPERROR2D "<< CGAL::squared_distance( t[0],origin_tri_2d) << endl;
                cout  <<"REMAPERROR "<< CGAL::squared_distance(tmp[0],origin_face) << endl;
                cout  <<"REMAPERRORPLANE "<< CGAL::squared_distance(tmp[0],origin_face.supporting_plane()) << endl;
            }
            else
                cout <<"ccok"<<endl;
            if(CGAL::squared_distance(tmp[1],origin_face) !=CGAL::Epeck::FT(0)){
                cout  <<"REMAPERROR2D "<< CGAL::squared_distance( t[1],origin_tri_2d) << endl;
                cout  <<"REMAPERROR "<< CGAL::squared_distance(tmp[1],origin_face) << endl;
                cout  <<"REMAPERRORPLANE "<< CGAL::squared_distance(tmp[1],origin_face.supporting_plane()) << endl;
            }
            else
                cout <<"ccok"<<endl;
            if(CGAL::squared_distance(tmp[2],origin_face) !=CGAL::Epeck::FT(0)){
                cout  <<"REMAPERROR2D "<< CGAL::squared_distance( t[2],origin_tri_2d) << endl;
                cout  <<"REMAPERROR "<< CGAL::squared_distance(tmp[2],origin_face) << endl;
                cout  <<"REMAPERRORPLANE "<< CGAL::squared_distance(tmp[2],origin_face.supporting_plane()) << endl;
            }
            else
                cout <<"ccok"<<endl;


            if(vec * vec1 < CGAL::Epeck::FT(0))
                swap(tmp[1],tmp[2]);

            ret.push_back(tmp);
        }
    }


    return ret;
}



int get_project_plane(const K2::Triangle_3& origin_face){
    K2::Vector_3 normal = CGAL::cross_product(origin_face.vertex(1) - origin_face.vertex(0), origin_face.vertex(2) - origin_face.vertex(0));
    if (normal.x()*normal.x() >= normal.y()*normal.y() && normal.x()*normal.x() >= normal.z()*normal.z()) {
        return 0;
    }
    else if (normal.y()*normal.y() >= normal.x()*normal.x() && normal.y()*normal.y() >= normal.z()*normal.z()) {
        return 1;
    }
    else
        return 2;
}

K2::Point_2 down_to_2d(const K2::Point_3& p,int project){
    if(project == 0){
        return K2::Point_2(p.y(),p.z());
    }
    else if(project == 1){
        return K2::Point_2(p.z(),p.x());
    }
    else
        return K2::Point_2(p.x(),p.y());

}
K2::Point_3 lift_to_3d(const K2::Point_2& p,int project,K2::Triangle_3 tri){
    if(project == 0){
        K2::Point_3 p0(0,p.x(),p.y());
        K2::Point_3 p1(1,p.x(),p.y());
        K2::Line_3 line(p0,p1);
        CGAL::cpp11::result_of<K2::Intersect_3(K2::Line_3 , K2::Triangle_3)>::type
                res_lt = intersection(line, tri);
        if (res_lt) {
            if (const K2::Point_3 *x = boost::get<K2::Point_3>(&*res_lt)) {
                return *x;
            }
        }
        cout <<"error"<<endl;
        exit(0);
        return K2::Point_3(0,0,0);
    }
    else if(project == 1){
        K2::Point_3 p0(p.y(),0,p.x());
        K2::Point_3 p1(p.y(),1,p.x());
        K2::Line_3 line(p0,p1);
        CGAL::cpp11::result_of<K2::Intersect_3(K2::Line_3 , K2::Triangle_3)>::type
                res_lt = intersection(line, tri);
        if (res_lt) {
            if (const K2::Point_3 *x = boost::get<K2::Point_3>(&*res_lt)) {
                return *x;
            }
        }
        cout <<"error"<<endl;
        exit(0);
        return K2::Point_3(0,0,0);
    }
    else {
        K2::Point_3 p0(p.x(),p.y(),0);
        K2::Point_3 p1(p.x(),p.y(),1);
        K2::Line_3 line(p0,p1);
        CGAL::cpp11::result_of<K2::Intersect_3(K2::Line_3 , K2::Triangle_3)>::type
                res_lt = intersection(line, tri);
        if (res_lt) {
            if (const K2::Point_3 *x = boost::get<K2::Point_3>(&*res_lt)) {
                return *x;
            }
        }
        cout <<"error"<<endl;
        exit(0);
        return K2::Point_3(0,0,0);
    }
}

vector<vector<K2::Point_3> > CGAL_CDT_NEW2(vector<K2::Point_3> sorted_bound_vertex, vector<K2::Segment_3> cs,K2::Triangle_3 origin_face) {

    //CGAL::make_conforming_Delaunay_2();
    // cout << sorted_bound_vertex.size() <<" "<< cs.size() << endl;
    int project_type = get_project_plane(origin_face);

    CDT cdt;// 为什么不使用射线法投影回去？

    K2::Plane_3 pla = origin_face.supporting_plane();

   // cout <<"cs.size()" <<cs.size() << endl;
   // cout << "origin_face"<<origin_face <<" "<<project_type<<" "<<origin_face.supporting_plane().orthogonal_vector()<< endl;



    std::vector<K2::Point_2> polygon_vertices;
    Arrangement_2 arrangement;
    for(int i=0;i<sorted_bound_vertex.size();i++){
        K2::Point_2 to2d = down_to_2d(sorted_bound_vertex[i],project_type);
        cdt.insert(to2d);
    }
    for(int i=0;i<3;i++){
        K2::Point_2 to2d = down_to_2d(origin_face.vertex(i),project_type);
        polygon_vertices.push_back(to2d);
    }

    Polygon_2 bounded_polygon(polygon_vertices.begin(),polygon_vertices.end());

    //cout <<"succ st"<<endl;
    for(int i=0;i<3;i++){
        K2::Segment_2 seg(polygon_vertices[i] ,polygon_vertices[(i+1)%3]);
        CGAL::insert(arrangement,seg);
        //cdt.insert_constraint(K2::Point_2(x0,y0) ,K2::Point_2(x1,y1));
    }

    for(auto s: cs){

        K2::Segment_2 seg(down_to_2d(s.vertex(0),project_type),down_to_2d(s.vertex(1),project_type));
        CGAL::insert(arrangement,seg);
        //cdt.insert_constraint(K2::Point_2(x0,y0) ,K2::Point_2(x1,y1));
    }
    //cout <<"succ end"<<endl;
    for (auto eit = arrangement.edges_begin(); eit != arrangement.edges_end(); ++eit) {
        K2::Segment_2 seg = eit->curve();
        K2::Point_2 source = seg.source();
        K2::Point_2 target = seg.target();
        cdt.insert_constraint(source,target);
//        std::cout << "Segment: ((" << source.x() << ", " << source.y() << "), ("
//                  << target.x() << ", " << target.y() << "))\n";
    }


    K2::Vector_3 vec = cross_product(origin_face.vertex(1)-origin_face.vertex(0),
                                     origin_face.vertex(2)-origin_face.vertex(0)
    );

    vector<vector<K2::Point_3> > ret;

    //cout <<"project_type:"<< project_type << endl;
//    K2::Point_2 p0_2d = down_to_2d(origin_face.vertex(0),project_type);
//    K2::Point_2 p1_2d = down_to_2d(origin_face.vertex(1),project_type);
//    K2::Point_2 p2_2d = down_to_2d(origin_face.vertex(2),project_type);
//
//    K2::Point_3 p0_l3d = lift_to_3d(p0_2d,project_type,origin_face);
//    K2::Point_3 p1_l3d = lift_to_3d(p1_2d,project_type,origin_face);
//    K2::Point_3 p2_l3d = lift_to_3d(p2_2d,project_type,origin_face);
//    cout <<p0_l3d<<" dist " << CGAL::squared_distance(p0_l3d,origin_face) << endl;
//    cout <<p1_l3d<<" dist " << CGAL::squared_distance(p1_l3d,origin_face) << endl;
//    cout <<p2_l3d<<" dist " << CGAL::squared_distance(p2_l3d,origin_face) << endl;
//    cout <<(p0_l3d == origin_face.vertex(0))<<endl;
//    cout <<(p1_l3d == origin_face.vertex(1))<<endl;
//    cout <<(p2_l3d == origin_face.vertex(2))<<endl;
//    for(int i=0;i<3;i++){
//        K2::Point_3 newp = base_point_3d + origin_tri_2d.vertex(i).x()*X_axis +  origin_tri_2d.vertex(i).y()*Y_axis;
//        if(CGAL::squared_distance(newp,origin_face)!=CGAL::Epeck::FT(0) ){
//            cout <<i<<"newp_dist_err"<< CGAL::squared_distance(newp,origin_face)<<endl;
//        }
//    }

   // exit(0);

    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        K2::Triangle_2 t = cdt.triangle(fit);
        if(bounded_polygon.has_on_bounded_side(CGAL::centroid(t))) {
            vector<K2::Point_3>tmp(3);
            for(int j=0;j<3;j++){
                tmp[j] = lift_to_3d(t[j],project_type,origin_face);
            }
            K2::Vector_3 vec1 = cross_product(tmp[1] - tmp[0],tmp[2] - tmp[0]);
//            if(CGAL::squared_distance(tmp[0],origin_face) !=CGAL::Epeck::FT(0)){
//                cout  <<"REMAPERROR2D "<< CGAL::squared_distance( t[0],origin_tri_2d) << endl;
//                cout  <<"REMAPERROR "<< CGAL::squared_distance(tmp[0],origin_face) << endl;
//                cout  <<"REMAPERRORPLANE "<< CGAL::squared_distance(tmp[0],origin_face.supporting_plane()) << endl;
//            }
//            else
//                cout <<"ccok"<<endl;
//            if(CGAL::squared_distance(tmp[1],origin_face) !=CGAL::Epeck::FT(0)){
//                cout  <<"REMAPERROR2D "<< CGAL::squared_distance( t[1],origin_tri_2d) << endl;
//                cout  <<"REMAPERROR "<< CGAL::squared_distance(tmp[1],origin_face) << endl;
//                cout  <<"REMAPERRORPLANE "<< CGAL::squared_distance(tmp[1],origin_face.supporting_plane()) << endl;
//            }
//            else
//                cout <<"ccok"<<endl;
//            if(CGAL::squared_distance(tmp[2],origin_face) !=CGAL::Epeck::FT(0)){
//                cout  <<"REMAPERROR2D "<< CGAL::squared_distance( t[2],origin_tri_2d) << endl;
//                cout  <<"REMAPERROR "<< CGAL::squared_distance(tmp[2],origin_face) << endl;
//                cout  <<"REMAPERRORPLANE "<< CGAL::squared_distance(tmp[2],origin_face.supporting_plane()) << endl;
//            }
//            else
//                cout <<"ccok"<<endl;
//

            if(vec * vec1 < CGAL::Epeck::FT(0))
                swap(tmp[1],tmp[2]);

            ret.push_back(tmp);
        }
    }


    return ret;
}



#endif //THICKEN2_CDT_H
