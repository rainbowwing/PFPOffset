//
// Created by rainbowwing on 2023/11/3.
//

#ifndef PFPOFFSET_DO_WINDING_NUMBER_H
#define PFPOFFSET_DO_WINDING_NUMBER_H

#include <utility>

#include "igl/winding_number.h"

void test_wining_num(){
    auto mesh0 = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../occ2/462k.input.obj"));
    auto mesh1 = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../occ2/462k.need_raoshu.obj"));
    std::vector<std::vector<double>> vertex;
    std::vector<std::vector<int>> face;
    for(int i=0;i<mesh0->VertexSize();i++){
        vertex.push_back({mesh0->fast_iGameVertex[i].x(),
                          mesh0->fast_iGameVertex[i].y(),
                          mesh0->fast_iGameVertex[i].z()}
                          );
    }
    for(int i=0;i<mesh0->FaceSize();i++){
        face.push_back({mesh0->fast_iGameFace[i].vh(0),
                          mesh0->fast_iGameFace[i].vh(1),
                          mesh0->fast_iGameFace[i].vh(2)}
                          );
    }
    //std::vector<std::vector<double>> vertex = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} };
    //std::vector<std::vector<int>> face = { {0, 1, 2} };

    // Convert to Eigen matrices:
    Eigen::MatrixXd V(vertex.size(), 3);
    Eigen::MatrixXi F(face.size(), 3);

    for (size_t i = 0; i < vertex.size(); ++i) {
        V(i, 0) = vertex[i][0];
        V(i, 1) = vertex[i][1];
        V(i, 2) = vertex[i][2];
    }

    for (size_t i = 0; i < face.size(); ++i) {
        F(i, 0) = face[i][0];
        F(i, 1) = face[i][1];
        F(i, 2) = face[i][2];
    }

    // Assuming the query points are in this structure:
    std::vector<std::vector<double>> query_points;

    for(int i=0;i<mesh1->FaceSize();i++){
        int vh0 = mesh1->fast_iGameFace[i].vh(0);
        int vh1 = mesh1->fast_iGameFace[i].vh(1);
        int vh2 = mesh1->fast_iGameFace[i].vh(2);
        MeshKernel::iGameVertex center = (mesh1->fast_iGameVertex[vh0] + mesh1->fast_iGameVertex[vh1] + mesh1->fast_iGameVertex[vh2]) / 3;
        query_points.push_back({center.x(), center.y(), center.z()});
    }

    // Convert query points to an Eigen matrix:
    Eigen::MatrixXd P(query_points.size(), 3);
//    std::ofstream debug_mesh_in("../occ2/304k.input2.winding_debugin.obj");
//    std::ofstream debug_mesh_out("../occ2/304k.input2.winding_debugout.obj");
    for (size_t i = 0; i < query_points.size(); ++i) {
        P(i, 0) = query_points[i][0];
        P(i, 1) = query_points[i][1];
        P(i, 2) = query_points[i][2];
    }

    // Compute the winding numbers:
    Eigen::VectorXd W;
    igl::winding_number(V, F, P, W);
    std::ofstream update_mesh_in("../occ2/462k.input.windingin.obj");
    std::ofstream update_mesh_out("../occ2/462k.input.windingout.obj");
    for(int i=0;i<mesh1->VertexSize();i++){
        update_mesh_in<<"v "<<mesh1->fast_iGameVertex[i].x()<<" "<<
                          mesh1->fast_iGameVertex[i].y()<<" "<<
                          mesh1->fast_iGameVertex[i].z()<<endl;
        update_mesh_out<<"v "<<mesh1->fast_iGameVertex[i].x()<<" "<<
                       mesh1->fast_iGameVertex[i].y()<<" "<<
                       mesh1->fast_iGameVertex[i].z()<<endl;
    }
    for(int i=0;i<mesh1->FaceSize();i++){
        std::cout << "Winding number for point " << i << ": " << W(i) << std::endl;
        if(W(i) >= 0.5){
            update_mesh_in<<"f "<<mesh1->fast_iGameFace[i].vh(0) +1
            <<" "<<mesh1->fast_iGameFace[i].vh(1) +1
            <<" "<<mesh1->fast_iGameFace[i].vh(2) +1 <<endl;
        }
        else{
            update_mesh_out<<"f "<<mesh1->fast_iGameFace[i].vh(0) +1
                           <<" "<<mesh1->fast_iGameFace[i].vh(1) +1
                           <<" "<<mesh1->fast_iGameFace[i].vh(2) +1 <<endl;
        }
    }
    update_mesh_in.close();
    update_mesh_out.close();
    // Output the results:
//    for (int i = 0; i < W.size(); ++i) {
//        if()
//        std::cout << "Winding number for point " << i << ": " << W(i) << std::endl;
//    }
    exit(0);
}

struct FinalFace{
    int f0,f1,f2;
    K2::Point_3 center;
    bool flag;
    FinalFace(){}
    FinalFace(int f0,int f1,int f2,K2::Point_3 center,bool flag){
        this->f0 = f0;
        this->f1 = f1;
        this->f2 = f2;
        this->center = std::move(center);
        this->flag = flag;
    }
    friend bool operator < (const FinalFace &a, const FinalFace& b){
        vector<int>va{a.f0,a.f1,a.f2};
        vector<int>vb{b.f0,b.f1,b.f2};
        sort(va.begin(),va.end());
        sort(vb.begin(),vb.end());
        if(va[0]!=vb[0])return va[0]<vb[0];
        else if(va[1]!=vb[1])return va[1]<vb[1];
        else return va[2]<vb[2];
    }
};




void winding_num(vector<FinalFace>& v){
    cout <<"st winding_num" << v.size() << endl;
    std::vector<std::vector<double>> vertex;
    std::vector<std::vector<int>> face;
    for(int i=0;i<mesh->VertexSize();i++){
        vertex.push_back({mesh->fast_iGameVertex[i].x(),
                          mesh->fast_iGameVertex[i].y(),
                          mesh->fast_iGameVertex[i].z()}
        );
    }
    for(int i=0;i<mesh->FaceSize();i++){
        face.push_back({mesh->fast_iGameFace[i].vh(0),
                        mesh->fast_iGameFace[i].vh(1),
                        mesh->fast_iGameFace[i].vh(2)}
        );
    }
    //std::vector<std::vector<double>> vertex = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} };
    //std::vector<std::vector<int>> face = { {0, 1, 2} };

    // Convert to Eigen matrices:
    Eigen::MatrixXd V(vertex.size(), 3);
    Eigen::MatrixXi F(face.size(), 3);

    for (size_t i = 0; i < vertex.size(); ++i) {
        V(i, 0) = vertex[i][0];
        V(i, 1) = vertex[i][1];
        V(i, 2) = vertex[i][2];
    }

    for (size_t i = 0; i < face.size(); ++i) {
        F(i, 0) = face[i][0];
        F(i, 1) = face[i][1];
        F(i, 2) = face[i][2];
    }

    // Assuming the query points are in this structure:
    std::vector<std::vector<double>> query_points;

    for(int i=0;i<v.size();i++){
        query_points.push_back({CGAL::to_double(v[i].center.x()),
                                CGAL::to_double(v[i].center.y()),
                                CGAL::to_double(v[i].center.z())});
    }
    Eigen::MatrixXd P(query_points.size(), 3);
    for (size_t i = 0; i < query_points.size(); ++i) {
        P(i, 0) = query_points[i][0];
        P(i, 1) = query_points[i][1];
        P(i, 2) = query_points[i][2];
    }
    Eigen::VectorXd W;
    igl::winding_number(V, F, P, W);
    for(int i=0;i< v.size(); i++){
        if( (running_mode == 1 && W(i) <= 0.5) || (running_mode == 2 && W(i) >= 0.5)){
            v[i].flag = true;
        }
        else{
            v[i].flag = false;
        }
        cout <<"v[i].flag: "<< v[i].flag << endl;
    }

}


vector<int> winding_num(vector<K2::Point_3>& v){
    vector<int>ret;
    ret.resize(v.size());
    cout <<"st winding_num" << v.size() << endl;
    std::vector<std::vector<double>> vertex;
    std::vector<std::vector<int>> face;
    for(int i=0;i<mesh->VertexSize();i++){
        vertex.push_back({mesh->fast_iGameVertex[i].x(),
                          mesh->fast_iGameVertex[i].y(),
                          mesh->fast_iGameVertex[i].z()}
        );
    }
    for(int i=0;i<mesh->FaceSize();i++){
        face.push_back({mesh->fast_iGameFace[i].vh(0),
                        mesh->fast_iGameFace[i].vh(1),
                        mesh->fast_iGameFace[i].vh(2)}
        );
    }
    //std::vector<std::vector<double>> vertex = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} };
    //std::vector<std::vector<int>> face = { {0, 1, 2} };

    // Convert to Eigen matrices:
    Eigen::MatrixXd V(vertex.size(), 3);
    Eigen::MatrixXi F(face.size(), 3);

    for (size_t i = 0; i < vertex.size(); ++i) {
        V(i, 0) = vertex[i][0];
        V(i, 1) = vertex[i][1];
        V(i, 2) = vertex[i][2];
    }

    for (size_t i = 0; i < face.size(); ++i) {
        F(i, 0) = face[i][0];
        F(i, 1) = face[i][1];
        F(i, 2) = face[i][2];
    }

    // Assuming the query points are in this structure:
    std::vector<std::vector<double> > query_points;

    for(int i=0;i<v.size();i++){
        query_points.push_back({CGAL::to_double(v[i].x()),
                                CGAL::to_double(v[i].y()),
                                CGAL::to_double(v[i].z())});
    }
    Eigen::MatrixXd P(query_points.size(), 3);
    for (size_t i = 0; i < query_points.size(); ++i) {
        P(i, 0) = query_points[i][0];
        P(i, 1) = query_points[i][1];
        P(i, 2) = query_points[i][2];
    }
    Eigen::VectorXd W;
    igl::winding_number(V, F, P, W);
    for(int i=0;i< v.size(); i++){
        if( (running_mode == 1 && W(i) <= 0.5) || (running_mode == 2 && W(i) >= 0.5)){
            ret[i] = 1;
        }
        else{
            ret[i] = 0;
        }
        //cout <<"v[i].flag: "<< v[i].flag << endl;
    }
    return ret;
}



#endif //PFPOFFSET_DO_WINDING_NUMBER_H
