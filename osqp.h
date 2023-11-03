//
// Created by rainbowwing on 2023/8/25.
//
#include <CGAL/Gmpq.h>

#ifndef THICKEN2_OSQP_H
#define THICKEN2_OSQP_H
int running_mode = 2;
MeshKernel::iGameVertex do_quadratic_error_metric_check(MeshKernel::iGameVertexHandle vh,vector<MeshKernel::iGameFaceHandle> neighbor_face_list,bool &is_succ,double& dist){
    MeshKernel::iGameVertex v = mesh->fast_iGameVertex[vh];
//    if(v.dist(MeshKernel::iGameVertex(50753.437500,-9133.144531,38876.300781)) <10){
//        cout <<"debug:" << vh << " "<<t << endl;
//    }
    int m = neighbor_face_list.size();
    if(m == 1){
        MeshKernel::iGameFaceHandle f = *neighbor_face_list.begin();
        K2::Vector_3 normal_k2 = K2::Triangle_3 (K2::Point_3(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])),
                                                 K2::Point_3(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)])),
                                                 K2::Point_3(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]))
        ).supporting_plane().orthogonal_vector();

        normal_k2 /= CGAL::approximate_sqrt(normal_k2.squared_length());

        MeshKernel::iGameVertex normal = MeshKernel::iGameVertex(CGAL::to_double(normal_k2.x()),
                                                                 CGAL::to_double(normal_k2.y()),
                                                                 CGAL::to_double(normal_k2.z()));
        is_succ = true;

        return v + normal * mesh->fast_iGameFace[f].move_dist;
    }


    Eigen::SparseMatrix<double> hessian(3, 3);     //P: n*n正定矩阵,必须为稀疏矩阵SparseMatrix
    hessian.setZero();
    Eigen::VectorXd gradient(3);                  //Q: n*1向量
    gradient.setZero();
    Eigen::SparseMatrix<double> linearMatrix(m, 3); //A: m*n矩阵,必须为稀疏矩阵SparseMatrix
    linearMatrix.setZero();
    Eigen::VectorXd lowerBound(m);                  //L: m*1下限向量
    lowerBound.setZero();
    Eigen::VectorXd upperBound(m);                  //U: m*1上限向量
    upperBound.setZero();
    int cnt = 0;
    dist = 0;

    double avg_move_dist=0;
    double max_move_dist = 0;

    for (auto f : neighbor_face_list){
        avg_move_dist += mesh->fast_iGameFace[f].move_dist;
        max_move_dist += mesh->fast_iGameFace[f].move_dist;
    }


    for (auto f : neighbor_face_list) {
        K2::Vector_3 normal_k2 = K2::Triangle_3 (K2::Point_3(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])),
                                                 K2::Point_3(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)])),
                                                 K2::Point_3(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]))
        ).supporting_plane().orthogonal_vector();




        normal_k2 /= CGAL::approximate_sqrt(normal_k2.squared_length());


        MeshKernel::iGameVertex normal = MeshKernel::iGameVertex(CGAL::to_double(normal_k2.x()),
                                                                 CGAL::to_double(normal_k2.y()),
                                                                 CGAL::to_double(normal_k2.z()));


        //cout <<"normal "<< normal.x() <<" "<< normal.y() <<" "<< normal.z() << endl;
        if (running_mode == 2)
            normal = normal * -1;
//        MeshKernel::iGameVertex new_v = v + normal * mesh->fast_iGameFace[f].move_dist;
//
//        MeshKernel::iGameVertex move_max_v = v + normal * mesh->fast_iGameFace[f].move_dist*1.35;
//        MeshKernel::iGameVertex move_min_v = v + normal * mesh->fast_iGameFace[f].move_dist*0.8;



        // MeshKernel::iGameVertex move_max_v = v + normal * avg_move_dist * (1.25+0.1*depth);
        double move_max_v =   mesh->fast_iGameFace[f].move_dist * max_distance_limit;
        double move_min_v =   mesh->fast_iGameFace[f].move_dist * min_distance_limit;


        Eigen::Vector3d p(normal.x(),  normal.y(), normal.z());
        Eigen::Matrix3d A = p * p.transpose();
        Eigen::Vector3d LimitV(normal.x(),normal.y(),normal.z());

        double lower =  move_min_v + (v.x() * normal.x()) + (v.y() * normal.y()) + (v.z() * normal.z());//normal.x()*move_min_v.x() + normal.y()*move_min_v.y() + normal.z()*move_min_v.z();
        double upper =  move_max_v + (v.x() * normal.x()) + (v.y() * normal.y()) + (v.z() * normal.z());//normal.x()*move_max_v.x() + normal.y()*move_max_v.y() + normal.z()*move_max_v.z();



//        for(int i=0;i<3;i++)
//            for(int j=0;j<3;j++)
//                hessian.coeffRef(i,j)+=A.coeff(i,j)*2;
//        for(int i=0;i<3;i++)
//            gradient.coeffRef(i) +=  D2.coeffRef(i);


        for(int i=0;i<3;i++)
            linearMatrix.coeffRef(cnt,i) = LimitV[i];
        lowerBound.coeffRef(cnt) = lower;
        upperBound.coeffRef(cnt) = upper;

        cnt++;
    }

    hessian.coeffRef(0,0) += (2.0);
    hessian.coeffRef(1,1) += (2.0);
    hessian.coeffRef(2,2) += (2.0);

    gradient.coeffRef(0) -= (2.0) * v.x();
    gradient.coeffRef(1) -= (2.0) * v.y();
    gradient.coeffRef(2) -= (2.0) * v.z();

    OsqpEigen::Solver solver;
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables(3);   //变量数n
    solver.data()->setNumberOfConstraints(m); //约束数m
    solver.data()->setHessianMatrix(hessian);
    solver.data()->setGradient(gradient);
    solver.data()->setLinearConstraintsMatrix(linearMatrix);
    solver.data()->setLowerBound(lowerBound);
    solver.data()->setUpperBound(upperBound);
    solver.initSolver();


    if(solver.solve()) {
        // if 长度过短
        is_succ = true;
        Eigen::VectorXd QPSolution = solver.getSolution();
        //dist
        MeshKernel::iGameVertex ret {QPSolution.coeffRef(0), QPSolution.coeffRef(1), QPSolution.coeffRef(2)};
        dist = 0;


        dist += abs((ret - v) .norm());
        return ret;
    }
    else{
        is_succ = false;
        return {0,0,0};
    }
}


#endif //THICKEN2_OSQP_H
