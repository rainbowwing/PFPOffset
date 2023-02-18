//
// Created by rainbowwing on 2023/2/18.
//

#ifndef THICKEN2OUT_DO_QUADRATIC_ERROR_METRIC_H
#define THICKEN2OUT_DO_QUADRATIC_ERROR_METRIC_H

vector<MeshKernel::iGameVertex> min_move_g;
vector<MeshKernel::iGameVertex> max_move_g;

MeshKernel::iGameVertex do_quadratic_error_metric(shared_ptr <MeshKernel::SurfaceMesh> mesh,MeshKernel::iGameVertexHandle vh,bool &is_succ,int depth=0){

    MeshKernel::iGameVertex v = mesh->fast_iGameVertex[vh];
    int m = mesh->FastNeighborFhOfVertex_[vh].size();
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

    double avg_move_dist=0;
    double max_move_dist = 0;
    MeshKernel::iGameVertex avg_move_vertex(0,0,0);
    for (auto f : mesh->FastNeighborFhOfVertex_[vh]){
        avg_move_dist += mesh->fast_iGameFace[f].move_dist;
        max_move_dist = max(max_move_dist, mesh->fast_iGameFace[f].move_dist);

    }
    avg_move_dist /= (1.0*mesh->FastNeighborFhOfVertex_[vh].size());

    for (auto f : mesh->FastNeighborFhOfVertex_[vh]) {
        MeshKernel::iGameVertex normal
                = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]
                    - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]) %
                   (mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]
                    - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])).normalize();


//        MeshKernel::iGameVertex new_v = v + normal * mesh->fast_iGameFace[f].move_dist;
//
//        MeshKernel::iGameVertex move_max_v = v + normal * mesh->fast_iGameFace[f].move_dist*1.35;
//        MeshKernel::iGameVertex move_min_v = v + normal * mesh->fast_iGameFace[f].move_dist*0.8;
#ifdef MODEIN
        normal = normal * -1;
#endif

        MeshKernel::iGameVertex new_v = v + normal * avg_move_dist;
        avg_move_vertex += normal * avg_move_dist;

        MeshKernel::iGameVertex move_max_v = v + normal * avg_move_dist * (1.35+0.15*depth);
        MeshKernel::iGameVertex move_min_v = v + normal * avg_move_dist * (0.80-0.15*depth);

        double d = -(normal.x() * new_v.x() + normal.y() * new_v.y() +  normal.z() * new_v.z());

        Eigen::Vector3d p(normal.x(),  normal.y(), normal.z());
        Eigen::Matrix3d A = p * p.transpose();
        Eigen::Vector3d D2(2*d*normal.x(),2*d*normal.y(),2*d*normal.z());
        Eigen::Vector3d LimitV(normal.x(),normal.y(),normal.z());

        double lower = normal.x()*move_min_v.x() + normal.y()*move_min_v.y() + normal.z()*move_min_v.z();
        double upper = normal.x()*move_max_v.x() + normal.y()*move_max_v.y() + normal.z()*move_max_v.z();



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
    avg_move_vertex/= mesh->FastNeighborFhOfVertex_[vh].size();

    min_move_g[vh] =  v + avg_move_vertex;

    avg_move_vertex = v + avg_move_vertex;;
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
    // return avg_move_vertex;
    if(solver.solve()) {
        is_succ = true;
        Eigen::VectorXd QPSolution = solver.getSolution();
        MeshKernel::iGameVertex res(QPSolution.coeffRef(0), QPSolution.coeffRef(1), QPSolution.coeffRef(2));
        if( (res-v).norm() > 2*avg_move_dist){
            return v+ (res-v)/(res-v).norm()*(2*avg_move_dist);
            // return avg_move_vertex;
        }

        return {QPSolution.coeffRef(0), QPSolution.coeffRef(1), QPSolution.coeffRef(2)};
    }
    else{
        is_succ = false;
        return do_quadratic_error_metric(mesh,vh,is_succ,depth+1);
        /* cout <<"v "<< avg_move_vertex.x()<<" "<<avg_move_vertex.y()<<" "<<avg_move_vertex.z() << endl;
         puts("use avg_move_vertex instead");
         return avg_move_vertex;*/
    }
}

#endif //THICKEN2OUT_DO_QUADRATIC_ERROR_METRIC_H
