//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_DP_H
#define THICKEN2_DP_H
#include <random>
double avg_edge_limit = 1e-9;
vector<int> get_sub_state(int state,int face_list_size){
    vector<int> ret;
    function<void(int,int,int)> dfs = [&](int state,int pos,int change){
        if(pos == -1 && change) {
            ret.push_back(state);
            return;
        }
        for(int i=pos;i>=0;i--){
            if(state & (1<<i)){
                dfs(state,i-1,change);
                dfs(state ^ (1<<i),i-1,1);
                break;
            }
            else if(i == 0){
                dfs(state,i-1,change);
                break;
            }
        }
    };
    dfs(state,face_list_size,0);
    return ret;
}
vector<MeshKernel::iGameVertex> run(MeshKernel::iGameVertexHandle vh,vector<MeshKernel::iGameFaceHandle> neighbor_face_list){
//    if(neighbor_face_list.size()>=21) {
//        cout << "neighbor_face_list.size()=" <<neighbor_face_list.size() <<">21 : this mesh please do remeshing before" << endl;
//        exit(0);
//    }
//    if(neighbor_face_list.size() ==0)
//        return {};
    //cout <<"v " <<mesh->fast_iGameVertex[vh].x() <<" "<<  mesh->fast_iGameVertex[vh].y() <<" "<<  mesh->fast_iGameVertex[vh].z() << endl;

    vector<double>dp;
    vector<int>dp_source;
    vector<MeshKernel::iGameVertex>dp_osqp_answer;
    dp.resize(1<<neighbor_face_list.size());
    dp_source.resize(1<<neighbor_face_list.size());
    dp_osqp_answer.resize(1<<neighbor_face_list.size());
    fill(dp.begin(),dp.end(),-1);
    fill(dp_source.begin(),dp_source.end(),-1);
    function<double(int)> dfs = [&](int state){
        //cout <<vh<< "dfs_state:" << state <<" "<<dp[state]<< endl;
        if(dp[state] >0 ) return dp[state];
        vector<MeshKernel::iGameFaceHandle> local_neighbor_face_list;
        for(int i=0;i<neighbor_face_list.size();i++){
            if(state & (1<<i))
                local_neighbor_face_list.push_back(neighbor_face_list[i]);
        }
        bool succ = false;
        double exceed_dist = 0;
        double self_value = 1e100;
        for(int times=0;times<2;times++) { //避免osqp求解器的不稳定性,多尝试几次
            MeshKernel::iGameVertex v_new = do_quadratic_error_metric_check(vh, local_neighbor_face_list, succ,
                                                                            exceed_dist);
            if (succ) {
                if((v_new - mesh->fast_iGameVertex[vh]).norm() < self_value){
                    dp[state] = (v_new - mesh->fast_iGameVertex[vh]).norm();
                    dp_osqp_answer[state] = v_new;
                    dp_source[state] = 0;
                    self_value = dp[state];
                }
                //break;
            }
        }
        double minx = 1e100;
        int min_from = -1;
//        if(local_neighbor_face_list.size() <= 1 && self_value <0){
//            cout <<"error in dfs dp" << endl;
//            exit(0);
//        }

        if(bitset<22>(state).count() >=1) {
            vector<int> sub_state = std::move(get_sub_state(state, neighbor_face_list.size()));
            int from = -1;
            for (auto next: sub_state) {
                if (next == 0 || state - next == 0)continue;
                if (next > state - next)continue;
                double sub_ans = dfs(next) + dfs(state - next);
                if (sub_ans < minx) {
                    min_from = next;
                    minx = sub_ans;
                }
            }
        }
        if(self_value > minx) {
            dp[state] = minx;
            dp_source[state] = min_from;
            return minx;
        }
        double ret_value = min(minx,self_value);
//        if(ret_value <= 0){
//            cout <<"ret value error"<< endl;
//            exit(0);
//        }
//        else{
//            cout <<"ret value not error"<<" "<<(ret_value <= 0.0)<<" "<<ret_value<< endl;
//        }
//
//        if(dp[state] <= 0){
//            cout << minx <<" "<< self_value << " "<<ret_value<<" "<<dp[state]<< endl;
//            cout <<"dp value error"<< endl;
//            exit(0);
//        }
        return ret_value;
    };
    // cout <<"stdfs"<< endl;
    dfs((1<<neighbor_face_list.size())-1);
    queue<int>q;
    vector<MeshKernel::iGameVertex>ret;

    q.push((1<<neighbor_face_list.size())-1);
    //cout<<"dp="<< dp[(1<<neighbor_face_list.size())-1] << endl;
    // cout <<"**"<< endl;
    while(!q.empty()){
        int now = q.front();
        //cout << now <<" : "<< dp_source[now] << endl;
        //cout << now <<" "<<dp_source[now] << endl;

        q.pop();
        if(dp_source[now] == 0){
            //ret.push_back(now);
            ret.push_back(dp_osqp_answer[now]);
        }
        else{
            q.push(dp_source[now]);
            q.push(now - dp_source[now]);
        }
    }
    return ret;
}

std::mt19937 mt;


vector<MeshKernel::iGameVertex> solve_by_dp(MeshKernel::iGameVertexHandle vh,vector<MeshKernel::iGameFaceHandle> neighbor_face_list){
    vector<MeshKernel::iGameFaceHandle> neighbor_face_list_tmp = neighbor_face_list;
    //cout <<"determins " <<neighbor_face_list.size() <<":"<<avg_edge_limit / 1000<< endl;
    neighbor_face_list.clear();
    for(int i=0;i<neighbor_face_list_tmp.size();i++){
        //cout <<"************************************" << endl;
        MeshKernel::iGameFaceHandle f = neighbor_face_list[i];
        K2::Triangle_3 this_facet(K2::Point_3(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])),
                                  K2::Point_3(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)])),
                                  K2::Point_3(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]))
        );
        if(!this_facet.is_degenerate()){
            neighbor_face_list.push_back(f);
        }

//        bool flag = true;
//        vector<double>len_list;
//        for(int j=0;j<3;j++) {
//            len_list.push_back((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(j)]
//                                - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh((j+1)%3)]).norm());
//        }
//        sort(len_list.begin(),len_list.end());
//        if(len_list[0] + len_list[1]> len_list[2] + avg_edge_limit/100){
//            neighbor_face_list.push_back(f);
//        }

//
//
//        if(flag) {
//
//
//
////            cout <<"v "<< mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)].x()
////            <<" "<<mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)].y() <<" "<<
////            mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)].z() << endl;
////            cout <<"v "<< mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)].x()
////                 <<" "<<mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)].y() <<" "<<
////                 mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)].z() << endl;
////            cout <<"v "<< mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)].x()
////                 <<" "<<mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)].y() <<" "<<
////                 mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)].z() << endl;
////            cout <<"f "<<1 <<" "<<2<<" "<<3<< endl;
////            MeshKernel::iGameVertex normal
////                    = ((mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]
////                        - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]) %
////                       (mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)]
////                        - mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)])).normalize();
////            cout <<"normal "<< normal.x() <<" "<< normal.y() <<" "<< normal.z() << endl;
////            K2::Triangle_3 tri(iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(0)]),
////                                iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(1)]),
////                                iGameVertex_to_Point_K2(mesh->fast_iGameVertex[mesh->fast_iGameFace[f].vh(2)])
////                               );
////            double xx = CGAL::to_double(tri.supporting_plane().orthogonal_vector().x());
////            double yy = CGAL::to_double(tri.supporting_plane().orthogonal_vector().y());
////            double zz = CGAL::to_double(tri.supporting_plane().orthogonal_vector().z());
////            cout << xx <<" "<< yy <<" "<< zz << endl;
//            neighbor_face_list.push_back(f);
//        }
    }
    //exit(0);
    //cout <<"determine "<< neighbor_face_list.size() << endl;
    if(neighbor_face_list.empty()){
        return {};
    }

    //cout << vh << endl;
    if(neighbor_face_list.size()<=10) {
        //cout <<"?? now : " << neighbor_face_list.size() << endl;
        return run(vh,neighbor_face_list);
    }
    else{
        //cout<<"1-ring size: "<< neighbor_face_list.size() <<" meeting limit 16 do random subdivide"<< endl;
        int s = (int)neighbor_face_list.size();
        int cnt = ((int)neighbor_face_list.size() - 1)/10 + 1;
        int each = neighbor_face_list.size() / cnt + 1;
        //if(each > 16)exit(0);
        //cout <<"each "<< each<< endl;
        double ans = 1e100;
        vector<MeshKernel::iGameVertex> ret;
        for(int times=0;times<2;times++){
            vector<MeshKernel::iGameVertex> now;
            std::shuffle(neighbor_face_list.begin(), neighbor_face_list.end(), mt);
            vector<MeshKernel::iGameFaceHandle> que;
            for(int i=0;i<neighbor_face_list.size();i++){
                //cout<<"times: " << times <<" i:"<< i << endl;
                que.push_back(neighbor_face_list[i]);
                if(que.size() == each){
                    vector<MeshKernel::iGameVertex> tmp = run(vh, que);
                    for(auto item : tmp) now.push_back(item);
                    que.clear();
                }
            }
            if(que.size()){
                //cout << "last run:"<<que.size() << endl;
                vector<MeshKernel::iGameVertex> tmp = run(vh, que);
                //cout << "last run end:"<<tmp.size() << endl;
                for(auto item : tmp) now.push_back(item);
                que.clear();
            }
            double now_ans = 0;
            for(int i=0;i<now.size();i++){
                now_ans += now[i].dist(mesh->fast_iGameVertex[vh]);
            }
            if(now_ans < ans){
                ans = now_ans;
                ret = now;
            }
        }
        return ret;
    }
}



#endif //THICKEN2_DP_H
