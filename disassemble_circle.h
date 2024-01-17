//
// Created by rainbowwing on 2023/12/5.
//

#ifndef PFPOFFSET_DISASSEMBLE_CIRCLE_H
#define PFPOFFSET_DISASSEMBLE_CIRCLE_H
void disassemble_circle(std::map<int, vector<int> >hole_line, vector<vector<int> > &each_hole){
    map<int,int>is_visit_cnt;
    int all_cnt = 0;
    for(auto i : hole_line){
        all_cnt += i.second.size();
        is_visit_cnt[i.first] = i.second.size();
    }
    int now_cnt = 0;
    for(auto i: hole_line){
        if(is_visit_cnt[i.first] == 0){
            continue;
        }
        int now_id = i.first;
//        cout <<"now_cnt" << now_cnt <<"all_cnt"<< all_cnt <<endl;

        if(now_id != -1){
            map<int,int>pre_node;
            queue<int>q;
//            cout <<"now_id"<< now_id << endl;
            q.push(now_id);
            pre_node[now_id] = now_id;
            int end_node = -1;
            while(!q.empty()){
                int node = q.front();
                q.pop();
//                cout <<" now node : "<< node<< endl;
                for(auto next_node : hole_line[node]){
                    if(next_node == now_id){
                        end_node = node;
                        break;
                    }
//                    cout << "next_node" <<" "<< next_node <<" "<< is_visit_cnt[next_node]<<" "<< pre_node.count(next_node) << endl;

                    if(is_visit_cnt[next_node] > 0 && pre_node.count(next_node) == 0){
//                        cout << "chose "<<endl;
                        pre_node[next_node] = node;
                        q.push(next_node);
//                        cout <<" push : "<< next_node<< endl;
                    }
                }
            }
//            cout <<"end_node:"<<end_node<<endl;
            if(end_node != -1){
                vector<int>res;
                int id = end_node;
                while(id != now_id){
                    res.push_back(id);
                    id = pre_node[id];
                }
                res.push_back(now_id);
                reverse(res.begin(),res.end());
                each_hole.push_back(res);
                now_cnt += res.size();
                cout <<"nowid:"<< now_id <<" res size : "<< res.size() << endl;
                for(int j=0;j<res.size();j++){
                    cout <<"res: " <<res[j] << endl;
//                    cout << res[i] << ",";
                    is_visit_cnt[res[j]]--;
                }
            } else {
//                cout <<" now node : "<< now_id<< endl;
//                cout <<"disassemble_circle error "<< endl;
//                exit(0);
            }
        }
    }
}

void fix_hole(vector<int>hole_line,vector<K2::Point_3>position,Tree * tree,vector<vector<int> >& ret){
    std::function<void(vector<int>)> dfs = [&](vector<int> v){
        if(v.size() == 3){
            vector<int> ans;
            for(int i=0;i<v.size();i++){
                ans.push_back(hole_line[v[i]]);
            }
            ret.push_back(ans);
            return ;
        }
        int st = -1;
        int en = -1;
        for(int i=0;i<v.size();i++) {
            int other = (i + 2) % v.size();
            while((other + 1) % v.size() != i){
                K2::Segment_3 seg(position[v[i]],position[v[other]]);//这里退化了
                if(seg.is_degenerate()){
                    st = i;
                    en = other;
                    break;
                }
                std::list< Tree::Intersection_and_primitive_id<K2::Triangle_3>::Type> intersections;
                //tree->all_intersections(seg,std::back_inserter(intersections));
                try {
                    tree->all_intersections(seg, std::back_inserter(intersections));
                } catch (const CGAL::Assertion_exception& e) {
                    // 打印异常信息
                    std::cerr << "Caught a CGAL assertion exception: " << e.what() << std::endl;
                    st = i;
                    en = other;
                    break;
                    // handle exception here
                } catch (const std::exception& e) {
                    // 这将捕获任何派生自std::exception的异常
                    std::cerr << "Caught a standard exception: " << e.what() << std::endl;
                    st = i;
                    en = other;
                    break;
                } catch (...) {
                    // 这将捕获任何其他类型的异常
                    std::cerr << "Caught an unknown exception" << std::endl;
                    st = i;
                    en = other;
                    break;
                }



                bool flag = true;
                for(auto item : intersections) { //todo 这里改成批量插入
                    if (const K2::Segment_3 *s = boost::get<K2::Segment_3>(&(item.first))) {

                    } else if (const K2::Point_3 *p = boost::get<K2::Point_3>(&(item.first))) {
                        if(*p != seg.vertex(0) && *p != seg.vertex(1)){
                            flag = false;
                            break;
                        }
                    }
                }
                if(flag){
                    st = i;
                    en = other;
                    break;
                }
                other = (other + 1) % v.size();
            }
            if(st !=-1)break;
        }
        if(st != -1) {
            vector<int> left;
            vector<int> right;
            for(int i=st;i!=en;i=(i+1)%v.size()){
                left.push_back(v[i]);
            }
            left.push_back(v[en]);

            for(int i=en;i!=st;i=(i+1)%v.size()){
                right.push_back(v[i]);
            }
            right.push_back(v[st]);
            dfs(left);
            dfs(right);
        }
        if(st == -1) {
            vector<int> ans;
            for(int i=0;i<v.size();i++){
                ans.push_back(hole_line[v[i]]);
            }
            ret.push_back(ans);
            return ;
        }
    };
    vector<int> d;
    for(int i=0;i<hole_line.size();i++){
        d.push_back(i);
    }
    dfs(d);
}


#endif //PFPOFFSET_DISASSEMBLE_CIRCLE_H
