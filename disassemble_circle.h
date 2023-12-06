//
// Created by rainbowwing on 2023/12/5.
//

#ifndef PFPOFFSET_DISASSEMBLE_CIRCLE_H
#define PFPOFFSET_DISASSEMBLE_CIRCLE_H
void disassemble_circle(std::map<int,vector<int> >hole_line,vector<vector<int> > &each_hole){
    map<int,int>is_visit_cnt;
    int all_cnt = 0;
    for(auto i : hole_line){
        all_cnt += i.second.size();
        is_visit_cnt[i.first] = i.second.size();
    }
    int now_cnt = 0;
    while(now_cnt != all_cnt){
        int now_id = -1;
        cout <<"now_cnt" << now_cnt <<"all_cnt"<< all_cnt<<endl;
        for(auto i: hole_line){
            if(is_visit_cnt[i.first]!=0){
                now_id = i.first;
                break;
            }
        }
        if(now_id != -1){
            map<int,int>pre_node;
            queue<int>q;
            q.push(now_id);
            pre_node[now_id] = now_id;
            int end_node = -1;
            while(!q.empty()){
                int node = q.front();
                q.pop();
                for(auto next_node : hole_line[node]){
                    if(next_node == now_id){
                        end_node = node;
                        break;
                    }
                    if(is_visit_cnt[next_node] > 0 && pre_node.count(next_node) == 0){
                        pre_node[next_node] = node;
                        q.push(next_node);
                    }
                }
            }
            cout <<"end_node:"<<end_node<<endl;
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
                for(int i=0;i<res.size();i++){
                    cout <<"res: "<< res[i] << endl;
                    is_visit_cnt[res[i]]--;
                }
            }else{
                cout <<"disassemble_circle error "<< endl;
                exit(0);
            }
        }
    }
}

#endif //PFPOFFSET_DISASSEMBLE_CIRCLE_H
