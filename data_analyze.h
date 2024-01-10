//
// Created by rainbowwing on 2024/1/8.
//

#ifndef PFPOFFSET_DATA_ANALYZE_H
#define PFPOFFSET_DATA_ANALYZE_H
struct MeshAnalyse{
    int vs;
    int fs;
    double t;
};


vector<string> data_split(string s){
    string tmp;
    vector<string> res;
    for(int i=0;i<s.size();i++){
        if(s[i]== '_'){
            if(tmp.size()){
                res.push_back(tmp);
            }
            tmp.clear();
        }
        else{
            tmp.push_back(s[i]);
        }
    }
    if(tmp.size()){
        res.push_back(tmp);
    }
    return res;
}
void data_analyze(){
    string s;
    map<string,map<pair<string,string>, MeshAnalyse> >mp;
    while(cin>>s){
        cout << s << endl;
        if(s=="OFF")break;
        vector<string> res = data_split(s);
        if(res.size()>1){
            string file_name = res[0];
            string dir = res[2];
            string dist = res[3];
            string type = res[4];
            cout << file_name <<"*"<<dir <<"*"<<dist <<"*"<<type <<endl;
            if(type[0] == 'r'){
                auto mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../res_100/"+s));
                mp[file_name][{dir,dist}].vs = mesh->VertexSize();
                mp[file_name][{dir,dist}].fs = mesh->FaceSize();
            }
            else if(type[0] == 't'){
                FILE * f = fopen(("../res_100/"+s).c_str(),"r");
                double xx;
                fscanf(f,"%lf",&xx);
                mp[file_name][{dir,dist}].t = xx;
            }
        }
        else{
            auto mesh = make_shared<MeshKernel::SurfaceMesh>(ReadObjFile("../res_100/"+s));
            string file_name = s.substr(0,s.find('.'));
            mp[file_name][{"self","0"}] = {(int)mesh->fast_iGameVertex.size(),(int)mesh->fast_iGameFace.size(),0.0};
        }
    }
    cout <<"model name\t"<<"vertices\tfacets\tvertices\tfacets\ttime\tvertices\tfacets\ttime\tvertices\tfacets\ttime"<<endl;
    for(auto i : mp){
        auto file_name = i.first;
        cout << file_name <<"\t"<< mp[file_name][{"self","0"}].vs <<"\t"<< mp[file_name][{"self","0"}].fs <<"\t";
        cout << mp[file_name][{"inward","0.003000"}].vs <<"\t"<< mp[file_name][{"inward","0.003000"}].fs <<"\t"<< mp[file_name][{"inward","0.003000"}].t<<"\t";
        cout << mp[file_name][{"inward","0.001000"}].vs <<"\t"<< mp[file_name][{"inward","0.001000"}].fs <<"\t"<< mp[file_name][{"inward","0.001000"}].t<<"\t";
        cout << mp[file_name][{"inward","0.000500"}].vs <<"\t"<< mp[file_name][{"inward","0.000500"}].fs <<"\t"<< mp[file_name][{"inward","0.000500"}].t<<"\t";
        cout << mp[file_name][{"outward","0.003000"}].vs <<"\t"<< mp[file_name][{"outward","0.003000"}].fs <<"\t"<< mp[file_name][{"outward","0.003000"}].t<<"\t";
        cout << mp[file_name][{"outward","0.001000"}].vs <<"\t"<< mp[file_name][{"outward","0.001000"}].fs <<"\t"<< mp[file_name][{"outward","0.001000"}].t<<"\t";
        cout << mp[file_name][{"outward","0.000500"}].vs <<"\t"<< mp[file_name][{"outward","0.000500"}].fs <<"\t"<< mp[file_name][{"outward","0.000500"}].t<<endl;

    }


}


#endif //PFPOFFSET_DATA_ANALYZE_H
