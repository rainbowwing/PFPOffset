//
// Created by rainbowwing on 2023/10/1.
//

#ifndef THICKEN2_ROUGH_VERTEX_H
#define THICKEN2_ROUGH_VERTEX_H
K2::Point_3 exact_start(0,0,0);
CGAL::Epeck::FT resolution(0.01);

struct RoughVertex{
    long long x,y,z;
    RoughVertex(){}
    RoughVertex(long long x,long long y,long long z){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    K2::Point_3 exact(){
        return exact_start + K2::Vector_3 (resolution * CGAL::Epeck::FT(x), resolution * CGAL::Epeck::FT(y) ,resolution * CGAL::Epeck::FT(z));
    }
};
std::hash<long long>long_hash;
struct rough_vertex_hash{
    size_t operator () (RoughVertex x) const {
        return long_hash(x.x) ^ int_hash(x.y<<1) ^ int_hash(x.z<<2);
    }};
struct rough_vertex_equal
{
    bool operator() (RoughVertex a,  RoughVertex b) const {
        return a.x == b.x  &&  a.y == b.y &&  a.z == b.z;
    }
};



long long get_rough_value(CGAL::Epeck::FT x,CGAL::Epeck::FT o){
    long long left = (long long)(CGAL::to_double((x - o)/resolution));
    if(CGAL::Epeck::FT(left) *resolution + o <= x && CGAL::Epeck::FT(left+1)+ o>x ){
        return left;
    }
    long long right = left * 2;
    long long ans = left;
    while(left<=right){
        long long mid = (left + right)>>1;
        if(CGAL::Epeck::FT(mid) *resolution + o <= x){
            ans = mid;
            left = mid + 1;
        }
        else
            right = mid - 1;
    }
    if(CGAL::Epeck::FT(ans) *resolution + o <= x && CGAL::Epeck::FT(ans+1)*resolution + o>x ){
        return ans;
    }
    else {
        cout <<"error"<<endl;
        cout <<"o:" <<CGAL::to_double(o) <<" x:"<<CGAL::to_double(x) <<" ans:"<<ans <<" +:"<< CGAL::Epeck::FT(ans) *resolution + o <<" ++:"<< CGAL::Epeck::FT(ans+1)*resolution + o<<endl;
        exit(0);
    }
    return ans;
}

RoughVertex transfer(K2::Point_3 v){
    long long x = get_rough_value(v.x(), exact_start.x());
    long long y = get_rough_value(v.y(), exact_start.y());
    long long z = get_rough_value(v.z(), exact_start.z());
    RoughVertex ret(x,y,z);
    CGAL::Epeck::FT min_miss = squared_distance(ret.exact(),v);

    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<2;k++){
                RoughVertex new_ret(x+i,y+j,z+k);
                CGAL::Epeck::FT new_miss = squared_distance(new_ret.exact(),v);
                if(new_miss < min_miss){
                    min_miss = new_miss;
                    ret = new_ret;
                }
            }
        }
    }
    return ret;
}

RoughVertex transfer(double x,double y,double z){
    return transfer(K2::Point_3(x,y,z));
}



#endif //THICKEN2_ROUGH_VERTEX_H
