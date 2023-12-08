//
// Created by rainbowwing on 2023/8/25.
//

#ifndef THICKEN2_OBJ_INPUT_H
#define THICKEN2_OBJ_INPUT_H
double default_move = -1;
double min_near_limit = 1e-5;
int thread_num = 12;
double max_distance_limit = 1.30;
double min_distance_limit = 1.0;
MeshKernel::SurfaceMesh ReadObjFile(const std::string &_InputFile) {
    //std::ifstream inputfile(_InputFile, std::ios::in);
    std::vector<MeshKernel::iGameVertex> vertices;
    std::vector<std::vector<MeshKernel::iGameVertexHandle> > faces;
    std::vector<double> move_dist;
    std::vector<std::vector<double>> normals;
    std::vector<std::vector<double>> uvs;
    std::unordered_map<int, int> V2N;// vertex to normal
    std::unordered_map<int, int> V2T;// vertex to uv
    //std::cout << "Reading " << _InputFile << " File" << std::endl;
    FILE *inputfile = fopen(_InputFile.c_str(), "r");
    char inputs[100];
    while (fscanf(inputfile, "%[^\n]\n", inputs) != EOF) {
        string line(inputs);
        if (line[0] == '#') {
            continue;
        }
        std::stringstream linestream;
        linestream.str(line);
        std::string flag;
        linestream >> flag;
        if (flag == "v") {
            double x, y, z;
            linestream >> x >> y >> z;
            vertices.push_back(MeshKernel::iGameVertex(x, y, z));
        } else if (flag == "f") {
            std::vector<std::string> vex;
            std::string tmp;
            while (linestream >> tmp) vex.push_back(tmp);
            std::vector<MeshKernel::iGameVertexHandle> face(3);
            for (size_t i = 0; i < 3; i++) {
                size_t idx = 0;
                while (idx < vex[i].length() && std::isdigit(vex[i][idx])) idx++;
                int vh = std::stoi(vex[i].substr(0, idx)) - 1;//  obj start v from 1
                face[i] = (MeshKernel::iGameVertexHandle)(vh);
                if (vh >= vertices.size()) {
                    std::cerr << vh << " " << vertices.size() << std::endl;
                }
            }
            if (vex.size() >= 4)
                move_dist.push_back(std::stod(vex[3]));
            else
                move_dist.push_back(default_move);
            faces.push_back(face);
        }
    }
    if (!normals.empty()) {
        int ncnt = normals.size();
        for (int i = 0; i < vertices.size(); ++i) {
            int nidx = V2N[i];
            //if (nidx < 0 || nidx >= ncnt) printf("error: nidx = %d\n", nidx);// debug 用
            assert(nidx >= 0 && nidx < ncnt);
            vertices[i].setNormal(normals[nidx]);
        }
    }
    auto mesh = MeshKernel::SurfaceMesh(vertices, faces, move_dist);
    fclose(inputfile);
    return mesh;
}
shared_ptr <MeshKernel::SurfaceMesh> mesh;
string input_filename;

double start_x;
double start_y;
double start_z;
void set_start(){
    mesh->initBBox();
    start_x = mesh->BBoxMin.x() + (mesh->BBoxMax.x() - mesh->BBoxMin.x()) / 2;
    start_y = mesh->BBoxMin.y() + (mesh->BBoxMax.y() - mesh->BBoxMin.y()) / 2;
    start_z = mesh->BBoxMin.z() + (mesh->BBoxMax.z() - mesh->BBoxMin.z()) / 2;
}

void update_model(){
    FILE *file_update = fopen( (input_filename + "_update.obj").c_str(), "w");
    set_start();
    for(int i=0;i<mesh->VertexSize();i++){
        fprintf(file_update,"v %lf %lf %lf\n",
                mesh->fast_iGameVertex[i].x() - start_x,
                mesh->fast_iGameVertex[i].y() - start_y,
                mesh->fast_iGameVertex[i].z() - start_z);
    }
    for(int i=0;i<mesh->FaceSize();i++){
        fprintf(file_update,"f %d %d %d\n",
                mesh->fast_iGameFace[i].vh(0)+1,
                mesh->fast_iGameFace[i].vh(1)+1,
                mesh->fast_iGameFace[i].vh(2)+1);
    }
    fclose(file_update);

    exit(0);
}
#endif //THICKEN2_OBJ_INPUT_H
