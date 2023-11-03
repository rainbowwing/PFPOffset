//
// Created by rainbowwing on 2023/11/3.
//

#ifndef PFPOFFSET_MODEL_UPDATE_H
#define PFPOFFSET_MODEL_UPDATE_H

void update_model(){
    FILE *file_update = fopen( (input_filename + "_update.obj").c_str(), "w");
    mesh->initBBox();
    double stx = mesh->BBoxMin.x() + (mesh->BBoxMax.x() - mesh->BBoxMin.x())/2;
    double sty = mesh->BBoxMin.y() + (mesh->BBoxMax.y() - mesh->BBoxMin.y())/2;
    double stz = mesh->BBoxMin.z() + (mesh->BBoxMax.z() - mesh->BBoxMin.z())/2;

    for(int i=0;i<mesh->VertexSize();i++){
        fprintf(file_update,"v %lf %lf %lf\n",
                mesh->fast_iGameVertex[i].x()-stx,
                mesh->fast_iGameVertex[i].y()-sty,
                mesh->fast_iGameVertex[i].z()-stz);
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


#endif //PFPOFFSET_MODEL_UPDATE_H
