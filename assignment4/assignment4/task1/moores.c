#include "moores.h"


void put_queue(int **vqRef, int vertex, int *qSize){
    (*qSize)++;

    *vqRef = (int *)realloc((*vqRef), (*qSize) * sizeof(int));
    (*vqRef)[(*qSize)-1] = vertex;
}

int get_queue(int **vqRef, int *qSize){
    int ret, i;
    
    if((*qSize) > 0){
        ret = (*vqRef)[0];
        
        //remvoe it from queue after fetching it
        (*qSize)--;
        for(i=0; i< (*qSize); i++)
            (*vqRef)[i] = (*vqRef)[i+1];

        *vqRef = (int *)realloc((*vqRef), (*qSize) * sizeof(int));
        
    }else{
        ret = INF;
    }
    return ret;
}

void mooresLaw(int N, int source, GNodePtr **adj_listheadRef, int **dist){
    int i;
    int *vertexQ;//vertex array
    int vi,vj,newdist_vj;
    int qSize = INIT_Q_SIZE;
    GNodePtr vj_p;
    
    vertexQ = (int *) malloc(qSize * sizeof(int));

    //initialize
    for(i=0;i<=N;i++)
        (*dist)[i] = INF;

    (*dist)[source] = 0;
    put_queue(&vertexQ, source, &qSize);

    // Loop over entries in queue
    while((vi = get_queue(&vertexQ, &qSize)) > 0){//get head of queue

        vj_p = (*adj_listheadRef)[vi];

        while(vj_p){
            vj = vj_p->vertex; //get vertex number
            newdist_vj = (*dist)[vi] + vj_p->weight; //distance thru vi

            if ((newdist_vj < (*dist)[vj]) || ((*dist)[vj]==INF)) {
                (*dist)[vj] = newdist_vj; // Update best distance to vj
                put_queue(&vertexQ, vj, &qSize); // add vj to q
            }
            vj_p = vj_p->next;
        }
    }

}
