#include "moores.h"


void put_queue(int **vqRef, int vertex, int *qSize){
    int temp;//create a var on each thread's stack

    (*qSize)++;

    *vqRef = (int *)realloc((*vqRef), (*qSize) * sizeof(int));
    
    temp = (*qSize)-1;

    (*vqRef)[temp] = vertex;
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
    int nullCtr;
    GNodePtr vj_p;
    
    vertexQ = (int *) malloc(qSize * sizeof(int));

    //initialize
    #pragma omp parallel shared(nullCtr, dist, vertexQ, qSize, N, source, adj_listheadRef) private(i,vi,newdist_vj, vj_p, vj)
    {
    #pragma omp for
    for(i=0;i<=N;i++){
         (*dist)[i] = INF;
    }

    #pragma omp single    
    {
        (*dist)[source] = 0;
        nullCtr = 0;
        put_queue(&vertexQ, source, &qSize);
    }
    #pragma omp barrier
    
    // Loop over entries in queue
    while(nullCtr < omp_get_num_threads()){
        
        #pragma omp single    
        {
        nullCtr=0;
        }
        #pragma omp barrier

        #pragma omp critical (two)
        {
        vi = get_queue(&vertexQ, &qSize);//get head of queue
        }
        #pragma omp barrier
        
        if(vi > 0){ 
            //fetch all vj_p at once
            for(vj_p = (*adj_listheadRef)[vi]; vj_p; vj_p = vj_p->next){
                vj = vj_p->vertex; //get vertex number

                //read from shared variable, but different location,  thread safe
                //write to private variable, thread safe
                //FIXME (*dist)[vi] & (*dist)[vj) causes trouble
                newdist_vj = (*dist)[vi] + vj_p->weight; //distance thru vi
                
                #pragma omp critical (one)
                {
                //we need mutual exclusion here, each thread has a different vj
                if ((newdist_vj < (*dist)[vj]) || ((*dist)[vj]==INF)) {
                    (*dist)[vj] = newdist_vj; // Update best distance to vj

                    put_queue(&vertexQ, vj, &qSize); // add vj to q
                }
                }
            }
        }else{
            #pragma omp atomic
            nullCtr++;
        }
        #pragma omp barrier
    }//while(1)
    #pragma omp barrier
    }//parallel 

}
