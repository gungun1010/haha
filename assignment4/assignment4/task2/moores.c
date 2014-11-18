#include "moores.h"

//read lock function, taken from course slide
void read_lock(omp_lock_t* lock, int* readers) {
    omp_set_lock(lock);

    #pragma omp atomic
    (*readers)++;

    omp_unset_lock(lock);
}

void read_unlock(int* readers) {
    #pragma omp atomic
    (*readers)--; 
}

void write_lock(omp_lock_t* lock, int readers) {
     omp_set_lock(lock);
     while (readers > 0);
}

void write_unlock(omp_lock_t* lock) {
     omp_unset_lock(lock);
}

void put_queue(int **vqRef, int vertex, int *qTail, int qSize, int qReaders, omp_lock_t *qLock, int myNum){
    int i;
/*
    #pragma omp flush (qTail)
    i = (*qTail) + myNum;
    
    if(i<qSize){
        write_lock(qLock, qReaders);
        (*vqRef)[i] = vertex;
        write_unlock(qLock);
       
       #pragma omp atomic
       (*qTail)++;
    }
*/
}

int get_queue(int **vqRef, int *qHead, int myNum, int qSize, int *qReaders, omp_lock_t *qLock){
    int ret, i;
    
    #pragma omp flush (qHead)
    i = (*qHead) + myNum;
    printf("%d: my i is %d\n",myNum, i);
    if(i<qSize){ 
        read_lock(qLock, qReaders);
        ret = (*vqRef)[i];
        read_unlock(qReaders);

    }else{
        ret = 0;
    }

    if(!ret){
        ret = INF;
    }else{
        #pragma omp atomic
        (*qHead)++;    
    }
    return ret;
}

void mooresLaw(int N, int source, GNodePtr **adj_listheadRef, int **dist){
    int i;
    int *vertexQ;//vertex array
    int vi,vj,newdist_vj;
    int qSize = N;
    int qHead = INIT_Q_SIZE;
    int qTail = INIT_Q_SIZE+1;
    int nullCtr, myNum, thisDis;
    GNodePtr vj_p;
    GNodePtr thisAdj;
    omp_lock_t lock;
    omp_lock_t distLock;
    omp_lock_t qLock;
    int readers, distReaders, qReaders;
    
    vertexQ = (int *) malloc(qSize * sizeof(int));

    //initialize
    #pragma omp parallel shared(nullCtr, dist, vertexQ, qHead, qTail, qSize, N, source, adj_listheadRef,lock, readers, distLock, distReaders, qLock, qReaders) private(i,vi,newdist_vj, vj_p, vj, myNum, thisDis, thisAdj)
    {
    myNum = omp_get_thread_num();
    #pragma omp for
    for(i=0;i<=N;i++){
         (*dist)[i] = INF;
    }
    //printf("thread %d at checkpoint 1\n",omp_get_thread_num());
    #pragma omp single    
    {
        (*dist)[source] = 0;
        nullCtr = 0;
        readers = 0;
        qReaders = 0;
        omp_init_lock (&lock);
        distReaders = 0;
        omp_init_lock (&distLock);
        omp_init_lock (&qLock);
        put_queue(&vertexQ, source, &qTail, qSize, qReaders, &qLock, myNum);
    }
    #pragma omp barrier
    
    //printf("thread %d at checkpoint 2\n",omp_get_thread_num());
    // Loop over entries in queue
    while(nullCtr < omp_get_num_threads()){
        
        #pragma omp master 
        {
        nullCtr=0;
        }
        //#pragma omp barrier
        //printf("---------------%d-----------------------------------------------\n",myNum);

        vi = get_queue(&vertexQ, &qHead, myNum, qSize, &qReaders, &qLock );//get head of queue
        //printf("%d-----------------------------------------------\n",vi);
        //printf("-----------------------------------------------%d\n",vi);
        printf("%d gets %d\n", myNum, vi); 
        #pragma omp barrier
        if(myNum == 0){
            printf("-----------------------------------\n");
        }
        sleep(1);

        if(vi > 0){ 
            //fetch all vj_p at once
            #pragma omp flush (adj_listheadRef)
            read_lock(&lock, &readers);
            thisAdj = (*adj_listheadRef)[vi];
            read_unlock(&readers);

            for(vj_p = thisAdj; vj_p; vj_p = vj_p->next){
                vj = vj_p->vertex; //get vertex number
                //printf("%d: vi %d gets vj %d\n", omp_get_thread_num(), vi,vj);
                //read from shared variable, but different location,  thread safe
                //write to private variable, thread safe
                //FIXME (*dist)[vi] & (*dist)[vj) causes trouble
                #pragma omp flush (dist)
                read_lock(&distLock, &distReaders);
                thisDis = (*dist)[vi];
                read_unlock(&distReaders);

                newdist_vj = thisDis + vj_p->weight; //distance thru vi
                
                read_lock(&distLock, &distReaders);
                thisDis = (*dist)[vj];
                read_unlock(&distReaders);

                //we need mutual exclusion here, each thread has a different vj
                if ((newdist_vj < thisDis) || (thisDis==INF)) {

                    write_lock(&distLock, distReaders);
                    (*dist)[vj] = newdist_vj; // Update best distance to vj
                    write_unlock(&distLock);

                    //printf("vj %d ready to q\n",vj); 
                    //exit(1);
                    put_queue(&vertexQ, vj, &qTail, qSize, qReaders, &qLock, myNum); // add vj to q
                }
            }
        }else{
            #pragma omp atomic
            nullCtr++;
            
            //printf("%d\n", nullCtr);
        }
        //#pragma omp barrier
    }//while(1)
    //#pragma omp barrier
    }//parallel 

}
