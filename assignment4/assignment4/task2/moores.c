#include "moores.h"


void put_queue(int **vqRef, int vertex, int *qSize, omp_lock_t *lock, int *readers){
    int temp;//create a var on each thread's stack
    
    #pragma omp atomic
    (*qSize)++;
    
    write_lock(lock, *readers);
    *vqRef = (int *)realloc((*vqRef), (*qSize) * sizeof(int));
    
    temp = (*qSize)-1;

    (*vqRef)[temp] = vertex;
    write_unlock(lock);
}

int get_queue(int **vqRef, int *qSize, omp_lock_t *lock, int *readers){
    int ret, i;
    
    if((*qSize) > 0){
        //protected read from queue, only read when there is no write
        read_lock(lock, readers);
        ret = (*vqRef)[0];
        read_unlock(readers);

        //remvoe it from queue after fetching it
        #pragma omp atomic
        (*qSize)--;
        
        //write protection, one write at a time
        write_lock(lock, *readers);
        for(i=0; i< (*qSize); i++)
            (*vqRef)[i] = (*vqRef)[i+1];

        *vqRef = (int *)realloc((*vqRef), (*qSize) * sizeof(int));
        write_unlock(lock);
    }else{
        ret = INF;
    }

    return ret;
}

void read_lock(omp_lock_t* lock, int* readers) {
    omp_set_lock(lock);

    #pragma omp atomic
    readers++;

    omp_unset_lock(lock);
}

void read_unlock(int* readers) {
    #pragma omp atomic
    readers--; 
}

void write_lock(omp_lock_t* lock, int readers) {
       omp_set_lock(lock);
       while (readers > 0);
}

void write_unlock(omp_lock_t* lock) {
       omp_unset_lock(lock);
}

void mooresLaw(int N, int source, GNodePtr **adj_listheadRef, int **dist){
    int i;
    int *vertexQ;//vertex array
    int vi,vj,newdist_vj;
    int qSize = INIT_Q_SIZE;
    int nullCtr;
    GNodePtr vj_p;
    omp_lock_t lock;
    omp_lock_t distLock;
    int readers, distReaders;
    
    vertexQ = (int *) malloc(qSize * sizeof(int));

    //initialize
    #pragma omp parallel shared(nullCtr, dist, vertexQ, qSize, N, source, adj_listheadRef, lock, readers, distLock, distReaders) private(i,vi,newdist_vj, vj_p, vj)
    {
    #pragma omp for
    for(i=0;i<=N;i++){
         (*dist)[i] = INF;
    }

    #pragma omp single    
    {
        (*dist)[source] = 0;
        nullCtr = 0;
        readers = 0;
        omp_init_lock (&lock);
        distReaders = 0;
        omp_init_lock (&distLock);
        put_queue(&vertexQ, source, &qSize, &lock, &readers);
    }
    #pragma omp barrier
    
    // Loop over entries in queue
    while(nullCtr < omp_get_num_threads()){
        
        #pragma omp single    
        {
        nullCtr=0;
        }
        #pragma omp barrier
        vi = get_queue(&vertexQ, &qSize, &lock, &readers);//get head of queue
        printf("%d get %d\n", omp_get_thread_num(), vi); 
        if(vi > 0){ 
            //fetch all vj_p at once
            for(vj_p = (*adj_listheadRef)[vi]; vj_p; vj_p = vj_p->next){
                vj = vj_p->vertex; //get vertex number
                
                //read from shared variable, but different location,  thread safe
                //write to private variable, thread safe
                //FIXME (*dist)[vi] & (*dist)[vj) causes trouble
                //when vj from another thread happens to be vi.. shits happen

                read_lock(&distLock, &distReaders);
                newdist_vj = (*dist)[vi] + vj_p->weight; //distance thru vi
                read_unlock(&distReaders);

                //we need mutual exclusion here, each thread has a different vj
                if ((newdist_vj < (*dist)[vj]) || ((*dist)[vj]==INF)) {
                    
                    write_lock(&distLock, distReaders);
                    (*dist)[vj] = newdist_vj; // Update best distance to vj
                    write_unlock(&distLock);

                    put_queue(&vertexQ, vj, &qSize, &lock, &readers); // add vj to q
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
