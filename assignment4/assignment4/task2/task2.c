//          #include <math.h>
#include "moores.h"


void add_gnode(int vi, int vj, int weight, GNodePtr *adj_listhead) {
  GNodePtr gnodep = (GNodePtr) malloc(sizeof(GNode));
  gnodep->vertex = vj;
  gnodep->weight = weight;
  gnodep->next = adj_listhead[vi];
  adj_listhead[vi] = gnodep;
}

int main(int argc, char *argv[]) {

  /*

    This is a fragmentary sample program to illustrate how to read the data for CPSC424/524 Assignment #4.

    Author: Andrew Sherman, Yale University

    Date: 3/20/2013

  */

  
  int N, M; // N=# of vertices; M=# of edges
  double et0,et1, etime0, etime1, cptime;
  int foundp = FALSE; // To indicate that a "p" card has been found
  int args_assigned = 0; // Counts number of args assigned from sscanf calls
  int nsources; // nsources = # of source vertices to use for this graph
  int i, j, k;
  FILE *fp;
  char line[128];
  int source, dest, weight; // Info about an edge in the graph
  int *sources;
  int totaledges = 0; // Count of graph edges
  GNodePtr *adj_listhead; // Use your own pointer type, depending on the structures you are using.
  int *dist; //distance array

  // -------------------------------- Read in Graph File -------------------------------------------

  fp = fopen(argv[1], "r"); // Assumes that arv[1] is the name of the .gr file
  if (fp != NULL) {
    while (fgets (line, sizeof line, fp) != NULL) { /* read a line */
      //      fputs(line, stderr);
      if (!foundp) { // If we haven't found the "p" line yet...
        args_assigned = sscanf (line, "p sp %d %d", &N, &M); // Matches only the "p" line: N=# of vertices; M=# of edges
        if (args_assigned == 2){
          printf("N and M are %d and %d\n", N, M);
          foundp = TRUE;
	  // Need to use proper cast in the following line. You need to define according to your own structure types for graph storage 
	  adj_listhead = (GNodePtr*) malloc((N+1) * sizeof(GNodePtr)); // Use N+1 entries because data files use vertex numbers from 1 to N (ignoring 0)
          for(i = 0; i <= N; i++){
            adj_listhead[i] = NULL; //Your adj_listhead or equivalent should be of the proper type for your data structures
          }
        }
      }
      else { // After we've found the "p" line, we only want to consider "a" lines
        args_assigned = sscanf (line, "a %d %d %d", &source, &dest, &weight); // Matches only "a" lines: edge source vertex, edge dest vertex, edge weight (all ints) 
        if (args_assigned == 3){
          add_gnode(source, dest, weight, adj_listhead); // Replace this with a call or code to insert the new edge into your adjacency list data structure
	  // DO NOT FORGET TO FREE THE GRAPH AT THE END!!
	  totaledges++; // Counting the edges
        }
      }
    }
    fclose(fp);
    printf ("Graph should have %d vertices and %d edges.\n It actually has %d edges.\n",N,M,totaledges); // Only for sanity checking
  }
  else { // Problem opening the .gr input file
    perror (argv[1]); 
    return(1);
  }

  if(!foundp) { // Uh Oh!
    perror ("Never found a properly-formed \"p\" line for the graph file.\n");
    return(1);
  }

  // -------------------------------- Read in file of sources to test ---------------
  foundp = FALSE;
  fp = fopen(argv[2], "r"); // Assumes that argv[2] is the name of the .ss file
  if (fp != NULL) {
    while (fgets (line, sizeof line, fp) != NULL) { /* read a line */
      //      fputs(line, stderr);
      if (!foundp) { // If we haven't found the "p" line yet...
        args_assigned = sscanf (line, "p aux sp ss %d", &nsources); // Matches only the "p" line: nsources = # sources to use
        if (args_assigned == 1){
          printf("Number of sources = %d\n", nsources);
          foundp = TRUE;
          sources = malloc(sizeof(int) * nsources); // allocate space for an array holding the sources
          j = 0; //index of next source vertex
        }
      }
      else if (j<nsources) { // After we've found the "p" line, we only want to consider nsource "s" lines
        args_assigned = sscanf (line, "s %d", &k);  // Matches only "s" lines: number of source vertex to use
        if (args_assigned == 1){
          sources[j] = k; // Save the source vertex
	  printf ("sources[%d] = %d\n",j,k);
          j++;
        }
      }
      else break; // Stop after the first nsource source vertices have been read
    }

    fclose(fp);

  }
  else {
    perror (argv[1]);
    return(1);
  }

  if(!foundp) { // Uh Oh!
    perror ("Never found a properly-formed \"p\" line for the sources file.\n");
    return(1);
  }
  
  //-------------------------------- My own code starts here--------------------------
  omp_set_num_threads(NUM_THREADS); //setting number of threads

  timing(&etime0,&cptime);
  //for(i=0; i<nsources; i++){ 
      dist = (int *)malloc((N+1) * sizeof(int));

      //mooresLaw, starts 
      timing(&et0,&cptime);
      mooresLaw(N, sources[0], &adj_listhead, &dist);
      timing(&et1,&cptime);

      printf("\nsource: %d\n",sources[0]);
      
      for(j=1; j<N+1; j=j+1000){
          printf("%d, %d\n",j, dist[j]);
          if(j >= 10001){
              break;
          }
      }

      printf("%d, %d\n",N, dist[N]);

      printf("this source run time %.5f seconds\n",(et1-et0));
      printf("------------------------------------\n");
  //}
  timing(&etime1,&cptime);
  printf("total run time: %.5f seconds\n",(etime1-etime0));
}
