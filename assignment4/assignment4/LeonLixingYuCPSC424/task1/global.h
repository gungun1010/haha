#define FALSE 0
#define TRUE 1
#define INF -1
#define INIT_Q_SIZE 0

typedef struct gnode{
  int vertex;
  int weight;
  struct gnode *next;
} GNode, *GNodePtr;


