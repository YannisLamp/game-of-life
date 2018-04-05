void getInputData(char* text,int* result,int n);
void makeRandomWorld(int* a,int n);
void computeNewWorld(int *world,int *newWorld,int n);
void checkInside(int *a,int *b,int n);
void checkPerimeter(int *a,int n,int ul,int ur,int dl,int dr,int *u,int *d,int *l,int*r,int *b);
void swapArrays(int **a,int **b);
int world_changed(int *world,int *newWorld,int n);
int find_new_a(int x,int neighbors);
void printWorld(int *world,int n);
