void getInputData(char* text,int* result,int n);
void makeArrayInBlocks(int* in,int* out,int sq_of_P,int k);
void makeRandomWorld(int* a,int n);
void findNeighbors(int rank, int sqrt_P, int P, int* up, int* down, int* left,
    int* right, int* up_left, int* up_right, int* down_left, int* down_right);
void checkInside(int *a,int *b,int n);
int find_new_a(int x,int neighbors);
void checkPerimeter(int* a, int n, int up_left, int up_right, int down_left, int down_right,
		int* up, int* down, int* left, int* right, int* b);
void swapArrays(int **a,int **b);
int world_changed(int *world,int *newWorld,int N);
void printWorld(int *world,int sq,int k);
