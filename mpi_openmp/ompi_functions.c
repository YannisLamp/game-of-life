#include<stdlib.h>
#include<stdio.h>
#include<time.h>

#ifdef _OPENMP
#include <omp.h>
#else
int omp_get_thread_num(void) { return 0; }
int omp_get_num_threads(void) { return 1; }
#endif


void getInputData(char* text, int* result, int n) {
	FILE* fp;
	char ch;
	int i = 0;
	fp = fopen(text, "r");
	do {
		ch = fgetc(fp);
		if (ch == '0' || ch == '1') {
			result[i] = ch - '0';
			i++;
		}
	} while (i < n && ch != EOF);
	fclose(fp);
	return ;
}

void makeArrayInBlocks(int *in, int *out, int sq_of_P, int k) {
	int i, j, z, w;
	int n = sq_of_P * k, x = 0;
	// #   pragma omp parallel for private(i,j) collapse(2)
 #   pragma omp parallel for private(i,j,z,w) collapse(4)
	for (i = 0; i < sq_of_P; i++)
		for (j = 0; j < sq_of_P; j++)
			for (z = 0; z < k; z++)
				for (w = 0; w < k; w++) {
					out[x++] = in[(i*k + z)*n + (j*k + w)];
				}
	return;
}


void makeRandomWorld(int *a, int n) {
	int i;
	srand(time(NULL));
#   pragma omp parallel for private(i)
	for (i = 0; i < n; i++){
		if (!(rand() % 2))
			a[i] = 1;
		else
			a[i] = 0;
	}
	return ;
}

void findNeighbors(int rank, int sq_of_P, int P, int* U, int* D, int* L,
				   int* R, int* Ul, int* Ur, int* Dl, int* Dr) {
	int Pmod = rank%sq_of_P;
	int Pdiv = rank/sq_of_P;  //rank = div*sP+mod
	*D  = ((Pdiv + 1) % sq_of_P)*sq_of_P + Pmod;
	*R  = Pdiv*sq_of_P + (Pmod + 1)%sq_of_P;
	*Dr = ((Pdiv + 1)%sq_of_P)*sq_of_P + (Pmod + 1)%sq_of_P;

// These two could be divided into two sections, but as tested, the overhead is too much
//#	pragma omp parallel sections
//	{
//#		pragma omp section
//		{
			if (Pdiv == 0) {
				*U  = P - sq_of_P + Pmod;
				*Ur = P - sq_of_P + (Pmod + 1)%sq_of_P;
				if (Pmod == 0)
					*Ul = P - 1;
				else
					*Ul = P - sq_of_P + Pmod - 1;
			}
			else {
				*U  = (Pdiv - 1)*sq_of_P + Pmod;
				*Ur = (Pdiv - 1)*sq_of_P + (Pmod + 1)%sq_of_P;
        		if (Pmod == 0)
        			*Ul = Pdiv*sq_of_P - 1;
        		else
        			*Ul = (Pdiv - 1)*sq_of_P + Pmod - 1;
			}
//		}
//#		pragma omp section
//		{
			if (Pmod == 0) {
				*L  = (Pdiv + 1)*sq_of_P - 1;
				*Dl =((Pdiv + 1)%sq_of_P)*sq_of_P + sq_of_P - 1;
			}
			else {
				*L  = Pdiv*sq_of_P + Pmod - 1;
				*Dl = ((Pdiv + 1)%sq_of_P)*sq_of_P + Pmod - 1;
			}
//		}
//	}
	return ;
}

int find_new_a(int x, int neighbors) {
	if (x == 1) {
    	if ((neighbors == 2) || (neighbors == 3)) {
            return 1; }
        else {
            return 0; }
    }
    else { //x==0
        if (neighbors == 3) {
            return 1; }
        else {
            return 0; }
    }
}

void checkInside(int* a, int* b, int n) {
	int i, j, neighbor;
#   pragma omp parallel for private(i,j,neighbor) collapse(2)
	for (i = 1; i < n - 1; i++) {
		for (j = 1; j < n - 1; j++) {
			neighbor = a[(i - 1)*n + (j - 1)] + a[(i - 1)*n + j] + a[(i - 1)*n + (j + 1)] +
				   a[(i + 1)*n + (j - 1)] + a[(i + 1)*n + j] + a[(i + 1)*n + (j + 1)] +
				   a[i*n + (j - 1)] + a[i*n + (j + 1)];
			b[i*n +j] = find_new_a(a[i*n + j], neighbor);
		}
	}
	return ;
}

void checkPerimeter(int* a, int n, int ul, int ur, int dl, int dr,
					int* u, int* d, int* l, int* r, int* b) {
	int i = 0, neighbor = 0;

// These four could also be divided into four sections, but as tested,
// the overhead is also too much

	//for a[0]
	neighbor = ul + u[0] + u[1] + l[0] + l[1] + a[1] + a[n] + a[n + 1];
	b[0] = find_new_a(a[0],neighbor);

	//for a[n - 1]
	neighbor = ur + u[n - 1] + u[n - 2] + r[0] + r[1] + a[n - 2] + a[2*n - 1] + a[2*n - 2];
	b[n - 1] = find_new_a(a[n - 1],neighbor);

	//for a[n*(n - 1)]
	neighbor = dl + l[n - 2] + l[n - 1] + d[0] + d[1] + a[n*(n - 1) + 1] + a[n*(n - 2)] + a[n*(n - 2) + 1];
	b[n*(n - 1)] = find_new_a(a[n*(n - 1)],neighbor);

	//for a[n*n - 1]
	neighbor = dr + d[n - 1] + d[n - 2] + r[n - 1] + r[n - 2] + a[n*n - 2] + a[n*(n - 1) - 1] + a[n*(n - 1) - 2];
	b[n*n - 1] = find_new_a(a[n*n - 1],neighbor);

#	pragma omp parallel
	{
		//for up
#   	pragma omp for private(i, neighbor) nowait
		for (i = 1; i < n - 1; i++) {
			neighbor = u[i - 1] + u[i] + u[i + 1] + a[i - 1] + a[i + 1] + a[i - 1 + n] + a[i + n] + a[i + 1 + n];
			b[i] = find_new_a(a[i],neighbor);
		}
		//for down
#   	pragma omp for private(i, neighbor) nowait
		for (i = 1; i < n - 1; i++) {
			neighbor = d[i - 1] + d[i] + d[i + 1] + a[n*(n - 1) + i - 1] + a[n*(n - 1) + i + 1] +
			   		a[n*(n - 2) + i] + a[n*(n - 2) + i - 1] + a[n*(n - 2) + i + 1];
					b[n*(n - 1) + i] = find_new_a(a[n*(n - 1) + i], neighbor);
		}
		//for left
#   	pragma omp for private(i, neighbor) nowait
		for (i = 1; i < n - 1; i++) {
			neighbor = l[i - 1] + l[i] + l[i + 1] + a[(i - 1)*n] + a[(i - 1)*n + 1] + a[i*n + 1] +
			    	a[(i + 1)*n] + a[(i + 1)*n + 1];
					b[n*i] = find_new_a(a[n*i],neighbor);
		}
		//for right
#   	pragma omp for private(i, neighbor) nowait
		for (i = 1;i < n - 1; i++) {
			neighbor = r[i - 1] + r[i] + r[i + 1] + a[i*n - 1] + a[(i + 2)*n - 1] +
					a[i*n - 2] + a[(i + 2)*n - 2] + a[(i + 1)*n - 2];
					b[(i + 1)*n - 1] = find_new_a(a[(i + 1)*n - 1], neighbor);
		}
	}
		return;
}


void swapArrays(int** a, int** b) {
	int *temp = *a;
	*a = *b;
	*b = temp;
	return;
}

// Different world_changed than the mpi one, could be implemented like that
// but because it is more likely for the world to have changed, the mpi
// implementation is chosen, even in the openmp, simply because of the fact that
// it is far more likely for us to encounter a change in the world and simply
// break, rather than execute all of the loops even in parallel

/*

int world_changed(int *world,int *newWorld,int N){
	int i = 0, dif = 0;
#   pragma omp parallel for shared(dif) private(i)
	for(i = 0; i < N; i++) {
		if(world[i] != newWorld[i]) {
#			pragma omp critical
			{
				dif += 1;
			}
		}
	}
#	pragma omp barrier

	if(dif == 0)
		return 0;
	else
		return 1;
}

*/

int world_changed(int *world, int *newWorld, int N){
	int i;
	for(i = 0; i < N; i++)
		if(world[i] != newWorld[i])
			return 1;
	return 0;
}

void printWorld(int *world,int sq,int k){
    int i,j,x,y;
    for(i = 0; i < sq; i++) {
        for(j = 0; j < k; j++) {
            for(x = 0; x < sq; x++) {
                for(y = 0; y < k; y++) {
                    printf("%d ",world[y+k*k*x+k*j+sq*k*k*i]);
                }
            }
        putchar('\n');
        }
    }
    putchar('\n');
    return ;
}
