#include<stdlib.h>
#include<stdio.h>
#include<time.h>

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

void makeRandomWorld(int *a, int n) {
        int i;
        srand(time(NULL));
        for (i = 0; i < n; i++){
                if (!(rand() % 4))
                        a[i] = 1;
                else
                        a[i] = 0;
        }
        return ;
}

int find_new_a(int x, int neighbors) {
        if (x) { //x==1
       	    if (neighbors == 2 || neighbors == 3)
                return 1;
            else
                return 0;
        }
        else { //x==0
            if (neighbors == 3)
                return 1;
            else
                return 0;
        }
}


void checkInside(int* a, int* b, int n) {
        int i, j, neighbor;
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
        int i,neighbor;
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

        //for up
        for (i = 1; i < n - 1; i++) {
                neighbor = u[i - 1] + u[i] + u[i + 1] + a[i - 1] + a[i + 1] + a[i - 1 + n] + a[i + n] + a[i + 1 + n];
                b[i] = find_new_a(a[i],neighbor);
        }
        //for down
        for (i = 1; i < n - 1; i++) {
                neighbor = d[i - 1] + d[i] + d[i + 1] + a[n*(n - 1) + i - 1] + a[n*(n - 1) + i + 1] +
                                   a[n*(n - 2) + i] + a[n*(n - 2) + i - 1] + a[n*(n - 2) + i + 1];
                b[n*(n - 1) + i] = find_new_a(a[n*(n - 1) + i], neighbor);
        }

        //for left
        for (i = 1; i < n - 1; i++) {
                neighbor = l[i - 1] + l[i] + l[i + 1] + a[(i - 1)*n] + a[(i - 1)*n + 1] + a[i*n + 1] +
                                   a[(i + 1)*n] + a[(i + 1)*n + 1];
                b[n*i] = find_new_a(a[n*i],neighbor);
        }
        //for right
        for (i = 1;i < n - 1; i++){
                neighbor = r[i - 1] + r[i] + r[i + 1] + a[i*n - 1] + a[(i + 2)*n - 1] +
                                a[i*n - 2] + a[(i + 2)*n - 2] + a[(i + 1)*n - 2];
                b[(i + 1)*n - 1] = find_new_a(a[(i + 1)*n - 1], neighbor);
        }
        return;
}

void computeNewWorld(int *w,int *nw,int n){
        int z;
        int *up=NULL,*down=NULL,*right=NULL,*left=NULL;
        up = malloc(n*sizeof(int));
        down = malloc(n*sizeof(int));
        right = malloc(n*sizeof(int));
        left = malloc(n*sizeof(int));
        checkInside(w,nw,n);//using the mpi functs to compute
        for(z=0;z<n;z++){
                right[z]=w[z*n];
                left[z]=w[z*n + (n-1)];
                up[z]=w[n*n-n + z];
                down[z]=w[z];
        }
        //up-left -> down-right = w[n*n-1]
        //up-right -> down-left = w[n*(n-1)]
        //down-left -> up-right = w[n-1]
        //down-right -> up-left = w[0]
        //checkPerimeter(w,n,w[n*n-1],w[n*n-n],w[n-1],w[0],up,down,right,left,nw);
        checkPerimeter(w,n,w[n*n-1],w[n*n-n],w[n-1],w[0],down,up,left,right,nw);
        free(up);
        free(down);
        free(right);
        free(left);
        return ;
}

void swapArrays(int** a, int** b) {
        int *temp = *a;
        *a = *b;
        *b = temp;
        return;
}

int world_changed(int *world,int *newWorld,int N){
        int i;
        for(i=0;i<N;i++)
                if(world[i] != newWorld[i])
                        return 1;
        return 0;
}

void printWorld(int *world,int n){
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
			printf("|%d",world[i*n+j]);
		printf("|\n");
	}
	return ;
}

