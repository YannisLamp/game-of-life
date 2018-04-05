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

void makeArrayInBlocks(int *in, int *out, int sqrt_P, int k) {
	int i, j, z, w;
	int n = sqrt_P * k, x = 0;
	for (i = 0; i < sqrt_P; i++)
		for (j = 0; j < sqrt_P; j++)
			for (z = 0; z < k; z++)
				for (w = 0; w < k; w++) {
					printf("%d\n", (i*k + z) * n + (j*k + w));
					out[x++] = in[(i*k + z)*n + (j*k + w)];
				}
	return;
}


void makeRandomWorld(int *a, int n) {
	int i;
	srand(time(NULL));
	for (i = 0; i < n; i++){
		if (!(rand() % 2))
			a[i] = 1;
		else
			a[i] = 0;
	}
	return ;
}

void findNeighbors(int rank, int sqrt_P, int P, int* up, int* down, int* left,
				   int* right, int* up_left, int* up_right, int* down_left, int* down_right) {
	int Pmod = rank%sqrt_P;
	int Pdiv = rank/sqrt_P;  //rank = div*sP+mod
	*down  = ((Pdiv + 1) % sqrt_P)*sqrt_P + Pmod;
	*right  = Pdiv*sqrt_P + (Pmod + 1)%sqrt_P;
	*down_right = ((Pdiv + 1)%sqrt_P)*sqrt_P + (Pmod + 1)%sqrt_P;
	if (Pdiv == 0) {
		*up  = P - sqrt_P + Pmod;
		*up_right = P - sqrt_P + (Pmod + 1)%sqrt_P;
		if (Pmod == 0)
			*up_left = P - 1;
		else
			*up_left = P - sqrt_P + Pmod - 1;
	}
	else {
		*up  = (Pdiv - 1)*sqrt_P + Pmod;
		*up_right = (Pdiv - 1)*sqrt_P + (Pmod + 1)%sqrt_P;
                if (Pmod == 0)
                        *up_left = Pdiv*sqrt_P - 1;
                else
                        *up_left = (Pdiv - 1)*sqrt_P + Pmod - 1;
	}
	if (Pmod == 0) {
		*left  = (Pdiv + 1)*sqrt_P - 1;
		*down_left =((Pdiv + 1)%sqrt_P)*sqrt_P + sqrt_P - 1;
	}
	else {
		*left  = Pdiv*sqrt_P + Pmod - 1;
		*down_left = ((Pdiv + 1)%sqrt_P)*sqrt_P + Pmod - 1;
	}
	return ;
}

int find_new_a(int x, int neighbors) {
        if (x==1) {
        	if ((neighbors == 2) || (neighbors == 3)){
                return 1;}
            else{
                return 0;}
        }
        else { //x==0
            if (neighbors == 3){
                return 1;}
            else{
                return 0;}
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

void checkPerimeter(int* a, int n, int up_left, int up_right, int down_left, int down_right,
					int* up, int* down, int* left, int* right, int* b) {
	int i,neighbor;
	//for a[0]
	neighbor = up_left + up[0] + up[1] + left[0] + left[1] + a[1] + a[n] + a[n + 1];
	b[0] = find_new_a(a[0],neighbor);
	//for a[n - 1]
	neighbor = up_right + up[n - 1] + up[n - 2] + right[0] + right[1] + a[n - 2] + a[2*n - 1] + a[2*n - 2];
	b[n - 1] = find_new_a(a[n - 1],neighbor);
	//for a[n*(n - 1)]
	neighbor = down_left + left[n - 2] + left[n - 1] + down[0] + down[1] + a[n*(n - 1) + 1] + a[n*(n - 2)] + a[n*(n - 2) + 1];
	b[n*(n - 1)] = find_new_a(a[n*(n - 1)],neighbor);
	//for a[n*n - 1]
	neighbor = down_right + down[n - 1] + down[n - 2] + right[n - 1] + right[n - 2] + a[n*n - 2] + a[n*(n - 1) - 1] + a[n*(n - 1) - 2];
	b[n*n - 1] = find_new_a(a[n*n - 1],neighbor);

	//for up
	for (i = 1; i < n - 1; i++) {
		neighbor = up[i - 1] + up[i] + up[i + 1] + a[i - 1] + a[i + 1] + a[i - 1 + n] + a[i + n] + a[i + 1 + n];
		b[i] = find_new_a(a[i],neighbor);
	}
	//for down
	for (i = 1; i < n - 1; i++) {
		neighbor = down[i - 1] + down[i] + down[i + 1] + a[n*(n - 1) + i - 1] + a[n*(n - 1) + i + 1] +
				   a[n*(n - 2) + i] + a[n*(n - 2) + i - 1] + a[n*(n - 2) + i + 1];
		b[n*(n - 1) + i] = find_new_a(a[n*(n - 1) + i], neighbor);
	}

	//for left
	for (i = 1; i < n - 1; i++) {
		neighbor = left[i - 1] + left[i] + left[i + 1] + a[(i - 1)*n] + a[(i - 1)*n + 1] + a[i*n + 1] +
				   a[(i + 1)*n] + a[(i + 1)*n + 1];
		b[n*i] = find_new_a(a[n*i],neighbor);
	}
	//for right
	for (i = 1;i < n - 1; i++){
		neighbor = right[i - 1] + right[i] + right[i + 1] + a[i*n - 1] + a[(i + 2)*n - 1] +
				a[i*n - 2] + a[(i + 2)*n - 2] + a[(i + 1)*n - 2];
		b[(i + 1)*n - 1] = find_new_a(a[(i + 1)*n - 1], neighbor);
	}
	return;
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

void printWorld(int *world,int sq,int k){
    int i,j,x,y;
    for(i=0;i<sq;i++){
        for(j=0;j<k;j++){
            for(x=0;x<sq;x++){
                for(y=0;y<k;y++){
                    printf("%d ",world[y+k*k*x+k*j+sq*k*k*i]);
                }
            }
        putchar('\n');
        }
    }
    putchar('\n');
    return ;
}
