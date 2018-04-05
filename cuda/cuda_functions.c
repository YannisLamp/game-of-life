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

void makeArrayInBlocks(int *in, int *out, int sq_of_P, int k) {
	int i, j, z, w;
	int n = sq_of_P * k, x = 0;
	for (i = 0; i < sq_of_P; i++)
		for (j = 0; j < sq_of_P; j++)
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

void swapArrays(int** a, int** b) {
	int *temp = *a;
	*a = *b;
	*b = temp;
	return;
}

void printWorld(int *world,int sq,int k) {
    int i,j,x,y;
    for(i = 0; i < sq; i++){
        for(j = 0; j < k; j++){
            for(x = 0; x < sq; x++){
                for(y = 0; y < k; y++){
                    printf("%d ",world[y+k*k*x+k*j+sq*k*k*i]);
                }
            }
        	putchar('\n');
        }
    }
    putchar('\n');
    return ;
}
