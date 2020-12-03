#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <wait.h>
#include <math.h>

double **matd, **matgx, **matgy, **matg;
int c, f, sharedmatd, sharedmatgx, sharedmatgy, sharedmatg, n, ix = 0, jx = 0, iy = 0, jy = 0;

pthread_t *hilos;
pthread_mutex_t pmatgx = PTHREAD_MUTEX_INITIALIZER, pmatgy = PTHREAD_MUTEX_INITIALIZER;

unsigned int getSize(int f, int c, size_t sizeElement){
    size_t size = f * (sizeof(void *) + (c * sizeElement));
    return size;
}

void indexarMatriz(void **m, int f, int c, size_t sizeElement){
    size_t sizeRow = c * sizeElement;
    m[0] = m+f;
    for(int i = 1; i < c; i++){      
        m[i] = (m[i-1]+sizeRow);
    }
}

void imprimirMatriz(double **m, int f, int c){
    for(int i = 0; i < f; i++){
        for(int j = 0; j < c; j++)
            printf("%.1f\t",m[i][j]);
        printf("\n");
    }
}

void* calcularGx(void *arg){
    int ci, cj;
    do{
        pthread_mutex_lock(&pmatgx);
            if(ix < f){
                ci = ix;
                cj = jx;
                jx++;
                if(jx == c){
                    jx = 0;
                    ix++;
                }
            }
            
        pthread_mutex_unlock(&pmatgx);

        if(ix == f) return NULL;
        matgx[ci][cj] = ( ((ci+1 < f? matd[ci+1][cj] : 1)-(ci-1 >= 0? matd[ci-1][cj] : 1))/2);
    }while(ix < f);
}

void* calcularGy(void *arg){
    int ci, cj;
    do{
        pthread_mutex_lock(&pmatgy);
            if(iy < f){
                ci = iy;
                cj = jy;
                jy++;
                if(jy == c){
                    jy = 0;
                    iy++;
                }
            }
            
        pthread_mutex_unlock(&pmatgy);

        if(iy == f) return NULL;

        matgy[ci][cj] = ( ((cj+1 < c? matd[ci][cj+1] : 1)-(cj-1 >= 0? matd[ci][cj-1] : 1))/2);
    }while(iy < f);
}

void calcularG(){
    for(int i = 0; i < f; i++){
        for(int j = 0; j < c; j++){
            matg[i][j] = sqrt(pow(matgx[i][j],2) + pow(matgy[i][j], 2));
        }
    }
}

int main(){
    system("clear");
    FILE *archivo = fopen("datos.txt","r");
    
    if(!archivo){
        printf("error\n");
        return 0;
    }

    fscanf(archivo,"%d",&f);
    fscanf(archivo,"%d",&c);

    size_t sizeMatrix = getSize(f, c, sizeof(double));

    sharedmatd = shmget(IPC_PRIVATE, sizeMatrix, IPC_CREAT | S_IRUSR | S_IWUSR);
    sharedmatgx = shmget(IPC_PRIVATE, sizeMatrix, IPC_CREAT | S_IRUSR | S_IWUSR);
    sharedmatgy = shmget(IPC_PRIVATE, sizeMatrix, IPC_CREAT | S_IRUSR | S_IWUSR);
    sharedmatg = shmget(IPC_PRIVATE, sizeMatrix, IPC_CREAT | S_IRUSR | S_IWUSR);


    matd = shmat(sharedmatd, NULL, 0);
    matgx = shmat(sharedmatgx, NULL, 0);
    matgy = shmat(sharedmatgy, NULL, 0);
    matg = shmat(sharedmatg, NULL, 0);

    indexarMatriz((void*)matd, f, c, sizeof(double));
    indexarMatriz((void*)matg, f, c, sizeof(double));
    indexarMatriz((void*)matgx, f, c, sizeof(double));
    indexarMatriz((void*)matgy, f, c, sizeof(double));

    for(int i = 0; i < f; i++){
        for(int j = 0; j < c; j++){
            fscanf(archivo, "%lf", &matd[i][j]);
            matg[i][j] = -1;
            matgx[i][j] = -1;
            matgy[i][j] = -1;
        }
    }

    fclose(archivo);
    printf("\n\n\t\tMatriz de datos\n");
    imprimirMatriz(matd, f, c);

    printf("\n\nNumero par de hilos: ");
    leerCantidad:
    scanf("%d",&n);

    if(n < 2 || n % 2 != 0){
        printf("El cantidad de hilos debe ser un numero par positivo... ");
        goto leerCantidad;
    }

    int z = 0;
    hilos = (pthread_t*)calloc(n, sizeof(pthread_t));

    for(z = 0; z < (n/2); z++){
        pthread_create(&hilos[z], NULL, calcularGx, NULL);
    }

    for(z = 0 ;z < n; z++){
        pthread_create(&hilos[z], NULL, calcularGy, NULL);
    }

    for(z = 0; z < n; z++){
        pthread_join(hilos[z], NULL);
    }

    printf("\n\n\t\tGX\n");
    imprimirMatriz(matgx, f, c);

    printf("\n\n\t\tGY\n");
    imprimirMatriz(matgy, f, c);

    calcularG();
    printf("\n\n\t\tMatriz G\n");
    imprimirMatriz(matg, f, c);

    shmdt(matd);
    shmdt(matg);
    shmdt(matgx);
    shmdt(matgy);

    shmctl(sharedmatg, IPC_RMID, 0);
    shmctl(sharedmatgx, IPC_RMID, 0);
    shmctl(sharedmatgy, IPC_RMID, 0);
    shmctl(sharedmatd, IPC_RMID, 0);
    
    return 0;   
}
