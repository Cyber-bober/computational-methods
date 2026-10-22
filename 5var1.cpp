#include <cstdio>
#include <cmath>
#include <algorithm>

using namespace std;

int main(){
    const int N = 12;
    double A[N][N+1] = {
        {1,0,0,0, 0,0,0,0, 0,0,0,0, 3},
        {0,0,0,0, 1,0,0,0, 0,0,0,0, 10},
        {0,0,0,0, 0,0,0,0, 1,0,0,0, 29},
        {1,1,1,1, 0,0,0,0, 0,0,0,0, 10},
        {0,0,0,0, 1,1,1,1, 0,0,0,0, 29},
        {0,0,0,0, 0,0,0,0, 1,1,1,1, 66},
        {0,1,2,3,-1,0,0,0, 0,0,0,0, 0},
        {0,0,0,0, 0,1,2,3,-1,0,0,0, 0},
        {0,0,2,6, 0,0,-2,0, 0,0,0,0, 0},
        {0,0,0,0, 0,0,2,6, 0,0,-2,0, 0},
        {0,0,2,0, 0,0,0,0, 0,0,0,0, 0},
        {0,0,0,0, 0,0,0,0, 0,0,2,6, 0}
    };

    // Прямой ход Гаусса
    for(int i=0;i<N;i++){
        int p=i;
        for(int j=i+1;j<N;j++)
            if(fabs(A[j][i])>fabs(A[p][i])) p=j;
        swap(A[i],A[p]);
        for(int j=i+1;j<N;j++){
            double k=A[j][i]/A[i][i];
            for(int t=i;t<=N;t++) A[j][t]-=k*A[i][t];
        }
    }

    // Обратный ход
    double x[N];
    for(int i=N-1;i>=0;i--){
        double s=A[i][N];
        for(int j=i+1;j<N;j++) s-=A[i][j]*x[j];
        x[i]=s/A[i][i];
    }

    const char* names[]={"a","b","c","d"};
    for(int i=0;i<3;i++){
        printf("S%d: ",i+1);
        for(int j=0;j<4;j++)
            printf("%s%d=%g ",names[j],i+1,x[i*4+j]);
        printf("\n");
    }
}
