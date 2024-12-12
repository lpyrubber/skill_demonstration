#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(){
    int N=16;
    int NP=15;
    int NC=3;
    int *index, *tag;
    float *min, *temp;

    index=(int*)malloc(NP*sizeof(int));
    min=(float*)malloc(NP*sizeof(float));
    tag=(int*)malloc(NP*sizeof(int));
    temp=(float*)malloc(NP*sizeof(float));

    srand(time(NULL));
    for(int i=0; i<NP; i++){
        index[i]=rand()%NC;
        min[i]=(float)(rand())/RAND_MAX;
        printf("%d: %d %f\n",i, index[i], min[i]);
    }

    for(int ic=0; ic<NC; ic++){
        int tid=NC+1;
        float mind=1e16;
        for(int i=0; i<NP; i++){
            if((index[i]==ic)&&(min[i]<mind)){
                tid=i;
                mind=min[i];
            }
        }
        printf("NC=%d: id=%d, min=%f\n",ic, tid, mind);
    }

    for(int ic=0; ic<NC; ic++){
        for(int i=0; i<NP; i++){
            if(index[i]==ic){
                tag[i]=i;
                temp[i]=min[i];
            }else{
                tag[i]=NC+2;
                temp[i]=1e16;
            }
        }

        for(int stride=N>>1; stride>0; stride>>=1){
            for(int i=0; i<stride; i++){
                if(i+stride<NP){
                    if((tag[i]!=NC+2)&&(tag[i+stride]!=NC+2)){
                    if(temp[i+stride]<temp[i]){
                            temp[i]=temp[i+stride];
                            tag[i]=tag[i+stride];
                        }
                    }else if(tag[i+stride]!=NC+2){
                        temp[i]=temp[i+stride];
                        tag[i]=tag[i+stride];
                    }
                }
            }
        }
        printf("NC=%d\n",ic);
        for(int i=0; i<NP; i++){
            printf("%d: %d %f\n",i, tag[i], temp[i]);
        }
    }
    free(temp);
    free(tag);
    free(index);
    free(min);
    return 0;

}