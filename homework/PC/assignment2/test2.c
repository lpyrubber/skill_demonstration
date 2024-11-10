#include <stdio.h>
#include <stdlib.h>
#include <time.h>



int main(){
    int array[5]={0,2,6,16,32};
    int n_list[4]={0,0,0,0};
    int offset[4]={0,0,0,0};
    int sum, temp, num, start, red, index,np, npf=2, i, j, end;
    np=4;

    index=0;
 
    start=0;
    end=0;
    temp=array[np]/npf;
    for(i=0; i<npf; i++){
        for(j=0; j<np; j++){
            n_list[j]=0;
        }
        num=(i<(array[np]-temp*np))?temp+1:temp;
        end+=num;
        while(end>(array[index+1])){
            n_list[index]=array[index+1]-start;
            start=array[index+1];
            index++;
        }
        n_list[index]=end-start; 
        start=end;
        for(j=0; j<np; j++){
            offset[j]=(j==0)?0:offset[j-1]+n_list[j-1];
            if(offset[j]==num){
                offset[j]--;
            }
        }
        printf("%d: ",i);
        for(j=0; j<np; j++){
            printf("%d ",n_list[j]);
        }
        printf("\n");
        printf("%d: ",i);
        for(j=0; j<np; j++){
            printf("%d ",offset[j]);
        }
        printf("\n");
    }
}