#include "main.h"
__global__ void GPUDOT(float *d_a, float *d_b, float *d_c, int ele);
__global__ void GPUSUM(float *d_a, float *d_b, int ele);
__global__ void GPUMDOT(float *d_a, float *d_b, float *d_c, int ele);
__global__ void GPU_UP(float *d_a, float *d_b, float *d_c, float *c_d, int ele);
__global__ void GPU_US(float *d_a, float *d_b, float *d_c, float *c_d, int ele);
__global__ void GPU_UX(float *d_a, float *d_b, float *d_c, float *c_d, int ele);
__global__ void GPU_UR(float *d_a, float *d_b, float *d_c, float *c_d, int ele);

__global__ void GPUDOT(float *d_a, float *d_b, float *d_c, int ele){
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	int I = threadIdx.x;
	__shared__ float temp[TPB];
	temp[I]=0;
	if(i<ele) temp[I]=d_a[i]*d_b[i];
	__syncthreads();
	for(int stride = blockDim.x/2; stride > 0; stride = stride/2){
		if(I<stride){
			temp[I]+=temp[I+stride];
		}
		__syncthreads();
	} 
	if(I==0) d_c[blockIdx.x]=temp[0];
	__syncthreads();
	if(i == 0 ){
		float sum = 0.0;
		for(int i1 = 0 ; i1 < gridDim.x; i1 += 2){
			sum += d_c[i1];
		}
		d_c[0] = sum;
	}else if(i == 1){
		float sum = 0.0;
		for(int i1 = 0 ; i1 < gridDim.x; i1 += 2){
			sum += d_c[i1];
		}
		d_c[1] = sum;
	}
	__syncthreads();
	if( i == 0 ) d_c[0] += d_c[1];
}

__global__ void GPUSUM(float *d_a, float *d_b, int ele){
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	float sum =0;
	if(i==0){
		for(int I=0; I<ele; I++){
			sum+=d_a[I];
		}
		*d_b = sum;
	}
}

__global__ void GPUMDOT(float *d_a, float *d_b, float *d_c, int ele){
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	float sum = 0;
	if(i<ele){
		for(int k=0; k<ele; k++){
			sum += d_a[i*N+k]*d_b[k];
		}
		d_c[i] = sum;
	}
}

__global__ void GPU_UP(float *d_a, float *d_b, float *d_c, float *d_d, int ele){
	int i=threadIdx.x+blockDim.x*blockIdx.x;
	if(i<ele){
		if(i==0) d_d[3]=(d_d[0]/d_d[1])*(d_d[2]/d_d[5]);
		__syncthreads(); 
		d_a[i]=d_b[i]+d_d[3]*(d_a[i]-d_c[i]*d_d[5]);
	}
}

__global__ void GPU_US(float *d_a, float *d_b, float *d_c, float *d_d, int ele){
	int i=threadIdx.x+blockDim.x*blockIdx.x;
	if(i<ele){
		if(i==0) d_d[2]=d_d[0]/d_d[6];
		__syncthreads(); 
		d_a[i]=d_b[i]-d_d[2]*d_c[i];
	}
}

__global__ void GPU_UX(float *d_a, float *d_b, float *d_c, float *d_d, int ele){
	int i=threadIdx.x+blockDim.x*blockIdx.x;
	if(i<ele){
		if(i==0){ 
			d_d[4]=d_d[6]/d_d[7];
			d_d[5]=d_d[4];
		}
		__syncthreads(); 
		d_a[i]=d_a[i]+d_d[2]*d_b[i]+d_d[4]*d_c[i];
	}
}

__global__ void GPU_UR(float *d_a, float *d_b, float *d_c, float *d_d, int ele){
	int i=threadIdx.x+blockDim.x*blockIdx.x;
	if(i<ele){
		if(i==0) d_d[1]=d_d[0];
		__syncthreads(); 
		d_a[i]=d_b[i]-d_d[4]*d_c[i];
	}
}

void Allocate_Memory(float **h_a, float **h_b, float **h_c){
	size_t size;
	size = N*N*sizeof(float);
	*h_a = (float*)malloc(size);
	size = N*sizeof(float);
	*h_b = (float*)malloc(size);
	*h_c = (float*)malloc(size);
}

void Allocate_Device_Memory(float **d_a, float **d_b, float **d_c, float **d_d, float **d_e, float **d_f, float **d_g, float **d_h, float **d_i, float **d_j){
	size_t size;
	cudaError_t Error;
	size = N*N*sizeof(float);
	Error = cudaMalloc((void**)d_a, size);
	if(DEBUG) printf("CUDA Error (malloc d_a)=%s\n",cudaGetErrorString(Error));
	size = N*sizeof(float);
	Error = cudaMalloc((void**)d_b, size);
	if(DEBUG) printf("CUDA Error (malloc d_b)=%s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)d_c, size);
	if(DEBUG) printf("CUDA Error (malloc d_c)=%s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)d_d, size);
	if(DEBUG) printf("CUDA Error (malloc d_d)=%s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)d_e, size);
	if(DEBUG) printf("CUDA Error (malloc d_e)=%s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)d_f, size);
	if(DEBUG) printf("CUDA Error (malloc d_f)=%s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)d_g, size);
	if(DEBUG) printf("CUDA Error (malloc d_g)=%s\n",cudaGetErrorString(Error));
	Error = cudaMalloc((void**)d_h, size);
	if(DEBUG) printf("CUDA Error (malloc d_h)=%s\n",cudaGetErrorString(Error));
	size = ((int)((N+TPB-1)/TPB))*sizeof(float);
	Error = cudaMalloc((void**)d_i, size);
	if(DEBUG) printf("CUDA Error (malloc d_i)=%s\n",cudaGetErrorString(Error));
	size =  8*sizeof(float);
	Error = cudaMalloc((void**)d_j, size);
	if(DEBUG) printf("CUDA Error (malloc d_j)=%s\n",cudaGetErrorString(Error));
}

void Free_Memory(float **h_a, float **h_b, float **h_c){
	free(*h_a);	
	free(*h_b);	
	free(*h_c);	
	
}

void Free_Device_Memory(float **d_a, float **d_b, float **d_c, float **d_d, float **d_e, float **d_f, float **d_g, float **d_h, float **d_i, float **d_j){
	cudaError_t Error;
	Error = cudaFree(*d_a);
	if(DEBUG) printf("CUDA Error (free d_a)=%s\n", cudaGetErrorString(Error));
	Error = cudaFree(*d_b);
	if(DEBUG) printf("CUDA Error (free d_b)=%s\n", cudaGetErrorString(Error));
	Error = cudaFree(*d_c);
	if(DEBUG) printf("CUDA Error (free d_c)=%s\n", cudaGetErrorString(Error));
	Error = cudaFree(*d_d);
	if(DEBUG) printf("CUDA Error (free d_d)=%s\n", cudaGetErrorString(Error));
	Error = cudaFree(*d_e);
	if(DEBUG) printf("CUDA Error (free d_e)=%s\n", cudaGetErrorString(Error));
	Error = cudaFree(*d_f);
	if(DEBUG) printf("CUDA Error (free d_f)=%s\n", cudaGetErrorString(Error));
	Error = cudaFree(*d_g);
	if(DEBUG) printf("CUDA Error (free d_g)=%s\n", cudaGetErrorString(Error));
	Error = cudaFree(*d_h);
	if(DEBUG) printf("CUDA Error (free d_h)=%s\n", cudaGetErrorString(Error));
	Error = cudaFree(*d_i);
	if(DEBUG) printf("CUDA Error (free d_i)=%s\n", cudaGetErrorString(Error));
	Error = cudaFree(*d_j);
	if(DEBUG) printf("CUDA Error (free d_j)=%s\n", cudaGetErrorString(Error));
}

void Copy_To_Device(float *h_a, float *h_b, float *h_c, float *d_a, float *d_b, float *d_c, float *d_d, float *d_e, float *d_f, float *d_g){
	cudaError_t Error;
	size_t size;
	float temp=1;
	size = N*N*sizeof(float);
	Error = cudaMemcpy(d_a, h_a, size, cudaMemcpyHostToDevice);
	if(DEBUG) printf("CUDA Error (h_a -> d_a)=%s\n", cudaGetErrorString(Error));
	size = N*sizeof(float);
	Error = cudaMemcpy(d_b, h_b, size, cudaMemcpyHostToDevice);
	if(DEBUG) printf("CUDA Error (h_b -> d_b)=%s\n", cudaGetErrorString(Error));
	Error = cudaMemcpy(d_c, h_b, size, cudaMemcpyHostToDevice);
	if(DEBUG) printf("CUDA Error (h_b -> d_c)=%s\n", cudaGetErrorString(Error));
	Error = cudaMemcpy(d_d, h_b, size, cudaMemcpyHostToDevice);
	if(DEBUG) printf("CUDA Error (h_b -> d_d)=%s\n", cudaGetErrorString(Error));
	Error = cudaMemcpy(d_e, h_c, size, cudaMemcpyHostToDevice);
	if(DEBUG) printf("CUDA Error (h_c -> d_e)=%s\n", cudaGetErrorString(Error));
	Error = cudaMemcpy(d_f, h_c, size, cudaMemcpyHostToDevice);
	if(DEBUG) printf("CUDA Error (h_c -> d_f)=%s\n", cudaGetErrorString(Error));
	size = sizeof(float);	
	Error = cudaMemcpy(d_g+1, &temp, size, cudaMemcpyHostToDevice);
	if(DEBUG) printf("CUDA Error (1 -> rho_old)=%s\n", cudaGetErrorString(Error));
	Error = cudaMemcpy(d_g+2, &temp, size, cudaMemcpyHostToDevice);
	if(DEBUG) printf("CUDA Error (1 -> alpha)=%s\n", cudaGetErrorString(Error));
	Error = cudaMemcpy(d_g+5, &temp, size, cudaMemcpyHostToDevice);
	if(DEBUG) printf("CUDA Error (1 -> omega_old)=%s\n", cudaGetErrorString(Error));
} 

void Copy_To_Host(float *h_a, float *d_a){
	cudaError_t Error;
	size_t size;
	size = N*sizeof(float);
	Error = cudaMemcpy(h_a, d_a, size, cudaMemcpyDeviceToHost);
	if(DEBUG) printf("CUDA Error (d_a -> h_a)=%s\n", cudaGetErrorString(Error));

}

void Initial(float *h_a, float *h_b, float *h_c){
	int i,j;
	for(j=0; j<N; j++){
		for(i=0; i<N; i++){
			h_b[j]=0;
			h_c[j]=(j==0||j==N-1)?0:1*dx*dx;
			if(i==j){
				h_a[i+N*j]=2.0;
			}else if((j==i-1)&&(j>0)){
				h_a[i+N*j]=-1.0;
			}else if((j==i+1)&&(j<N-1)){
				h_a[i+N*j]=-1.0;
			}else{
				h_a[i+N*j]=0;
			}
		}
	}
}

void Save_Result(float *h_a){
	int i;
	FILE *fp;
	fp = fopen("data.txt", "w");
	for(i=0; i<N; i++){
		fprintf(fp, "%f\n", h_a[i]);
	}
	fclose(fp);

}

void GPU_VV_Dot(float *d_a, float *d_b, float *d_c, float *d_d, int ele){
	int threadsperblock = TPB;
	int blockspergrid = (ele+threadsperblock-1)/threadsperblock;
	GPUDOT<<<blockspergrid, threadsperblock>>>(d_a, d_b, d_c, ele);
	GPUSUM<<<1,1>>>(d_c,d_d, blockspergrid);
}

void GPU_MV_Dot(float *d_a, float *d_b, float *d_c, int ele){
	int threadsperblock = TPB;
	int blockspergrid = (ele+threadsperblock-1)/threadsperblock;
	GPUMDOT<<<blockspergrid, threadsperblock>>>(d_a, d_b, d_c, ele);	
}

void GPU_Update_P(float *d_a, float *d_b, float *d_c, float *d_d, int ele){
	int threadsperblock = TPB;
	int blockspergrid = (ele+threadsperblock-1)/threadsperblock;
	GPU_UP<<<blockspergrid, threadsperblock>>>(d_a, d_b, d_c, d_d, ele);
}	

void GPU_Update_S(float *d_a, float *d_b, float *d_c, float *d_d, int ele){
	int threadsperblock = TPB;
	int blockspergrid = (ele+threadsperblock-1)/threadsperblock;
	GPU_US<<<blockspergrid, threadsperblock>>>(d_a, d_b, d_c, d_d, ele);
}

void GPU_Update_X(float *d_a, float *d_b, float *d_c, float *d_d, int ele){
	int threadsperblock = TPB;
	int blockspergrid = (ele+threadsperblock-1)/threadsperblock;
	GPU_UX<<<blockspergrid, threadsperblock>>>(d_a, d_b, d_c, d_d, ele);
}
void GPU_Update_R(float *d_a, float *d_b, float *d_c, float *d_d, int ele){
	int threadsperblock = TPB;
	int blockspergrid = (ele+threadsperblock-1)/threadsperblock;
	GPU_UR<<<blockspergrid, threadsperblock>>>(d_a, d_b, d_c, d_d, ele);
}

void GPU_Update_temp(float *h_a, float *d_a){
	size_t size = sizeof(float);
	cudaError_t Error;
	Error = cudaMemcpy(h_a, d_a, size, cudaMemcpyDeviceToHost);
	if(DEBUG) printf("CUDA Error (Memcpy temp)=%s\n", cudaGetErrorString(Error));
}


