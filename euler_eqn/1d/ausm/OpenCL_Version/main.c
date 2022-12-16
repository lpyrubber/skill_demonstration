#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#ifdef __APPLE__
    #include <OpenCL/cl.h>
#else
    #include <CL/cl.h>
#endif

#define L 1.0
#define N 20000
#define DX (L/N)
#define DT (0.01*DX)
#define Z (DT/DX)
#define NO_STEP 320000
#define gamma 1.4
#define R 1.0
#define CV (R/(gamma-1))
#define CP (CV+R)

#define LOCAL 256
#define GLOBAL (((N+LOCAL-1)/LOCAL)*LOCAL)
#define PLATFORM 0
#define DEVICE 0

#define CaseReturnString(x) case x: return #x;

char* clErrorStr(cl_int err);

int main(int argc, char* argv[]){
    float *rho, *v, *T, Z_cl;
    float *u;
    FILE *fp;
    char *source_str;
    size_t source_size, program_size, array_size, global_size, local_size;
    int i, j, N_cl;
    clock_t start_t, end_t;
    double total_t;

    start_t=clock();
    Z_cl = Z;
    N_cl = N;
    global_size = GLOBAL;
    local_size = LOCAL;

    //variable for OpenCL

    cl_platform_id *platforms = NULL;
    cl_device_id *device_list = NULL;
    cl_context context;
    cl_command_queue command_queue;
    cl_mem U_clmem, RHO_clmem, V_clmem, T_clmem, ML_clmem;
    cl_mem MR_clmem, PL_clmem, PR_clmem, F_clmem;
    cl_program program;
    cl_kernel kernel;
    cl_uint num_platforms, num_devices;
    cl_int clStatus;

    //load program
    fp = fopen( "kernel.cl", "rb");
	if( !fp ){
		printf( "Error: failed to open cl file\n");
	}
	fseek( fp , 0 , SEEK_END );
	program_size = ftell( fp );
	rewind( fp );
	source_size = program_size / sizeof( char ) + 1;
	source_str = ( char* )malloc( source_size * sizeof( char ) );
	source_str[ program_size ] = '\0';
	fread( source_str , sizeof(char) , source_size - 1 , fp );
	fclose(fp);

	//allocate memory
	array_size = N*sizeof(float);
	rho = (float*)malloc(array_size);
	v = (float*)malloc(array_size);
	T = (float*)malloc(array_size);
	u = (float*)malloc(3*array_size);

	//initial
	for(i=0; i<N; ++i){
		if((i-1)<0.5*N){
			rho[i] = 10.0;
			v[i] = 0.0;
			T[i] = 1.0;
		}else{
			rho[i] = 1.0;
			v[i] = 0.0;
			T[i] = 1.0;
		}
		u[i]=rho[i];
		u[i+N]=rho[i]*v[i];
		u[i+N*2]=rho[i]*(CV*T[i]+0.5*v[i]*v[i]);
	}

	//get platform and device information
	clStatus = clGetPlatformIDs(0, NULL, &num_platforms);
	platforms = (cl_platform_id*)malloc(num_platforms*sizeof(cl_platform_id));
	clStatus = clGetPlatformIDs(num_platforms, platforms, NULL);
	printf("num_platforms=%d\n",(unsigned int)num_platforms);

	//get the devices list
	clStatus = clGetDeviceIDs(platforms[PLATFORM], CL_DEVICE_TYPE_GPU, 1, NULL, &num_devices);
	device_list = (cl_device_id*)malloc(num_devices*sizeof(cl_device_id));
	clStatus = clGetDeviceIDs(platforms[PLATFORM], CL_DEVICE_TYPE_GPU, num_devices, device_list, NULL);
	printf("num_devices=%d, choose device #%d\n", num_devices, DEVICE);
	//create context and command queue
	context = clCreateContext(NULL, num_devices, device_list, NULL, NULL, &clStatus );
#ifdef __APPLE__
	command_queue = clCreateCommandQueue(context, device_list[DEVICE], 0, &clStatus);
#else
	command_queue = clCreateCommandQueueWithProperties(context, device_list[DEVICE], 0, &clStatus);
#endif
	if(clStatus != CL_SUCCESS){
		printf("Error: Failed to create command queue!\n%s\n", clErrorStr(clStatus));
	}
	//create memory buffers on the device and copy data to the device
	RHO_clmem = clCreateBuffer(context, CL_MEM_READ_WRITE, array_size, NULL, &clStatus);
	V_clmem   = clCreateBuffer(context, CL_MEM_READ_WRITE, array_size, NULL, &clStatus);
	T_clmem   = clCreateBuffer(context, CL_MEM_READ_WRITE, array_size, NULL, &clStatus);
	U_clmem   = clCreateBuffer(context, CL_MEM_READ_WRITE, 3*array_size, NULL, &clStatus);
    F_clmem   = clCreateBuffer(context, CL_MEM_READ_WRITE, 3*array_size, NULL, &clStatus);
    PL_clmem   = clCreateBuffer(context, CL_MEM_READ_WRITE, array_size, NULL, &clStatus);
    PR_clmem   = clCreateBuffer(context, CL_MEM_READ_WRITE, array_size, NULL, &clStatus);
    ML_clmem   = clCreateBuffer(context, CL_MEM_READ_WRITE, array_size, NULL, &clStatus);
    MR_clmem   = clCreateBuffer(context, CL_MEM_READ_WRITE, array_size, NULL, &clStatus);
	if(!RHO_clmem||!V_clmem||!T_clmem||!U_clmem||!F_clmem||!PL_clmem||!PR_clmem||!ML_clmem||!MR_clmem){
		printf("Error: Failed to allocate device memory!\n%s\n", clErrorStr(clStatus));
	}
	//copy to the device
	clStatus  = clEnqueueWriteBuffer(command_queue, RHO_clmem, CL_TRUE, 0, array_size, rho, 0, NULL, NULL);
	clStatus |= clEnqueueWriteBuffer(command_queue, V_clmem, CL_TRUE, 0, array_size, v, 0, NULL, NULL);
	clStatus |= clEnqueueWriteBuffer(command_queue, T_clmem, CL_TRUE, 0, array_size, T, 0, NULL, NULL);
	clStatus |= clEnqueueWriteBuffer(command_queue, U_clmem, CL_TRUE, 0, 3*array_size, u, 0, NULL, NULL);
	if(clStatus!=CL_SUCCESS){
		printf("Error: Failed to sent host memory to device array!\n%s\n", clErrorStr(clStatus));
	}
	//create a program from the kernel source and debugging
	program = clCreateProgramWithSource(context, 1, (const char**)&source_str, NULL, &clStatus);
	if(clStatus!=CL_SUCCESS){
		printf("Error: Failed to create program form the kernel source!\n%s\n", clErrorStr(clStatus));
		exit(1);
	}
	clStatus=clBuildProgram(program, 1, device_list+DEVICE , NULL, NULL, NULL);
	if(clStatus!=CL_SUCCESS){
		size_t len;
		char buffer[2048];
		printf("Error:Failed to build program executable!\n%s\n", clErrorStr(clStatus));
		clGetProgramBuildInfo(program, device_list[DEVICE], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		exit(1);
	}
	kernel=clCreateKernel(program, "GPU_Calc", &clStatus);
	if(!kernel||clStatus!=CL_SUCCESS){
		printf("Error: Failed to create compute kernel!\n%s\n", clErrorStr(clStatus));
		exit(1);
	}
	for( i = 0 ; i < NO_STEP; ++i){
		//Input argument of the kernel
		clStatus  = clSetKernelArg( kernel, 0, sizeof(int), (void*) &N_cl );
		clStatus |= clSetKernelArg( kernel, 1, sizeof(float), (void*) &Z_cl);
		clStatus |= clSetKernelArg( kernel, 2, sizeof(cl_mem) , (void*) &U_clmem );
		clStatus |= clSetKernelArg( kernel, 3, sizeof(cl_mem) , (void*) &RHO_clmem );
		clStatus |= clSetKernelArg( kernel, 4, sizeof(cl_mem) , (void*) &V_clmem );
		clStatus |= clSetKernelArg( kernel, 5, sizeof(cl_mem) , (void*) &T_clmem );
		clStatus |= clSetKernelArg( kernel, 6, sizeof(cl_mem) , (void*) &F_clmem );
		clStatus |= clSetKernelArg( kernel, 7, sizeof(cl_mem) , (void*) &PL_clmem );
        clStatus |= clSetKernelArg( kernel, 8, sizeof(cl_mem) , (void*) &PR_clmem );
        clStatus |= clSetKernelArg( kernel, 9, sizeof(cl_mem) , (void*) &ML_clmem );
        clStatus |= clSetKernelArg( kernel,10, sizeof(cl_mem) , (void*) &MR_clmem );
		if(clStatus != CL_SUCCESS){
			printf("Error: Failed to set argument for kernels!\n%s\n", clErrorStr(clStatus));
			exit(1);
		}
		//Execute the OpenCL Kernel on the queue
		clStatus = clEnqueueNDRangeKernel( command_queue, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
		if(clStatus != CL_SUCCESS){
			printf("Error: Failed to excute kernel!\n%s\n", clErrorStr(clStatus));
			exit(1);
		}
	}
	//Send the device memory to the host
	clStatus  = clEnqueueReadBuffer( command_queue, U_clmem, CL_TRUE, 0, 3*array_size, u, 0, NULL, NULL);
	clStatus |= clEnqueueReadBuffer( command_queue, RHO_clmem, CL_TRUE, 0, array_size, rho, 0, NULL, NULL);
	clStatus |= clEnqueueReadBuffer( command_queue, V_clmem, CL_TRUE, 0, array_size, v, 0, NULL, NULL);
	clStatus |= clEnqueueReadBuffer( command_queue, T_clmem, CL_TRUE, 0, array_size, T, 0, NULL, NULL);
	if(clStatus != CL_SUCCESS){
		printf("Error: Failed to send device memory to host memory!\n%s\n", clErrorStr(clStatus));
		exit(1);
	}
	//Save Result
	fp = fopen("data.txt","w");
	for(i=0; i<N; ++i){
		fprintf(fp, "%e %e %e %e\n", i*DX, rho[i], v[i], T[i] );
	}
	fclose(fp);
	clStatus = clReleaseKernel( kernel );
	clStatus = clReleaseProgram( program );
	clStatus = clReleaseMemObject( U_clmem );
	clStatus = clReleaseMemObject( F_clmem );
	clStatus = clReleaseMemObject( PL_clmem );
    clStatus = clReleaseMemObject( PR_clmem );
    clStatus = clReleaseMemObject( ML_clmem );
    clStatus = clReleaseMemObject( MR_clmem );
	clStatus = clReleaseMemObject( RHO_clmem );
	clStatus = clReleaseMemObject( V_clmem );
	clStatus = clReleaseMemObject( T_clmem );
	clStatus = clReleaseCommandQueue( command_queue );
	clStatus = clReleaseContext( context );

	free( source_str );
	free( u );
	free( rho );
	free( v );
	free( T );
	free( platforms );
	free( device_list );
	
	end_t = clock();
	total_t = (double)(end_t-start_t)/CLOCKS_PER_SEC;
	printf("CPU runtime = %lf sec\n",total_t);
	return 0;

}

char *clErrorStr(cl_int err)
{
    switch (err)
    {
        CaseReturnString(CL_SUCCESS                        )
        CaseReturnString(CL_DEVICE_NOT_FOUND               )
        CaseReturnString(CL_DEVICE_NOT_AVAILABLE           )
        CaseReturnString(CL_COMPILER_NOT_AVAILABLE         )
        CaseReturnString(CL_MEM_OBJECT_ALLOCATION_FAILURE  )
        CaseReturnString(CL_OUT_OF_RESOURCES               )
        CaseReturnString(CL_OUT_OF_HOST_MEMORY             )
        CaseReturnString(CL_PROFILING_INFO_NOT_AVAILABLE   )
        CaseReturnString(CL_MEM_COPY_OVERLAP               )
        CaseReturnString(CL_IMAGE_FORMAT_MISMATCH          )
        CaseReturnString(CL_IMAGE_FORMAT_NOT_SUPPORTED     )
        CaseReturnString(CL_BUILD_PROGRAM_FAILURE          )
        CaseReturnString(CL_MAP_FAILURE                    )
        CaseReturnString(CL_MISALIGNED_SUB_BUFFER_OFFSET   )
        CaseReturnString(CL_COMPILE_PROGRAM_FAILURE        )
        CaseReturnString(CL_LINKER_NOT_AVAILABLE           )
        CaseReturnString(CL_LINK_PROGRAM_FAILURE           )
        CaseReturnString(CL_DEVICE_PARTITION_FAILED        )
        CaseReturnString(CL_KERNEL_ARG_INFO_NOT_AVAILABLE  )
        CaseReturnString(CL_INVALID_VALUE                  )
        CaseReturnString(CL_INVALID_DEVICE_TYPE            )
        CaseReturnString(CL_INVALID_PLATFORM               )
        CaseReturnString(CL_INVALID_DEVICE                 )
        CaseReturnString(CL_INVALID_CONTEXT                )
        CaseReturnString(CL_INVALID_QUEUE_PROPERTIES       )
        CaseReturnString(CL_INVALID_COMMAND_QUEUE          )
        CaseReturnString(CL_INVALID_HOST_PTR               )
        CaseReturnString(CL_INVALID_MEM_OBJECT             )
        CaseReturnString(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR)
        CaseReturnString(CL_INVALID_IMAGE_SIZE             )
        CaseReturnString(CL_INVALID_SAMPLER                )
        CaseReturnString(CL_INVALID_BINARY                 )
        CaseReturnString(CL_INVALID_BUILD_OPTIONS          )
        CaseReturnString(CL_INVALID_PROGRAM                )
        CaseReturnString(CL_INVALID_PROGRAM_EXECUTABLE     )
        CaseReturnString(CL_INVALID_KERNEL_NAME            )
        CaseReturnString(CL_INVALID_KERNEL_DEFINITION      )
        CaseReturnString(CL_INVALID_KERNEL                 )
        CaseReturnString(CL_INVALID_ARG_INDEX              )
        CaseReturnString(CL_INVALID_ARG_VALUE              )
        CaseReturnString(CL_INVALID_ARG_SIZE               )
        CaseReturnString(CL_INVALID_KERNEL_ARGS            )
        CaseReturnString(CL_INVALID_WORK_DIMENSION         )
        CaseReturnString(CL_INVALID_WORK_GROUP_SIZE        )
        CaseReturnString(CL_INVALID_WORK_ITEM_SIZE         )
        CaseReturnString(CL_INVALID_GLOBAL_OFFSET          )
        CaseReturnString(CL_INVALID_EVENT_WAIT_LIST        )
        CaseReturnString(CL_INVALID_EVENT                  )
        CaseReturnString(CL_INVALID_OPERATION              )
        CaseReturnString(CL_INVALID_GL_OBJECT              )
        CaseReturnString(CL_INVALID_BUFFER_SIZE            )
        CaseReturnString(CL_INVALID_MIP_LEVEL              )
        CaseReturnString(CL_INVALID_GLOBAL_WORK_SIZE       )
        CaseReturnString(CL_INVALID_PROPERTY               )
        CaseReturnString(CL_INVALID_IMAGE_DESCRIPTOR       )
        CaseReturnString(CL_INVALID_COMPILER_OPTIONS       )
        CaseReturnString(CL_INVALID_LINKER_OPTIONS         )
        CaseReturnString(CL_INVALID_DEVICE_PARTITION_COUNT )
        default: return "Unknown OpenCL error code";
    }
}
