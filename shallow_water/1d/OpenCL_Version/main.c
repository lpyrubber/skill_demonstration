#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef __APPLE__
    #include <OpenCL/cl.h>
#else
    #include <CL/cl.h>
#endif

#define CaseReturnString(x) case x: return #x;
#define L       100.0
#define N       200
#define DX      (L/N)
#define DT      (0.01*DX)

#define NO_STEP 80000
#define G       9.81

#define DEVICE 0
char* clErrorStr(cl_int err);

int main(int argc, const char * argv[]) {
    
    // Load Kernel
    FILE *fp;
    char *source_str;
    size_t source_size , program_size;
    fp = fopen( "kernel.cl" , "rb" );
    if( !fp ){
        printf( "Error: Failed to open cl file!\n" );
        return 1;
    }
    fseek( fp , 0 , SEEK_END );
    program_size = ftell( fp );
    rewind( fp );
    source_size = program_size/sizeof(char) + 1;
    source_str = ( char* )malloc( source_size * sizeof( char ) );
    source_str[ program_size ] = '\0';
    fread( source_str , sizeof(char) , source_size - 1 , fp );
    fclose( fp );
    
    float *u, *s;
    int i;
    int N_cl = N + 2 ;
    float Z = DT/DX;
    
    // Allocate Memory
    size_t array_size = ( N + 2 ) * sizeof( float );
    u = ( float* )malloc( array_size );
    s = ( float* )malloc( array_size );
    
    // Initiate
    for( i = 0 ; i < N + 2 ; ++i ){
        if( ( i - 1 ) < 0.5 * N ){
            u[ i ] = 10.0;
            s[ i ] = 0.0;
        }else{
            u[ i ] = 1.0;
            s[ i ] = 0.0;
        }
    }
    
    // Get platform and device information
    cl_platform_id *platforms = NULL;
    cl_uint num_platforms;
    // Set up the platform
    cl_int clStatus;
    clStatus = clGetPlatformIDs( 0 , NULL, &num_platforms );
    platforms = ( cl_platform_id* )malloc( sizeof( cl_platform_id ) * num_platforms );
    clStatus = clGetPlatformIDs( num_platforms , platforms , NULL );
    printf( " num_platforms = %d\n", ( unsigned int )num_platforms );
    
    // Get the devices list and choose the device you want to run on
    cl_device_id *device_list = NULL;
    cl_uint num_devices;
    
    clStatus = clGetDeviceIDs( platforms[ 0 ], CL_DEVICE_TYPE_GPU , 1 , NULL , &num_devices );
    device_list = ( cl_device_id *)malloc( sizeof( cl_device_id ) * num_devices );
    clStatus = clGetDeviceIDs( platforms[ 0 ], CL_DEVICE_TYPE_GPU , num_devices , device_list , NULL);
    
    // Create one OpenCL context for each device in the platform
    cl_context context;
    context = clCreateContext( NULL , num_devices , device_list , NULL , NULL , &clStatus );
    if(clStatus != CL_SUCCESS){
	    printf("Error: Failed to create context!\n%s\n", clErrorStr(clStatus));
    }
    // Create a command queue
    cl_command_queue command_queue;
#ifdef __APPLE__
    command_queue = clCreateCommandQueue( context , device_list[ DEVICE ] , 0 , &clStatus );
#else
    command_queue = clCreateCommandQueueWithProperties( context , device_list[ DEVICE ] , 0 , &clStatus );
#endif
    if(clStatus != CL_SUCCESS){
	    printf("Error: Failed to create command queue!\n%s\n", clErrorStr(clStatus));
    }

    // Create memory buffers on the device for each array
    cl_mem U_clmem = clCreateBuffer( context , CL_MEM_READ_WRITE , 2 * array_size * sizeof( float ) , NULL , &clStatus );
    cl_mem S_clmem = clCreateBuffer( context , CL_MEM_READ_WRITE , 2 * array_size * sizeof( float ) , NULL , &clStatus );
    cl_mem F_clmem = clCreateBuffer( context , CL_MEM_READ_WRITE , 4 * array_size * sizeof( float ) , NULL , &clStatus );
    if( !U_clmem || !S_clmem || !F_clmem ){
        printf( "Error: Failed to allocate device memory!\n%s\n", clErrorStr(clStatus));
        exit(1);
    }
    
    // Copy the Buffer u and f to the device
    clStatus = clEnqueueWriteBuffer( command_queue , U_clmem , CL_TRUE , 0 , array_size , u , 0 , NULL , NULL );
    if( clStatus != CL_SUCCESS ){
        printf( "Error: Failed to write source array u!\n%s\n", clErrorStr(clStatus));
    }
    clStatus = clEnqueueWriteBuffer( command_queue , S_clmem , CL_TRUE , 0 , array_size , s , 0 , NULL , NULL );
    if( clStatus != CL_SUCCESS ){
        printf( "Error: Failed to write source array s!\n%s\n", clErrorStr(clStatus));
    }
    
    // Create a program from the kernel source
    cl_program program = clCreateProgramWithSource( context , 1 , (const char**)&source_str , NULL , &clStatus );
    if( clStatus != CL_SUCCESS ){
        printf( "Error: Failed to create program from the kernel source!\n%s\n", clErrorStr(clStatus));
        exit(1);
    }
    
    //Build the program
    clStatus = clBuildProgram( program , 1 , &device_list[ DEVICE ] , NULL , NULL , NULL );
    //DEBUG for device
    if(clStatus != CL_SUCCESS){
        size_t len;
        char buffer[2048];
        printf("Error: Failed to build program executable!\n%s\n", clErrorStr(clStatus));
        clGetProgramBuildInfo(program, device_list[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        exit(1);
    }
    cl_kernel kernel = clCreateKernel( program , "GPU_Calc" , &clStatus );
    if( !kernel || clStatus != CL_SUCCESS ){
        printf( "Error: Failed to Create Compute kernel!\n%s\n", clErrorStr(clStatus));
        exit(1);
    }
    size_t global_size = 2;
    char flag = 1;
    while(flag){
	global_size *= 2;
	if( global_size > N ){
	    flag = 0;
	}
    }
    size_t local_size = 64;
    
    for( i = 0; i < NO_STEP; ++i){
        // Input argument of the kernel
        clStatus = clSetKernelArg( kernel , 0 , sizeof(int) , (void*) &N_cl );
        clStatus = clSetKernelArg( kernel , 1 , sizeof(float) , (void*) &Z );
        clStatus = clSetKernelArg( kernel , 2 , sizeof(cl_mem) , (void*) &U_clmem );
        clStatus = clSetKernelArg( kernel , 3 , sizeof(cl_mem) , (void*) &S_clmem );
        clStatus = clSetKernelArg( kernel , 4 , sizeof(cl_mem) , (void*) &F_clmem );
        
        if(clStatus != CL_SUCCESS){
            printf("Error: Failed to set argument for kernels!\n%s\n", clErrorStr(clStatus));
            exit(1);
        }
        
        // Execute the OpenCL kernel on the list
        clStatus = clEnqueueNDRangeKernel( command_queue, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        if(clStatus != CL_SUCCESS){
            printf("Error: Failed to excute kernl!\n%s\n", clErrorStr(clStatus));
            exit(1);
        }
        
        // Update the value of u and s
        clStatus = clEnqueueCopyBuffer( command_queue , U_clmem , U_clmem , N_cl * sizeof(float) , 0 , N_cl * sizeof(float) , 0 , NULL , NULL );
        clStatus = clEnqueueCopyBuffer( command_queue , S_clmem , S_clmem , N_cl * sizeof(float) , 0 , N_cl * sizeof(float) , 0 , NULL , NULL );
        if(clStatus != CL_SUCCESS){
            printf( "Error: Failed to copy buffer to buffer!\n%s\n", clErrorStr(clStatus));
        }
    }
    
    // Read the cl memory C_clmem on the device to the host variable C
    clStatus = clEnqueueReadBuffer(command_queue, U_clmem, CL_TRUE, 0, ( N + 2 ) * sizeof(float), u, 0, NULL, NULL);
    clStatus = clEnqueueReadBuffer(command_queue, S_clmem, CL_TRUE, 0, ( N + 2 ) * sizeof(float), s, 0, NULL, NULL);
    if(clStatus != CL_SUCCESS){
        printf("Error: Faild to write sourec array  U and S !\n%s\n", clErrorStr(clStatus));
        exit(1);
    }    clStatus = clReleaseKernel( kernel );
    
    // Save Result
    fp = fopen( "data_0.txt", "w");
    for( i = 0 ; i < N + 2 ; ++i ){
        fprintf( fp, "%e %e %e\n", i*DX , u[ i ] , s[ i ] / u[ i ] );
    }
    fclose(fp);
    clStatus = clReleaseProgram( program );
    clStatus = clReleaseMemObject( U_clmem );
    clStatus = clReleaseMemObject( S_clmem );
    clStatus = clReleaseMemObject( F_clmem );
    clStatus = clReleaseCommandQueue( command_queue );
    clStatus = clReleaseContext( context );
    
    free( source_str );
    free( u );
    free( s );
    free( platforms );
    free( device_list );
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
