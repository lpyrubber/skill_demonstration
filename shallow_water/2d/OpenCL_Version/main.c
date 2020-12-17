#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef __APPLE__
	#include <OpenCL/cl.h>
#else
	#include <CL/cl.h>
#endif

#define Lx      100.0
#define Ly      100.0
#define Nx      2000
#define Ny      2000
#define DX      (Lx/Nx)
#define DY      (Ly/Ny)
#define DT      (0.01*DX)
#define Zx      (DT/DX)
#define Zy      (DT/DY)
#define NO_STEP 800
#define LOCAL   256
#define GLOBAL  ((((Nx+2)*(Ny+2)+LOCAL-1)/LOCAL)*LOCAL)

#define DEVICE     0
#define PLATFORM   0
#define LOCAL_SIZE 256

#define CaseReturnString(x) case x: return #x;

char* clErrorStr(cl_int err);

int main(int argc, const char *argv[]){
	float *u , Zx_cl , Zy_cl;
	FILE *fp;
	char *source_str;
	size_t source_size , program_size , array_size , global_size , local_size;
	int i , j , Nx_cl , Ny_cl;
	
	Zx_cl = Zx;
	Zy_cl = Zy;
	Nx_cl = Nx;
	Ny_cl = Ny;
	
	// Variables for OpenCL
	
	cl_platform_id *platforms = NULL;
	cl_device_id *device_list = NULL;
	cl_context context;
	cl_command_queue command_queue;
	cl_mem U_clmem , FM_clmem , FP_clmem;
	cl_program program;
	cl_kernel kernel;
	cl_uint num_platforms , num_devices;
	cl_int clStatus;

	// Load Program

	fp = fopen( "kernel.cl" , "rb" );
	if( !fp ){
		printf( "Error: failed to open cl file!\n" );
	}
	fseek( fp , 0 , SEEK_END );
	program_size = ftell( fp );
	rewind( fp );
	source_size = program_size / sizeof( char ) + 1 ;
	source_str = ( char* )malloc( source_size * sizeof( char ) );
	source_str[ program_size ] = '\0' ;
	fread( source_str , sizeof( char ) , source_size - 1 , fp );
	fclose( fp );
	
	// Allcate Memory

	array_size = 3 * ( Nx + 2 ) * ( Ny + 2 );
	u = ( float* )malloc( 3* GLOBAL * sizeof( float ) );
	
	// Initial
	
	for( j = 0 ; j < Ny + 2 ; ++j ){
		for( i = 0 ; i < Nx + 2 ; ++i ){
			if( ( i - 1 ) < 0.5 * Nx ){
				u[ i + ( Nx + 2 ) * j                               ] = 10.0;
				u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ] = 0.0;
				u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] = 0.0;
			}else{
				u[ i + ( Nx + 2 ) * j                               ] = 1.0;
				u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ] = 0.0;
				u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] = 0.0;
			}
		}
	}

	// Get Platform and Device Information
	// * Set up the plateform
	
	clStatus = clGetPlatformIDs( 0 , NULL , &num_platforms );
	platforms = ( cl_platform_id* )malloc( num_platforms * sizeof( cl_platform_id ) );
	clStatus = clGetPlatformIDs( num_platforms , platforms , NULL);
	printf( "num_platforms = %d\n " , ( unsigned int )num_platforms );
	if(clStatus != CL_SUCCESS) printf("Error: Failed to get platform!\n%s\n", clErrorStr(clStatus));

	// * Get the devices list and choose the device you wnat to run on
	
	clStatus = clGetDeviceIDs( platforms[ PLATFORM ] , CL_DEVICE_TYPE_GPU , 1 , NULL , &num_devices );
	device_list = ( cl_device_id* )malloc( num_devices * sizeof( cl_device_id ) );
	clStatus = clGetDeviceIDs( platforms[ PLATFORM ] , CL_DEVICE_TYPE_GPU , num_devices , device_list , NULL );
	if(clStatus != CL_SUCCESS) printf("Error: Failed to get devices!\n%s\n", clErrorStr(clStatus));

	// Create one OpenCL contex for each device in the choosed platform and a command queue for the specific device
	// * Set up the context
	
	context = clCreateContext( NULL , num_devices , device_list , NULL , NULL , &clStatus );
	if(clStatus != CL_SUCCESS) printf("Error: Failed to create context!\n%s\n", clErrorStr(clStatus));
	// * Set up the command queue
	
#ifdef __APPLE__
	command_queue = clCreateCommandQueue( context , device_list[ DEVICE ] , 0 , &clStatus );
#else
	command_queue = clCreateCommandQueueWithProperties( context , device_list[ DEVICE ] , 0 , &clStatus );
#endif
	if(clStatus != CL_SUCCESS) printf("Error: failed to create queue!\n%s\n", clErrorStr(clStatus));
	// Create memory buffers on the device and copy data to the device
	// * Set up memory bUffers on the device for each arrary
		
	array_size = 3 * GLOBAL * sizeof(float);
	U_clmem = clCreateBuffer( context , CL_MEM_READ_WRITE , array_size , NULL , &clStatus );
	FM_clmem = clCreateBuffer( context , CL_MEM_READ_WRITE , 2 * array_size , NULL , &clStatus );
	FP_clmem = clCreateBuffer( context , CL_MEM_READ_WRITE , 2 * array_size , NULL , &clStatus );
	if( !U_clmem || !FM_clmem || !FP_clmem ) printf( "Error: Failed to allocate device memory!\n%s\n",clErrorStr(clStatus) );
	// * Copy to the device
	
	array_size = 3 * ( Nx + 2 ) * ( Ny + 2 ) * sizeof(float);
	clStatus = clEnqueueWriteBuffer( command_queue , U_clmem , CL_TRUE , 0 , array_size , u , 0 , NULL , NULL );
	if( clStatus != CL_SUCCESS ) printf( "Error: Failed to sent host u to device!\n%s\n", clErrorStr(clStatus));
	
	// Create a program from the kernel source and Debugging
	// * Set up the program
	
	program = clCreateProgramWithSource( context , 1 , ( const char**) &source_str , NULL , &clStatus );
	if( clStatus != CL_SUCCESS ){
		printf( "Error: Failed to create program form the kernel source!\n%s\n", clErrorStr(clStatus) );
		exit( 1 );
	}

	// * Build the program

	clStatus = clBuildProgram( program , 1 , &device_list[ DEVICE ] , NULL , NULL , NULL );

	// * Get debug information

	if( clStatus != CL_SUCCESS ){
		size_t len;
		char buffer[ 2048 ];
		printf( "Error: Failed to build program executable!\n" );
		clGetProgramBuildInfo( program, device_list[ DEVICE ] , CL_PROGRAM_BUILD_LOG , sizeof(buffer) , buffer , &len );
		printf( "%s\n" , buffer );
		exit( 1 );
	}
	
	// * Set up Kernel
	
	kernel = clCreateKernel( program , "GPU_Calc" , &clStatus );
	if( !kernel || clStatus != CL_SUCCESS ){
		printf( "Error: Failed to create Compute kernel!\n%s\n", clErrorStr(clStatus) );
		exit( 1 );
	}

	// Set up local & global size for the device
	
	local_size = LOCAL_SIZE;
	global_size = ( ( ( Nx + 2 ) * ( Ny + 2 ) + local_size - 1 ) / local_size ) * local_size;
	printf("global_size=%d\n",(int)global_size );

	for( i = 0 ; i < NO_STEP ; ++i ){

		// Input the argument of the kernel

		clStatus  = clSetKernelArg( kernel , 0 , sizeof( int ) , ( void* ) &Nx_cl );
		clStatus |= clSetKernelArg( kernel , 1 , sizeof( int ) , ( void* ) &Ny_cl );
		clStatus |= clSetKernelArg( kernel , 2 , sizeof( float ) , ( void* ) &Zx_cl );
		clStatus |= clSetKernelArg( kernel , 3 , sizeof( float ) , ( void* ) &Zy_cl );
		clStatus |= clSetKernelArg( kernel , 4 , sizeof( cl_mem ) , ( void* ) &U_clmem );
		clStatus |= clSetKernelArg( kernel , 5 , sizeof( cl_mem ) , ( void* ) &FM_clmem );
		clStatus |= clSetKernelArg( kernel , 6 , sizeof( cl_mem ) , ( void* ) &FP_clmem );

		if(clStatus != CL_SUCCESS ){
			printf( "Error: Failed to set argument for kernels!\n%s\n", clErrorStr(clStatus) );
			exit(1);
		}

		// Execute the OpenCL kernel on the list
	
		clStatus = clEnqueueNDRangeKernel( command_queue , kernel , 1 , NULL , &global_size , &local_size , 0 , NULL , NULL );
		if( clStatus != CL_SUCCESS ){
			printf( "Error: Failed to excute kernel!\n%s\n" , clErrorStr(clStatus));
			exit( 1 );
		}
	}
		
	// Read the cl memory on the device to the host variable
	
	array_size = 3 * ( Nx + 2 ) * ( Ny + 2 ) * sizeof( float );
	clStatus = clEnqueueReadBuffer( command_queue , U_clmem , CL_TRUE , 0 , array_size , u , 0 , NULL , NULL );
	if( clStatus != CL_SUCCESS ){
		printf( "Error: Failded to send device array u to host!\n%s\n" , clErrorStr(clStatus));
		exit( 1 );
	}

	fp = fopen( "data.txt" , "w" );
	for( j = 1 ; j < Ny + 1 ; ++j ){
		for ( i = 1 ; i < Nx + 1 ; ++i ){
			fprintf( fp , "%e %e %e %e %e\n"\
				, i * DX\
				, j * DY\
				, u[ i + ( Nx + 2 ) * j                               ]\
				, u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ]\
				, u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]);
		}
	}

	fclose(fp);

	clStatus = clReleaseKernel( kernel );
	clStatus = clReleaseProgram( program );
	clStatus = clReleaseMemObject( U_clmem );
	clStatus = clReleaseMemObject( FM_clmem );
	clStatus = clReleaseMemObject( FP_clmem );
	clStatus = clReleaseCommandQueue( command_queue );
	clStatus = clReleaseContext( context );
    
	free( source_str );
	free( u );
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
