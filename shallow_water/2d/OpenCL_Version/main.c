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
#define Nx      200
#define Ny      200
#define DX      (Lx/Nx)
#define DY      (Ly/Ny)
#define DT      (0.01*DX)
#define Zx      (DT/DX)
#define Zy      (DT/DY)
#define NO_STEP 800

#define DEVICE     0
#define PLATFORM   0
#define LOCAL_SIZE 64

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
	printf("%d\n ", (int)array_size);
	u = ( float* )malloc( array_size * sizeof( float ) );
	
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

	fp = fopen( "data_0.txt" , "w" );
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


	// Get Platform and Device Information
	// * Set up the plateform
	
	clStatus = clGetPlatformIDs( 0 , NULL , &num_platforms );
	platforms = ( cl_platform_id* )malloc( num_platforms * sizeof( cl_platform_id ) );
	clStatus = clGetPlatformIDs( num_platforms , platforms , NULL);
	printf( "num_platforms = %d\n " , ( unsigned int )num_platforms );

	// * Get the devices list and choose the device you wnat to run on
	
	clStatus = clGetDeviceIDs( platforms[ PLATFORM ] , CL_DEVICE_TYPE_GPU , 1 , NULL , &num_devices );
	device_list = ( cl_device_id* )malloc( num_devices * sizeof( cl_device_id ) );
	clStatus = clGetDeviceIDs( platforms[ PLATFORM ] , CL_DEVICE_TYPE_GPU , num_devices , device_list , NULL );

	// Create one OpenCL contex for each device in the choosed platform and a command queue for the specific device
	// * Set up the context
	
	context = clCreateContext( NULL , num_devices , device_list , NULL , NULL , &clStatus );

	// * Set up the command queue
	
#ifdef __APPLE__
	command_queue = clCreateCommandQueue( context , device_list[ DEVICE ] , 0 , &clStatus );
#else
	command_queue = clCreateCommandQueueWithProperties( context , device_list[ DEVICE ] , 0 , &clStatus );
#endif

	// Create memory buffers on the device and copy data to the device
	// * Set up memory bUffers on the device for each arrary
	
	U_clmem = clCreateBuffer( context , CL_MEM_READ_WRITE , array_size * sizeof( float ) , NULL , &clStatus );
	FM_clmem = clCreateBuffer( context , CL_MEM_READ_WRITE , 2 * array_size * sizeof( float ) , NULL , &clStatus );
	FP_clmem = clCreateBuffer( context , CL_MEM_READ_WRITE , 2 * array_size * sizeof( float ) , NULL , &clStatus );
	if( !U_clmem || !FM_clmem || !FP_clmem ){
		printf( "Error: Failed to allocate device memory!\n" );
	}
	// * Copy to the device
	
	clStatus = clEnqueueWriteBuffer( command_queue , U_clmem , CL_TRUE , 0 , array_size * sizeof(float)  , u , 0 , NULL , NULL );
	if( clStatus != CL_SUCCESS ){
		printf( "Error: Flaid to write source array u!\n" );
	}
	
	// Create a program from the kernel source and Debugging
	// * Set up the program
	
	program = clCreateProgramWithSource ( context , 1 , ( const char**) &source_str , NULL , &clStatus );
	if( clStatus != CL_SUCCESS ){
		printf( "Error: Failed to create program form the kernel source!\n" );
		exit( 1 );
	}

	// * Build the program

	clStatus = clBuildProgram( program , 1 , device_list , NULL , NULL , NULL );

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
		printf( "Error: failed to create Compute kernel!\n" );
		exit( 1 );
	}

	// Set up local & global size for the device
	
	global_size = 2;
	while( global_size < ( Nx + 2) * ( Ny + 2 ) ){
		global_size *= 2;
	}
	local_size = LOCAL_SIZE;

	for( i = 0 ; i < NO_STEP ; ++i ){

		// Input the argument of the kernel

		clStatus  = clSetKernelArg( kernel , 0 , sizeof( int ) , ( void* ) &Nx_cl );
		clStatus &= clSetKernelArg( kernel , 1 , sizeof( int ) , ( void* ) &Ny_cl );
		clStatus &= clSetKernelArg( kernel , 2 , sizeof( float ) , ( void* ) &Zx_cl );
		clStatus &= clSetKernelArg( kernel , 3 , sizeof( float ) , ( void* ) &Zy_cl );
		clStatus &= clSetKernelArg( kernel , 4 , sizeof( cl_mem ) , ( void* ) &U_clmem );
		clStatus &= clSetKernelArg( kernel , 5 , sizeof( cl_mem ) , ( void* ) &FM_clmem );
		clStatus &= clSetKernelArg( kernel , 6 , sizeof( cl_mem ) , ( void* ) &FP_clmem );

		if(clStatus != CL_SUCCESS ){
			printf( "Error: Failed to set argument for kernels!\n" );
			exit(1);
		}

		// Execute the OpenCL kernel on the list
	
		clStatus = clEnqueueNDRangeKernel( command_queue , kernel , 1 , NULL , &global_size , &local_size , 0 , NULL , NULL );
		if( clStatus != CL_SUCCESS ){
			printf( "Error: Failed to excute kernel!\n" );
			exit( 1 );
		}
	}
		
	// Read the cl memory on the device to the host variable
	
	clStatus = clEnqueueReadBuffer( command_queue , U_clmem , CL_TRUE , 0 , array_size * sizeof( float ) , u , 0 , NULL , NULL );
	if( clStatus != CL_SUCCESS ){
		printf( "Error: Failded to write source array u!\n" );
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
