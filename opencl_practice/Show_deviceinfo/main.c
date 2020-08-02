#include <stdio.h>
#include <stdlib.h>
#include <CL/cl.h>

#define NUM_VALUES 1024

static int validate(cl_float *input , cl_float *output){
	int i;
	for( i = 0 ; i < NUM_VALUES ; ++i ){
		if( output[ i ] != input[ i ] * input[ i ] ){
//			fprintf(stdout,
//				"Error: Element %d did not match expected output.\n", i );
//			fprintf(stdout,
//				"Error: saw %1.4f, expected %1.4f\n", output[i],
//				input[i] * input[i]);
//			fflush(studout);
			return 0;
		}
	}
	return 1;
}

int main( int argc , const char *argv[] ){
	int i;
	char name[128];
	cl_int err;
	cl_uint numPlatforms;
	cl_platform_id *platforms;
	
	err = clGetPlatformIDs( 0 , NULL , &numPlatforms );
	if( err == CL_SUCCESS )
		printf( "\nDetected OpenCL platforms: %d\n", numPlatforms );
	else
		printf( "\nError calling clGetPlatformIDs. Error code: %d\n", err );
	platforms = ( cl_platform_id* )malloc( sizeof( cl_platform_id ) * numPlatforms);
	err = clGetPlatformIDs( numPlatforms , &platforms[0] , &numPlatforms );
	if( err != CL_SUCCESS ){
		printf( "Unable to get Platforms\n" );
	}

	cl_uint numDevices = 0;
	cl_device_id *devices;
	err = clGetDeviceIDs( platforms[ 0 ] , CL_DEVICE_TYPE_GPU , 0 , NULL , &numDevices);
	if( numDevices == 0 ){
		printf( "NO GPU device aavailable.\n");
	}else{
		printf( "Number of GPU device is %d\n", numDevices );
		devices = ( cl_device_id* )malloc( sizeof( cl_device_id ) * numDevices );
		err = clGetDeviceIDs( platforms[0] , CL_DEVICE_TYPE_GPU, numDevices , devices, NULL);

		clGetDeviceInfo(devices[0], CL_DEVICE_NAME, 128, name, NULL);
		printf("Using device: %s\n", name);
	}

	free(platforms);
	free(devices);
	return 0;
}
