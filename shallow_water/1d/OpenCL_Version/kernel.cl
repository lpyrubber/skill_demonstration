__kernel void GPU_Calc(
                       int N,
                       float Z,
                       __global float *A,
                       __global float *B,
                       __global float *C)
{
    float vel , acc, F1 , F2 , Fr , FL1 , FL2 , FR1 , FR2;
    int i = get_global_id(0);
    float g = 9.81;
    //compute flux
    if( i < N  ){
        vel = B[ i ] / A[ i ];
        acc = sqrt( g * A[ i ] );
        Fr = vel / acc;
        if ( Fr > 1 ) Fr = 1;
        if ( Fr < -1 ) Fr = -1;
        F1 = A[ i ] * vel;
        F2 = A[ i ] * vel * vel + 0.5 * g * A[ i ] * A[ i ];
        C[ i         ] =   0.5 * ( F1 * ( Fr + 1 ) + A[ i ] * acc * ( 1 - Fr * Fr ) );
        C[ i + N     ] =   0.5 * ( F2 * ( Fr + 1 ) + B[ i ] * acc * ( 1 - Fr * Fr ) );
        C[ i + N * 2 ] = - 0.5 * ( F1 * ( Fr - 1 ) + A[ i ] * acc * ( 1 - Fr * Fr ) );
        C[ i + N * 3 ] = - 0.5 * ( F2 * ( Fr - 1 ) + B[ i ] * acc * ( 1 - Fr * Fr ) );
    }
    barrier( CLK_LOCAL_MEM_FENCE );
    
    //compute U
    if( ( i > 0 ) && ( i < N - 1 ) ){
        FL1 = C[ i - 1 ] + C[ i +     N * 2 ];
        FR1 = C[ i     ] + C[ i + 1 + N * 2 ];
        FL2 = C[ i - 1 + N ] + C[ i +     N * 3 ];
        FR2 = C[ i     + N ] + C[ i + 1 + N * 3 ];
        A[ i + N ] = A[ i ] - Z * ( FR1 - FL1 );
        B[ i + N ] = B[ i ] - Z * ( FR2 - FL2 );
    }
    barrier( CLK_LOCAL_MEM_FENCE );

    //set reflective boundary
    if( i == 0 ){
        A[ i + N ] =   A[ i + 1 ];
        B[ i + N ] = - B[ i + 1 ];
    }
    if( i == N - 1 ){
        A[ i + N ] =   A[ i - 1 ];
        B[ i + N ] = - B[ i - 1 ];
    }
}
