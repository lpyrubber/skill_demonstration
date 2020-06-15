#include <stdio.h>

void Hello( int a\
	    , int b\
	  );
int main(){
	Hello\
	( 10 \
	  , 55);
	return 0;
}

void Hello( int a, int b\
		){
	printf("%d %d\n" , a ,b);
}
