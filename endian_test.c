#include <stdio.h>
#include <stdint.h>

int main( int argc, char **argv )
{
	uint16_t word = 0xff00;
	char *byte1 = (char *) &word;
	char *byte2 = byte1 + 1;
	printf( "byte1=%x\nbyte2=%x\n", *byte1 & 0xff, *byte2 & 0xff );
	if( *byte1 ) printf( "This system is big-endian.\n" );
	else printf( "This system is little-endian.\n" );
	

	union
	{
		float f[2];
		double d;
	} f2d;
	f2d.f[0]=1.;
	f2d.f[1]=1.;
	printf("Float contents: %9.1E %9.1E. Double: %9.1E\n",f2d.f[0],f2d.f[1],f2d.d);

	char* words[]={"first","second","third"};
	for (int i=0;i<3;i++)
		printf("%s\n",words[i]);
	// The size print works only because the array was declared right here
	printf("Number of elements in words: %d\n",(int)sizeof(words)/sizeof(words[1]));

}

