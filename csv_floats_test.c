#include <stdio.h>
#include "csvwriter.h"    
    
int main()
{    
	char *header_values[] = {"first", "second", "third"};
	int row_ints[3] = {1,2,3};
	float row_floats[3] = {1./3.,2./7.,3./1100000000000000000000000000000000000000.};
	char row_ints_printed[3][17]; // 16 recommended + 1 terminating
	char row_floats_printed[3][25]; // 24 recommended + 1 terminating

	printf("Length: %lld\n", sizeof(row_ints) / sizeof(row_ints[0]));

	CsvWriter *csvWriter = CsvWriter_new("test.csv", ";", 0);
	for (int i = 0 ; i < 3 ; i++)
	{    
		CsvWriter_writeField(csvWriter, header_values[i]); 
	}
	CsvWriter_nextRow(csvWriter);

	char int_cell_buf[16];
	char float_cell_buf[24];
	for (int i = 0 ; i < 3 ; i++)
	{
		sprintf(row_ints_printed[i],"%d",row_ints[i]);
		sprintf(row_floats_printed[i],"%E",row_floats[i]);
	}

	for (int i = 0 ; i < 3 ; i++)
	{
		CsvWriter_writeField(csvWriter, row_ints_printed[i]); 
	}
	CsvWriter_nextRow(csvWriter);

	for (int i = 0 ; i < 3 ; i++)
	{
		CsvWriter_writeField(csvWriter, row_floats_printed[i]);
	}
	CsvWriter_destroy(csvWriter);
	return 0;
}