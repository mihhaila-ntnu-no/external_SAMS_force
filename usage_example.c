#include <stdio.h>
    
#include "csvwriter.h"    
    
int main() {    
	char *firstLine[] = {"this", "is the first line", "that ends here"};    
	char *secLine[] = {"this field, contains a comma", "this field, \" contains a quote", "this field\ncontains a newline"};    
	char *thirdLine[] = {"this line is crazy", "12345\"\"sdgdsag,adad\"\"\"\nabcdefg", "\","};    
	char **data[] = {firstLine, secLine, thirdLine};    
	
	CsvWriter *csvWriter = CsvWriter_new("Book1.csv", ",", 0);	    
	int i, j;    
	for (i = 0 ; i < 3 ; i++) {    
		for (j = 0 ; j < 3 ; j++) {    
			if (CsvWriter_writeField(csvWriter, data[i][j])) {    
				printf("Error: %s\n", CsvWriter_getErrorMessage(csvWriter));    
				return 1;    
			}    
		}    
		CsvWriter_nextRow(csvWriter);    
	}    
	CsvWriter_destroy(csvWriter);    
	    
	return 0;    
}