#include <stdio.h>
#include <windows.h>   
    
int main(void)
{    
	printf("Hello FileAPI\n");

	HANDLE hFile, hFile1;
	LPCWSTR fname1 = L"C:\\Users\\mihha\\source\\repos\\external_SAMS_force\\filesize_test.txt";
	DWORD dwFileSize;
	DWORD dwFileType;

	hFile1=CreateFile(fname1,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,NULL);
	if (hFile1==INVALID_HANDLE_VALUE)
	{
		printf("Could not open file %S, error %d\n",fname1,GetLastError());
		return 1;
	}
	dwFileType=GetFileType(hFile1);
	dwFileSize=GetFileSize(hFile1,NULL);
	printf("%S size is %d bytes and file type is %d\n",fname1, dwFileSize, dwFileType);
	CloseHandle(hFile1);
	return 0;
}