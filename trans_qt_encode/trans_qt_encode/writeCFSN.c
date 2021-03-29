#define _CRT_SECURE_NO_WARNINGS
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"
void writeCFSN(float*cf0,int lg)
{
	int*cf0Int32 = calloc(lg, sizeof(int));
	//printf("64位int字节大小 = %d", sizeof(int));
	for (int i = 0; i < lg; i++)
	{
		cf0Int32[i] = cf0[i];
	}

	FILE* fp = fopen("CF0SN.txt", "wb");
	fwrite(cf0Int32, sizeof(int), lg, fp);
	fclose(fp);

}