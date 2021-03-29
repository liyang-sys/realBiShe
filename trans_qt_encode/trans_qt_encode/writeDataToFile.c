#define _CRT_SECURE_NO_WARNINGS
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"
void writeDataToFile(uint *r, int lenr)
{
	FILE* fp = fopen("R.txt", "wb");
	fwrite(r, sizeof(unsigned int), lenr, fp);
	fclose(fp);

}