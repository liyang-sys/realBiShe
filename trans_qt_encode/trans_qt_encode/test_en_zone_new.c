#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

int ptr;
uchar *bin;
void test_en_zone_new()
{
	ptr = 0;
	bin = calloc(4096000, sizeof(uchar));
	Uint32_Dat nc,r;
	nc.len = 126176;
	r.len = 723;
	FILE *fp = fopen("Nc.txt", "rb");
	nc.dat = (unsigned int *)calloc(nc.len, sizeof(unsigned int));
	fread(nc.dat, sizeof(unsigned int), nc.len, fp);
	fclose(fp);

	fp = fopen("newR.txt", "rb");
	r.dat = (unsigned int *)calloc(r.len, sizeof(unsigned int));
	fread(r.dat, sizeof(unsigned int), r.len, fp);
	fclose(fp);
	en_zone_new(&r, &nc);
	writeBinToFile(bin, ptr);
	int x = 0;
}