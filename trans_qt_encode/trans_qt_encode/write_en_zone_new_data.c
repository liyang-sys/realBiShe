#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

void write_en_zone_new_data(Uint32_Dat*r, Uint32_Dat*nc)
{
	FILE*fp = fopen("newR.txt", "wb");
	fwrite(r->dat, sizeof(unsigned int), r->len, fp);
	fclose(fp);

	fp = fopen("Nc.txt", "wb");
	fwrite(nc->dat, sizeof(unsigned int), nc->len , fp);
	fclose(fp);
}