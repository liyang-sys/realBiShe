#ifndef _QUANT_INV_H
#define _QUANT_INV_H

#define uchar unsigned char
#define uint unsigned int

float* rst_sub(int* qcf, float delta, int lg, uchar qnt, int qctr, uchar trim, uint ctr1);
float* rstTHDctr(int *qcf, char *sgn, float T, float delta, int* nc, int maxcf, int lg);
float* rstEVENctr(int* qcf, char *sgn, float delta, int* nc, int maxcf, float ctr1, int lg);
void scan_inv(int a, int b, int c, int lg, int *m, int ptv, int harr, int dc, float *qcf);
void scan_inv_L2(int H, int W, int a, int b, int c, int lg, int *m, int harr, float *qcf);
int find_qctr(uchar *bin, uchar *qnt, uchar *trim, uint *ctr1);

#endif