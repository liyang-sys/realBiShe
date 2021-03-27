#ifndef _MODIFY_HIBERT2D_H_
#define _MODIFY_HIBERT2D_H_
#include<assert.h>

/*find_block.h*/
typedef struct
{
   int block[200];
   int order[200];
   int num;
}block;

typedef struct
{
   block basic_block;
   block h_order;
}total;


int Sgn(int value);										//�ж�input������
int power_of_2(int num);								//2��������
block find_basic_block(int length);						//��2���ݴη��ָ�ɴӴ�С������
total judge_order(total basic, int w, int h);			//�����жϳ��ȺͿ�ȵļ���˳��
block find_2s_inc(block bb, int dH);					//��������Ҫ��������鳤��
total judge_order_2(total basic, int w, int h);    //
block find_2s_inc_2(int H0, int dH, block inc0);



/*hilbert_wiki.h*/
int xy_point2d(int n, int x, int y);						//������ض���Сͼ���и�xy�����ж�Ӧ��hilbert��
void rot(int n, int *x, int *y, int rx, int ry);
void d2xy(int n, int d, int *x, int *y);					//���hilbert�����ض���С��ͼ���ж�Ӧ��xy����




/*draw_hilbert.h*/
int **create_arr(int width, int height);
int **draw_hibert(int W, int H);
int **HorizontalMirror(int **arr, int n, int m);
int *row(int n);
int **rotate_arr(int w, int h, int **old);




/*add.h*/
int **add(int **orginal_arr, int width, int height, int num);    //add
int **sub(int **orginal_arr, int width, int height, int num);    //sub

int *modify_Hilbert2D(int W, int H);


#endif
