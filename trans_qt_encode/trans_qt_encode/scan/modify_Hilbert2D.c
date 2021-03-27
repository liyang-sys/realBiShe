/*
*�޸����ڣ�2021-1-13
*bug��1,calculation ��������ֿ�ʱ���Լ����������ֳ����ӿ�����������
*     2,judge_order ���¼���ֿ�ʱ�����·ֿ�����Լ������������
*
*�޸����ڣ�2021-1-14
*���ƣ�·�����к��е������ʾ������ת������
*
*�޸����ڣ�2021-1-15---2021-1-18
*�߼��ϵ��޸ģ����»����ӿ飬����дlfind_2s_inc����������find_2s_inc_2
*              ʹ����ɨ������MATLAB��һ��
*              ��draw_hibert��calculation������Ӧ�Ķ�
*
*----��ϼ
*/


#include "modify_Hilbert2D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*���ڶ�̬�������ⳤ�ȶ�ά���鲢����*/
int **create_arr(int width, int height)
{
	int **new_arr = (int **)malloc(width*sizeof(int *));										//������
	for (int i = 0; i<width; i++)
		new_arr[i] = (int *)malloc(height*sizeof(int));											//������
	return new_arr;
}

//���ݻ���square�ĳ���width,�����Ҫ�ӳɵĸ߶�height������hҪ�ӵ����У�����һ������Ķ�ά����
int **bc_operation(int width, int height, block h)
{
	int i, j;
	int len = width;
	int times = h.num;
	int **block_one = create_arr(len, height);
	int **square = create_arr(len, len);
	int idh = width;
	int w = len; 
	for (i = 0; i < len; i++){
		for (j = 0; j < len; j++)
			square[i][j] = xy_point2d(len, i, j);
	}

	for (i = 0; i < times; i++){
		if (h.block[i] > 0){
			square = add(square, len, idh, h.block[i]);
			idh += h.block[i];
		}
		else if (h.block[i] < 0){
			square = sub(square, len, idh, abs(h.block[i]));
			idh += h.block[i];
		}
	}
	for (i = 0; i < height; i++){
		for (j = 0; j < len; j++)
			block_one[j][i] = square[j][i];
	}
	return block_one;
}


total calculation(int w, int h)
{
	total basic;
   int MaxDel = 0, ia = 0;
   block th_block;

	basic.basic_block = find_basic_block(w);													//��������ڿ�ȵ�����
	basic.h_order = find_basic_block(h);														//��������ڸ߶ȵ�����
   basic = judge_order_2(basic, w, h);
	//basic = judge_order(basic, w, h);															//��ϸ߶ȺͿ�ȵ������жϲ��޸���������
	//basic.h_order = find_2s_inc(basic.h_order, h - basic.basic_block.block[0]);					//����߶ȵ������Ƿ��ܹ�����
      
   th_block = find_2s_inc_2(basic.basic_block.block[0], h - basic.basic_block.block[0], basic.h_order);
   
   while (ia<th_block.num&&th_block.block[ia]>0){
      MaxDel += th_block.block[ia];
      ia++;
   }
   if (basic.basic_block.block[0] == MaxDel){
      return basic;
   }
   else{
      int i = 0;
      while (th_block.block[i] != 0){
         basic.h_order.block[i] = th_block.block[i];
         i++;
      }
      while (basic.h_order.block[i] != 0){
         basic.h_order.block[i] = 0;
         i++;
      }
      basic.h_order.num = th_block.num;
   }

	return basic;
}

int **draw_hibert(int W, int H)
{
	int flip = 0;;
	int i, j, k;
	total chunk;
   block h_new, tem_h1, tem_w1;
	int w , h;																					//���ڱ�ʾÿ��Ҫ�����w��h
	int W1, H1, W2, H2;
	int w1, h1, w2, h2;																			//�߶ȺͿ�ȵ�תʱʹ��
	int flag = 0;																				//���ڱ�־�Ƿ������ֻ���һ��
   int ptx = 0, pty = 0, ptw = 0, ptz = 0;																//������Ϊ���Ҫ��Ӹ������λ��ָ��
	int offset = 0;																				//��ʾ��ʱ����ӵ���
   int rotation = 0, flg1 = 1, flg2 = 1, f1 = 0, f2 = 0;			//�����ж������Ƿ���Ҫ��ת��0��ʾ��ת0�ȣ�1��ʾ˳ʱ����ת270��,flgһ���ж��Ƿ�ı䷽��,f������������ӵĴ���
	int **img;																					//����������ܵ�ͼ������
	int **part, **part_rot, **part_culumn;														//ÿ�����������������
	int *part_row;			
	int **img_flip;

	if (W < H){																					//�����ʼ���С�ڸ߶ȣ��ѿ�͸ߵ�������
		flip = 1;
		w = H;
		h = W;
		W2 = H;
		H2 = W;
	}
	else{
		flip = 0;
		w = W;
		h = H;
		W1 = W;
		H1 = H;
	}
	img = create_arr(w, h);
	w1 = w;
	h1 = h;
	while (flag == 0){
		if ((rotation == 0||flg1==0)&&flg2!=0){																		//���ж��Ƿ���Ҫ��ת
         f2 = 0;
         if (f1 > 0)pty -= h1;
         if (w1 < h1&&flg1==1){
				rotation = 1;
				h2 = w1;
				w2 = h1;
			}
			else if (w1 >= h1 && h1 == 1){														//�ж��Ƿ�������߶�ֻ��һ�е����
				part_row = row(w1);
				for (i = 0; i < w1; i++)
					img[ptx + i][pty] = part_row[i] + offset;
				ptx += w1;
				pty += 1;
            ptz = pty;
			}
			else if ( h1 > 1){										//w1 >= h1 &&	//�߶ȴ���1
				chunk = calculation(w1, h1);	
				for (i = 0; i < chunk.basic_block.num; i++){
					part = bc_operation(chunk.basic_block.block[i], h1, chunk.h_order);
					for (j = 0; j < chunk.basic_block.block[i]; j++){
						for (k = 0; k < h1; k++)
							img[ptx + j][pty + k] = part[j][k] + offset;
					}
               ptz = pty + h1;

					for (int ii = 0; ii < chunk.basic_block.block[i]; ii++)					
						free(part[ii]);														
					free(part);																	
					offset += chunk.basic_block.block[i] * h1;
					ptx += chunk.basic_block.block[i];               
					if (chunk.basic_block.block[i + 1] != chunk.basic_block.block[i]){
						i += 1;
						break;
					}
				}
            flg1 = 1;
				rotation = 1;																	//����˳ʱ��ı�270��
				w -= i*chunk.basic_block.block[0];												//���ڵĿ�ȱ仯
				w2 = h1;																		//��һ�����̵Ŀ��
				h2 = w;																			//��һ�����̵ĸ߶�
            tem_h1 = find_basic_block(h1);
            if ((chunk.basic_block.block[0] << 1) > h1&&tem_h1.block[0] <= w){
               flg1 = 0;
               w1 = w;
               f1++;
               pty += h1;
               ptz = pty;
            }

			}
			if (flip == 0){
				if (ptx == W && ptz == H)														//�����ж��Ƿ�������
					flag = 1;
				else
					flag = 0;
			}
			else{
				if (ptx == H && ptz == W)														//�����ж��Ƿ�������
					flag = 1;
				else
					flag = 0;
			}
		}
		else {    //if (rotation == 1)
         if (f2 > 0)ptx -= w1;;
         f1 = 0;
			if (h2 == 1){
				part_culumn = create_arr(h2, w2);
				for (i = 0; i < w2; i++)
					part_culumn[0][i] = i;
				for (i = 0; i < w2; i++)
					img[ptx][pty + i] = part_culumn[0][i] + offset;
				pty += w2;
				ptx += 1;
				ptw = ptx;
			}
			else if (h2 > 1){
				h_new = find_basic_block(w2);
				w2 = h_new.block[0];
				chunk = calculation(w2, h2);
				for (i = 0; i < chunk.basic_block.num; i++){
					part = bc_operation(chunk.basic_block.block[i], h2, chunk.h_order);
					part_rot = rotate_arr(chunk.basic_block.block[i], h2, part);
					for (j = 0; j < chunk.basic_block.block[i]; j++){
						for (k = 0; k < h2; k++)
							img[ptx + k][pty + j] = part_rot[k][j] + offset;
					}
					ptw = ptx + h2;
					//�ͷ�������ά����
					for (int ii = 0; ii < chunk.basic_block.block[i]; ii++)						//�ͷ��ڴ�
						free(part[ii]);															//�ͷ���
					free(part);																	//�ͷ���
					for (int ii = 0; ii < h2; ii++)												//�ͷ��ڴ�
						free(part_rot[ii]);														//�ͷ���
					free(part_rot);																//�ͷ���

					offset += chunk.basic_block.block[i] * h2;
					pty += chunk.basic_block.block[i];
					if (chunk.basic_block.block[i + 1] != chunk.basic_block.block[i])
						break;
				}
				rotation = 0;																	//����˳ʱ��ı�90��
            flg2 = 1;
            h -= w2;
				w1 = h2;																		//��һ�����
				h1 = h;																			//��һ���߶�
            tem_h1 = find_basic_block(h1);
            tem_w1 = find_basic_block(w1);
            if ((tem_h1.block[0] << 1) > h2&&tem_w1.block[0] <= h1){
               flg2 = 0;
               w2 = h;
               f2++;
               ptx += w1;
               ptw = ptx;
            }

			}
			else if (h2 == 0){
				flag = 1;
				break;
			}
			if (flip == 0){
				if (ptw == W && pty == H)														//�����ж��Ƿ�������
					flag = 1;
				else
					flag = 0;
			}
			else if (flip == 1){
				if (ptw == H && pty == W)														//�����ж��Ƿ�������
					flag = 1;
				else
					flag = 0;
			}
		}      

	}
	if (flip == 1){
		img_flip = rotate_arr(H, W, img);
		return img_flip;
	}
	else if (flip == 0)
		return img;
}

//����һά���һ������
int *row(int n)
{
	int *arr = (int *)malloc(n*sizeof(int *));
	for (int i = 0; i < n; i++)
		arr[i] = i;
	return arr;
}


//////������ת��ά���飬wΪԭ�����hΪԭ�����
int **rotate_arr(int w, int h, int **old)
{
	int i, j;
	int **new, **mirror;
	new = create_arr(h, w);
	mirror = HorizontalMirror(old, w, h);				//����Գƺ�����ת
	////////270,old[w][h],new[h][w]
	for (i = 0; i<h; i++)
	{
		for (j = 0; j<w; j++)
			new[i][j] = mirror[w - 1 - j][i];
	}
	return new;
}

//��ά����ˮƽ����
int **HorizontalMirror(int **arr, int n, int m)
{
	int r, c;
	for (r = 0; r < m; r++){
		for (c = 0; c < n / 2; c++)						//ÿһ��ת����n/2����
		{												//�����������ԳƵ�Ԫ�ؽ��н���	
			int temp = arr[c][r];
			arr[c][r] = arr[n-c-1][r];
			arr[n-c-1][r] = temp;
		}
	}
	return arr;
}

int *modify_Hilbert2D(int W, int H)
{
	int i, j, k, lg;
	int **m0, **img, *m;
	lg = H*W;

	img = draw_hibert(W, H);

   m0 = create_arr(3, lg);     //��·�����к��е������ʾ;m0[0]:��,m0[1]:��,m0[3]:���,m0[3]�ڶ�ά��ȫΪ0��
   for (i = 0; i < lg; i++)
      m0[2][i] = 0;

   int x = 0, y = 0;
   m0[0][0] = 0, m0[1][0] = 0;
   for (int i0 = 1; i0 < lg; i0++){
      if (x>0 && img[x - 1][y] == img[x][y] + 1){
         x = x - 1;
         m0[0][i0] = y;
         m0[1][i0] = x;
      }
      else if (x<W - 1 && img[x + 1][y] == img[x][y] + 1){
         x = x + 1;
         m0[0][i0] = y;
         m0[1][i0] = x;
      }
      else if (y>0 && img[x][y - 1] == img[x][y] + 1){
         y = y - 1;
         m0[0][i0] = y;
         m0[1][i0] = x;
      }
      else if (y<H - 1 && img[x][y + 1] == img[x][y] + 1){
         y = y + 1;
         m0[0][i0] = y;
         m0[1][i0] = x;
      }

   }
   for (i = 0; i < W; i++)
      free(img[i]);
   free(img);

   /*************************/
	//for (int i = 5524; i < 5534; i++)
	//	printf("%d ", m0[0][i]);
	//printf("\n");
	//for (int i = 5524; i < 5534; i++)
	//	printf("%d ", m0[1][i]);

	//FILE *fp2;
	//fopen_s(&fp2, "m0.dat", "wb");
	//for (i = 0; i<2; i++)
	//{
	//	fwrite(&m0[i][0], sizeof(int), lg, fp2);
	//}
	//fclose(fp2);
	/*************/

   m = (int*)calloc(lg, sizeof(int));
   if (!m){
      printf("��������mʧ�ܣ�\n");
      exit(1);
   }
	for (i = 0; i < lg - 1; i++){   /*����ά��m��Ϊһά��m1*/
		if (m0[0][i + 1] - m0[0][i] == 1)/*��*/
			m[i] = 2;
		else if (m0[0][i + 1] - m0[0][i] == -1)
			m[i] = -2;
		else if (m0[1][i + 1] - m0[1][i] == 1)/*��*/
			m[i] = 3;
		else if (m0[1][i + 1] - m0[1][i] == -1)
			m[i] = -3;
	}
	free(m0[0]); free(m0[1]); free(m0[2]); free(m0);

	return m;
}

