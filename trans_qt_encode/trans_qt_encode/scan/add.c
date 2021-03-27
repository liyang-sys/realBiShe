#include "modify_Hilbert2D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//�ú�������ʵ�����¼ӵĹ��ܣ��������������Լ����ĳ��Ϳ��Լ�Ҫ��ӵĳ��ȣ���������ϵ�����
int **add(int **arr, int width, int height, int num)
{
	int i, j, k;
	int counter = height*width;
	int flag = 0;
	int offset;
	int portion;
	int idy, idx,pt,pt_new;
	int unit = num * 2;
	int new_unit = num * 2 + num;
	int orginal = height*width;
	int sum = width*(height+num);
	int pointer, pointer0;
	int subnum = width / (num*2);

	int *st = (int *)malloc((subnum+1)*sizeof(int *));
	int *orginal_row = (int *)malloc(orginal*sizeof(int *));								//���ڴ洢��������ĺ�����
	int *orginal_culumn = (int *)malloc(orginal*sizeof(int *));								//���ڴ洢���������������
	int *sum_row = (int *)malloc(sum*sizeof(int *));										//���ڴ洢��ɼӷ��������ĺ�����
	int *sum_culumn = (int *)malloc(sum*sizeof(int *));										//���ڴ洢��ɼӷ��������ĺ�����
	int *transit_row = (int *)malloc(sum*sizeof(int *));									//��������,��sum�������
	int *transit_culumn = (int *)malloc(sum*sizeof(int *));									//��������,��sum�������
		
	int *row = (int *)malloc((unit*unit)*sizeof(int *));									//���ڴ洢�����ģ��ĺ�����
	int *culumn = (int *)malloc((unit*unit)*sizeof(int *));									//���ڴ洢�����ģ���������
	int *new_row = (int *)malloc((new_unit*unit)*sizeof(int *));							//���ڴ洢��Ӻ�ģ��ĺ�����
	int *new_culumn = (int *)malloc((new_unit*unit)*sizeof(int *));							//���ڴ洢��Ӻ�ģ���������

	int **orginal_arr = (int **)malloc(width*sizeof(int *));								//���ڱ�ʾ��������
	int **new_arr = (int **)malloc(width*sizeof(int *));									//���ڱ�ʾ�������
	for (int i = 0; i<width; i++)
		new_arr[i] = (int *)malloc((height+num)*sizeof(int));								
	for (int i = 0; i<width; i++)
		orginal_arr[i] = (int *)malloc(height*sizeof(int));									
///////////////////////////////////////////////////////////////////////////////////////////////��ʼ�������
	//*****��������Ӹ߶���ԭ�߶���ȣ��������һ�ּӷ�
	if (width == height && num == height){													
		idx = 0;
		idy = 0;
		offset = num*num / 4;
		portion = num*num / 4;
		for (j = 0; j < height/2; j++){														//��һ�׶Σ���ԭ��һ����һ��
			for (i = 0; i < width/2; i++)
				new_arr[i][j] = arr[i][j];
		}
		idy = height / 2;
		for (j = height/2; j < height; j++){												//�ڶ��׶Σ���ԭ��һ����һ�¼���offset
			for (i = 0; i < width / 2; i++)
				new_arr[i][j] = arr[i][j - idy] + offset;
		}
		idy = height;
		offset = offset + portion;
		for (j = height; j < height + num; j++){											//�����׶Σ���ԭ��������һֱ
			for (i = 0; i < width; i++)
				new_arr[i][j] = arr[i][j - idy] + offset;
		}
		offset = offset + portion;
		idy = height / 2;
		for (j = height / 2; j < height; j++){												//���Ľ׶Σ���ԭ���Ĳ���һ�¼���offset
			for (i = width / 2; i < width; i++)
				new_arr[i][j] = arr[i][j - idy] + offset;
		}
		offset = offset + portion;
		for (j = 0; j < height / 2; j++){													//����׶Σ���ԭ���Ĳ���һ�¼���offset
			for (i = width / 2; i < width; i++)
				new_arr[i][j] = arr[i][j] + offset;
		}
	}
	//********�����Ҫ���ӵĲ��ֱ�ԭ�߶�С����ʹ�ó���ӷ�
	else{
		//���ڵõ����������ò�������
		for (i = 0; i < width; i++){
			for (j = 0; j < height; j++)
				orginal_arr[i][j] = arr[i][j];
		}
		idy = 0;
		idx = 0;
		for (i = 0; i < unit*unit; i++)															//�����ҵ�Ҫ������ӵĵ�Ԫ
			d2xy(unit, i, &row[i], &culumn[i]);
		//////////////��ʼƴ����һ���µ�
		offset = unit*unit / 4;
		portion = unit*unit / 4;
		pt = 0;
		pt_new = 0;
		for (i = pt_new; i < offset; i++){														//��һ�׶Σ���ԭ��һ��
			new_row[i] = row[i];
			new_culumn[i] = culumn[i];
		}
		idx += 0;
		idy += num;
		pt_new = offset;
		offset += portion;
		for (i = pt_new; i < offset; i++){														//�ڶ��׶Σ�����һ�׶�һ����x��ͬ��y+num
			new_row[i] = new_row[i - num*num];
			new_culumn[i] = new_culumn[i - num*num] + num;
		}
		pt_new = offset;
		pt += portion;
		offset += 2 * portion;
		idx += 0;
		idy += num;
		for (i = pt_new, j = pt; i < offset; i++, j++){											//�����׶Σ�����ԭ��Ԫ�ĵڶ�����
			new_row[i] = row[j];
			new_culumn[i] = culumn[j] + num;
		}
		pt_new = offset;
		pt += 2 * portion;
		idx += num;
		idy -= num;
		offset += portion;
		for (i = pt_new, j = pt; i < offset; i++, j++) {										//���Ľ׶Σ���ԭ��Ԫ���������һ�£�y+num;
			d2xy(unit, j, &new_row[i], &new_culumn[i]);
			new_culumn[i] += num;
		}
		pt_new = offset;
		idx += 0;
		idy -= num;
		offset += portion;
		for (i = pt_new, j = pt; i < offset; i++, j++) {										//����׶Σ���ԭ��Ԫ���������һ��
			new_row[i] = row[j];
			new_culumn[i] = culumn[j];
		}
		/////////��Ҫ�����ÿ��Ҫ�ӵĽ���λ������
		idx = 0;
		idy = height - num * 2;
		for (i = 0; i < subnum; i++, idx += unit)
			st[i] = orginal_arr[idx][idy];
		//�˴��ǰѴ��������Ǹ���ά����ת�����и��е�������Ӧx��y�����һά����
		for (i = 0; i < width; i++){
			for (j = 0; j < height; j++){
				orginal_culumn[orginal_arr[i][j]] = j;
				orginal_row[orginal_arr[i][j]] = i;
			}
		}
		//ʹ���ɺ�������ԭʼ����
		for (i = 0; i < width*height; i++){
			transit_row[i] = orginal_row[i];
			transit_culumn[i] = orginal_culumn[i];
		}
		idx = 0;
		idy = height - num * 2;
		/////////������Ҫ��������ȫ���ı�һ�Σ������ǵ�Ԫ���Լӵ���len
		for (i = 0; i < subnum; i++){															//ѭ�����	
			pointer = st[i];
			st[i + 1] = st[i + 1] + unit*num*(i + 1);
			counter += num*unit;
			for (j = 0; j < pointer; j++){														//��һ������ԭ����һ��
				sum_row[j] = transit_row[j];
				sum_culumn[j] = transit_culumn[j];
			}
			pointer += new_unit*unit;
			for (j = st[i], k = 0; j < pointer; j++, k++){										//�ڶ����ְ����޸Ĳ�����ӽ�ȥ
				sum_row[j] = new_row[k] + idx;
				sum_culumn[j] = new_culumn[k] + idy;
			}
			pointer0 = st[i] + unit*unit;
			for (j = pointer, k = pointer0; j < counter; j++, k++){								//����������ԭ����һ��
				sum_row[j] = transit_row[k];
				sum_culumn[j] = transit_culumn[k];
			}
			idx += unit;
			for (j = 0; j < sum; j++){															//ʹ����������ڵ�ǰ����������һ�μ���
				transit_row[j] = sum_row[j];
				transit_culumn[j] = sum_culumn[j];
			}
		}
		//������һά�����Ӧ��x��y������ָ��ɶ�ά����
		for (k = 0; k < width*(height + num); k++)
			new_arr[sum_row[k]][sum_culumn[k]] = k;
	}
	return new_arr;
}

//�ú�������ʵ�����ϼ��Ĺ��ܣ������ά�����Լ����ĳ��Ϳ��Լ�Ҫ��ȥ�ĳ���
int **sub(int **arr, int width, int height, int num)
{
	int i, j, k;
	int unit = num * 2;
	int unit_size = unit*unit;
	int new_unit = num * 2 - num;
	int new_unit_size = new_unit*unit;
	int subnum = width / unit;
	int pt, pt_new, pointer, pointer0;
	int counter = width*height;
	int idx, idy;
	int orginal = width*height;
	int sum = width*(height - num);

	int *st = (int *)malloc((subnum + 1)*sizeof(int *));
	int *orginal_row = (int *)malloc(orginal*sizeof(int *));									//ԭ�����ܳ�Ŷ
	int *orginal_culumn = (int *)malloc(orginal*sizeof(int *));
	int *transit_row = (int *)malloc(orginal*sizeof(int *));									//��������,��orginalһ����
	int *transit_culumn = (int *)malloc(orginal*sizeof(int *));
	int *transit_row1 = (int *)malloc(orginal*sizeof(int *));									//��������,��orginalһ����
	int *transit_culumn1 = (int *)malloc(orginal*sizeof(int *));
	int *sum_row = (int *)malloc(sum*sizeof(int *));											//�����ܵ��µ�������ܳ�������
	int *sum_culumn = (int *)malloc(sum*sizeof(int *));

	int *row = (int *)malloc((unit_size)*sizeof(int *));										//����Ҫsub���Ǹ�unit��x��������飿
	int *culumn = (int *)malloc((unit_size)*sizeof(int *));										//����Ҫsub���Ǹ�unit��y��������飿
	int *new_row = (int *)malloc((new_unit_size)*sizeof(int *));								//����ԭ�����Ǹ�unit�����ϲ��ֵ�x��������飿
	int *new_culumn = (int *)malloc((new_unit_size)*sizeof(int *));								//����ԭ���Ǹ�unit���ϲ��ֵ�y��������飿

	int **orginal_arr = (int **)malloc(width*sizeof(int *));									//������
	for (int i = 0; i<width; i++)
		orginal_arr[i] = (int *)malloc(height*sizeof(int));										//������

	//���ڵõ����������ò�������
	for (i = 0; i < width; i++){
		for (j = 0; j < height; j++)
			orginal_arr[i][j] = arr[i][j];
	}
	for (i = 0; i < unit_size; i++)																//�����ҵ�Ҫ������ӵĵ�Ԫ
		d2xy(unit, i, &row[i], &culumn[i]);
	idx = 0;
	idy = -1*(num);
	//��ʼƴ����һ���µģ���1��3����ɾ��
	pt = unit_size / 4; 
	pt_new = pt*3;
	for (i = pt, j = 0; i < pt_new; i++, j++){
		new_row[j] = row[i] + idx;
		new_culumn[j] = culumn[i] + idy;
	}
	//�Ȱ�Ҫ���м����λ���ҳ���
	idx = 0;
	idy = height - num * 2;
	for (i = 0; i < subnum; i++, idx += unit)
		st[i] = orginal_arr[idx][idy];
	//�˴��ǰѴ��������Ǹ���ά����ת�����и��е�������Ӧx��y�����һά����
	for (i = 0; i < width; i++){
		for (j = 0; j < height; j++){
			orginal_culumn[orginal_arr[i][j]] = j;
			orginal_row[orginal_arr[i][j]] = i;
		}
	}
	for (i = 0; i < width*height; i++){															//��ʼ���������飬������ԭ������ͬ
		transit_row[i] = orginal_row[i];
		transit_culumn[i] = orginal_culumn[i];
		transit_row1[i] = orginal_row[i];
		transit_culumn1[i] = orginal_culumn[i];
	}
	//��ʼ���µ�Ԫװ�ص�ԭ���Ĵ�������
	idx = 0;
	idy = height - unit;
	/////////������Ҫ��������ȫ���ı�һ�Σ������ǵ�Ԫ���Լӵ���len
	for (i = 0; i < subnum; i++){																//ѭ����
		pointer = st[i];
		st[i + 1] = st[i + 1] - new_unit_size*(i + 1);
		counter -= new_unit_size;
		for (j = 0; j < pointer; j++){															//��ԭ�����һ����һ��
			transit_row1[j] = transit_row[j];
			transit_culumn1[j] = transit_culumn[j];
		}
		pointer += new_unit_size;
		for (j = st[i], k = 0; j < pointer; j++, k++){											//�ڶ����������޸ĵ�Ԫ����
			transit_row1[j] = new_row[k] + idx;
			transit_culumn1[j] = new_culumn[k] + idy;
		}
		pointer0 = st[i] + unit_size;
		for (j = pointer, k = pointer0; j < counter; j++, k++){									//����������ԭ��������һ��
			transit_row1[j] = transit_row[k];
			transit_culumn1[j] = transit_culumn[k];
		}
		idx += unit;
		for (j = 0; j < orginal; j++){																
			transit_row[j] = transit_row1[j];
			transit_culumn[j] = transit_culumn1[j];
		}
	}
	for (i = 0; i < sum; i++){																	//�����������
		sum_row[i] = transit_row1[i];
		sum_culumn[i] = transit_culumn1[i];
	}
	//������һά�����Ӧ��x��y������ָ��ɶ�ά����
	for (k = 0; k < sum; k++)
		arr[sum_row[k]][sum_culumn[k]] = k;
	return arr;
}

