#include "Fusian.h"
#include "armadillo"
#include <map>
#define uchar unsigned char
using namespace arma;
using namespace std;
fcube MaxMinStretch2fcube(int rowHeightY, int colWidthX, int bandNum, fcube resFcube)
{
	fcube outCube(rowHeightY, colWidthX, bandNum);
	//qDebug() << outCube(0, 0, 0);

	//ֱ��ͼ����
	for (int band_n = 0; band_n < bandNum; band_n++)
	{
		//��ÿ�����ν��д���
		//��ò������ֵ
		float max = resFcube(0, 0, band_n);
		for (int i = 0; i < rowHeightY; i++) {
			for (int j = 0; j < colWidthX; j++) {
				if (resFcube(i, j, band_n) > max)
					max = resFcube(i, j, band_n);
			}
		}
		//��ò�����Сֵ
		float min = resFcube(0, 0, band_n);//������̨
		for (int i = 0; i < rowHeightY; i++) {
			for (int j = 0; j < colWidthX; j++) {
				if (resFcube(i, j, band_n) < min)
					min = resFcube(i, j, band_n);
			}
		}
		for (int i = 0; i < rowHeightY; i++)
		{
			for (int j = 0; j < colWidthX; j++)
			{
				outCube(i, j, band_n) = 255 * ((resFcube(i, j, band_n) - min) / (max - min));
			}
		}
	}
	return outCube;
}
uchar_cube PR(int rowHeightY, int colWidthX, int bandNum, fcube resFcube, fcube resFcube2)
{
	//ֱ��ͼ���⻯

	uchar_cube outCube1(rowHeightY, colWidthX, bandNum);
	for (int t = 0; t < resFcube.n_slices; ++t) {
		for (int i = 0; i < rowHeightY; i++)
		{
			for (int j = 0; j < colWidthX; j++)
			{
				outCube1(i, j, t) = (unsigned char)resFcube(i, j, t);
			}
		}
	}
	uchar_cube outCube2(rowHeightY, colWidthX, bandNum);
	for (int t = 0; t < resFcube2.n_slices; ++t)
	{
		for (int i = 0; i < rowHeightY; i++)
		{
			for (int j = 0; j < colWidthX; j++)
			{
				outCube2(i, j, t) = (unsigned char)resFcube2(i, j, t);
			}
		}
	}
	int gray[256] = { 0 };
	double gray_prob[256] = { 0 };
	double gray_distribution[256] = { 0 };
	int gray_equal[256] = { 0 };

	int gray2[256] = { 0 };
	double gray_prob2[256] = { 0 };
	double gray_distribution2[256] = { 0 };
	int gray_equal2[256] = { 0 };

	int gray_sum = 0;
	map<int, int> S2Z;
	map<int, int> r2z;
	for (int t = 0; t < bandNum; ++t)
	{
		gray_sum = colWidthX * rowHeightY;

		//ͳ��ÿ���Ҷ��µ����ظ���
		for (int i = 0; i < rowHeightY; i++)
			for (int j = 0; j < colWidthX; j++)
			{
				int vaule = outCube1(i, j, t);
				gray[vaule]++;
			}

		//ͳ�ƻҶ�Ƶ��
		for (int i = 0; i < 256; i++)
		{
			gray_prob[i] = ((double)gray[i] / gray_sum);
		}

		//�����ۼ��ܶ�
		gray_distribution[0] = gray_prob[0];
		for (int i = 1; i < 256; i++)
		{
			gray_distribution[i] = gray_distribution[i - 1] + gray_prob[i];
		}

		//���¼�����⻯��ĻҶ�ֵ
		for (int i = 0; i < 256; i++)
		{
			gray_equal[i] = (uchar)(255 * gray_distribution[i] + 0.5);
		}

		//ͳ��ÿ���Ҷ��µ����ظ���
		for (int i = 0; i < rowHeightY; i++)
			for (int j = 0; j < colWidthX; j++)
			{
				int vaule = outCube2(i, j, t);
				gray2[vaule]++;
			}

		//ͳ�ƻҶ�Ƶ��
		for (int i = 0; i < 256; i++)
		{
			gray_prob2[i] = ((double)gray2[i] / gray_sum);
		}

		//�����ۼ��ܶ�
		gray_distribution2[0] = gray_prob2[0];
		for (int i = 1; i < 256; i++)
		{
			gray_distribution2[i] = gray_distribution2[i - 1] + gray_prob2[i];
		}
		//���¼�����⻯��ĻҶ�ֵ���������롣�ο���ʽ��(N-1)*T+0.5
		for (int i = 0; i < 256; i++)
		{
			gray_equal2[i] = (uchar)(255 * gray_distribution2[i] + 0.5);
		}
		//��G(Z)=S ����S->Z��ӳ���
		for (int i = 0; i < 256; i++)
		{
			for (int j = 1; j < 256; j++)
			{
				if (abs(gray_equal[i] - gray_equal2[j - 1]) < abs(gray_equal[i] - gray_equal2[j]))
				{
					S2Z[gray_equal[i]] = j - 1;
					break;
				}
			}
		}
		S2Z[gray_equal[255]] = gray_equal2[255];

		for (int i = 0; i < 256; i++)
			r2z[i] = S2Z[gray_equal[i]];

		for (int i = 0; i < rowHeightY; i++)
		{
			for (int j = 0; j < colWidthX; j++)
			{
				outCube1(i, j, t) = r2z[outCube1(i, j, t)];
			}
		}
	}

	//�Ҷ�����
	for (int band_n = 0; band_n < bandNum; band_n++)
	{
		//��ÿ�����ν��д���
		//��ò������ֵ
		float max = outCube1(0, 0, band_n);
		for (int i = 0; i < rowHeightY; i++)
		{
			for (int j = 0; j < colWidthX; j++)
			{
				if (outCube1(i, j, band_n) > max)
					max = outCube1(i, j, band_n);
			}
		}
		//��ò�����Сֵ
		float min = outCube1(0, 0, band_n);
		for (int i = 0; i < rowHeightY; i++)
		{
			for (int j = 0; j < colWidthX; j++)
			{
				if (outCube1(i, j, band_n) < min)
					min = outCube1(i, j, band_n);
			}
		}
		for (int i = 0; i < rowHeightY; i++)
		{
			for (int j = 0; j < colWidthX; j++)
			{
				outCube1(i, j, band_n) = (unsigned char)255 * ((outCube1(i, j, band_n) - min) / (max - min));
			}
		}
	}
	return outCube1;
}


