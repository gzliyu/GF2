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

	//直方图拉伸
	for (int band_n = 0; band_n < bandNum; band_n++)
	{
		//对每个波段进行处理
		//获得波段最大值
		float max = resFcube(0, 0, band_n);
		for (int i = 0; i < rowHeightY; i++) {
			for (int j = 0; j < colWidthX; j++) {
				if (resFcube(i, j, band_n) > max)
					max = resFcube(i, j, band_n);
			}
		}
		//获得波段最小值
		float min = resFcube(0, 0, band_n);//擂主上台
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
	//直方图均衡化

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

		//统计每个灰度下的像素个数
		for (int i = 0; i < rowHeightY; i++)
			for (int j = 0; j < colWidthX; j++)
			{
				int vaule = outCube1(i, j, t);
				gray[vaule]++;
			}

		//统计灰度频率
		for (int i = 0; i < 256; i++)
		{
			gray_prob[i] = ((double)gray[i] / gray_sum);
		}

		//计算累计密度
		gray_distribution[0] = gray_prob[0];
		for (int i = 1; i < 256; i++)
		{
			gray_distribution[i] = gray_distribution[i - 1] + gray_prob[i];
		}

		//重新计算均衡化后的灰度值
		for (int i = 0; i < 256; i++)
		{
			gray_equal[i] = (uchar)(255 * gray_distribution[i] + 0.5);
		}

		//统计每个灰度下的像素个数
		for (int i = 0; i < rowHeightY; i++)
			for (int j = 0; j < colWidthX; j++)
			{
				int vaule = outCube2(i, j, t);
				gray2[vaule]++;
			}

		//统计灰度频率
		for (int i = 0; i < 256; i++)
		{
			gray_prob2[i] = ((double)gray2[i] / gray_sum);
		}

		//计算累计密度
		gray_distribution2[0] = gray_prob2[0];
		for (int i = 1; i < 256; i++)
		{
			gray_distribution2[i] = gray_distribution2[i - 1] + gray_prob2[i];
		}
		//重新计算均衡化后的灰度值，四舍五入。参考公式：(N-1)*T+0.5
		for (int i = 0; i < 256; i++)
		{
			gray_equal2[i] = (uchar)(255 * gray_distribution2[i] + 0.5);
		}
		//令G(Z)=S 建立S->Z的映射表
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

	//灰度拉伸
	for (int band_n = 0; band_n < bandNum; band_n++)
	{
		//对每个波段进行处理
		//获得波段最大值
		float max = outCube1(0, 0, band_n);
		for (int i = 0; i < rowHeightY; i++)
		{
			for (int j = 0; j < colWidthX; j++)
			{
				if (outCube1(i, j, band_n) > max)
					max = outCube1(i, j, band_n);
			}
		}
		//获得波段最小值
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


