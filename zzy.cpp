#include "zzy.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QString>
#include <QDebug>


#include "gdal_priv.h"
#include "cpl_conv.h" 
#include "armadillo"

using namespace arma;
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
		gray_sum = colWidthX* rowHeightY;

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
zzy::zzy(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);
	ui.lineEdit_input->setEnabled(false);
	ui.lineEdit_output->setEnabled(false);
	//�źŲ�����
	connect(ui.pushButton_inputSelect, SIGNAL(clicked()), this, SLOT(on_inputSelect_clicked()));
	connect(ui.pushButton_outputSelect, SIGNAL(clicked()), this, SLOT(on_outputSelect_clicked()));
	connect(ui.pushButton_process, SIGNAL(clicked()), this, SLOT(on_process_clicked()));
}

void zzy::on_inputSelect_clicked()
{
	this->inputFilePath = QFileDialog::getOpenFileName(this, QStringLiteral("ѡ������ң��Ӱ��"), "",
		tr("TIFF(*.tif *.tiff);;JPEG(*.jpg *.jpeg);;" "All files(*)"));
	ui.lineEdit_input->setText(this->inputFilePath);
}

void zzy::on_outputSelect_clicked()
{
	this->outputFilePath = QFileDialog::getSaveFileName(this, QStringLiteral("ѡ�����ң��Ӱ��"), "",
		tr("TIFF(*.tif *.tiff);;JPEG(*.jpg *.jpeg);;" "All files(*)"));
	ui.lineEdit_output->setText(this->outputFilePath);
}

void zzy::on_process_clicked()
{
	//**************1-��GDAL�ж�ȡӰ�����ص�����************************//
	//��Ӱ�񣬻�ȡ���ݼ�
	GDALDataset  *poDataset;
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");    //����GDAL֧������·��
	qDebug() << this->inputFilePath.toStdString().c_str();
	poDataset = (GDALDataset *)GDALOpen(this->inputFilePath.toStdString().c_str(), GA_ReadOnly);
	if (poDataset == NULL)
	{
		QMessageBox::warning(this, "Warning", "Invaild input file name.");
		return;
	}
	//��ȡӰ��ߴ���Ϣ
	const int rowHeightY = poDataset->GetRasterYSize();;
	const int colWidthX = poDataset->GetRasterXSize();
	const int bandNum = poDataset->GetRasterCount();
	long bandSize = rowHeightY*colWidthX;
	//��ȡӰ���������ݵ�����
	float *pInputData = new float[rowHeightY * colWidthX * bandNum];
	int *panBandMap = new int[bandNum];
	for (int i = 0; i < bandNum; ++i)
		panBandMap[i] = i + 1;
	if (poDataset->RasterIO(GF_Read, 0, 0, colWidthX, rowHeightY, pInputData, colWidthX, rowHeightY,
		GDT_Float32, bandNum, panBandMap, 0, 0, 0) != CE_None)
		return;
	//�ر����ݼ�
	GDALClose((GDALDatasetH)poDataset);


	//Ӱ����������洢������--������
	float *pProcessData = new float[rowHeightY * colWidthX * bandNum];
	for (int b = 0; b < bandNum; ++b)
		for (int i = 0; i < rowHeightY; ++i)
			for (int j = 0; j < colWidthX; ++j)
				pProcessData[b*bandSize + j*rowHeightY + i] = pInputData[b*bandSize + i*colWidthX + j];


	//�����������ʼ��Ϊarmadillo�еĶ�ά�������ݽṹ
	fcube myFcube(pProcessData, rowHeightY, colWidthX, bandNum);





	//**************2-����Ӱ��************************//
	
	//�������Ӱ��չ��,�ع�����
	//cubeת���ɶ�ά���飬��ԭ�ȵĶ�ά����ת����������
	fmat mss_vec(rowHeightY_mss* colWidthX_mss, bandNum_mss);
	for (int t = 0; t < bandNum_mss; t++)
	{
		int count = 0;
		for (int i = 0; i < rowHeightY_mss; i++)
		{
			for (int j = 0; j < colWidthX_mss; j++)
			{
				mss_vec(j+count*colWidthX_mss,t)=mss(i, j, t);
			}
			count += 1;
		}
	}
	//�������ϵ������
	fmat corMat=cor(mss_vec);
	//�������ϵ�����������ֵ����������

	fvec eigval;//����ֵ
	fmat eigvec;//������������
	arma::eig_sym(eigval, eigvec, corMat);//��ֵ
	
	fvec characteristic_vec(bandNum_mss);//ȡ���������ֵ��Ӧ����������
	for (int t = 0; t < bandNum_mss; t++)
	{
		characteristic_vec(t)=eigvec(t, bandNum_mss-1);
	}
	//����õ���һ���ɷ�ͼ��
	fmat first=mss_vec*characteristic_vec;//��һ���ɷ�ͼ��

	//ֱ��ͼ�涨��
	fmat transs(rowHeightY_mss, colWidthX_mss);//ת�����飬�õ���һ���ɷ�Ӱ�����
	int count = 0;
	for (int i = 0; i < rowHeightY_mss; i++)
	{
		for (int j = 0; j < colWidthX_mss; j++)
		{
			transs(i, j) = first(j + count*colWidthX_mss, 0);
		}
		count++;
	}

	//�����һ���ɷ�ͼ��
	fcube mss_pca_1(rowHeightY_mss, colWidthX_mss, 1);
	mss_pca_1.slice(0)= transs;

	//�Ƚ������죬��ƥ��

	//����ȫɫӰ��
	fcube pan_forPR= MaxMinStretch2fcube(rowHeightY_pan, colWidthX_pan, 1, pan);
	//�������ɷ�Ӱ��
	fcube mss_pca_1_forPR=MaxMinStretch2fcube(rowHeightY_mss, colWidthX_mss, 1, mss_pca_1);
	//ֱ��ͼ�涨
	uchar_cube outCube_afterPR = PR(rowHeightY_pan, colWidthX_pan, 1, pan_forPR, mss_pca_1_forPR);//�涨���������
	//������ɷ�Ӱ��
	fmat mss_pca_vec(rowHeightY_mss* colWidthX_mss, bandNum_mss);//
	mss_pca_vec = mss_vec*eigvec;//���ɷ�Ӱ�� 
	//���涨�����������
	fmat pan_vec(rowHeightY_pan* colWidthX_pan, 1);
	int count2 = 0;
	for (int i = 0; i < rowHeightY_pan; i++)
	{
		for (int j = 0; j < colWidthX_pan; j++)
		{
			pan_vec(j + count2*colWidthX_pan, 0) = outCube_afterPR(i, j, 0);//ֱ����ȫɫ�滻��������ƥ��
		}
		count2 += 1;
	}
	//�涨�����ȫɫӰ���滻���ɷ�
	for (int i = 0; i < rowHeightY_pan* colWidthX_pan; i++)//480643*1
	{
		mss_pca_vec(i,2)= pan_vec(i, 0);
	}
	//���任
	fmat mid_mss_pca_vec(rowHeightY_pan* colWidthX_pan, bandNum_mss);//�м�480643*7
	mid_mss_pca_vec = mss_pca_vec*trans(eigvec);
	fcube outTrans(rowHeightY_pan, colWidthX_pan, bandNum_mss);
	for (int t = 0; t < bandNum_mss; t++)
	{
		int count3 = 0;
		for (int i = 0; i < rowHeightY_pan; i++)
		{
			for (int j = 0; j < colWidthX_pan; j++)
			{
				outTrans(i, j, t)=mid_mss_pca_vec(j + count3*colWidthX_pan, t);
			}
			count3 += 1;
		}

	}
	

	//��ȡ��ά�����������ָ��
	float *pResData = resFcube.memptr();

	//Ӱ����������洢������--������
	float *pOutputData = new float[rowHeightY * colWidthX * bandNum];
	for (int b = 0; b < bandNum; ++b)
		for (int i = 0; i < rowHeightY; ++i)
			for (int j = 0; j < colWidthX; ++j)
				pOutputData[b*bandSize + i * colWidthX + j] = pResData[b*bandSize + j * rowHeightY + i];

	//**************3-д��Ӱ��************************//
	//�������ݼ�
	const char* pszDstFilename = this->outputFilePath.toStdString().c_str();
	const char *pszFormat = "GTiff";
	GDALDriver *poDriver;
	char **papszMetadata;
	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	GDALDataset *poDstDS;
	char **papszOptions = NULL;
	poDstDS = poDriver->Create(pszDstFilename, colWidthX, rowHeightY, bandNum, GDT_Float32,
		papszOptions);
	//д������
	if (poDstDS->RasterIO(GF_Write, 0, 0, colWidthX, rowHeightY, pOutputData, colWidthX, rowHeightY, GDT_Float32, bandNum, panBandMap, 0, 0, 0) != CE_None)
	{
		QMessageBox::warning(this, "Warning", "Write raster failed.");
		return;
	}
	//�ر����ݼ�
	GDALClose((GDALDatasetH)poDstDS);

	//�����ڴ�
	delete[] pInputData;
	delete[] panBandMap;
	delete[] pProcessData;
	//delete[] pResData;
	delete[] pOutputData;

	//������ʾ
	QMessageBox::information(this, "Information", "Successfully!");

}
