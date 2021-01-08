#include "ImageProcessing.h"
#include<QMessageBox>
#include <QDebug>
#include<gdal.h>
#include <gdal_priv.h>
#include "cpl_conv.h" 
#include <fstream>
#include <string>
#include <QTextCodec>
#include "BWM.h"
#include "Filter.h"
using namespace std;
using namespace arma;
ImageProcessing::ImageProcessing(QWidget *parent):QMainWindow(parent)
{
    ui.setupUi(this); 
    connect(ui.pushButton, SIGNAL(clicked()), this, SLOT(on_loadButton_clicked()));
    connect(ui.pushButton_2, SIGNAL(clicked()), this, SLOT(on_selectInputPathButton_clicked()));
    connect(ui.pushButton_3, SIGNAL(clicked()), this, SLOT(on_selectOutputPathButton_clicked()));
	ui.lineEditInput->setText(inputpath);
	ui.lineEditOutput->setText(outputpath);
	connect(ui.actionGF, &QAction::triggered, this, &ImageProcessing::on_GF_selected);
	connect(ui.actionBWM, &QAction::triggered, this, &ImageProcessing::on_BWM_selected);
	connect(ui.actionMean_2, &QAction::triggered, this, &ImageProcessing::on_meanfilter_selected);
	connect(ui.actionMid_2, &QAction::triggered, this, &ImageProcessing::on_midfilter_selected);
	connect(ui.actionGau_2, &QAction::triggered, this, &ImageProcessing::on_gaufilter_selected);
	connect(ui.actionRoberts_2, &QAction::triggered, this, &ImageProcessing::on_roberts_seletected);
	connect(ui.actionPrewitt_2, &QAction::triggered, this, &ImageProcessing::on_prewitt_selected);
	connect(ui.actionSobel_2, &QAction::triggered, this, &ImageProcessing::on_sobel_selected);
	connect(ui.actionLaplacian_2, &QAction::triggered, this, &ImageProcessing::on_laplacian_selected);
	connect(ui.actionLoG, &QAction::triggered, this, &ImageProcessing::on_log_selected);
	connect(ui.action_Laplacian, &QAction::triggered, this, &ImageProcessing::on_lap2_selected);
	connect(ui.action_fusion, &QAction::triggered, this, &ImageProcessing::on_fusion_selected);


}
void ImageProcessing::on_loadButton_clicked()
{
	ui.errMsgEdit->textCursor().insertText("Loading data...\n");
	//**************1-从GDAL中读取影像像素到数组************************//
//打开影像，获取数据集
	GDALDataset* poDataset;
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "YES");    //设置GDAL支持中文路径
	qDebug() << this->inputpath.toStdString().c_str();
	poDataset = (GDALDataset*)GDALOpen(this->inputpath.toStdString().c_str(), GA_ReadOnly);
	
	if (poDataset == NULL)
	{
		QMessageBox::warning(this, "Warning", "Invaild input file name.");
		return;
	}
	//获取影像尺寸信息
	const int rowHeightY = poDataset->GetRasterYSize();;
	const int colWidthX = poDataset->GetRasterXSize();
	const int bandNum = poDataset->GetRasterCount();
	long bandSize = rowHeightY * colWidthX;

	//读取影像像素数据到数组
	float* pInputData = new float[rowHeightY * colWidthX * bandNum];
	unsigned short* tmpdata = new unsigned short[rowHeightY * colWidthX * bandNum];
	int* panBandMap = new int[bandNum];
	for (int i = 0; i < bandNum; ++i)
		panBandMap[i] = i + 1;
	if (poDataset->RasterIO(GF_Read, 0, 0, colWidthX, rowHeightY, tmpdata, colWidthX, rowHeightY,
		GDT_Int16, bandNum, panBandMap, 0, 0, 0) != CE_None)
		return;
	for (int i = 0; i < rowHeightY * colWidthX * bandNum; i++) {
		pInputData[i] = tmpdata[i];
	}
	delete[] tmpdata;
	//关闭数据集
	GDALClose((GDALDatasetH)poDataset);

	//影像像素数组存储行主序--列主序
	float* pProcessData = new float[rowHeightY * colWidthX * bandNum];
	for (int b = 0; b < bandNum; ++b)
		for (int i = 0; i < rowHeightY; ++i)
			for (int j = 0; j < colWidthX; ++j)
				pProcessData[b * bandSize + j * rowHeightY + i] = pInputData[b * bandSize + i * colWidthX + j];

	fcube myFcube(pProcessData, rowHeightY, colWidthX, bandNum);
	//**************2-图像处理************************//
	fcube rescube(rowHeightY, colWidthX, bandNum);
	switch (this->tool) {
		case 1: {
			ui.errMsgEdit->textCursor().insertText("Running Gaussian filter...\n");
			Filter f;
			f.addData(&myFcube);
			f.execute("hpfilter", "Gaussian");
			rescube.slice(0) = f.result.slice(0);
			rescube.slice(1) = myFcube.slice(1);
			rescube.slice(2) = myFcube.slice(2);
			break;
		}
		case 2: {
			ui.errMsgEdit->textCursor().insertText("Running image blind watermark...\n");
			BWM b;
			arma::fmat wm;
			wm.load("watermark.csv", csv_ascii);
			rescube = b.encode(myFcube, wm);
			break;
		}
		case 3: {
			ui.errMsgEdit->textCursor().insertText("Mean\n");
			rescube = myFcube;
			for (int s = 0; s < bandNum; s++) {
				fmat A = rescube.slice(s);
				for (int j = 1; j < rowHeightY - 1; j++)
				{
					for (int i = 1; i < colWidthX - 1; i++)
					{
						A(j, i) = A(j - 1, i - 1) + A(j, i - 1) + A(j + 1, i - 1) + A(j - 1, i) + A(j, i) + A(j + 1, i) + A(j - 1, i + 1) + A(j, i + 1) + A(j + 1, i + 1);
						A(j, i) = A(j, i) / 9;
					}
				}
				rescube.slice(s) = A;
			}

			break;
		}
		case 4: {
			ui.errMsgEdit->textCursor().insertText("Mid\n");
			rescube = myFcube;
			for (int s = 0; s < bandNum; s++) {
				fmat A = rescube.slice(s);
				fmat T = A;
				for (int j = 1; j < rowHeightY - 1; j++)
				{
					for (int i = 1; i < colWidthX - 1; i++)
					{
						fmat B(3, 3);
						for (int k = -1; k < 2; k++)
						{
							for (int l = -1; l < 2; l++)
							{
								B(k + 1, l + 1) = A(j + k, i + l);
							}
						}
						fvec v = B.as_col();
						T(j, i) = median(v);
					}
					A = T;
				}
				rescube.slice(s) = A;
			}
			break;
		}
		case 5: {
			ui.errMsgEdit->textCursor().insertText("高斯滤波\n");
			rescube = myFcube;
			for (int s = 0; s < bandNum; s++) {
				fmat A = rescube.slice(s);
				for (int j = 1; j < rowHeightY - 1; j++)
				{
					for (int i = 1; i < colWidthX - 1; i++)
					{
						A(j, i) = A(j - 1, i - 1) + 2 * A(j, i - 1) + A(j + 1, i - 1) + 2 * A(j - 1, i) + 4 * A(j, i) + 2 * A(j + 1, i) + A(j - 1, i + 1) + 2 * A(j, i + 1) + A(j + 1, i + 1);
						A(j, i) = A(j, i) / 16;
					}
				}
				rescube.slice(s) = A;
			}
			break;
		}
		case 6: {
			ui.errMsgEdit->textCursor().insertText("Roberts\n");
			rescube = myFcube;
			for (int s = 0; s < bandNum; s++) {
				fmat A = rescube.slice(s);
				fmat T = A;
				for (int j = 1; j < rowHeightY - 1; j++)
				{
					for (int i = 1; i < colWidthX - 1; i++)
					{
						T(j, i) = abs(A(j + 1, i + 1) - A(j, i)) + abs(A(j, i + 1) - A(j + 1, i));
					}
				}
				rescube.slice(s) = T;
			}
			break;
		}
		case 7: {
			ui.errMsgEdit->textCursor().insertText("Prewitt\n");
			rescube = myFcube;
			for (int s = 0; s < bandNum; s++) {
				fmat A = rescube.slice(s);
				fmat T = A;
				for (int j = 1; j < rowHeightY - 1; j++)
				{
					for (int i = 1; i < colWidthX - 1; i++)
					{
						T(j, i) = abs(A(j - 1, i - 1) + A(j - 1, i) + A(j - 1, i + 1) - A(j + 1, i - 1) - A(j + 1, i) - A(j + 1, i + 1))
							+ abs(A(j - 1, i - 1) + A(j, i - 1) + A(j + 1, i - 1) - A(j - 1, i + 1) - A(j, i + 1) - A(j + 1, i + 1));
					}
				}
				rescube.slice(s) = T;
			}
			break;
		}
		case 8: {
			ui.errMsgEdit->textCursor().insertText("Sobel\n");
			rescube = myFcube;
			for (int s = 0; s < bandNum; s++) {
				fmat A = rescube.slice(s);
				fmat T = A;
				for (int j = 1; j < rowHeightY - 1; j++)
				{
					for (int i = 1; i < colWidthX - 1; i++)
					{
						T(j, i) = abs(A(j - 1, i - 1) + 2 * A(j - 1, i) + A(j - 1, i + 1) - A(j + 1, i - 1) - 2 * A(j + 1, i) - A(j + 1, i + 1))
							+ abs(A(j - 1, i - 1) + 2 * A(j, i - 1) + A(j + 1, i - 1) - A(j - 1, i + 1) - 2 * A(j, i + 1) - A(j + 1, i + 1));
					}
				}
				rescube.slice(s) = T;
			}
			break;
		}
		case 9: {
			ui.errMsgEdit->textCursor().insertText("Laplacian\n");
			rescube = myFcube;
			for (int s = 0; s < bandNum; s++) {
				fmat A = rescube.slice(s);
				fmat T = A;
				for (int j = 1; j < rowHeightY - 1; j++)
				{
					for (int i = 1; i < colWidthX - 1; i++)
					{
						T(j, i) = abs(A(j - 1, i - 1) + 4 * A(j, i - 1) + A(j + 1, i - 1)
							+ 4 * A(j - 1, i) - 20 * A(j, i) + 4 * A(j + 1, i)
							+ A(j - 1, i + 1) + 4 * A(j, i + 1) + A(j + 1, i + 1));
					}
				}
				rescube.slice(s) = T;
			}
			break;
		}
		case 10: {
			ui.errMsgEdit->textCursor().insertText("LoG\n");
			rescube = myFcube;
			for (int s = 0; s < bandNum; s++) {
				fmat A = rescube.slice(s);
				fmat T = A;
				for (int j = 2; j < rowHeightY - 2; j++)
				{
					for (int i = 2; i < colWidthX - 2; i++)
					{
						T(j, i) = abs(-A(j - 2, i)
							- A(j - 1, i - 1) - 2 * A(j - 1, i) - A(j - 1, i + 1)
							- A(j, i - 2) - 2 * A(j, i - 1) + 16 * A(j, i) - 2 * A(j, i + 1) - A(j, i + 2)
							- A(j + 1, i - 1) - 2 * A(j + 1, i) - A(j + 1, i + 1)
							- A(j + 2, i)
						);
					}
				}
				rescube.slice(s) = T;
			}	
			break;
		}
		case 11: {
			ui.errMsgEdit->textCursor().insertText("Improved Laplacian\n");
			rescube = myFcube;
			for (int s = 0; s < bandNum; s++) {
				fmat A = rescube.slice(s);
				fmat T = A;
				for (int j = 1; j < rowHeightY - 1; j++)
				{
					for (int i = 1; i < colWidthX - 1; i++)
					{
						int a1 = A(j - 1, i - 1) + 4 * A(j, i - 1) + A(j + 1, i - 1)
							+ 4 * A(j - 1, i) - 20 * A(j, i) + 4 * A(j + 1, i)
							+ A(j - 1, i + 1) + 4 * A(j, i + 1) + A(j + 1, i + 1);
						T(j, i) = 1.8 * T(j, i) - a1;
					}
				}
				rescube.slice(s) = T;
			}
			break;
		}
	}
	//获取多维矩阵的数据裸指针
	const float* pResData = rescube.memptr();

	float* pOutputData = new float[bandSize * bandNum];

	//影像像素数组存储列主序--行主序
	for (int b = 0; b < bandNum; ++b)
		for (int i = 0; i < rowHeightY; ++i)
			for (int j = 0; j < colWidthX; ++j)
				pOutputData[b * bandSize + i * colWidthX + j] = pResData[b * bandSize + j * rowHeightY + i];
	//**************3-写出影像************************//
	//创建数据集
	const char* pszDstFilename = this->outputpath.toStdString().c_str();
	const char* pszFormat = "GTiff";
	GDALDriver* poDriver;
	char** papszMetadata;
	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	GDALDataset* poDstDS;
	char** papszOptions = NULL;
	poDstDS = poDriver->Create(pszDstFilename, colWidthX, rowHeightY, bandNum, GDT_Int16,
		papszOptions);
	//写出数据
	
	if (poDstDS->RasterIO(GF_Write, 0, 0, colWidthX, rowHeightY, pOutputData, colWidthX, rowHeightY, GDT_Float32, bandNum, panBandMap, 0, 0, 0) != CE_None)
	{
		QMessageBox::warning(this, "Warning", "Write raster failed.");
		return;
	}

	//关闭数据集
	GDALClose((GDALDatasetH)poDstDS);

	////清理内存
	delete[] pInputData;
	delete[] panBandMap;
	delete[] pProcessData;
	delete[] pOutputData;
	ui.errMsgEdit->textCursor().insertText("Finished...\n");
}
void ImageProcessing::on_selectInputPathButton_clicked()
{
	QString inputpath = QFileDialog::getOpenFileName(this, tr(""), "", tr("Tif(*.tif)"));
	ui.lineEditInput->setText(inputpath);
	ui.errMsgEdit->textCursor().insertText("Added input file path!\n");

}
void ImageProcessing::on_selectOutputPathButton_clicked()
{
	QString outputpath = QFileDialog::getOpenFileName(this, tr(""), "", tr("Tif(*.tif)"));
	ui.lineEditOutput->setText(outputpath);
	ui.errMsgEdit->textCursor().insertText("Added input file path!\n");
}
void ImageProcessing::on_GF_selected()
{
	this->tool = 1;
	ui.errMsgEdit->textCursor().insertText("High-pass filter(Gaussian)\n");
}
void ImageProcessing::on_BWM_selected()
{
	this->tool = 2;	
	ui.errMsgEdit->textCursor().insertText("Blind water mark! Please make sure the size of input image is 400*400.\n");

}
void ImageProcessing::on_meanfilter_selected()              //均值滤波
{
	this->tool = 3;
	ui.errMsgEdit->textCursor().insertText("Mean filter!\n");

}
void ImageProcessing::on_midfilter_selected()               //中值滤波
{
	this->tool = 4;
	ui.errMsgEdit->textCursor().insertText("Mid filter!\n");
}
void ImageProcessing::on_gaufilter_selected()               //高斯滤波
{
	this->tool = 5;
	ui.errMsgEdit->textCursor().insertText("Gaussian filter!\n");
}
void ImageProcessing::on_roberts_seletected()
{
	this->tool = 6;
	ui.errMsgEdit->textCursor().insertText("Image sharpening based on Roberts operator!\n");
}
void ImageProcessing::on_prewitt_selected()
{
	this->tool = 7;
	ui.errMsgEdit->textCursor().insertText("Image sharpening based on Prewitt operator!\n");
}
void ImageProcessing::on_sobel_selected()
{
	this->tool = 8;
	ui.errMsgEdit->textCursor().insertText("Image sharpening based on Sobel operator!\n");
}
void ImageProcessing::on_laplacian_selected()
{
	this->tool = 9;
	ui.errMsgEdit->textCursor().insertText("Image sharpening based on Laplacian operator!\n");
}
void ImageProcessing::on_log_selected()
{
	this->tool = 10;
	ui.errMsgEdit->textCursor().insertText("Edge detection based on LoG operator!\n");
}
void ImageProcessing::on_lap2_selected()
{
	this->tool = 11;
	ui.errMsgEdit->textCursor().insertText("Improved image sharpening method based on Laplacian operator!\n");
}

void ImageProcessing::on_fusion_selected()
{
	this->tool = 12;
	ui.errMsgEdit->textCursor().insertText("Image fusion!\n");

}
