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
		gray_sum = colWidthX* rowHeightY;

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
zzy::zzy(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);
	ui.lineEdit_input->setEnabled(false);
	ui.lineEdit_output->setEnabled(false);
	//信号槽连接
	connect(ui.pushButton_inputSelect, SIGNAL(clicked()), this, SLOT(on_inputSelect_clicked()));
	connect(ui.pushButton_outputSelect, SIGNAL(clicked()), this, SLOT(on_outputSelect_clicked()));
	connect(ui.pushButton_process, SIGNAL(clicked()), this, SLOT(on_process_clicked()));
}

void zzy::on_inputSelect_clicked()
{
	this->inputFilePath = QFileDialog::getOpenFileName(this, QStringLiteral("选择输入遥感影像"), "",
		tr("TIFF(*.tif *.tiff);;JPEG(*.jpg *.jpeg);;" "All files(*)"));
	ui.lineEdit_input->setText(this->inputFilePath);
}

void zzy::on_outputSelect_clicked()
{
	this->outputFilePath = QFileDialog::getSaveFileName(this, QStringLiteral("选择输出遥感影像"), "",
		tr("TIFF(*.tif *.tiff);;JPEG(*.jpg *.jpeg);;" "All files(*)"));
	ui.lineEdit_output->setText(this->outputFilePath);
}

void zzy::on_process_clicked()
{
	//**************1-从GDAL中读取影像像素到数组************************//
	//打开影像，获取数据集
	GDALDataset  *poDataset;
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");    //设置GDAL支持中文路径
	qDebug() << this->inputFilePath.toStdString().c_str();
	poDataset = (GDALDataset *)GDALOpen(this->inputFilePath.toStdString().c_str(), GA_ReadOnly);
	if (poDataset == NULL)
	{
		QMessageBox::warning(this, "Warning", "Invaild input file name.");
		return;
	}
	//获取影像尺寸信息
	const int rowHeightY = poDataset->GetRasterYSize();;
	const int colWidthX = poDataset->GetRasterXSize();
	const int bandNum = poDataset->GetRasterCount();
	long bandSize = rowHeightY*colWidthX;
	//读取影像像素数据到数组
	float *pInputData = new float[rowHeightY * colWidthX * bandNum];
	int *panBandMap = new int[bandNum];
	for (int i = 0; i < bandNum; ++i)
		panBandMap[i] = i + 1;
	if (poDataset->RasterIO(GF_Read, 0, 0, colWidthX, rowHeightY, pInputData, colWidthX, rowHeightY,
		GDT_Float32, bandNum, panBandMap, 0, 0, 0) != CE_None)
		return;
	//关闭数据集
	GDALClose((GDALDatasetH)poDataset);


	//影像像素数组存储行主序--列主序
	float *pProcessData = new float[rowHeightY * colWidthX * bandNum];
	for (int b = 0; b < bandNum; ++b)
		for (int i = 0; i < rowHeightY; ++i)
			for (int j = 0; j < colWidthX; ++j)
				pProcessData[b*bandSize + j*rowHeightY + i] = pInputData[b*bandSize + i*colWidthX + j];


	//将像素数组初始化为armadillo中的多维矩阵数据结构
	fcube myFcube(pProcessData, rowHeightY, colWidthX, bandNum);





	//**************2-处理影像************************//
	
	//将多光谱影像展开,重构数组
	//cube转换成二维数组，将原先的二维数组转换成列向量
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
	//计算相关系数矩阵
	fmat corMat=cor(mss_vec);
	//计算相关系数矩阵的特征值和特征向量

	fvec eigval;//特征值
	fmat eigvec;//特征向量矩阵
	arma::eig_sym(eigval, eigvec, corMat);//求值
	
	fvec characteristic_vec(bandNum_mss);//取出最大特征值对应的特征向量
	for (int t = 0; t < bandNum_mss; t++)
	{
		characteristic_vec(t)=eigvec(t, bandNum_mss-1);
	}
	//计算得到第一主成分图像：
	fmat first=mss_vec*characteristic_vec;//第一主成分图像

	//直方图规定化
	fmat transs(rowHeightY_mss, colWidthX_mss);//转换数组，得到第一主成分影像矩阵
	int count = 0;
	for (int i = 0; i < rowHeightY_mss; i++)
	{
		for (int j = 0; j < colWidthX_mss; j++)
		{
			transs(i, j) = first(j + count*colWidthX_mss, 0);
		}
		count++;
	}

	//输出第一主成分图像：
	fcube mss_pca_1(rowHeightY_mss, colWidthX_mss, 1);
	mss_pca_1.slice(0)= transs;

	//先进行拉伸，再匹配

	//拉伸全色影像
	fcube pan_forPR= MaxMinStretch2fcube(rowHeightY_pan, colWidthX_pan, 1, pan);
	//拉伸主成分影像
	fcube mss_pca_1_forPR=MaxMinStretch2fcube(rowHeightY_mss, colWidthX_mss, 1, mss_pca_1);
	//直方图规定
	uchar_cube outCube_afterPR = PR(rowHeightY_pan, colWidthX_pan, 1, pan_forPR, mss_pca_1_forPR);//规定化结果矩阵
	//获得主成分影像
	fmat mss_pca_vec(rowHeightY_mss* colWidthX_mss, bandNum_mss);//
	mss_pca_vec = mss_vec*eigvec;//主成分影像 
	//将规定化结果向量化
	fmat pan_vec(rowHeightY_pan* colWidthX_pan, 1);
	int count2 = 0;
	for (int i = 0; i < rowHeightY_pan; i++)
	{
		for (int j = 0; j < colWidthX_pan; j++)
		{
			pan_vec(j + count2*colWidthX_pan, 0) = outCube_afterPR(i, j, 0);//直接用全色替换，不进行匹配
		}
		count2 += 1;
	}
	//规定化后的全色影像替换主成分
	for (int i = 0; i < rowHeightY_pan* colWidthX_pan; i++)//480643*1
	{
		mss_pca_vec(i,2)= pan_vec(i, 0);
	}
	//反变换
	fmat mid_mss_pca_vec(rowHeightY_pan* colWidthX_pan, bandNum_mss);//中间480643*7
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
	

	//获取多维矩阵的数据裸指针
	float *pResData = resFcube.memptr();

	//影像像素数组存储列主序--行主序
	float *pOutputData = new float[rowHeightY * colWidthX * bandNum];
	for (int b = 0; b < bandNum; ++b)
		for (int i = 0; i < rowHeightY; ++i)
			for (int j = 0; j < colWidthX; ++j)
				pOutputData[b*bandSize + i * colWidthX + j] = pResData[b*bandSize + j * rowHeightY + i];

	//**************3-写出影像************************//
	//创建数据集
	const char* pszDstFilename = this->outputFilePath.toStdString().c_str();
	const char *pszFormat = "GTiff";
	GDALDriver *poDriver;
	char **papszMetadata;
	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	GDALDataset *poDstDS;
	char **papszOptions = NULL;
	poDstDS = poDriver->Create(pszDstFilename, colWidthX, rowHeightY, bandNum, GDT_Float32,
		papszOptions);
	//写出数据
	if (poDstDS->RasterIO(GF_Write, 0, 0, colWidthX, rowHeightY, pOutputData, colWidthX, rowHeightY, GDT_Float32, bandNum, panBandMap, 0, 0, 0) != CE_None)
	{
		QMessageBox::warning(this, "Warning", "Write raster failed.");
		return;
	}
	//关闭数据集
	GDALClose((GDALDatasetH)poDstDS);

	//清理内存
	delete[] pInputData;
	delete[] panBandMap;
	delete[] pProcessData;
	//delete[] pResData;
	delete[] pOutputData;

	//处理提示
	QMessageBox::information(this, "Information", "Successfully!");

}
