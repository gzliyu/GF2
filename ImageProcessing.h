#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_ImageProcessing.h"
#include<qstring.h>
#include "armadillo"
#include <QFileDialog>
#include <QLineEdit>
#include <QTextEdit>
#include <QAction>
class ImageProcessing : public QMainWindow
{
    Q_OBJECT

public:
    ImageProcessing(QWidget *parent = Q_NULLPTR);
private:
    Ui::ImageProcessingClass ui;
    QString inputpath = "D:\\data\\RS_TEST_DATA\\FUSION\\GF2-MSS.tif";
    QString outputpath="D:\\data\\RS_TEST_DATA\\FUSION\\test.tif";
    int tool = -1;
    
private slots:
    void on_loadButton_clicked();
    void on_selectInputPathButton_clicked();
    void on_selectOutputPathButton_clicked();
    void on_GF_selected();                      //频率域滤波（高斯滤波器）
    void on_BWM_selected();                     //盲水印
    void on_meanfilter_selected();              //均值滤波
    void on_midfilter_selected();               //中值滤波
    void on_gaufilter_selected();               //高斯滤波
    void on_roberts_seletected();               //Robert 交叉锐化
    void on_prewitt_selected();                 //Prewitt锐化
    void on_sobel_selected();                   //sobel锐化
    void on_laplacian_selected();               //laplacian锐化
    void on_log_selected();                     //LOG边缘检测
    void on_lap2_selected();                    //更厉害的拉普拉斯滤波
    void on_fusion_selected();
};
