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
    void on_GF_selected();                      //Ƶ�����˲�����˹�˲�����
    void on_BWM_selected();                     //äˮӡ
    void on_meanfilter_selected();              //��ֵ�˲�
    void on_midfilter_selected();               //��ֵ�˲�
    void on_gaufilter_selected();               //��˹�˲�
    void on_roberts_seletected();               //Robert ������
    void on_prewitt_selected();                 //Prewitt��
    void on_sobel_selected();                   //sobel��
    void on_laplacian_selected();               //laplacian��
    void on_log_selected();                     //LOG��Ե���
    void on_lap2_selected();                    //��������������˹�˲�
    void on_fusion_selected();
};
