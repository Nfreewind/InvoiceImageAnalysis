#include "stdafx.h"
#include <windows.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <io.h>
#include <direct.h>

#include <opencv2/opencv.hpp>

#include "QczEllipsePkg.h"
#include "QczUtils.h"
#include "QczTemplateMaker.h"
#include "QczPreProcess.h"
#include <FastPCCR.h>

using namespace cv;
using namespace std;

void main(){

	string SrcImgFold = "F:\\qcz_pro\\实验室项目\\格式文档\\飞机票\\data\\Original_sub";
	string DstImgFold = "Img_Out_Rotate";
	if (_access(DstImgFold.c_str(), 0)){ _mkdir(DstImgFold.c_str()); }

	vector<string> ImgPaths;
	QczFile::getFiles_(SrcImgFold, ImgPaths);

	int a = InitializeRecognizer(".\\dictionary");
	if (a != 0){ cout << "init failed!---" << endl; return; }



	int minNum = 300;

	for (int i = 0; i < min((int)ImgPaths.size(), minNum); i++){
		string Name = QczFile::splitFileName(ImgPaths[i]);
		string savePath = DstImgFold + "\\" + Name;
		Mat im = imread(ImgPaths[i]);

		cout << "dealing with " << i << " th img : " << Name << endl;
		////////////////////////////////////////////////////////////////////////

		//vector<struct Ellipse> ellipses;
		//struct Ellipse eTemplate;
		//int idx = 0;
		//int Flag = 0;
		//double scale = 0.0;
		//Flag = EllipseDetector(im, ellipses,idx, eTemplate, scale, savePath);
		//if (Flag != 0){
		//	cout << "Detect " << i << " th img : " << Name << " failure!---" << endl;
		//}
		//else{
		//	cout << "Detect " << i << " th img : " << Name << " success!---" << endl;
		//}

		////////////////////////////////////////////////////////////////////////
		//if (i == 0){
		//	string TemplatePath = "template.txt";
		//	TemplateMaker(im, TemplatePath);
		//}
		////////////////////////////////////////////////////////////////////////

		//struct Ellipse TopEllipse, BottomEllipse;
		//Barcode barcode;
		//TemplateLoader(TemplatePath, TopEllipse, BottomEllipse, barcode);

		////////////////////////////////////////////////////////////////////////

		string TemplatePath = "template.txt";
		Mat out;
		Rect rect;
		string BarcodeStr;
		PreProcess(im, TemplatePath, out, rect, BarcodeStr, savePath, DstImgFold);



	}


}
