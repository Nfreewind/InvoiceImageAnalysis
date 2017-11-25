#include "stdafx.h"

#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cv.h>
#include <highgui.h>
#include <float.h>
#include <fstream>


#include "recogBarcodeDLL.h"
#include "QczPreProcess.h"
#include <FastPCCR.h>


using namespace std;
using namespace cv;

int recogBarcode(std::string PlanePath, std::string& BarcodeStr){

	int a = InitializeRecognizer(".\\dictionary");
	if (a != 0){ cout << "init failed!---" << endl; return -4; }
	else{ cout << "init success!---" << endl; }
	Mat im = imread(PlanePath);
	if (!im.data){
		return -3;
	}
	else{
		cout << "read image success!----" << endl;
	}

	Mat out;
	Rect rect;
	string TemplatePath = "Template/template.txt";

	int flag = PreProcess(im, TemplatePath, out, rect, BarcodeStr);
	return flag;
}
