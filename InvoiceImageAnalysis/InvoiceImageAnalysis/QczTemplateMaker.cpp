#include "stdafx.h"

#include "QczTemplateMaker.h"
#include "AntigerResize.h"
#include "AntigerColor.h"
#include <string>
#include <iostream>
#include <fstream>
#include "QczUtils.h"

void TemplateMaker(Mat im, string TemplatePath){

	Mat InputImg, RedImg;

	vector<Barcode> barcodes;
	//CBarcodeReader barcodeReader;
	HANDLE hBarcode;
	InitBarcodeReader(hBarcode);

	im.copyTo(InputImg);

	BarcodeRead(InputImg, hBarcode, barcodes);


	double ratio;
	Mat NormImg;
	ratio = ImgResizeWithFixedLongside(InputImg, NormImg, 640);
	GetRednessMap(NormImg, RedImg);

	Size sz = RedImg.size();
	int		iThLength = 16; //16
	float	fThObb = 3.0f;
	float	fThPos = 1.0f;
	float	fTaoCenters = 0.05f;
	int 	iNs = 16;
	float	fMaxCenterDistance = sqrt(float(sz.width*sz.width + sz.height*sz.height)) * fTaoCenters;

	float	fThScoreScore = 0.5f;//0.5f

	// Other constant parameters settings. 

	// Gaussian filter parameters, in pre-processing
	Size	szPreProcessingGaussKernelSize = Size(5, 5);
	double	dPreProcessingGaussSigma = 1.0;

	float	fDistanceToEllipseContour = 0.1f;	// (Sect. 3.3.1 - Validation)
	float	fMinReliability = 0.7f;	// Const parameters to discard bad ellipses 0.54£¬ 0.7


	// Initialize Detector with selected parameters
	CEllipseDetectorYaed yaed;
	yaed.SetParameters(szPreProcessingGaussKernelSize,
		dPreProcessingGaussSigma,
		fThPos,
		fMaxCenterDistance,
		iThLength,
		fThObb,
		fDistanceToEllipseContour,
		fThScoreScore,
		fMinReliability,
		iNs
		);
	// Detect
	vector<struct Ellipse> ellsYaed;
	Mat1b gray;
	gray = RedImg;
	yaed.Detect(gray, ellsYaed);

	int idx = 0;
	struct Ellipse ellipse;
	double maxRatio = 0.0;
	if (ellsYaed.size() == 1){
		ellipse = ellsYaed[0];
	}
	else{
		for (int i = 0; i < ellsYaed.size(); i++){
			double tmpRatio = max(ellsYaed[i]._a, ellsYaed[i]._b) / min(ellsYaed[i]._a, ellsYaed[i]._b);
			if (tmpRatio > maxRatio){ maxRatio = tmpRatio; idx = i; }
		}
	}
	ellipse = ellsYaed[idx];

	// resize the ellipse to the src img
	for (int i = 0; i < ellsYaed.size(); i++){
		ellsYaed[i]._xc = (int)(ellsYaed[i]._xc * ratio); ellsYaed[i]._yc = (int)(ellsYaed[i]._yc * ratio);
		ellsYaed[i]._a = (int)(ellsYaed[i]._a * ratio); ellsYaed[i]._b = (int)(ellsYaed[i]._b * ratio);
	}


	Mat3b resultImage = InputImg.clone();
	yaed.DrawDetectedEllipses(resultImage, ellsYaed, 0, 2);

	imwrite("template_show.jpg", resultImage);
	ofstream outfile(TemplatePath.c_str());

	
	outfile << ellsYaed[0]._xc << endl << ellsYaed[0]._yc << endl << ellsYaed[0]._a << endl << ellsYaed[0]._b << endl << ellsYaed[0]._rad << endl;
	outfile << ellsYaed[1]._xc << endl << ellsYaed[1]._yc << endl << ellsYaed[1]._a << endl << ellsYaed[1]._b << endl << ellsYaed[1]._rad << endl;
	
	string src(barcodes[0].sBarCodeData);
	string re = QczStr::replace_all_distinct(src, " ", "_");
	outfile << barcodes[0].iAngle << endl << barcodes[0].iType << endl << re << endl;
	outfile << barcodes[0].rcBarcodeRegion.x << endl << barcodes[0].rcBarcodeRegion.y << endl << barcodes[0].rcBarcodeRegion.width << endl << barcodes[0].rcBarcodeRegion.height << endl;

	/*outfile << ellsYaed[0]._xc << "\t" << ellsYaed[0]._yc << "\t" << ellsYaed[0]._a << "\t" << ellsYaed[0]._b << "\t" << ellsYaed[0]._rad << endl;
	outfile << ellsYaed[1]._xc << "\t" << ellsYaed[1]._yc << "\t" << ellsYaed[1]._a << "\t" << ellsYaed[1]._b << "\t" << ellsYaed[1]._rad << endl;
	outfile << barcodes[0].iAngle << "\t" << barcodes[0].iType << "\t" << barcodes[0].sBarCodeData << "\t";
	outfile << barcodes[0].rcBarcodeRegion.x << "\t" << barcodes[0].rcBarcodeRegion.y << "\t" << barcodes[0].rcBarcodeRegion.width << "\t" << barcodes[0].rcBarcodeRegion.height << endl;*/
	outfile.close();
}

void TemplateLoader(string TemplatePath, struct Ellipse& TopSeal, struct Ellipse& BottomSeal, Barcode& barcode, Rect& BarcodeStrRect){


	ifstream infile(TemplatePath);
	infile >> TopSeal._xc >> TopSeal._yc >> TopSeal._a >> TopSeal._b >> TopSeal._rad;
	infile >> BottomSeal._xc >> BottomSeal._yc >> BottomSeal._a >> BottomSeal._b >> BottomSeal._rad;
	infile >> barcode.iAngle >> barcode.iType;
	string str,re;
	infile >> str;
	re = QczStr::replace_all_distinct(str, "_", " ");
	strcmp(barcode.sBarCodeData, re.c_str());
	//re.copy(barcode.sBarCodeData, re.length(), 0);
	//*(barcode.sBarCodeData + re.length()) = '/0';
	infile >> barcode.rcBarcodeRegion.x >> barcode.rcBarcodeRegion.y >> barcode.rcBarcodeRegion.width >> barcode.rcBarcodeRegion.height;
	infile >> BarcodeStrRect.x >> BarcodeStrRect.y >> BarcodeStrRect.width >> BarcodeStrRect.height;
	//file["TopSeal_xc"] >> TopSeal._xc; file["TopSeal_yc"] >> TopSeal._yc;
	//file["TopSeal_a"] >> TopSeal._a; file["TopSeal_b"] >> TopSeal._b;
	//file["TopSeal_rad"] >> TopSeal._rad;

	//file["BottomSeal_xc"] >> BottomSeal._xc; file["BottomSeal_yc"] >> BottomSeal._yc;
	//file["BottomSeal_a"] >> BottomSeal._a; file["BottomSeal_b"] >> BottomSeal._b;
	//file["BottomSeal_rad"] >> BottomSeal._rad;

	//file["Barcode_iType"] >> barcode.iType;
	//file["Barcode_rcBarcodeRegion"] >> barcode.rcBarcodeRegion;
	//string data;
	//file["Barcode_sBarCodeData"] >> data;

	//int i = 0;
	//for (i = 0; i < data.length(); i++){ barcode.sBarCodeData[i] = data[i]; }
	//barcode.sBarCodeData[i] = '\0';
	//data.copy(barcode.sBarCodeData, (int)data.length(), 0);
	//*(barcode.sBarCodeData + (int)data.length()) = '/0';

}