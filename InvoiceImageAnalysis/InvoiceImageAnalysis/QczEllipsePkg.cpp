#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <io.h>
#include <direct.h>


#include <opencv2/opencv.hpp>

#include "QczEllipsePkg.h"
#include "AntigerResize.h"
#include "AntigerColor.h"

using namespace cv;
using namespace std;


int EllipseDetector(Mat InputImg, vector<struct Ellipse>& ellipses, int& idx, struct Ellipse Template, double& scale, string savePath){

	//int w, h;
	//w = InputImg.size().width;
	//h = InputImg.size().height;
	//if (w > h)
	//{
	//	transpose(InputImg, InputImg);
	//	flip(InputImg, InputImg, 0);
	//}
	double ratio;
	Mat NormImg, RedImg;
	ratio = ImgResizeWithFixedLongside(InputImg, NormImg, 640);
	
	GetRednessMap(NormImg, RedImg);

	// Parameters Settings (Sect. 4.2)
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
	float	fMinReliability = 0.7f;	// Const parameters to discard bad ellipses 0.54�� 0.7


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
	//vector<struct Ellipse> ellsYaed;
	Mat1b gray;
	gray = RedImg;
	yaed.Detect(gray, ellipses);

	if (ellipses.size() == 0){
		return -1;
	}

	// Recoginze the top seal
	idx = 0;
	struct Ellipse ellipse;
	double maxRatio = 0.0;
	if (ellipses.size() == 1){
		ellipse = ellipses[0];
	}
	else{
		for (int i = 0; i < ellipses.size(); i++){
			double tmpRatio = max(ellipses[i]._a, ellipses[i]._b) / min(ellipses[i]._a, ellipses[i]._b);
			if (tmpRatio > maxRatio){ maxRatio = tmpRatio; idx = i; }
		}
	}
	ellipse = ellipses[idx];
	// Compare with the template and generate the output scale
	scale = (ellipse._a / Template._a + ellipse._b / Template._b) / 2;
	scale *= ratio;

	// resize the ellipse to the src img
	for (int i = 0; i < ellipses.size(); i++){
		ellipses[i]._xc = (int)(ellipses[i]._xc * ratio); ellipses[i]._yc = (int)(ellipses[i]._yc * ratio);
		ellipses[i]._a = (int)(ellipses[i]._a * ratio); ellipses[i]._b = (int)(ellipses[i]._b * ratio);
	}
	vector<struct Ellipse> out;
	out.push_back(ellipse);

	Mat3b resultImage = InputImg.clone();
	yaed.DrawDetectedEllipses(resultImage, out, 0, 3);

	if (savePath != ""){
		imwrite(savePath, resultImage);
	}

	return 0;
}

