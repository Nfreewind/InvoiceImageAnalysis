#include "stdafx.h"

#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cv.h>
#include <highgui.h>
#include <float.h>
#include <fstream>


#include <opencv2/opencv.hpp>

#include "AntigerBarcode.h"
#include "AntigerEllipse.h"
#include "QczUtils.h"
#include "QczPreProcess.h"
#include "QczTemplateMaker.h"
#include "QczEllipsePkg.h"
#include <FastPCCR.h>

#define DEBUG_SEAL 0
#define DEBUG_RECOG 0

int PreProcess(Mat im, string TemplatePath, Mat& out, Rect& rect, string& BarcodeStr, string savePath, string fold){

	int step = 0;

	//MatchPoints matchPoints, matchPointsTem;
	vector<struct Ellipse> ellipses;
	int idx = 0;
	//	load template
	struct Ellipse TemTopEllipse, TemBottomEllipse;
	Barcode barcode;
	Rect BarcodeStrRect;
	TemplateLoader(TemplatePath, TemTopEllipse, TemBottomEllipse, barcode, BarcodeStrRect);

	if (TemTopEllipse._a == 0 || TemTopEllipse._b == 0 || TemBottomEllipse._a == 0 || TemBottomEllipse._b == 0){
		return -5;
	}

	//matchPointsTem.Barcode = Point(barcode.rcBarcodeRegion.x + barcode.rcBarcodeRegion.width / 2, barcode.rcBarcodeRegion.y + barcode.rcBarcodeRegion.height / 2);
	//matchPointsTem.BarcodeStr = Point(BarcodeStrRect.x + BarcodeStrRect.width / 2, BarcodeStrRect.y + BarcodeStrRect.height / 2);
	//matchPointsTem.TopSeal = Point(TemTopEllipse._xc, TemTopEllipse._yc);
	//matchPointsTem.BottomSeal = Point(TemBottomEllipse._xc, TemBottomEllipse._yc);


	Mat InputImg;
	im.copyTo(InputImg);
	int w, h;
	w = InputImg.size().width;
	h = InputImg.size().height;
	if (w > h)
	{
		transpose(InputImg, InputImg);
		flip(InputImg, InputImg, 1);
	}
	/*if (index.size() == 1){
		matchPoints.TopSeal = Point(ellipses[index[0]]._xc, ellipses[index[0]]._yc);
	}
	else if (index.size() == 2){
		matchPoints.TopSeal = Point(ellipses[index[0]]._xc, ellipses[index[0]]._yc);
		matchPoints.BottomSeal = Point(ellipses[index[1]]._xc, ellipses[index[1]]._yc);
	}*/
	//detect barcode

	vector<Barcode> barcodes;
	//CBarcodeReader barcodeReader;
	HANDLE hBarcode;
	InitBarcodeReader(hBarcode);


	BarcodeRead(InputImg, hBarcode, barcodes);
	GlobalFree(hBarcode);
	// recognition barcode
	if (barcodes.size() != 0){
		BarcodeStr = barcodes[0].sBarCodeData;
		// delete blank
		int pos = BarcodeStr.find_last_of(" ");
		BarcodeStr.erase(BarcodeStr.begin() + pos);

		Barcode detectedBarcode = barcodes[0];
		return 0;
	}


	// detect seal
	double scale = 0.0;
	vector<int> index;
	//EllipseDetector(InputImg, ellipses, index, TemTopEllipse, scale, savePath, fold);
	int ellipseFlag = EllipseDetector(InputImg, ellipses, index, TemTopEllipse, scale);
	if (ellipseFlag < 0){
		return ellipseFlag;
	}
	struct Ellipse TopEllipse = ellipses[index[0]];

	//rotate img 
	//get dst center point 
	Mat M, ImgRotate;
	double degree = TopEllipse._rad - TemTopEllipse._rad;
	QczVision::rotateImg(InputImg, ImgRotate, M, degree);

	double axis[3];
	axis[0] = TopEllipse._xc; axis[1] = TopEllipse._yc;
	axis[2] = 1;
	Mat srcAxis(3, 1, CV_64F, axis);
	Mat dstAxis = M * srcAxis;
	Point2f sealCenterAfterRotate(dstAxis.at<double>(0, 0), dstAxis.at<double>(1, 0));

	Rect TemRectFusion = BarcodeStrRect | barcode.rcBarcodeRegion;
	//Rect TemRectFusion = BarcodeStrRect;

	Point2f TemRectCentral(TemRectFusion.x + TemRectFusion.width / 2, TemRectFusion.y + TemRectFusion.height / 2);
	Point2f TemSealCentral(TemTopEllipse._xc, TemTopEllipse._yc);
	double xDis = TemRectCentral.x - TemSealCentral.x;
	double yDis = TemRectCentral.y - TemSealCentral.y;
	Point2f RectCenterAfterRotate(sealCenterAfterRotate.x + scale*xDis, sealCenterAfterRotate.y + scale*yDis);
	Point2f lt(RectCenterAfterRotate.x - TemRectFusion.width*scale / 2,RectCenterAfterRotate.y - TemRectFusion.height*scale/2);
	int RotateRectX = lt.x - TemRectFusion.width * scale / 1.5;
	Rect RotateRect(RotateRectX, lt.y, TemRectFusion.width * scale * 1.6, TemRectFusion.height * scale * 1.2);

	//Mat ImgRotate;
	//warpAffine(InputImg, ImgRotate, M, Size(width_rotate,height_rotate), CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, BORDER_CONSTANT,cvScalarAll(0));
	//rectangle(ImgRotate, RotateRect, Scalar(0, 0, 255), 3);
	
	// threshold  + find character
	Rect Total(0, 0, ImgRotate.cols, ImgRotate.rows);
	RotateRect = RotateRect & Total;

	Mat wordPart;
	ImgRotate(RotateRect).copyTo(wordPart);
	if (wordPart.channels() == 3){
		cvtColor(wordPart, wordPart, CV_BGR2GRAY);
	}
	wordPart.convertTo(wordPart,CV_8UC1);

	Mat binWordPart;
	//threshold(wordPart, binWordPart, 0, 255, CV_THRESH_BINARY_INV | CV_THRESH_OTSU);
	int blocksize = 25;
	int constValue = 10;
	//Mat local;
	adaptiveThreshold(wordPart, binWordPart, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY_INV, blocksize, constValue);

	////////////////////////////////////////////////////////////////////////////
	// find contours 
	Mat binResult; binWordPart.copyTo(binResult);
	vector<Mat> contours(10000);
	Mat hierarchy;
	cv::findContours(binWordPart, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);

	vector<Rect> rects;
	for (int i = 0; i < contours.size(); i++){
		Rect tmp = boundingRect(contours[i]);
		rects.push_back(tmp);
	}
	////////////////////////////////////////////////////////////////////////////
	//  remove noises
	double ratioThreshold = 3;
	for (int i = rects.size() - 1; i >= 0; i--){
		double ratio = (double)max(rects[i].width, rects[i].height) / min(rects[i].width, rects[i].height);
		double areaRatio = (double)rects[i].area() / (binResult.cols * binResult.rows);
		if (ratio > ratioThreshold || areaRatio < 0.0005 || areaRatio > 0.02){
			rects.erase(rects.begin() + i);
		}
	}
	
	Mat binThreshold;
	cvtColor(binResult, binThreshold, CV_GRAY2BGR);
	for (int i = 0; i < rects.size(); i++){
		rectangle(binThreshold, rects[i], Scalar(0, 0, 255), 1);
	}
	///////////////////////////////////////////////////////////////////////////
	// remove rects contained
	sort(rects.begin(), rects.end(), QczSort::compareRectArea);
	for (int i = 0; i < rects.size() - 1; i++){
		for (int j = i + 1; j < rects.size(); j++){
			if ((rects[i] & rects[j]) == rects[j]){
				rects.erase(rects.begin() + j);
			}
		}
	}

	// remove
	for (int i = rects.size() - 1; i >= 0; i--){
		double ratio = (double)max(rects[i].width, rects[i].height) / min(rects[i].width, rects[i].height);
		double areaRatio = (double)rects[i].area() / (binResult.cols * binResult.rows);
		double h = rects[i].height;
		if (ratio > ratioThreshold || areaRatio < 0.0005 || areaRatio > 0.02 || h > binResult.rows / 10){
			rects.erase(rects.begin() + i);
		}
	}
	///////////////////////////////////////////////////////////////////////
	// merge horizontal
	for (int i = 0; i < rects.size() - 1; i++){
		Rect ri = rects[i];
		for (int j = i + 1; j < rects.size(); j++){
			Point ci(ri.x + ri.width / 2, ri.y + ri.height / 2);
			Point cj(rects[j].x + rects[j].width / 2, rects[j].y + rects[j].height / 2);
			if (abs(ci.y - cj.y) <= (ri.height / 2 + rects[j].height / 2) && abs(ci.x - cj.x) - (ri.width / 2 + rects[j].width / 2) <= (min(ri.height, rects[j].height)) / 5){
				ri = ri | rects[j];
				rects[i] = ri;
				rects.erase(rects.begin() + j);
				j--;
			}
		}
	}
	for (int i = rects.size() - 1; i >= 0; i--){
		double ratio = (double)max(rects[i].width, rects[i].height) / min(rects[i].width, rects[i].height);
		double areaRatio = (double)rects[i].area() / (binResult.cols * binResult.rows);
		double h = rects[i].height;
		if (ratio > ratioThreshold || areaRatio < 0.001 || areaRatio > 0.02 || h > binResult.rows / 10){
			rects.erase(rects.begin() + i);
		}
	}
	///////////////////////////////////////////////////////////////////////
	// project for cluster
	sort(rects.begin(), rects.end(), QczSort::compareRectWidth);
	vector<bool> cover(rects.size());
	for (int i = 0; i < cover.size(); i++){ cover[i] = true; }
	vector<vector<int> > cluster;
	for (int i = 0; i < rects.size() - 1; i++){
		if (cover[i] == false){ continue; }
		vector<int> ICover; 
		ICover.push_back(i);
		cover[i] = false;
		for (int j = i + 1; j < rects.size(); j++){
			Rect TopBottom(rects[i].x, 0, rects[i].width, binResult.rows);
			Rect tmp = TopBottom & rects[j];
			if (tmp.area() != 0){
				ICover.push_back(j);
				cover[j] = false;
			}
		}
		cluster.push_back(ICover);
	}
	// get the word cluster
	int cluster_index = 0;
	double cluster_variance = DBL_MAX;
	for (int i = 0; i < cluster.size(); i++){
		if (cluster[i].size() < 9){ continue; }
		else{
			vector<double> areas(cluster[i].size());
			double meanArea = 0.0;
			for (int j = 0; j < cluster[i].size(); j++){
				//areas[j] = rects[cluster[i][j]].area(); // adopt area
				areas[j] = rects[cluster[i][j]].width;	// adopt width
				//areas[j] = rects[cluster[i][j]].x;	// adopt x
				meanArea += areas[j];
			}
			meanArea /= (double)cluster[i].size();
			double variance = 0.0;
			for (int j = 0; j < areas.size(); j++){
				variance += (meanArea - areas[j])*(meanArea - areas[j]);
			}
			variance /= (double)(areas.size() * meanArea);
			if (variance <= cluster_variance){
				cluster_variance = variance;
				cluster_index = i;
			}
		}
	}

	// fusion per cluster
	vector<Rect> clusterRect(cluster.size());
	for (int i = 0; i < cluster.size(); i++){
		Rect tmp = rects[cluster[i][0]];
		for (int j = 0; j < cluster[i].size(); j++){
			tmp |= rects[cluster[i][j]];
		}
		clusterRect[i] = tmp;
	}
	///////////////////////////////////////////////////////////
	// remove noise rect and keep the top 11
	vector<int> wordCluster = cluster[cluster_index];

	if (wordCluster.size() > 11){
		double meanWidth = 0;
		for (int i = 0; i < wordCluster.size(); i++){
			meanWidth += rects[wordCluster[i]].width;
		}
		meanWidth = meanWidth / wordCluster.size();
		vector<pair<double, int>> disIndex;
		for (int i = 0; i < wordCluster.size(); i++){
			disIndex.push_back(make_pair(abs(meanWidth - rects[wordCluster[i]].width), wordCluster[i]));
		}
		sort(disIndex.begin(), disIndex.end(), QczSort::compareIndexD<double,int>);
		wordCluster.clear();
		for (int i = 0; i < 11; i++){
			wordCluster.push_back(disIndex[i].second);
		}
	}
	///////////////////////////////////////////////////////////
	// draw the rotateRect
	vector<Point> corner;
	for (int i = 0; i < wordCluster.size(); i++){
		Rect tmp = rects[wordCluster[i]];
		corner.push_back(Point(tmp.x, tmp.y));
		corner.push_back(Point(tmp.x, tmp.y + tmp.height -1));
		corner.push_back(Point(tmp.x + tmp.width -1, tmp.y));
		corner.push_back(Point(tmp.x + tmp.width -1, tmp.y + tmp.height - 1));
	}

	RotatedRect rRect = minAreaRect(corner);

	///////////////////////////////////////////////////////////
	double angle_roRect = fabs(rRect.angle) > 45 ? (180 + rRect.angle) : (90 + rRect.angle);
	double degree_roRect = angle_roRect * CV_PI / 180;
	//double cos_theta_roRect = cos(degree_roRect); double sin_theta_roRect = sin(degree_roRect);
	Mat roMat;
	wordPart.copyTo(roMat);
	Mat M_roRect, ImgRotate_roRect;
	QczVision::rotateImg(roMat, ImgRotate_roRect, M_roRect, degree_roRect);
	int roWidth = roMat.cols, roHeight = roMat.rows;

	Point2f vertices[4];
	rRect.points(vertices);
	Mat srcWordAxis = Mat::zeros(3, 4, CV_64F);
	for (int i = 0; i < 4; i++){
		srcWordAxis.at<double>(0, i) = vertices[i].x;
		srcWordAxis.at<double>(1, i) = vertices[i].y;
		srcWordAxis.at<double>(2, i) = 1;
	}
	Mat dstWordAxis = M_roRect * srcWordAxis;
	double minX = 0.0, maxX = 0.0, minY = 0.0, maxY = 0.0;
	Point minXP = 0, maxXP = 0, minYP = 0, maxYP = 0;
	minMaxLoc(dstWordAxis(Range(0, 1), Range::all()), &minX, &maxX, &minXP, &maxXP);
	minMaxLoc(dstWordAxis(Range(1, 2), Range::all()), &minY, &maxY, &minYP, &maxYP);

	Mat recogWordMat;
	Rect recogWordRect((int)minX, (int)minY, (int)(maxX - minX + 1), (int)(maxY - minY + 1));
	ImgRotate_roRect(recogWordRect).copyTo(recogWordMat);

	///////////////////////////////////////////////////////////
	// recongize the word
	Mat recogWordBin, recogWordMean;
	int meanValue = mean(recogWordMat)[0];
	threshold(recogWordMat, recogWordBin, meanValue, 1, CV_THRESH_BINARY_INV);
	//threshold(recogWordMat, recogWordBin, 0, 255, CV_THRESH_BINARY_INV | CV_THRESH_OTSU);

	Mat segWordPos;
	reduce(recogWordBin, segWordPos, 0, CV_REDUCE_SUM,CV_32S);
	vector<int> segPos;
	bool flag = false;
	for(int i = 0; i < segWordPos.cols; i++){
		if (i == 0 && segWordPos.at<int>(0, i) > 0){
			segPos.push_back(i);
			flag = !flag;
		}
		if ((segWordPos.at<int>(0, i) > 0) != flag){
			segPos.push_back(i);
			flag = !flag;
		}
	}
	if ((int)segPos.size() % 2 != 0){
		segPos.push_back(segWordPos.cols - 1);
	}
	vector<Mat> characters;
	for (int i = 0; i < segPos.size(); i += 2){
		Mat tmp;
		recogWordBin(Range::all(), Range(segPos[i], segPos[i + 1] + 1)).copyTo(tmp);
		characters.push_back(tmp);
	}
	//for (int i = 0; i < characters.size(); i++){
	//	string name = QczStr::int2string(i) + ".jpg";
	//	imwrite(name, characters[i]);
	//}
	string recogResult;
	for (int i = 0; i < characters.size(); i++){
		int normSize = 64;
		int col = characters[i].cols, row = characters[i].rows;
		Mat tmp = Mat::zeros(normSize, normSize*col / row, IPL_DEPTH_8U);
		resize(characters[i], tmp, tmp.size());
		//tmp = tmp / 255;
		CharacterArray array;
		array.pArrayData = (unsigned char*)malloc(tmp.cols*tmp.rows*sizeof(unsigned char) * 2);
		array.width = tmp.cols;
		array.height = tmp.rows;

		for (int j = 0; j<array.height; j++)
		{
			for (int k = 0; k<array.width; k++)
			{
				array.pArrayData[j*array.width + k] = tmp.data[j*tmp.step + k];
			}
		}

		unsigned short LanguageSetOption = 7;
		int Recog_Candidate_Num = 1;
		char* InnerCode = (char *)calloc(Recog_Candidate_Num, 2 * sizeof(char));
		int *Dis = (int *)calloc(Recog_Candidate_Num, sizeof(int));
		int num = Recognize(&array, InnerCode, Dis, Recog_Candidate_Num, LanguageSetOption);
		recogResult.push_back(InnerCode[1]);
		free(array.pArrayData);
		free(InnerCode);
		free(Dis);
	}
	if (recogResult.length() != 11){
		return -2;
	}
	BarcodeStr = recogResult;

	//////////////////////////////////////////////////////////
	return 0;
}