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

int PreProcess(Mat im, string TemplatePath, Mat& out, Rect& rect, string& BarcodeStr, string savePath, string fold){

	int step = 0;

	vector<struct Ellipse> ellipses;
	int idx = 0;
	//	load template
	struct Ellipse TemTopEllipse, TemBottomEllipse;
	Barcode barcode;
	Rect BarcodeStrRect;
	TemplateLoader(TemplatePath, TemTopEllipse, TemBottomEllipse, barcode, BarcodeStrRect);


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

	//if (savePath != ""){
	//	imwrite(savePath, InputImg);
	//}

	// detect seal
	double scale = 0.0;
	EllipseDetector(InputImg, ellipses, idx, TemTopEllipse, scale);

	//detect barcode

	vector<Barcode> barcodes;
	CBarcodeReader barcodeReader;
	HANDLE hBarcode;
	InitBarcodeReader(barcodeReader, hBarcode);


	BarcodeRead(InputImg, barcodeReader, hBarcode, barcodes);

	// recognition barcode
	if (barcodes.size() != 0){
		BarcodeStr = barcodes[0].sBarCodeData;
		//return 0;
	}
	else{

	}

	//rotate img 
	struct Ellipse TopEllipse = ellipses[idx];
	int width = InputImg.cols, height = InputImg.rows;

	Point2f center;
	center.x = width / 2.0 + 0.5;
	center.y = height / 2.0 + 0.5;

	double degree = TopEllipse._rad - TemTopEllipse._rad;
	double angle = degree * 180.0 / CV_PI;
	double sin_theta = sin(degree), cos_theta = cos(degree);

	////////////////////////////////////////////////////////////////////////////
	// tag the five points
	/*
	//for (int i = 0; i < 5; i++){
	//	char num[100];
	//	sprintf(num, "%d", i);
	//	string num_s(num);
	//	putText(InputImg, num_s, points[i], FONT_HERSHEY_PLAIN, 3, Scalar(0, 0, 255), 3);
	//}
	//if (savePath != ""){
	//	imwrite(savePath, InputImg);
	//}
	*/
	////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	// get affine matrix to detect the word
	/*
	//double degree = TopEllipse._rad - 0.0174539089;
	//double angle = degree * 180.0 / CV_PI;
	//double sin_theta = sin(degree), cos_theta = cos(degree);

	//Point2f points[5];
	//int x = TopEllipse._xc, y = TopEllipse._yc;
	//int l_len = TopEllipse._a, s_len = TopEllipse._b;
	//points[0] = Point2f(x + s_len * sin_theta, y + s_len*cos_theta);
	//points[1] = Point2f(x - l_len * cos_theta, y - l_len*sin_theta);
	//points[2] = Point2f(x - s_len * sin_theta, y - s_len*cos_theta);
	//points[3] = Point2f(x + l_len * cos_theta, y + l_len*sin_theta);
	//points[4] = Point2f(x, y);
	//Point2f TemPoints[5];
	//int x_t = TemTopEllipse._xc, y_t = TemTopEllipse._yc;
	//int l_t_len = TemTopEllipse._a, s_t_len = TemTopEllipse._b;
	//double degree_t = TemTopEllipse._rad - 0.0174539089;
	//double angle_t = degree_t * 180.0 / CV_PI;
	//double sin_t_theta = sin(degree_t), cos_t_theta = cos(degree_t);
	//TemPoints[0] = Point2f(x_t + s_t_len * sin_t_theta, y_t + s_t_len*cos_t_theta);
	//TemPoints[1] = Point2f(x_t - l_t_len * cos_t_theta, y_t - l_t_len*sin_t_theta);
	//TemPoints[2] = Point2f(x_t - s_t_len * sin_t_theta, y_t - s_t_len*cos_t_theta);
	//TemPoints[3] = Point2f(x_t + l_t_len * cos_t_theta, y_t + l_t_len*sin_t_theta);
	//TemPoints[4] = Point2f(x_t, y_t);

	//Mat M1 = getAffineTransform(TemPoints, points);

	//cout << M1.at<double>(0, 0) << endl;
	//cout << M1.at<double>(0, 1) << endl;
	//cout << M1.at<double>(0, 2) << endl;
	//cout << M1.at<double>(1, 0) << endl;
	//cout << M1.at<double>(1, 1) << endl;
	//cout << M1.at<double>(1, 2) << endl;

	//Point2f TemStrPoints[4];
	//x = BarcodeStrRect.x, y = BarcodeStrRect.y;
	//int w = BarcodeStrRect.width, h = BarcodeStrRect.height;

	//double axis[12];
	//axis[0] = x; axis[1] = x + w - 1; axis[2] = x; axis[3] = x + w - 1;
	//axis[4] = y; axis[5] = y; axis[6] = y + h - 1; axis[7] = y + h - 1;
	//axis[8] = 1; axis[9] = 1; axis[10] = 1; axis[11] = 1;
	//Mat TemStrPointsMat(3, 4, CV_64F, axis);
	//Mat StrPoints = M1 * TemStrPointsMat;
	//Point2f pointLT(StrPoints.at<double>(0, 0), StrPoints.at<double>(1, 0));
	//Point2f pointRT(StrPoints.at<double>(0, 1), StrPoints.at<double>(1, 1));
	//Point2f pointLB(StrPoints.at<double>(0, 2), StrPoints.at<double>(1, 2));
	//Point2f pointRB(StrPoints.at<double>(0, 3), StrPoints.at<double>(1, 3));

	//rectangle(InputImg, pointLT, pointRB, Scalar(0, 0, 255), 3);
	//if (savePath != ""){
	//	imwrite(savePath, InputImg);
	//}
	*/
	//////////////////////////////////////////////////////////////////////////

	int width_rotate = width * fabs(cos_theta) + height * fabs(sin_theta);
	int height_rotate = width * fabs(sin_theta) + height *fabs(cos_theta);

	float m[6];
	Mat M(2, 3, CV_32F, m);
	M = getRotationMatrix2D(center, angle, 1);
	//cout << M.at<double>(0, 0) << endl;
	//cout << M.at<double>(0, 1) << endl;
	//cout << M.at<double>(0, 2) << endl;
	//cout << M.at<double>(1, 0) << endl;
	//cout << M.at<double>(1, 1) << endl;
	//cout << M.at<double>(1, 2) << endl;
	M.at<double>(0, 2) += (width_rotate - width) / 2;
	M.at<double>(1, 2) += (height_rotate - height) / 2;


	//get dst center point 
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

	Mat ImgRotate;
	warpAffine(InputImg, ImgRotate, M, Size(width_rotate,height_rotate), CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, BORDER_CONSTANT,cvScalarAll(0));
	rectangle(ImgRotate, RotateRect, Scalar(0, 0, 255), 3);
	
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

	//////////////////////////////////////////////////////////////////////////
	// erode + dilate
	//Mat erodeImg, dilateImg;
	//int erodeSize = 1;
	//int dilateSize = 1;
	//Mat eEllement = getStructuringElement(MORPH_RECT, Size(2 * erodeSize + 1, 2 * erodeSize + 1), Point(erodeSize, erodeSize));
	//Mat dEllement = getStructuringElement(MORPH_RECT, Size(2 * dilateSize + 1, 2 * dilateSize + 1), Point(dilateSize, dilateSize));
	//dilate(binWordPart, dilateImg, dEllement);
	//erode(binWordPart, erodeImg, eEllement);
	//if (savePath != ""){
	//	string name = QczFile::splitFileName(savePath);
	//	string subName, ext;
	//	QczFile::splitExt(name, subName, ext);
	//	string path = fold + "\\" + subName + "_dilate" + ext;
	//	imwrite(path, dilateImg);
	//}
	/////////////////////////////////////////////////////////////////////////

	if (savePath != ""){
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_src" + ext;
		imwrite(path, wordPart);
		step++;
	}

	if (savePath != ""){
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_bin" + ext;
		imwrite(path, binWordPart);
		step++;
	}

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
	if (savePath != ""){
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_bin_threshold" + ext;
		step++;
		imwrite(path, binThreshold);
	}

	///////////////////////////////////////////////////////////////////////////
	// merge rects with intersect
	//for (int i = 0; i < rects.size() - 1; i++){
	//	for (int j = i + 1; j < rects.size(); j++){
	//		Rect tmp = rects[i] & rects[j];
	//		if (tmp.area() != 0){
	//			rects[j] = rects[i] | rects[j];
	//			rects.erase(rects.begin() + i);
	//			i--;
	//			break;
	//		}
	//	}
	//}

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

	if (savePath != ""){
		Mat binResultColor;
		cvtColor(binResult, binResultColor, CV_GRAY2BGR);
		for (int i = 0; i < rects.size(); i++){
			rectangle(binResultColor, rects[i], Scalar(0, 0, 255), 1);
		}
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_bin_remove_contain" + ext;
		step++;
		imwrite(path, binResultColor);
	}
	///////////////////////////////////////////////////////////////////////

	// fit line
	/*
	std::vector<cv::Point> points;
	for (int i = 0; i < rects.size(); i++){
		Point tmp(rects[i].x + rects[i].width / 2, rects[i].y + rects[i].height / 2);
		points.push_back(tmp);
	}
	cv::Vec4f line;
	fitLine(points, line, CV_DIST_HUBER, 0, 0.01, 0.01);
	double cos_theta_l = line[0], sin_theta_l = line[1], x0 = line[2], y0 = line[3];
	double phi = atan2(sin_theta_l, cos_theta_l) + CV_PI / 2.0;
	double rho = y0*cos_theta_l - x0*sin_theta_l;
	double k = sin_theta_l / cos_theta_l;
	double b = y0 - k*x0;
	double x = 0;
	double y = k*x + b;
	cv::line(binResultColor, Point(x0, y0), Point(x, y), Scalar(0, 0, 255), 2);

	if (savePath != ""){
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_line" + ext;
		imwrite(path, binResultColor);
	}
	*/

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

	if (savePath != ""){
		Mat binResultColor;
		cvtColor(binResult, binResultColor, CV_GRAY2BGR);
		for (int i = 0; i < rects.size(); i++){
			rectangle(binResultColor, rects[i], Scalar(0, 0, 255), 1);
		}
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_bin_merge_horizontal" + ext;
		step++;
		imwrite(path, binResultColor);
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

	if (savePath != ""){
		Mat binClusterColor;
		cvtColor(binResult, binClusterColor, CV_GRAY2BGR);
		for (int i = 0; i < clusterRect.size(); i++){
			rectangle(binClusterColor, clusterRect[i], Scalar(0, 0, 255), 2);
		}
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_cluster" + ext;
		step++;
		imwrite(path, binClusterColor);
	}

	///////////////////////////////////////////////////////////
	// remove noise rect and keep the top 11
	vector<int> wordCluster = cluster[cluster_index];

	if (savePath != ""){
		Mat binClusterColor;
		cvtColor(binResult, binClusterColor, CV_GRAY2BGR);
		for (int i = 0; i < wordCluster.size(); i++){
			rectangle(binClusterColor, rects[wordCluster[i]], Scalar(0, 0, 255), 2);
		}
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_cluster_sub_rect" + ext;
		step++;
		imwrite(path, binClusterColor);
	}



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

	if (savePath != ""){
		Mat binClusterColor;
		cvtColor(binResult, binClusterColor, CV_GRAY2BGR);
		for (int i = 0; i < wordCluster.size(); i++){
			rectangle(binClusterColor, rects[wordCluster[i]], Scalar(0, 0, 255), 2);
		}
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_cluster_sub_rect_after" + ext;
		step++;
		imwrite(path, binClusterColor);
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
	if (savePath != ""){
		Mat binRotateRect;
		cvtColor(binResult, binRotateRect, CV_GRAY2BGR);
		Point2f vertices[4];
		rRect.points(vertices);
		for (int i = 0; i < 4; i++){
			line(binRotateRect, vertices[i], vertices[(i + 1) % 4], Scalar(0, 0, 255));
		}
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_rotate_rect" + ext;
		step++;
		imwrite(path, binRotateRect);
	}
	double angle_roRect = fabs(rRect.angle) > 45 ? (180 + rRect.angle) : (90 + rRect.angle);
	double degree_roRect = angle_roRect * CV_PI / 180;
	double cos_theta_roRect = cos(degree_roRect); double sin_theta_roRect = sin(degree_roRect);
	Mat roMat;
	wordPart.copyTo(roMat);
	int roWidth = roMat.cols, roHeight = roMat.rows;

	int width_after_rotate = roWidth * fabs(cos_theta_roRect) + roHeight * fabs(sin_theta_roRect);
	int height_after_rotate = roWidth * fabs(sin_theta_roRect) + roHeight *fabs(cos_theta_roRect);
	float m_roRect[6];
	Mat M_roRect(2, 3, CV_32F, m_roRect);
	Point2f center_roRect;
	center_roRect.x = roWidth / 2.0 + 0.5;
	center_roRect.y = roHeight / 2.0 + 0.5;
	M_roRect = getRotationMatrix2D(center_roRect, angle_roRect, 1);
	M_roRect.at<double>(0, 2) += (width_after_rotate - roWidth) / 2;
	M_roRect.at<double>(1, 2) += (height_after_rotate - roHeight) / 2;
	Mat ImgRotate_roRect;
	warpAffine(wordPart, ImgRotate_roRect, M_roRect, Size(width_after_rotate, height_after_rotate), CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, BORDER_CONSTANT, cvScalarAll(0));

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


	if (savePath != ""){
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_word_for_recog" + ext;
		step++;
		imwrite(path, recogWordMat);
	}
	///////////////////////////////////////////////////////////
	// recongize the word
	Mat recogWordBin;
	threshold(recogWordMat, recogWordBin, 0, 255, CV_THRESH_BINARY_INV | CV_THRESH_OTSU);
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
		//imshow("bin", tmp);
		//waitKey(0);
	}
	string recogResult;
	for (int i = 0; i < characters.size(); i++){
		int normSize = 64;
		int col = characters[i].cols, row = characters[i].rows;
		Mat tmp = Mat::zeros(normSize, normSize*col / row, IPL_DEPTH_8U);
		resize(characters[i], tmp, tmp.size());
		tmp = tmp / 255;
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
	}
	if (savePath != ""){
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_word_for_recog_bin" + ext;
		step++;
		imwrite(path, recogWordBin);
	}
	if (savePath != ""){
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_word_for_recog_bin_string.txt";
		step++;
		ofstream outfile(path.c_str());
		outfile << recogResult << endl;
		outfile.close();
	}


	///////////////////////////////////////////////////////////
	Rect wordRect;

	for (int i = 0; i < wordCluster.size(); i++){
		if (i == 0){ wordRect = rects[wordCluster[0]]; }
		else{ wordRect |= rects[wordCluster[i]]; }
	}

	if (savePath != ""){
		Mat binClusterWordColor;
		cvtColor(binResult, binClusterWordColor, CV_GRAY2BGR);
		rectangle(binClusterWordColor, wordRect, Scalar(0, 0, 255), 2);
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_word" + ext;
		step++;
		imwrite(path, binClusterWordColor);
	}

	// re binary and character recognition
	Mat wordMat;
	wordPart(wordRect).copyTo(wordMat);
	// rotate 90
	transpose(wordMat, wordMat);
	flip(wordMat, wordMat, 0);
	if (savePath != ""){
		string name = QczFile::splitFileName(savePath);
		string subName, ext;
		QczFile::splitExt(name, subName, ext);
		string path = fold + "\\" + subName + "_" + QczStr::int2string(step) + "_word_flip" + ext;
		step++;
		imwrite(path, wordMat);
	}
	// re threshold


	return 0;
}