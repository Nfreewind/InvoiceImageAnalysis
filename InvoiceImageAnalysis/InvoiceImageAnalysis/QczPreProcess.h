#pragma once

#include <stdio.h>

#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

//struct MatchPoints{
//	Point TopSeal;
//	Point BottomSeal;
//	Point Barcode;
//	Point BarcodeStr;
//
//	MatchPoints() : TopSeal(Point(0, 0)), BottomSeal(Point(0, 0)), Barcode(Point(0, 0)), BarcodeStr(Point(0, 0)) {};
//};


/*
Input:
	TemplatePath	: template for barcodes, barcodeStr and ellipse.
	im				: image read from disk
	savePath		: default
	fold			: default
Output:
	out				: output image							not implement
	rect			: output barcodeStr rect				not implement
	BarcodeStr		: output barcode recognition result
Return:		int
	-1				: ellipse not detected
	0				: success
	-2				: recognization output are not 11 chars
*/

int PreProcess(Mat im, string TemplatePath, Mat& out, Rect& rect, string& BarcodeStr, string savePath = "", string fold = "");