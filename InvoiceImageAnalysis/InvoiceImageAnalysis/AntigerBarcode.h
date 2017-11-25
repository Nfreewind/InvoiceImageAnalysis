#pragma once
#include <cv.h>
#include <opencv2/opencv.hpp>
#include <highgui.h>
//#include "DynamsoftBarcodeReader.h"
#include <windows.h>
#include "SoftekBarcodeDLL.h"

using namespace cv;

typedef HANDLE  HDIB;

typedef struct {
	int             iType;              // recognized barcode type
	char			sBarCodeData[1000]; // recognized barcode string, maxlen = 999
	int             iAngle;             // recognized barcode direction
	Rect            rcBarcodeRegion;    // recognized barcode location on the bitmap
} Barcode;

#define WIDTHBYTES(width, bitCount)    (((width)*(bitCount) + 31) / 32 * 4)

HDIB Mat2Dib(Mat &Img);

void InitBarcodeReader(HANDLE &hBarcode);
void BarcodeRead(Mat &image, HANDLE &hBarcode, vector<Barcode>& barcodes);