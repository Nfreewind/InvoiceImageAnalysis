#pragma once

#include <stdio.h>

#include <opencv2/opencv.hpp>

#include "AntigerEllipse.h"
#include "AntigerBarcode.h"


using namespace std;
using namespace cv;



void TemplateMaker(Mat im, string TemplatePath);
void TemplateLoader(string TemplatePath, struct Ellipse& TopSeal, struct Ellipse& BottomSeal, Barcode& barcode, Rect& BarcodeStrRect);