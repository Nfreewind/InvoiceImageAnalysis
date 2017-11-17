#pragma once

#include <stdio.h>

#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

int PreProcess(Mat im, string TemplatePath, Mat& out, Rect& rect, string& BarcodeStr, string savePath = "", string fold = "");