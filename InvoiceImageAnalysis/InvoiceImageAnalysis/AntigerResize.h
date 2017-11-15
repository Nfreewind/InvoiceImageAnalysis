#pragma once

#include <cv.h>
#include <opencv2/opencv.hpp>
#include <highgui.h>

using namespace cv;

double ImgResizeWithFixedLongside(Mat &inputImage, Mat &outputImage, int longSize = 320);
int ImgResizeWithFixedShortside(Mat &inputImage, Mat &outputImage, int shortSize = 240);
int ImgResizeWithFixedSize(Mat &inputImage, Mat &outputImage, int iSize = 76800); //76800=320*240