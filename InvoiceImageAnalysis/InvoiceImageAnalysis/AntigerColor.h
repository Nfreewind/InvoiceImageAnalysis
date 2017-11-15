#pragma once
#include <cv.h>
#include <opencv2/opencv.hpp>
#include <highgui.h>

using namespace cv;

void Decolorization(const cv::Mat &im, cv::Mat &out);
void FastDecolorization(const cv::Mat &im, cv::Mat &out);
void GetRednessMap(const Mat& BgrImage, Mat& RedImage);