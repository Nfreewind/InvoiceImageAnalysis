#pragma once
#include <cv.h>
#include <opencv2/opencv.hpp>
#include <highgui.h>

using namespace cv;

void BgEstimation(const Mat &InputImage, Mat &BgImage, int iBlockSize = 3);
void FastBgEstimation(const Mat &InputImage, Mat &BgImage, int iBlockSize = 15);
void BgCompensation(const Mat &InImg, Mat &BgImg, Mat &EnhancedImg);

int  getHistGram(const Mat& InputGrayImage, int *HistGram, int &total, InputArray mask = noArray());
void getOtsuThreshold(int *HistGram, int total, int* pThreshold);
void getEntropyThreshold(int *HistGram, int total, int* pThreshold);
void getFuzzyThreshold(int *HistGram, int total, int* pThreshold);
void getCorrelationThreshold(int *HistGram, int total, int* pThreshold);
void getNiblackThreshold(const Mat& InputGrayImage, int* pThreshold, bool bInvertFlag = 0);
void Nick_BinarizeImage(Mat &InputImage, Mat &OutputBinarizedImage, int iBlockSize);

void DeSpeckles(Mat &BinarizedImage, int aperture_size = 3, int threshold = 3);

