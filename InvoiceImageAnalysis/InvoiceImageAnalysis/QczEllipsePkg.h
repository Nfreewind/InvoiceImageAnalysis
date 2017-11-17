#pragma once

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include "AntigerEllipse.h"

using namespace std;
using namespace cv;

/*
Output: 
	0	: normal
	-1	: none seal detected
Input:
	InputImg	: read from disk
	ellipses	: output detected Ellipses
	Template	: ellipse for template
	scale		: output scale for word matching / scale = ellipse / template 
	saveName	: for debug save the result to disk
*/


int EllipseDetector(Mat InputImg, vector<struct Ellipse>& ellipses, int& idx, struct Ellipse Template, double& scale, string savePath = "");
