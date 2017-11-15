#pragma once

#include <cv.h>
#include <opencv2/opencv.hpp>
#include <highgui.h>

using namespace cv;
/********************* Defines *********************/

#ifndef CV_SQR
#  define CV_SQR(x)  ((x)*(x))
#endif

#ifndef CV_CUBE
#  define CV_CUBE(x)  ((x)*(x)*(x))
#endif

#ifndef CV_INIT_VECTOR
#  define CV_INIT_VECTOR(vname, type, ...) \
    static const type vname##_a[] = __VA_ARGS__; \
    std::vector <type> vname(vname##_a, \
    vname##_a + sizeof(vname##_a) / sizeof(*vname##_a))
#endif
/*! fictitious type to highlight that function
*  can process n-channels arguments */
typedef cv::Mat NChannelsMat;

/*! random forest used to detect edges */
struct RandomForest
{
	/*! random forest options, e.g. number of trees */
	struct RandomForestOptions
	{
		// model params

		int numberOfOutputChannels; /*!< number of edge orientation bins for output */

		int patchSize;              /*!< width of image patches */
		int patchInnerSize;         /*!< width of predicted part inside patch*/

									// feature params

		int regFeatureSmoothingRadius;    /*!< radius for smoothing of regular features
										  *   (using convolution with triangle filter) */

		int ssFeatureSmoothingRadius;     /*!< radius for smoothing of additional features
										  *   (using convolution with triangle filter) */

		int shrinkNumber;                 /*!< amount to shrink channels */

		int numberOfGradientOrientations; /*!< number of orientations per gradient scale */

		int gradientSmoothingRadius;      /*!< radius for smoothing of gradients
										  *   (using convolution with triangle filter) */

		int gradientNormalizationRadius;  /*!< gradient normalization radius */
		int selfsimilarityGridSize;       /*!< number of self similarity cells */

										  // detection params
		int numberOfTrees;            /*!< number of trees in forest to train */
		int numberOfTreesToEvaluate;  /*!< number of trees to evaluate per location */

		int stride;                   /*!< stride at which to compute edges */

	} options;

	int numberOfTreeNodes;

	std::vector <int> featureIds;     /*!< feature coordinate thresholded at k-th node */
	std::vector <float> thresholds;   /*!< threshold applied to featureIds[k] at k-th node */
	std::vector <int> childs;         /*!< k --> child[k] - 1, child[k] */

	std::vector <int> edgeBoundaries; /*!< ... */
	std::vector <int> edgeBins;       /*!< ... */
};

void CreateRandomForest(RandomForest& __rf, const std::string &filename);
void RandomForestEdges(const cv::Mat &src, CV_OUT cv::Mat &dst, RandomForest& __rf);
void DeleteRandomForest(RandomForest& __rf);

//Use canny algorithm to binarize and thin the RandomForestEdges
void CannyEdges(const Mat &srcImg, const Mat &rfImg, CV_OUT Mat &edgeImg);
