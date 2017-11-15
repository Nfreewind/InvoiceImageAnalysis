#include "stdafx.h"
#include "AntigerEdge.h"
#include "AntigerBin.h"

using namespace std;

/*!
* Lightweight wrapper over cv::resize
*
* \param src : source image to resize
* \param dst : destination image size
* \return resized image
*/
static cv::Mat imresize(const cv::Mat &src, const cv::Size &nSize)
{
	cv::Mat dst;
	if (nSize.width < src.size().width
		&&  nSize.height < src.size().height)
		cv::resize(src, dst, nSize, 0.0, 0.0, cv::INTER_AREA);
	else
		cv::resize(src, dst, nSize, 0.0, 0.0, cv::INTER_LINEAR);

	return dst;
}

/*!
* The function filters src with triangle filter with radius equal rad
*
* \param src : source image to filter
* \param rad : radius of filtering kernel
* \return filtering result
*/
static cv::Mat imsmooth(const cv::Mat &src, const int rad)
{
	if (rad == 0)
		return src;
	else
	{
		const float p = 12.0f / rad / (rad + 2) - 2;
		cv::Mat dst;

		if (rad <= 1)
		{
			CV_INIT_VECTOR(kernelXY, float, { 1 / (p + 2), p / (p + 2), 1 / (p + 2) });
			cv::sepFilter2D(src, dst, -1, kernelXY, kernelXY);
		}
		else
		{
			float nrml = CV_SQR(rad + 1.0f);

			std::vector <float> kernelXY(2 * rad + 1);
			for (int i = 0; i <= rad; ++i)
			{
				kernelXY[2 * rad - i] = (i + 1) / nrml;
				kernelXY[i] = (i + 1) / nrml;
			}
			sepFilter2D(src, dst, -1, kernelXY, kernelXY);
		}

		return dst;
	}
}

/*!
*  The function implements rgb to luv conversion in a way similar
*  to UCSD computer vision toolbox
*
* \param src : source image (RGB, float, in [0;1]) to convert
* \return converted image in luv colorspace
*/
static cv::Mat rgb2luv(const cv::Mat &src)
{
	cv::Mat dst(src.size(), src.type());

	const float a = CV_CUBE(29.0f) / 27;
	const float y0 = 8.0f / a;

	const float mX[] = { 0.430574f, 0.341550f, 0.178325f };
	const float mY[] = { 0.222015f, 0.706655f, 0.071330f };
	const float mZ[] = { 0.020183f, 0.129553f, 0.939180f };

	const float maxi = 1.0f / 270;
	const float minu = -88 * maxi;
	const float minv = -134 * maxi;

	const float un = 0.197833f;
	const float vn = 0.468331f;

	// build (padded) lookup table for y->l conversion assuming y in [0,1]
	std::vector <float> lTable(1024);
	for (int i = 0; i < 1024; ++i)
	{
		float y = i / 1024.0f;
		float l = y > y0 ? 116 * powf(y, 1.0f / 3.0f) - 16 : y*a;

		lTable[i] = l*maxi;
	}
	for (int i = 0; i < 40; ++i)
		lTable.push_back(*--lTable.end());

	const int nchannels = 3;

	for (int i = 0; i < src.rows; ++i)
	{
		const float *pSrc = src.ptr<float>(i);
		float *pDst = dst.ptr<float>(i);

		for (int j = 0; j < src.cols*nchannels; j += nchannels)
		{
			const float rgb[] = { pSrc[j + 0], pSrc[j + 1], pSrc[j + 2] };

			const float xyz[] = { mX[0] * rgb[0] + mX[1] * rgb[1] + mX[2] * rgb[2],
				mY[0] * rgb[0] + mY[1] * rgb[1] + mY[2] * rgb[2],
				mZ[0] * rgb[0] + mZ[1] * rgb[1] + mZ[2] * rgb[2] };
			const float nz = 1.0f / float(xyz[0] + 15 * xyz[1] + 3 * xyz[2] + 1e-35);

			const float l = pDst[j] = lTable[cvFloor(1024 * xyz[1])];

			pDst[j + 1] = l * (13 * 4 * xyz[0] * nz - 13 * un) - minu;;
			pDst[j + 2] = l * (13 * 9 * xyz[1] * nz - 13 * vn) - minv;
		}
	}

	return dst;
}

/*!
* The function computes gradient magnitude and weighted (with magnitude)
* orientation histogram. Magnitude is additionally normalized
* by dividing on imsmooth(M, gnrmRad) + 0.01;
*
* \param src : source image
* \param magnitude : gradient magnitude
* \param histogram : gradient orientation nBins-channels histogram
* \param nBins : number of gradient orientations
* \param pSize : factor to downscale histogram
* \param gnrmRad : radius for magnitude normalization
*/
static void gradientHist(const cv::Mat &src, cv::Mat &magnitude, cv::Mat &histogram,
	const int nBins, const int pSize, const int gnrmRad)
{
	cv::Mat phase, Dx, Dy;

	magnitude.create(src.size(), cv::DataType<float>::type);
	phase.create(src.size(), cv::DataType<float>::type);
	histogram.create(cv::Size(cvCeil(src.size().width / float(pSize)),
		cvCeil(src.size().height / float(pSize))),
		CV_MAKETYPE(cv::DataType<float>::type, nBins));

	histogram.setTo(0);

	cv::Sobel(src, Dx, cv::DataType<float>::type,
		1, 0, 1, 1.0, 0.0, cv::BORDER_REFLECT);
	cv::Sobel(src, Dy, cv::DataType<float>::type,
		0, 1, 1, 1.0, 0.0, cv::BORDER_REFLECT);

	int nchannels = src.channels();

	for (int i = 0; i < src.rows; ++i)
	{
		const float *pDx = Dx.ptr<float>(i);
		const float *pDy = Dy.ptr<float>(i);

		float *pMagnitude = magnitude.ptr<float>(i);
		float *pPhase = phase.ptr<float>(i);

		for (int j = 0; j < src.cols*nchannels; j += nchannels)
		{
			float fMagn = float(-1e-5), fdx = 0, fdy = 0;
			for (int k = 0; k < nchannels; ++k)
			{
				float cMagn = CV_SQR(pDx[j + k]) + CV_SQR(pDy[j + k]);
				if (cMagn > fMagn)
				{
					fMagn = cMagn;
					fdx = pDx[j + k];
					fdy = pDy[j + k];
				}
			}

			pMagnitude[j / nchannels] = sqrtf(fMagn);

			float angle = cv::fastAtan2(fdy, fdx) / 180.0f - 1.0f * (fdy < 0);
			if (std::fabs(fdx) + std::fabs(fdy) < 1e-5)
				angle = 0.5f;
			pPhase[j / nchannels] = angle;
		}
	}

	magnitude /= imsmooth(magnitude, gnrmRad)
		+ 0.01*cv::Mat::ones(magnitude.size(), magnitude.type());

	int pHistSize = histogram.cols*histogram.channels() - 1;
	for (int i = 0; i < phase.rows; ++i)
	{
		const float *pPhase = phase.ptr<float>(i);
		const float *pMagn = magnitude.ptr<float>(i);

		float *pHist = histogram.ptr<float>(i / pSize);

		for (int j = 0; j < phase.cols; ++j)
		{
			int index = cvRound((j / pSize + pPhase[j])*nBins);
			index = std::max(0, std::min(index, pHistSize));
			pHist[index] += pMagn[j] / CV_SQR(pSize);
		}
	}
}
/*!
* The method extracts features from img and store them to features.
* Extracted features are appropriate for StructuredEdgeDetection::predictEdges.
*
* \param src : source image (RGB, float, in [0;1]) to extract features
* \param features : destination feature image
*
* \param gnrmRad : __rf.options.gradientNormalizationRadius
* \param gsmthRad : __rf.options.gradientSmoothingRadius
* \param shrink : __rf.options.shrinkNumber
* \param outNum : __rf.options.numberOfOutputChannels
* \param gradNum : __rf.options.numberOfGradientOrientations
*/
void getFeatures(const Mat &src, Mat &features, const int gnrmRad, const int gsmthRad, const int shrink, const int outNum, const int gradNum)
{
	cv::Mat luvImg = rgb2luv(src);

	std::vector <cv::Mat> featureArray;

	cv::Size nSize = Size(src.size().width / float(shrink), src.size().height / float(shrink));

	split(imresize(luvImg, nSize), featureArray);

	CV_INIT_VECTOR(scales, float, { 1.0f, 0.5f });

	for (size_t i = 0; i < scales.size(); ++i)
	{
		int pSize = std::max(1, int(shrink*scales[i]));

		cv::Mat magnitude, histogram;
		gradientHist(/**/ imsmooth(imresize(luvImg, Size(src.size().width * scales[i], src.size().height * scales[i])), gsmthRad),
			magnitude, histogram, gradNum, pSize, gnrmRad /**/);

		featureArray.push_back(/**/ imresize(magnitude, nSize).clone() /**/);
		featureArray.push_back(/**/ imresize(histogram, nSize).clone() /**/);
	}

	// Mixing
	int resType = CV_MAKETYPE(cv::DataType<float>::type, outNum);
	features.create(nSize, resType);

	std::vector <int> fromTo;
	for (int i = 0; i < 2 * outNum; ++i)
		fromTo.push_back(i / 2);

	mixChannels(featureArray, features, fromTo);
}

/*!
* This constructor loads __rf model from filename
*
* \param filename : name of the file where the model is stored
*/
void CreateRandomForest(RandomForest& __rf, const std::string &filename)
{
	cv::FileStorage modelFile(filename, FileStorage::READ);
	CV_Assert(modelFile.isOpened());

	__rf.options.stride
		= modelFile["options"]["stride"];
	__rf.options.shrinkNumber
		= modelFile["options"]["shrinkNumber"];
	__rf.options.patchSize
		= modelFile["options"]["patchSize"];
	__rf.options.patchInnerSize
		= modelFile["options"]["patchInnerSize"];

	__rf.options.numberOfGradientOrientations
		= modelFile["options"]["numberOfGradientOrientations"];
	__rf.options.gradientSmoothingRadius
		= modelFile["options"]["gradientSmoothingRadius"];
	__rf.options.regFeatureSmoothingRadius
		= modelFile["options"]["regFeatureSmoothingRadius"];
	__rf.options.ssFeatureSmoothingRadius
		= modelFile["options"]["ssFeatureSmoothingRadius"];
	__rf.options.gradientNormalizationRadius
		= modelFile["options"]["gradientNormalizationRadius"];

	__rf.options.selfsimilarityGridSize
		= modelFile["options"]["selfsimilarityGridSize"];

	__rf.options.numberOfTrees
		= modelFile["options"]["numberOfTrees"];
	__rf.options.numberOfTreesToEvaluate
		= modelFile["options"]["numberOfTreesToEvaluate"];

	__rf.options.numberOfOutputChannels =
		2 * (__rf.options.numberOfGradientOrientations + 1) + 3;
	//--------------------------------------------

	cv::FileNode childs = modelFile["childs"];
	cv::FileNode featureIds = modelFile["featureIds"];

	std::vector <int> currentTree;

	for (cv::FileNodeIterator it = childs.begin();
		it != childs.end(); ++it)
	{
		(*it) >> currentTree;
		std::copy(currentTree.begin(), currentTree.end(),
			std::back_inserter(__rf.childs));
	}

	for (cv::FileNodeIterator it = featureIds.begin();
		it != featureIds.end(); ++it)
	{
		(*it) >> currentTree;
		std::copy(currentTree.begin(), currentTree.end(),
			std::back_inserter(__rf.featureIds));
	}

	cv::FileNode thresholds = modelFile["thresholds"];
	std::vector <float> fcurrentTree;

	for (cv::FileNodeIterator it = thresholds.begin();
		it != thresholds.end(); ++it)
	{
		(*it) >> fcurrentTree;
		std::copy(fcurrentTree.begin(), fcurrentTree.end(),
			std::back_inserter(__rf.thresholds));
	}

	cv::FileNode edgeBoundaries = modelFile["edgeBoundaries"];
	cv::FileNode edgeBins = modelFile["edgeBins"];

	for (cv::FileNodeIterator it = edgeBoundaries.begin();
		it != edgeBoundaries.end(); ++it)
	{
		(*it) >> currentTree;
		std::copy(currentTree.begin(), currentTree.end(),
			std::back_inserter(__rf.edgeBoundaries));
	}

	for (cv::FileNodeIterator it = edgeBins.begin();
		it != edgeBins.end(); ++it)
	{
		(*it) >> currentTree;
		std::copy(currentTree.begin(), currentTree.end(),
			std::back_inserter(__rf.edgeBins));
	}

	__rf.numberOfTreeNodes = int(__rf.childs.size()) / __rf.options.numberOfTrees;
}

/*!
* Private method used by process method. The function
* predict edges in n-channel feature image and store them to dst.
*
* \param features : source image (n-channels, float) to detect edges
* \param dst : destination image (grayscale, float, in [0;1]) where edges are drawn
*/
void predictEdges(const NChannelsMat &features, cv::Mat &dst, RandomForest& __rf)
{
	int shrink = __rf.options.shrinkNumber;
	int rfs = __rf.options.regFeatureSmoothingRadius;
	int sfs = __rf.options.ssFeatureSmoothingRadius;

	int nTreesEval = __rf.options.numberOfTreesToEvaluate;
	int nTrees = __rf.options.numberOfTrees;
	int nTreesNodes = __rf.numberOfTreeNodes;

	const int nchannels = features.channels();
	int pSize = __rf.options.patchSize;

	int nFeatures = CV_SQR(pSize / shrink)*nchannels;
	int outNum = __rf.options.numberOfOutputChannels;

	int stride = __rf.options.stride;
	int ipSize = __rf.options.patchInnerSize;
	int gridSize = __rf.options.selfsimilarityGridSize;

	const int height = cvCeil(double(features.rows*shrink - pSize) / stride);
	const int width = cvCeil(double(features.cols*shrink - pSize) / stride);
	// image size in patches with overlapping

	//-------------------------------------------------------------------------

	NChannelsMat regFeatures = imsmooth(features, cvRound(rfs / float(shrink)));
	NChannelsMat  ssFeatures = imsmooth(features, cvRound(sfs / float(shrink)));

	NChannelsMat indexes(height, width, CV_MAKETYPE(DataType<int>::type, nTreesEval));

	std::vector <int> offsetI(/**/ CV_SQR(pSize / shrink)*nchannels, 0);
	for (int i = 0; i < CV_SQR(pSize / shrink)*nchannels; ++i)
	{
		int z = i / CV_SQR(pSize / shrink);
		int y = (i % CV_SQR(pSize / shrink)) / (pSize / shrink);
		int x = (i % CV_SQR(pSize / shrink)) % (pSize / shrink);

		offsetI[i] = x*features.cols*nchannels + y*nchannels + z;
	}
	// lookup table for mapping linear index to offsets

	std::vector <int> offsetE(/**/ CV_SQR(ipSize)*outNum, 0);
	for (int i = 0; i < CV_SQR(ipSize)*outNum; ++i)
	{
		int z = i / CV_SQR(ipSize);
		int y = (i % CV_SQR(ipSize)) / ipSize;
		int x = (i % CV_SQR(ipSize)) % ipSize;

		offsetE[i] = x*dst.cols*outNum + y*outNum + z;
	}
	// lookup table for mapping linear index to offsets

	std::vector <int> offsetX(CV_SQR(gridSize)*(CV_SQR(gridSize) - 1) / 2 * nchannels, 0);
	std::vector <int> offsetY(CV_SQR(gridSize)*(CV_SQR(gridSize) - 1) / 2 * nchannels, 0);

	int hc = cvRound((pSize / shrink) / (2.0*gridSize));
	// half of cell
	std::vector <int> gridPositions;
	for (int i = 0; i < gridSize; i++)
		gridPositions.push_back(int((i + 1)*(pSize / shrink + 2 * hc - 1) / (gridSize + 1.0) - hc + 0.5f));

	for (int i = 0, n = 0; i < CV_SQR(gridSize)*nchannels; ++i)
		for (int j = (i%CV_SQR(gridSize)) + 1; j < CV_SQR(gridSize); ++j, ++n)
		{
			int z = i / CV_SQR(gridSize);

			int x1 = gridPositions[i%CV_SQR(gridSize) % gridSize];
			int y1 = gridPositions[i%CV_SQR(gridSize) / gridSize];

			int x2 = gridPositions[j%gridSize];
			int y2 = gridPositions[j / gridSize];

			offsetX[n] = x1*features.cols*nchannels + y1*nchannels + z;
			offsetY[n] = x2*features.cols*nchannels + y2*nchannels + z;
		}
	// lookup tables for mapping linear index to offset pairs
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < height; ++i)
	{
		float *regFeaturesPtr = regFeatures.ptr<float>(i*stride / shrink);
		float  *ssFeaturesPtr = ssFeatures.ptr<float>(i*stride / shrink);

		int *indexPtr = indexes.ptr<int>(i);

		for (int j = 0, k = 0; j < width; ++k, j += !(k %= nTreesEval))
			// for j,k in [0;width)x[0;nTreesEval)
		{
			int baseNode = (((i + j) % (2 * nTreesEval) + k) % nTrees)*nTreesNodes;
			int currentNode = baseNode;
			// select root node of the tree to evaluate

			int offset = (j*stride / shrink)*nchannels;
			while (__rf.childs[currentNode] != 0)
			{
				int currentId = __rf.featureIds[currentNode];
				float currentFeature;

				if (currentId >= nFeatures)
				{
					int xIndex = offsetX[currentId - nFeatures];
					float A = ssFeaturesPtr[offset + xIndex];

					int yIndex = offsetY[currentId - nFeatures];
					float B = ssFeaturesPtr[offset + yIndex];

					currentFeature = A - B;
				}
				else
					currentFeature = regFeaturesPtr[offset + offsetI[currentId]];

				// compare feature to threshold and move left or right accordingly
				if (currentFeature < __rf.thresholds[currentNode])
					currentNode = baseNode + __rf.childs[currentNode] - 1;
				else
					currentNode = baseNode + __rf.childs[currentNode];
			}

			indexPtr[j*nTreesEval + k] = currentNode;
		}
	}

	NChannelsMat dstM(dst.size(),
		CV_MAKETYPE(DataType<float>::type, outNum));
	dstM.setTo(0);

	float step = 2.0f * CV_SQR(stride) / CV_SQR(ipSize) / nTreesEval;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < height; ++i)
	{
		int *pIndex = indexes.ptr<int>(i);
		float *pDst = dstM.ptr<float>(i*stride);

		for (int j = 0, k = 0; j < width; ++k, j += !(k %= nTreesEval))
		{// for j,k in [0;width)x[0;nTreesEval)

			int currentNode = pIndex[j*nTreesEval + k];

			int start = __rf.edgeBoundaries[currentNode];
			int finish = __rf.edgeBoundaries[currentNode + 1];

			if (start == finish)
				continue;

			int offset = j*stride*outNum;
			for (int p = start; p < finish; ++p)
				pDst[offset + offsetE[__rf.edgeBins[p]]] += step;
		}
	}

	cv::reduce(dstM.reshape(1, int(dstM.total())), dstM, 2, CV_REDUCE_SUM);
	imsmooth(dstM.reshape(1, dst.rows), 1).copyTo(dst);
	dst.convertTo(dst, cv::DataType<unsigned char>::type, 255.0);
}

/*!
* The function detects edges in src and draw them to dst
*
* \param src : source image (RGB, float, in [0;1]) to detect edges
* \param dst : destination image (grayscale, in [0;255])
*              where edges are drawn
*/
void RandomForestEdges(const cv::Mat &src, CV_OUT cv::Mat &dst, RandomForest& __rf)
{
	cv::Mat fsrc;
	cv::Mat rgbsrc;

	cv::cvtColor(src, rgbsrc, CV_BGR2RGB);
	rgbsrc.convertTo(fsrc, cv::DataType<float>::type, 1 / 255.0);

	dst.create(fsrc.size(), cv::DataType<float>::type);

	int padding = (__rf.options.patchSize
		- __rf.options.patchInnerSize) / 2;

	cv::Mat nSrc;
	copyMakeBorder(fsrc, nSrc, padding, padding,
		padding, padding, BORDER_REFLECT);

	NChannelsMat features;
	getFeatures(nSrc, features,
		__rf.options.gradientNormalizationRadius,
		__rf.options.gradientSmoothingRadius,
		__rf.options.shrinkNumber,
		__rf.options.numberOfOutputChannels,
		__rf.options.numberOfGradientOrientations);
	predictEdges(features, dst, __rf);
}

void DeleteRandomForest(RandomForest& __rf)
{
	__rf.featureIds.clear();
	__rf.thresholds.clear();
	__rf.childs.clear();
	__rf.edgeBoundaries.clear();
	__rf.edgeBins.clear();
}

void CannyEdges(const Mat &srcImg, const Mat &rfImg, CV_OUT Mat &edgeImg)
{
	Mat src;
	cv::cvtColor(srcImg, src, CV_BGR2GRAY); 
	CV_Assert(src.depth() == CV_8U);

	edgeImg.create(src.size(), CV_8UC1);
	Mat dst = edgeImg;

	const int cn = src.channels();
	Mat dx(src.rows, src.cols, CV_16SC(cn));
	Mat dy(src.rows, src.cols, CV_16SC(cn));

	Sobel(src, dx, CV_16S, 1, 0, 3, 1, 0, cv::BORDER_REPLICATE);
	Sobel(src, dy, CV_16S, 0, 1, 3, 1, 0, cv::BORDER_REPLICATE);

	Mat mag = rfImg;
	int low,high;
	int HistGram[256];
	int iTotal;
	getHistGram(rfImg, HistGram, iTotal);
	getOtsuThreshold(HistGram, iTotal, &high);
	getFuzzyThreshold(HistGram, iTotal, &low);
	//printf("µÍ %d£¬ ¸ß %d!\n\n", low,high);
	high = high * 4 / 5;
	if (low > high)
		std::swap(low, high);

	ptrdiff_t mapstep = src.cols + 2;
	AutoBuffer<uchar> buffer((src.cols + 2)*(src.rows + 2) + cn * mapstep * 3 * sizeof(int));

	int* mag_buf[3];
	mag_buf[0] = (int*)(uchar*)buffer;
	mag_buf[1] = mag_buf[0] + mapstep*cn;
	mag_buf[2] = mag_buf[1] + mapstep*cn;
	memset(mag_buf[0], 0, /* cn* */mapstep * sizeof(int));

	uchar* map = (uchar*)(mag_buf[2] + mapstep*cn);
	memset(map, 1, mapstep);
	memset(map + mapstep*(src.rows + 1), 1, mapstep);

	int maxsize = std::max(1 << 10, src.cols * src.rows / 10);
	std::vector<uchar*> stack(maxsize);
	uchar **stack_top = &stack[0];
	uchar **stack_bottom = &stack[0];

	/* sector numbers
	(Top-Left Origin)

	1   2   3
	*  *  *
	* * *
	0*******0
	* * *
	*  *  *
	3   2   1
	*/

#define CANNY_PUSH(d)    *(d) = uchar(2), *stack_top++ = (d)
#define CANNY_POP(d)     (d) = *--stack_top

	// calculate magnitude and angle of gradient, perform non-maxima suppression.
	// fill the map with one of the following values:
	//   0 - the pixel might belong to an edge
	//   1 - the pixel can not belong to an edge
	//   2 - the pixel does belong to an edge
	for (int i = 0; i <= src.rows; i++)
	{
		int* _norm = mag_buf[(i > 0) + 1] + 1;
		if (i < src.rows)
		{
			short* _dx = dx.ptr<short>(i);
			short* _dy = dy.ptr<short>(i);
			unsigned char *_rfmag = mag.ptr<unsigned char>(i);
			for (int j = 0; j < src.cols*cn; j++)
				_norm[j] = _rfmag[j];
			_norm[-1] = _norm[src.cols] = 0;
		}
		else
			memset(_norm - 1, 0, /* cn* */mapstep * sizeof(int));

		// at the very beginning we do not have a complete ring
		// buffer of 3 magnitude rows for non-maxima suppression
		if (i == 0)
			continue;

		uchar* _map = map + mapstep*i + 1;
		_map[-1] = _map[src.cols] = 1;

		int* _mag = mag_buf[1] + 1; // take the central row
		ptrdiff_t magstep1 = mag_buf[2] - mag_buf[1];
		ptrdiff_t magstep2 = mag_buf[0] - mag_buf[1];

		const short* _x = dx.ptr<short>(i - 1);
		const short* _y = dy.ptr<short>(i - 1);

		if ((stack_top - stack_bottom) + src.cols > maxsize)
		{
			int sz = (int)(stack_top - stack_bottom);
			maxsize = std::max(sz + src.cols, maxsize * 3 / 2);
			stack.resize(maxsize);
			stack_bottom = &stack[0];
			stack_top = stack_bottom + sz;
		}

		int prev_flag = 0;
		for (int j = 0; j < src.cols; j++)
		{
#define CANNY_SHIFT 15
			const int TG22 = (int)(0.4142135623730950488016887242097*(1 << CANNY_SHIFT) + 0.5);

			int m = _mag[j];

			if (m > low)
			{
				int xs = _x[j];
				int ys = _y[j];
				int x = std::abs(xs);
				int y = std::abs(ys) << CANNY_SHIFT;

				int tg22x = x * TG22;

				if (y < tg22x)
				{
					if (m > _mag[j - 1] && m >= _mag[j + 1]) goto __ocv_canny_push;
				}
				else
				{
					int tg67x = tg22x + (x << (CANNY_SHIFT + 1));
					if (y > tg67x)
					{
						if (m > _mag[j + magstep2] && m >= _mag[j + magstep1]) goto __ocv_canny_push;
					}
					else
					{
						int s = (xs ^ ys) < 0 ? -1 : 1;
						if (m > _mag[j + magstep2 - s] && m > _mag[j + magstep1 + s]) goto __ocv_canny_push;
					}
				}
			}
			prev_flag = 0;
			_map[j] = uchar(1);
			continue;
		__ocv_canny_push:
			if (!prev_flag && m > high && _map[j - mapstep] != 2)
			{
				CANNY_PUSH(_map + j);
				prev_flag = 1;
			}
			else
				_map[j] = 0;
		}

		// scroll the ring buffer
		_mag = mag_buf[0];
		mag_buf[0] = mag_buf[1];
		mag_buf[1] = mag_buf[2];
		mag_buf[2] = _mag;
	}

	// now track the edges (hysteresis thresholding)
	while (stack_top > stack_bottom)
	{
		uchar* m;
		if ((stack_top - stack_bottom) + 8 > maxsize)
		{
			int sz = (int)(stack_top - stack_bottom);
			maxsize = maxsize * 3 / 2;
			stack.resize(maxsize);
			stack_bottom = &stack[0];
			stack_top = stack_bottom + sz;
		}

		CANNY_POP(m);

		if (!m[-1])         CANNY_PUSH(m - 1);
		if (!m[1])          CANNY_PUSH(m + 1);
		if (!m[-mapstep - 1]) CANNY_PUSH(m - mapstep - 1);
		if (!m[-mapstep])   CANNY_PUSH(m - mapstep);
		if (!m[-mapstep + 1]) CANNY_PUSH(m - mapstep + 1);
		if (!m[mapstep - 1])  CANNY_PUSH(m + mapstep - 1);
		if (!m[mapstep])    CANNY_PUSH(m + mapstep);
		if (!m[mapstep + 1])  CANNY_PUSH(m + mapstep + 1);
	}

	// the final pass, form the final image
	const uchar* pmap = map + mapstep + 1;
	uchar* pdst = dst.ptr();
	for (int i = 0; i < src.rows; i++, pmap += mapstep, pdst += dst.step)
	{
		for (int j = 0; j < src.cols; j++)
			pdst[j] = (uchar)-(pmap[j] >> 1);
	}
}