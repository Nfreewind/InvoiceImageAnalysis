#include "stdafx.h"
#include "AntigerResize.h"


double ImgResizeWithFixedLongside(Mat &inputImage, Mat &outputImage, int longSize)
{
	double scaleRatio;
	int inW, inH, outW, outH;
	inW = inputImage.size().width;
	inH = inputImage.size().height;

	if (inW < inH) // height > width
	{
		scaleRatio = inH / (double)longSize;
		if (inH == longSize)
		{
			outputImage = inputImage;
			return scaleRatio;
		}
		outH = longSize;
		outW = cvRound(inW*longSize / (float)inH);
		if (inH < longSize) //need upscaling
		{
			cv::resize(inputImage, outputImage, Size(outW, outH), 0.0, 0.0, cv::INTER_LINEAR);
		}
		else //need downscaling
		{
			cv::resize(inputImage, outputImage, Size(outW, outH), 0.0, 0.0, cv::INTER_AREA);
		}
	}
	else// height <= width
	{
		scaleRatio = inW / (double)longSize;
		if (inW == longSize)
		{
			outputImage = inputImage;
			return scaleRatio;
		}
		outH = cvRound(inH*longSize / (float)inW);
		outW = longSize;
		if (inW < longSize) //need upscaling
		{
			cv::resize(inputImage, outputImage, Size(outW, outH), 0.0, 0.0, cv::INTER_LINEAR);
		}
		else //need downscaling
		{
			cv::resize(inputImage, outputImage, Size(outW, outH), 0.0, 0.0, cv::INTER_AREA);
		}
	}
	return scaleRatio;
}

int ImgResizeWithFixedShortside(Mat &inputImage, Mat &outputImage, int shortSize)
{
	int inW, inH, outW, outH;
	inW = inputImage.size().width;
	inH = inputImage.size().height;

	if (inW > inH) // height < width
	{
		if (inH == shortSize)
		{
			outputImage = inputImage;
			return 0;
		}
		outH = shortSize;
		outW = cvRound(inW*shortSize / (float)inH);
		if (inH < shortSize) //need upscaling
		{
			cv::resize(inputImage, outputImage, Size(outW, outH), 0.0, 0.0, cv::INTER_LINEAR);
		}
		else //need downscaling
		{
			cv::resize(inputImage, outputImage, Size(outW, outH), 0.0, 0.0, cv::INTER_AREA);
		}
	}
	else// height >= width
	{
		if (inW == shortSize)
		{
			outputImage = inputImage;
			return 0;
		}
		outH = cvRound(inH*shortSize / (float)inW);
		outW = shortSize;
		if (inW < shortSize) //need upscaling
		{
			cv::resize(inputImage, outputImage, Size(outW, outH), 0.0, 0.0, cv::INTER_LINEAR);
		}
		else //need downscaling
		{
			cv::resize(inputImage, outputImage, Size(outW, outH), 0.0, 0.0, cv::INTER_AREA);
		}
	}
	return 0;
}

int ImgResizeWithFixedSize(Mat &inputImage, Mat &outputImage, int iSize) //35200=220*160
{
	int inW, inH, outW, outH;
	inW = inputImage.size().width;
	inH = inputImage.size().height;

	int inSize = inW * inH;
	double ratio;
	ratio = sqrt((double)iSize / inSize);

	outH = cvRound(inH*ratio);
	outW = cvRound(inW*ratio);
	if (ratio > 1.0) //need upscaling
	{
		cv::resize(inputImage, outputImage, Size(outW, outH), 0.0, 0.0, cv::INTER_LINEAR);
	}
	else //need downscaling
	{
		cv::resize(inputImage, outputImage, Size(outW, outH), 0.0, 0.0, cv::INTER_AREA);
	}
	return 0;
}