#include "stdafx.h"
#include "AntigerColor.h"

using namespace std;

static std::vector<cv::Vec3f> wei()
{
	std::vector<cv::Vec3f> W;

	for (int i = 0; i <= 10; ++i)
		for (int j = 0; j <= (10 - i); ++j)
		{
			int k = 10 - (i + j);
			W.push_back(cv::Vec3f(i / 10.f, j / 10.f, k / 10.f));
		}

	return W;
}

static void add(const cv::Vec3b &c0, const cv::Vec3b &c1, std::vector<cv::Vec3f> &P, std::vector<float> &det)
{
	float Rx = c0[2] / 255.f - c1[2] / 255.f;
	float Gx = c0[1] / 255.f - c1[1] / 255.f;
	float Bx = c0[0] / 255.f - c1[0] / 255.f;

	float d = sqrt(Rx * Rx + Gx * Gx + Bx * Bx) / 1.41f;

	if (d >= 0.05f)
	{
		P.push_back(cv::Vec3f(Rx, Gx, Bx));
		det.push_back(d);
	}
}

void Decolorization(const cv::Mat &im, cv::Mat &out)
{
	float s = 64.f / sqrt((float)(im.rows * im.cols));
	int cols = (s * im.cols) + 0.5f;
	int rows = (s * im.rows) + 0.5f;

	const float sigma = 0.05f;
	std::vector<cv::Vec3f> W = wei();

	std::vector<cv::Vec3f> P;
	std::vector<float> det;

	std::vector<cv::Point> pos0;
	for (int i = 0; i < cols; ++i)
	{
		for (int j = 0; j < rows; ++j)
		{
			int x = (i + 0.5f) * im.cols / cols;
			int y = (j + 0.5f) * im.rows / rows;
			pos0.push_back(cv::Point(x, y));
		}
	}

	std::vector<cv::Point> pos1 = pos0;
	std::random_shuffle(pos1.begin(), pos1.end());

	for (std::size_t i = 0; i < pos0.size(); ++i)
	{
		cv::Vec3b c0 = im.at<cv::Vec3b>(pos0[i].y, pos0[i].x);
		cv::Vec3b c1 = im.at<cv::Vec3b>(pos1[i].y, pos1[i].x);

		add(c0, c1, P, det);
	}

	cols /= 2;
	rows /= 2;

	for (int i = 0; i < cols - 1; ++i)
	{
		for (int j = 0; j < rows; ++j)
		{
			int x0 = (i + 0.5f) * im.cols / cols;
			int x1 = (i + 1.5f) * im.cols / cols;
			int y = (j + 0.5f) * im.rows / rows;

			cv::Vec3b c0 = im.at<cv::Vec3b>(y, x0);
			cv::Vec3b c1 = im.at<cv::Vec3b>(y, x1);

			add(c0, c1, P, det);
		}
	}

	for (int i = 0; i < cols; ++i)
	{
		for (int j = 0; j < rows - 1; ++j)
		{
			int x = (i + 0.5f) * im.cols / cols;
			int y0 = (j + 0.5f) * im.rows / rows;
			int y1 = (j + 1.5f) * im.rows / rows;

			cv::Vec3b c0 = im.at<cv::Vec3b>(y0, x);
			cv::Vec3b c1 = im.at<cv::Vec3b>(y1, x);

			add(c0, c1, P, det);
		}
	}

	float maxEs = -FLT_MAX;
	int bw;
	for (std::size_t i = 0; i < W.size(); ++i)
	{
		float Es = 0;
		for (std::size_t j = 0; j < P.size(); ++j)
		{
			float L = P[j][0] * W[i][0] + P[j][1] * W[i][1] + P[j][2] * W[i][2];
			float detM = det[j];

			float a = (L + detM) / sigma;
			float b = (L - detM) / sigma;

			Es += log(exp(-a * a) + exp(-b * b));
		}
		Es /= P.size();

		if (Es > maxEs)
		{
			maxEs = Es;
			bw = i;
		}
	}

	std::vector<cv::Mat> c;
	cv::split(im, c);

	out.release();
	cv::addWeighted(c[2], W[bw][0], c[1], W[bw][1], 0.0, out);
	cv::addWeighted(c[0], W[bw][2], out, 1.0, 0.0, out);

}

void FastDecolorization(const cv::Mat &im, cv::Mat &out)
{
	cv::cvtColor(im, out, CV_BGR2GRAY);
}

void GetRednessMap(const Mat& BgrImage, Mat& RedImage)
{
	Mat BgrImg;
	Size size;
	int i, j;
	unsigned char *rLine;
	cv::Vec3b* cLine;

	BgrImg = BgrImage;

	size = BgrImg.size();

	RedImage.create(size, CV_8UC1);

	if (BgrImg.isContinuous())
	{
		size.width *= size.height;
		size.height = 1;
	}

	for (i = 0; i < size.height; i++)
	{
		rLine = RedImage.ptr<unsigned char>(i);
		cLine = BgrImg.ptr<cv::Vec3b>(i);
		for (j = 0; j < size.width; j++)
		{
			int r,g,b,mingb;
			r = cLine[j][2];
			g = cLine[j][1];
			b = cLine[j][0];
			rLine[j] = (r > g) ? ((r - g)) : 0;
		}
	}
}

cv::Mat gammaTransform(cv::Mat& srcImage, float kFactor)
{
	// 建立查表文件LUT
	unsigned char LUT[256];
	for (int i = 0; i < 256; i++)
	{
		// Gamma变换表达式
		LUT[i] = saturate_cast<uchar>(pow((float)(
			i / 255.0), kFactor) * 255.0f);
	}
	cv::Mat resultImage = srcImage.clone();
	// 输入通道为单通道时 直接进行变换
	if (srcImage.channels() == 1)
	{
		cv::MatIterator_<uchar> iterator =
			resultImage.begin<uchar>();
		cv::MatIterator_<uchar> iteratorEnd =
			resultImage.end<uchar>();
		for (; iterator != iteratorEnd; iterator++)
			*iterator = LUT[(*iterator)];
	}
	else
	{
		// 输入通道为三通道时 需对每个通道分别进行变换
		cv::MatIterator_<cv::Vec3b> iterator =
			resultImage.begin<Vec3b>();
		cv::MatIterator_<cv::Vec3b> iteratorEnd =
			resultImage.end<Vec3b>();
		//  通过查找表进行转换
		for (; iterator != iteratorEnd; iterator++)
		{
			(*iterator)[0] = LUT[((*iterator)[0])];
			(*iterator)[1] = LUT[((*iterator)[1])];
			(*iterator)[2] = LUT[((*iterator)[2])];
		}
	}
	return resultImage;
}
