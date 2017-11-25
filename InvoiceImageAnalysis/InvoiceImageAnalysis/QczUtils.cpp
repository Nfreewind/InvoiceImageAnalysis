#include "stdafx.h"
#include "QczUtils.h"
#include <sstream>

namespace QczFile{
	void getFiles_(string path, vector<string>& files)
	{
		//文件句柄
		intptr_t   hFile = 0;
		//文件信息
		struct _finddata_t fileinfo;
		string p;
		if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
		{
			do
			{
				//如果是目录,迭代之
				//如果不是,加入列表
				if ((fileinfo.attrib &  _A_SUBDIR))
				{
					if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
						getFiles_(p.assign(path).append("\\").append(fileinfo.name), files);
				}
				else
				{
					files.push_back(p.assign(path).append("\\").append(fileinfo.name));
				}
			} while (_findnext(hFile, &fileinfo) == 0);
			_findclose(hFile);
		}
	}
	string splitFileName(string path) {
		size_t start = path.find_last_of("/\\") + 1;
		//size_t end = file.find_last_of(".");
		return path.substr(start);
	}
	void splitExt(string name, string& subName, string& ext){
		int pos = name.find_last_of(".");
		subName = name.substr(0, pos);
		ext = name.substr(pos, name.length());
	}
}

namespace QczStr{
	// recursive replace the substr
	string&   replace_all(string&   str, const   string&   old_value, const   string&   new_value)
	{
		while (true)   {
			string::size_type   pos(0);
			if ((pos = str.find(old_value)) != string::npos)
				str.replace(pos, old_value.length(), new_value);
			else   break;
		}
		return   str;
	}
	// replace the substr once
	string&   replace_all_distinct(string&   str, const   string&   old_value, const   string&   new_value)
	{
		for (string::size_type pos(0); pos != string::npos; pos += new_value.length())   {
			if ((pos = str.find(old_value, pos)) != string::npos)
				str.replace(pos, old_value.length(), new_value);
			else   break;
		}
		return   str;
	}
	// int to string
	string int2string(int input){
		stringstream ss;
		string str;
		ss << input;
		ss >> str;
		return str;
	}
}
namespace QczVision{

	void rotateImg(cv::Mat src, cv::Mat& dst, cv::Mat& M, double degree){

		int width = src.cols, height = src.rows;

		cv::Point2f center;
		center.x = width / 2.0 + 0.5;
		center.y = height / 2.0 + 0.5;

		double angle = degree * 180.0 / CV_PI;
		double sin_theta = sin(degree), cos_theta = cos(degree);

		int width_rotate = width * fabs(cos_theta) + height * fabs(sin_theta);
		int height_rotate = width * fabs(sin_theta) + height *fabs(cos_theta);

		float m[6];
		cv::Mat M_tmp(2, 3, CV_32F, m);
		M_tmp = getRotationMatrix2D(center, angle, 1);
		M_tmp.at<double>(0, 2) += (width_rotate - width) / 2;
		M_tmp.at<double>(1, 2) += (height_rotate - height) / 2;

		warpAffine(src, dst, M_tmp, cv::Size(width_rotate, height_rotate), CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, cv::BORDER_CONSTANT, cvScalarAll(0));
		M_tmp.copyTo(M);
	}

}