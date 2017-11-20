#include "stdafx.h"
#include "QczUtils.h"
#include <sstream>

namespace QczFile{
	void getFiles_(string path, vector<string>& files)
	{
		//文件句柄
		long   hFile = 0;
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