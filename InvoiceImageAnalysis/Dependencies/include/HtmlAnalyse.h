#include <iostream>  
#include <fstream> 
#include <string>
#include <cstring>
#if defined(_MSC_VER)
#define strcasecmp _stricmp
#endif
//////////////////////////////////////////////////////////////////////////////////////////////
//去掉字符串中的占位符和乱码
//void trim(string &s);
//对HTML文件进行解析并用字符串输出
//filename为HTML文件路径
//查验单号：GP170731221178
//印刷序号：65357262003
//ml27 : TAO / JIANHUA
//from : 北京
//from : STOARLANDAAPT
//from : 北京
//from : VOID
//flight : 国航
//flight : 国航
//flight : VOID
//flight_no : CA911
//flight_no : CA912
//class :V
//class :Z
//date : 08 - 20
//date : 08 - 24
//time : 13 : 50
//time : 19 : 10
//cny - total : CNY17493.0
//e - ticket - no : 9995735369492
//insurance : XXX
/////////////////////////////////////////////////////////////////////////
int HtmlAnalyse(std::string fileName);
/////////////////////////////////////////////////////////////////////////
//返回结果:(-10):输出："姓名信息错误"
//////////(-11):输出："无识别信息"
//////////(n)：正常，输出上列识别条目数量