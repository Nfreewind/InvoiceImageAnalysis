#pragma once

#include <stdio.h>
#include <string>

#define RECOGBAR_DLL _declspec(dllimport)


#ifdef __cplusplus
extern "C" {
#endif
	/*****************************************
	// 识别条形码
	输入:
	PlanePath		:	飞机票图片地址

	输出:
	BarcodeStr		:	识别结果
	函数返回值		:
						0	:	识别成功
						-1	:	未检测到红章
						-2	:	识别结果格式错误
						-3	:	图片打开失败
						-4	:	识别引擎加载失败
						-5	:	模板加载失败
	*****************************************/
	RECOGBAR_DLL int recogBarcode(std::string PlanePath, std::string& BarcodeStr);

#ifdef __cplusplus
}
#endif
