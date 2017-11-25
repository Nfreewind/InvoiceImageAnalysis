#ifndef _RECOGBARCODEDLL_H
#define _RECOGBARCODEDLL_H

#include <stdio.h>

#define RECOGBAR_DLL _declspec(dllexport)


#ifdef __cplusplus
extern "C" {
#endif
	/*****************************************
	// 识别条形码
	输入:
	PlanePath		:	飞机票图片地址
	UserName		:	用户姓名
	PlaneFlag		:	飞机票标志位
								1	:	国内机票
								0	:	国外或境外机票

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
	RECOGBAR_DLL int recogPlaneDLL(char* PlanePath, char* UserName, int PlaneFlag);

#ifdef __cplusplus
}
#endif


#endif