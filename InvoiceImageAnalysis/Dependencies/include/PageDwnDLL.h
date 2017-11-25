#ifndef _PAGEDWN_DLL_H
#define _PAGEDWN_DLL_H

#include <stdio.h>
#include <string>

#define PAGEDWN_DLL _declspec(dllimport)


#ifdef __cplusplus
extern "C" {
#endif
	/*****************************************
	// 下载网页
	输入:
		Barcode		:	条形码的识别结果
		UserName	:	飞机票用户姓名  中文
		Flag		:	标志位	
								1	:	国内机票
								0	:	国外或境外机票
	输出:
		PagePath	:	下载网页 保存地址
		函数返回值	:	
						0	:	下载成功
						1	:	网络异常
						2	:	输入错误
	*****************************************/
	PAGEDWN_DLL int PageDwn(std::string Barcode, std::string UserName, int Flag, std::string& PagePath);


#ifdef __cplusplus
}
#endif

#endif