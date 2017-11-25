#ifndef _PAGEDWN_DLL_H
#define _PAGEDWN_DLL_H

#include <stdio.h>
#include <string>

#define PAGEDWN_DLL _declspec(dllimport)


#ifdef __cplusplus
extern "C" {
#endif
	/*****************************************
	// ������ҳ
	����:
		Barcode		:	�������ʶ����
		UserName	:	�ɻ�Ʊ�û�����  ����
		Flag		:	��־λ	
								1	:	���ڻ�Ʊ
								0	:	��������Ʊ
	���:
		PagePath	:	������ҳ �����ַ
		��������ֵ	:	
						0	:	���سɹ�
						1	:	�����쳣
						2	:	�������
	*****************************************/
	PAGEDWN_DLL int PageDwn(std::string Barcode, std::string UserName, int Flag, std::string& PagePath);


#ifdef __cplusplus
}
#endif

#endif