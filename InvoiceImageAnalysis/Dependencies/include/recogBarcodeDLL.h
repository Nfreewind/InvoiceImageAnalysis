#pragma once

#include <stdio.h>
#include <string>

#define RECOGBAR_DLL _declspec(dllimport)


#ifdef __cplusplus
extern "C" {
#endif
	/*****************************************
	// ʶ��������
	����:
	PlanePath		:	�ɻ�ƱͼƬ��ַ

	���:
	BarcodeStr		:	ʶ����
	��������ֵ		:
						0	:	ʶ��ɹ�
						-1	:	δ��⵽����
						-2	:	ʶ������ʽ����
						-3	:	ͼƬ��ʧ��
						-4	:	ʶ���������ʧ��
						-5	:	ģ�����ʧ��
	*****************************************/
	RECOGBAR_DLL int recogBarcode(std::string PlanePath, std::string& BarcodeStr);

#ifdef __cplusplus
}
#endif
