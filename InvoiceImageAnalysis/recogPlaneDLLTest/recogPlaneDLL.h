#ifndef _RECOGBARCODEDLL_H
#define _RECOGBARCODEDLL_H

#include <stdio.h>

#define RECOGBAR_DLL _declspec(dllexport)


#ifdef __cplusplus
extern "C" {
#endif
	/*****************************************
	// ʶ��������
	����:
	PlanePath		:	�ɻ�ƱͼƬ��ַ
	UserName		:	�û�����
	PlaneFlag		:	�ɻ�Ʊ��־λ
								1	:	���ڻ�Ʊ
								0	:	��������Ʊ

	���:
	BarcodeStr		:	ʶ����
	��������ֵ		:
								n	:	��Ҫʶ�����Ŀ
								0	:	ʶ��ɹ�
								-1	:	δ��⵽����
								-2	:	ʶ������ʽ����
								-3	:	ͼƬ��ʧ��
								-4	:	ʶ���������ʧ��
								-5	:	ģ�����ʧ��
								-10	:	������Ϣ����
								-11	:	��ʶ����Ϣ
								-12	:	�����쳣
								-13	:	�������
	*****************************************/
	RECOGBAR_DLL int recogPlaneDLL(char* PlanePath, char* UserName, int PlaneFlag);

#ifdef __cplusplus
}
#endif


#endif