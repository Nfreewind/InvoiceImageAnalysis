#include <iostream>  
#include <fstream> 
#include <string>
#include <cstring>
#if defined(_MSC_VER)
#define strcasecmp _stricmp
#endif
//////////////////////////////////////////////////////////////////////////////////////////////
//ȥ���ַ����е�ռλ��������
//void trim(string &s);
//��HTML�ļ����н��������ַ������
//filenameΪHTML�ļ�·��
//���鵥�ţ�GP170731221178
//ӡˢ��ţ�65357262003
//ml27 : TAO / JIANHUA
//from : ����
//from : STOARLANDAAPT
//from : ����
//from : VOID
//flight : ����
//flight : ����
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
//���ؽ��:(-10):�����"������Ϣ����"
//////////(-11):�����"��ʶ����Ϣ"
//////////(n)���������������ʶ����Ŀ����