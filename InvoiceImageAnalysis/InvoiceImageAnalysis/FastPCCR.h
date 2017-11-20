#ifndef _FAST_PCCR_H
#define _FAST_PCCR_H

#include <Common.h>
#ifdef __cplusplus
extern "C"
{
#endif

int  InitializeRecognizer(char* sDicPath);
int  Recognize(CharacterArray *pArray, char *InnerCode, int *Dis, int CandidateNum, unsigned short LanguageSetOption);
void CloseRecognizer();
void ReturnFeature(unsigned char *pFeature);

#ifdef __cplusplus
}
#endif

#endif 