#include <stdio.h>
#include <string>
#include <iostream>


#include "recogPlaneDLL.h"
#include "HtmlAnalyse.h"
#include "PageDwnDll.h"
#include "recogBarcodeDLL.h"

using namespace std;

int recogPlaneDLL(char* PlanePath, char* UserName, int PlaneFlag){




	string BarcodeStr;
	int recogBarcodeFlag = recogBarcode(string(PlanePath), BarcodeStr);
	if (recogBarcodeFlag != 0){
		return recogBarcodeFlag;
	}
	cout << BarcodeStr << endl;

	string PagePath;
	int pageDwnFlag = PageDwn(BarcodeStr, string(UserName), PlaneFlag, PagePath);
	if (pageDwnFlag != 0){
		return pageDwnFlag;
	}

	int htmlAnalyseFlag = HtmlAnalyse(PagePath);
	
	return htmlAnalyseFlag;

}