#include <stdio.h>
#include <string>
#include <iostream>

#include "recogBarcodeDLL.h"

using namespace std;

void main(){
	string PlanePath = "img00001.jpg";
	string BarcodeStr;

	int flag = recogBarcode(PlanePath, BarcodeStr);
	
	cout << flag << endl;
	cout << BarcodeStr << endl;
	system("pause");
}