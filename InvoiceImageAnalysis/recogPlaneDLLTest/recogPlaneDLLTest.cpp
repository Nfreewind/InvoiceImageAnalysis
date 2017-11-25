#include <stdio.h>
#include <string>
#include <iostream>
#include "recogPlaneDLL.h"

using namespace std;

void main(int argc,char**argv){

	string PlanePath = argv[1];
	string UserName = argv[2];
	int PlaneFlag = atoi(argv[3]);

	int recogPlaneFlag = recogPlaneDLL((char*)PlanePath.c_str(), (char*)UserName.c_str(), PlaneFlag);
	cout << recogPlaneFlag << endl;
}