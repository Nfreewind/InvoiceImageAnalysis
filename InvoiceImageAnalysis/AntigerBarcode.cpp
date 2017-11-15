#include "stdafx.h"
#include "AntigerBarcode.h"
#ifdef _WIN64
#pragma comment(lib, "DBRx64.lib")
#pragma comment(lib, "SoftekBarcode64DLL.lib")
#else
#pragma comment(lib, "DBRx86.lib")
#pragma comment(lib, "SoftekBarcodeDLL.lib")
#endif


HDIB Mat2Dib(Mat &Img)
{
	HDIB hImage;
	UINT8* pDib;                        //pointer to the DIB
	BITMAPINFOHEADER *pHead;			//pointer to the bitmap header
	UINT8*         pImageData;          //pointer to the bitmap bits
	RGBQUAD *pPalette;			        //pointer to the color Palette

	UINT32 dSize;
	INT32  colormapEntries, BufWidth, nWidth, nHeight;
	UINT8 *pSource, *pTarget;
	INT32 i;

	hImage = NULL;

	nWidth = Img.cols;
	nHeight = Img.rows;

	if (nWidth <= 0 || nHeight <= 0 || Img.data == NULL || Img.depth() != CV_8U || (Img.channels() != 1 && Img.channels() != 3))
	{
		return NULL;
	}
	//
	if (Img.channels() == 1)
	{
		colormapEntries = 256;
		BufWidth = WIDTHBYTES(nWidth, 8);
	}
	else
	{
		colormapEntries = 0;
		BufWidth = WIDTHBYTES(nWidth, 24);
	}

	dSize = sizeof(BITMAPINFOHEADER) + colormapEntries * sizeof(RGBQUAD) + BufWidth*nHeight;
	hImage = GlobalAlloc(GMEM_MOVEABLE | GMEM_ZEROINIT, dSize);
	pDib = (UINT8*)GlobalLock(hImage);
	pHead = (BITMAPINFOHEADER*)pDib;
	pPalette = (RGBQUAD*)(pDib + sizeof(BITMAPINFOHEADER));
	pImageData = pDib + sizeof(BITMAPINFOHEADER) + colormapEntries * sizeof(RGBQUAD);
	pHead->biSize = sizeof(BITMAPINFOHEADER);
	pHead->biWidth = nWidth;
	pHead->biHeight = nHeight;

	if (Img.channels() == 1)
	{
		pHead->biBitCount = 8;
		pHead->biClrUsed = 256;
		pHead->biSizeImage = nWidth*nHeight;
	}
	else
	{
		pHead->biBitCount = 24;
		pHead->biClrUsed = 0;
		pHead->biSizeImage = nWidth*nHeight * 3;
	}
	pHead->biClrImportant = 0;
	pHead->biCompression = BI_RGB;
	pHead->biPlanes = 1;
	pHead->biXPelsPerMeter = 0;
	pHead->biYPelsPerMeter = 0;

	for (i = 0; i < colormapEntries; i++)
	{
		pPalette[i].rgbBlue = i;
		pPalette[i].rgbGreen = i;
		pPalette[i].rgbRed = i;
		pPalette[i].rgbReserved = 0;
	}
	
	if (Img.channels() == 1)
	{
		pTarget = pImageData + (nHeight - 1)*BufWidth;
		for (i = 0; i < nHeight; i++)
		{
			pSource = Img.ptr<unsigned char>(i);
			memcpy(pTarget, pSource, nWidth);
			pTarget -= BufWidth;
		}
	}
	else
	{
		pTarget = pImageData + (nHeight - 1)*BufWidth;
		for (i = 0; i < nHeight; i++)
		{
			pSource = (unsigned char*)(Img.ptr<cv::Vec3b>(i));
			memcpy(pTarget, pSource, nWidth * 3);
			pTarget -= BufWidth;
		}
	}

	GlobalUnlock(hImage);
	return hImage;
}

void InitBarcodeReader(CBarcodeReader &reader, HANDLE &hBarcode)
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// DynamSoft
	reader.InitLicense("t0068MgAAAKOmEhKZl9J5EwQWUSBJNVAPFURAUnYk/qmSR4RSOwjyQ8nOBegubmtEZWmlqRXcCK4EhWkj6tABCCQhUZpvbqM=");

	//Set Property
	reader.SetBarcodeFormats(BF_CODE_39);
	reader.SetMaxBarcodesNumPerPage(1);
	reader.SetUseOneDDeblur(1);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Softek
	hBarcode = mtCreateBarcodeInstance();
	mtSetLicenseKey(hBarcode, "88H8V1J178ZZD5IKL7ABKN5SK451V2DP");

	// 设定识别条码类型
	mtSetReadCode128(hBarcode, 0);
	mtSetReadCode39(hBarcode, 1);
	mtSetReadCode25(hBarcode, 0);
	mtSetReadEAN13(hBarcode, 0);
	mtSetReadEAN8(hBarcode, 0);
	mtSetReadUPCA(hBarcode, 0);
	mtSetReadUPCE(hBarcode, 0);
	mtSetReadCodabar(hBarcode, 0);
	mtSetReadPDF417(hBarcode, 0);
	mtSetReadDataMatrix(hBarcode, 0);
	mtSetReadDatabar(hBarcode, 0);
	mtSetReadMicroPDF417(hBarcode, 0);
	mtSetReadQRCode(hBarcode, 0);
	
	// 强制设定只识别一个条码以加快速度
	mtSetMultipleRead(hBarcode, 0);
	mtSetMaxBarcodesPerPage(hBarcode, 1);
	mtSetUseFastScan(hBarcode, 1);

	// 设定条码周围空白区大小
	mtSetQuietZoneSize(hBarcode, 0);

	// 设定扫描线步长
	mtSetLineJump(hBarcode, 1);

	// 设定全方向扫描
	mtSetScanDirection(hBarcode, 15);

	// 设定扫描角度0-5，对应0-45度
	//mtSetSkewTolerance(hBarcode, 0);
	mtSetSkewedLinear(hBarcode, 1);

	// 设定彩色图像搜索次数
	mtSetColorProcessingLevel(hBarcode, 2);

	// 设定条码包含字符数.(机票专用)
	mtSetMinLength(hBarcode, 11);
	mtSetMaxLength(hBarcode, 12);

	// 设定中值滤波.
	mtSetMedianFilter(hBarcode, 0);

	// 设定显示检验码
	mtSetShowCheckDigit(hBarcode, 1);
}

void BarcodeRead(Mat &image, CBarcodeReader &reader , HANDLE &hBarcode, vector<Barcode>& barcodes)
{
	barcodes.clear();

	//Mat 2 DIB
	HDIB hbmp = Mat2Dib(image);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Softek
	// 打开图像句柄读取条码
	int nBarCodes = mtScanBarCodeFromDIB(hBarcode, (long)hbmp);
	if (nBarCodes <= 0)// 条码识别错误或者图像中没找到条码
	{
		//mtDestroyBarcodeInstance(hBarcode);
		//重试一次
		mtSetColorProcessingLevel(hBarcode, 3);
		int nBarCodes_Re = mtScanBarCodeFromDIB(hBarcode, (long)hbmp);
		mtSetColorProcessingLevel(hBarcode, 2);
		nBarCodes = nBarCodes_Re;
	}

	if (nBarCodes > 0)//Softek找到条码
	{
		// 将条码识别结果传回调用方
		for (int iIndex = 1; iIndex <= nBarCodes; iIndex++)
		{
			Barcode recogBarcode;
			recogBarcode.iType = BF_CODE_39;
			sprintf(recogBarcode.sBarCodeData, "%s", mtGetBarString(hBarcode, iIndex));
			int nDirection;
			nDirection = mtGetBarStringDirection(hBarcode, iIndex);
			switch (nDirection)
			{
			case 1:
				recogBarcode.iAngle = 0;
				break;
			case 2:
				recogBarcode.iAngle = 90;
				break;
			case 4:
				recogBarcode.iAngle = 180;
				break;
			case 8:
				recogBarcode.iAngle = -90;
				break;
			case 16:
				recogBarcode.iAngle = 0;
				break;
			case 32:
				recogBarcode.iAngle = 45;
				break;
			case 64:
				recogBarcode.iAngle = 135;
				break;
			case 128:
				recogBarcode.iAngle = -135;
				break;
			default:
				recogBarcode.iAngle = -45;
				break;
			}
			long topLeftX, topLeftY, bottomRightX, bottomRightY;
			mtGetBarStringPos(hBarcode, iIndex, &topLeftX, &topLeftY, &bottomRightX, &bottomRightY);
			recogBarcode.rcBarcodeRegion = Rect(Point(topLeftX, topLeftY), Point(bottomRightX, bottomRightY));
			barcodes.push_back(recogBarcode);
		}
	}
	else //Softek没找到条码，改换DynamSoft试一试
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// DynamSoft
		// Read barcode
		int iRet = reader.DecodeDIB(hbmp);
			
		//output barcode
		SBarcodeResultArray *paryResult = NULL;
		reader.GetBarcodes(&paryResult);
		
		for (int iIndex = 0; iIndex < paryResult->iBarcodeCount; iIndex++)
		{
			Barcode recogBarcode;
			recogBarcode.iType = BF_CODE_39;
			sprintf(recogBarcode.sBarCodeData, "%s", paryResult->ppBarcodes[iIndex]->pBarcodeData);
			recogBarcode.iAngle = paryResult->ppBarcodes[iIndex]->iAngle;
			recogBarcode.rcBarcodeRegion = Rect(paryResult->ppBarcodes[iIndex]->iLeft, paryResult->ppBarcodes[iIndex]->iTop, paryResult->ppBarcodes[iIndex]->iWidth, paryResult->ppBarcodes[iIndex]->iHeight);
			barcodes.push_back(recogBarcode);
		}
		
		reader.FreeBarcodeResults(&paryResult);
	}

	// 清理现场
	GlobalFree(hbmp);
}