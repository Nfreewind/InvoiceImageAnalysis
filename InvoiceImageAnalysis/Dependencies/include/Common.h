#ifndef COMMON_H_
#define COMMON_H_
#define  MAX(a,b)      (((a)>(b))?(a):(b))
#define  MIN(a,b)      (((a)<(b))?(a):(b))


#define D7(pArray,pos,i,j)	( ((i)>0 && (j)<((pArray)->width-1)) ? ((pArray)->pArrayData[(pos)-((pArray)->width)+1]) : 0 )
#define D6(pArray,pos,i,j)	( ((i)>0) ? ((pArray)->pArrayData[(pos)-((pArray)->width)]) : 0 )
#define D5(pArray,pos,i,j)	( ((i)>0&&(j)>0) ? ((pArray)->pArrayData[(pos)-((pArray)->width)-1]) : 0 )
#define D4(pArray,pos,i,j)	( ((j)>0) ? ((pArray)->pArrayData[(pos)-1]) : 0 )
#define D3(pArray,pos,i,j)	( ((i)<((pArray)->height-1)&&(j)>0) ? ((pArray)->pArrayData[(pos)+((pArray)->width)-1]) : 0 )
#define D2(pArray,pos,i,j)	( ((i)<((pArray)->height-1)) ? ((pArray)->pArrayData[(pos)+((pArray)->width)]) : 0 )
#define D1(pArray,pos,i,j)	( ((i)<((pArray)->height-1)&&(j)<((pArray)->width-1)) ? ((pArray)->pArrayData[(pos)+((pArray)->width)+1]) : 0)
#define	D0(pArray,pos,i,j)	( ((j)<((pArray)->width-1)) ? ((pArray)->pArrayData[(pos)+1]) : 0 )

#define N_BITS 32
#define MAX_BIT ((N_BITS + 1) / 2 - 1)


#define NORM_SIZE_C (64)
//#define NORM_SIZE_C_32 (32)

#define NORM_TOTALSIZE_C (4096)
#define NORM_TOTALSIZE_C_32 (32*32)


#define NORM_SIZE_SLANT (127)
#define MAX_EDGEPOINTS (1529)

#define FEATURE_SIZE_EPDC (896)
#define FEATURE_SIZE_EXTEND (32)

#define DIM (64)
#define LIMITEDCHINNESE_NUM (3755+20+3)

//#define LIMITEDITrainStation_NUM (1260-7)//////火车站站名

#define LIMITEDITrainStation_NUM (1260)//////火车站站名

#define LIMITEDITrainNumber_NUM (9)/////火车具体的车厢号
#define LIMITEDITrainDrivingRequest_NUM (23)//////火车乘车要求：限乘当日当次车
#define LIMITEDITrainPayType_NUM (21)//////火车付款方式
#define LIMITEDITrainSeatType_NUM (16)//////火车座位类型
#define LIMITEDSaleTicket_NUM (1260+1) /////售票点
#define LIMITEDSaleTicketSymbol_NUM (LIMITEDSaleTicket_NUM+311)///////// 311怎么来的 70 -6729
#define LIMITEDTrainNumberSymbol_NUM (LIMITEDITrainNumber_NUM+311)///////// 311怎么来的？？？

//#define CHINESE_NUM (3755+18)/////6729
#define CHINESE_NUM (6729)/////6729
#define SYMBOL_NUM  (179)
#define RADICAL_NUM  (96)
#define TOP_NONCHARACTER_NUM (7)
#define BOTTOM_NONCHARACTER_NUM (8)
#define BRACKET_NUM (14)
#define SIMILAR_NUM (239)///239, 相似字从limited里面的limitInx2个数读取
#define SIMILAR_TrainNumber_NUM (1)///239, 相似字从limited里面的limitInx2个数读取
#define SIMILAR_TrainSaleTicket_NUM (53)

#define LimitedNumandEng_NUM (53+2)/////纯数字和英文， 去掉了相似字，2016.05.25
#define LimitedNumandUpEng_NUM (35)////限制级大写英文和数字，去掉了英文O
#define LimitedNumandLowEng_NUM (34)////限制级小写英文和数字,去除i和l


#define NUMERAL_NUM (10)
#define ENGLISH_NUM  (45)
#define VERTICAL_NUM  (16)

#define Date_NUM (3)/////年月日, 2016年6月22日改
#define LIMITEDDate_NUM (Date_NUM+311)/////限制集日期，2016年6月22日修改
#define CSN_NUM    (CHINESE_NUM + 311)///////


#define TaixNumandUpEngSymbol_NUM (38) //////出租车
#define TaixNumandUpEngSymbolFea_NUM (38*4)//////////////////出租车限制集
#define TaixNum_NUM (10)
#define TaixUpEng_NUM (24)
#define TaixSymbol_NUM (4)



#define SRN_NUM    (457)
#define SN_NUM     (840)/////数字英文等
#define T_NUM      (496)
#define B_NUM      (504)
#define BN_NUM     (28)

#define LimitedNumandEngFea_NUM (55*4)
#define LimitedNumandUpEngFea_NUM (35*4)////限制级大写英文和数字，去掉了英文O
#define LimitedNumandLowEngFea_NUM (34*4)////限制级大写英文和数字，去掉了英文O
#define LimitedICBusinessFea_NUM (22*4)//////////////////工商银行第二列：业务产品种类

#define LIMITEDJSBorrowLoadFea_NUM (2)//////江苏银行借贷
#define LIMITEDCCBProofKindFea_NUM (27)////////建设银行凭证种类
#define LIMITEDChinaTradeTypeFea_NUM (32)/////////中国银行交易类型
#define LIMITEDGuangFaTradeType_NUM (31)////////广发银行限制集,交易类型


#define  TaixNumandSymbol_NUM (14)/////用于提取HOG的数字和符号
#define  TaixNumandSymbolFea_NUM (14*4)/////用于提取HOG的数字和符号和特征数目
#define  HOGDIM 576 ////////hog特征维度
//#define  HOGDIM 2304 ////////hog特征维度

typedef struct{
	short			width;
	short			height;
	unsigned char*  pArrayData;
} CharacterArray;

typedef struct{
	unsigned char left;
	unsigned char right;
	unsigned char up;
	unsigned char down;
	unsigned char leftup;
	unsigned char rightup;
	unsigned char leftdown;
	unsigned char rightdown;
}Point_8D_Runlen;


#endif