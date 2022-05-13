//+------------------------------------------------------------------+
//|                                                       Matrix.mq5 |
//|                                Copyright 2022, Omega Joctan Ltd. |
//|                        https://www.mql5.com/en/users/omegajoctan |
//+------------------------------------------------------------------+
#property copyright "Copyright 2022, MetaQuotes Ltd."
#property link      "https://www.mql5.com/en/users/omegajoctan"
#property version   "1.00"

#include "MatrixRegression.mqh";
#include "LinearRegressionLib.mqh";
CSimpleMatLinearRegression matlr;
CSimpleLinearRegression lr;
//+------------------------------------------------------------------+
//| Script program start function                                    |
//+------------------------------------------------------------------+
void OnStart()
  {
    //double x[] = {651,762,856,1063,1190,1298,1421,1440,1518}; //stands for sales
    //double y[] = {23,26,30,34,43,48,52,57,58}; //money spent on ads
//---
    double x[], y[];
    string file_name = "NASDAQ_DATA.csv", delimiter = ",";
    
    lr.GetDataToArray(x,file_name,delimiter,1);
    lr.GetDataToArray(y,file_name,delimiter,2);
    
    //ArrayPrint(x);
    //ArrayPrint(y);
    
    matlr.Init(x,y);
    matlr.LinearRegression();
  }
//+------------------------------------------------------------------+
