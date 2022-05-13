//+------------------------------------------------------------------+
//|                                     multipleMatRegTestScript.mq5 |
//|                                    Copyright 2022, Omega Joctan. |
//|                           https://www.mql5.com/users/omegajoctan |
//+------------------------------------------------------------------+
#property copyright "Copyright 2022, Omega Joctan."
#property link      "https://www.mql5.com/users/omegajoctan"
#property version   "1.00"
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
#include "multipleMatLinearReg.mqh";
CMultipleMatLinearReg matreg;
//+------------------------------------------------------------------+
//| Script program start function                                    |
//+------------------------------------------------------------------+
void OnStart()
  {
//---
      string filename= "NASDAQ_DATA.csv";
      matreg.Init(2,"1,3,4",filename);
      matreg.MultipleMatLinearRegMain();
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+

