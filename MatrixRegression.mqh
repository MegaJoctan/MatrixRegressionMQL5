//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
#property copyright "Copyright 2022, Omega Joctan."
#property link      "https://mql5.com/en/users/omegajoctan"
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+

#define DBL_MAX_MIN(val) if (val>DBL_MAX) Alert("Function ",__FUNCTION__,"\n Maximum Double value Allowed reached"); if (val<DBL_MIN && val>0) Alert("Function ",__FUNCTION__,"\n MInimum Double value Allowed reached") 

class CSimpleMatLinearRegression
  {
      private:
                           int m_rowsize;
                           double Betas[]; //vector for our model coefficient
                           double m_xvalues[];
                           double m_yvalues[];
                           bool   m_debug;
                          
      protected:      
                           void   MatrixInverse(double &Matrix[], double& output_mat[]);
                           int    MatrixtypeSquare(int sizearr);
                           void   MatrixDetectType(double &Matrix[], int rows, int &__r__,int &__c__);
                           void   MatrixTranspose(double& Matrix[],string mat_type = "4x4");
                           
                           void   MatrixPrint(double &Matrix[],int rows,int cols,int digits=0);
                           void   MatrixMultiply(double &A[],double &B[],double &output_arr[],int row1,int col1,int row2,int col2);
                           void   MatrixUnTranspose(double &Matrix[],int torows, int tocolumns);
                           
      public:
                           CSimpleMatLinearRegression(void);
                          ~CSimpleMatLinearRegression(void);
                           
                           void Init(double& x[], double& y[], bool debugmode=true);
                           void LinearRegression();
                                        
                          
  };
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
CSimpleMatLinearRegression::CSimpleMatLinearRegression(void) {};
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
CSimpleMatLinearRegression::~CSimpleMatLinearRegression(void) {};
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CSimpleMatLinearRegression::Init(double &x[],double &y[], bool debugmode=true)
 {
    ArrayResize(Betas,2); //since it is simple linear Regression we only have two variables x and y
        
    if (ArraySize(x) != ArraySize(y))
      Alert("There is variance in number of independent variables and dependent variables \n Calculations may fall short");
      
    m_rowsize = ArraySize(x);
    
    ArrayResize(m_xvalues,m_rowsize+m_rowsize); //add one row size space for the filled values 
    ArrayFill(m_xvalues,0,m_rowsize,1); //fill the first row with one(s)
    
    ArrayCopy(m_xvalues,x,m_rowsize,0,WHOLE_ARRAY); //add x values to the array starting where the filled values ended
    ArrayCopy(m_yvalues,y);
    
    //Print("Design matrix");
    //ArrayPrint(m_xvalues);
    
    m_debug=debugmode;
 }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CSimpleMatLinearRegression::LinearRegression(void)
 {
   /* To find the betas the formula is  B = (xT x) -1 (xT y)
      so let's first find the xT and x
   */
    int _digits = 7; 

//---
    
    double xTx[]; //x transpose matrix times x
    double xT[];
    ArrayCopy(xT,m_xvalues);
    
    int tr_rows = m_rowsize,
        tr_cols = 1+1; //since we have one independent variable we add one for the space created by those values of one
        
    //Print("Transposed Matrix");
    //MatrixPrint(m_xvalues,tr_rows,tr_cols);
    
    MatrixUnTranspose(m_xvalues,tr_cols,tr_rows);
    
    Print("UnTransposed Matrix");
    MatrixPrint(m_xvalues,tr_cols,tr_rows);
//---
    MatrixMultiply(xT,m_xvalues,xTx,tr_cols,tr_rows,tr_rows,tr_cols);
    
    Print("xTx");
    MatrixPrint(xTx,tr_cols,tr_cols,5);
//---

  
    double inverse_xTx[];
    MatrixInverse(xTx,inverse_xTx); //findind the inverse of matrix
    
    //Print("x values");
    //ArrayPrint(m_xvalues);
    //Print("y values");
    //ArrayPrint(m_yvalues);
    double xTy[];
    MatrixMultiply(xT,m_yvalues,xTy,tr_cols,tr_rows,tr_rows,1); //1 at the end is because the y values matrix will always have one column which is it
    
    
    //to find our coefficient 
    
   if (m_debug) 
    {
       Print("xTy");
       MatrixPrint(xTy,tr_rows,1,_digits);
       
       Print("inverse xtx");
       MatrixPrint(inverse_xTx,2,2,_digits); //inverse of simple lr will always be a 2x2 matrix
    }
    
    MatrixMultiply(inverse_xTx,xTy,Betas,2,2,2,1);
    
    if (m_debug)
     {
      Print("coefficients");
      MatrixPrint(Betas,2,1,5); // for simple lr our betas matrix will be a 2x1
     } 
 }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CSimpleMatLinearRegression::MatrixInverse(double &Matrix[],double &output_mat[])
 {
// According to Matrix Rules the Inverse of a matrix can only be found when the 
// Matrix is Identical Starting from a 2x2 matrix so this is our starting point
   
   int matrix_size = ArraySize(Matrix);

   if (matrix_size > 4)
     Print("Matrix allowed using this method is a 2x2 matrix Only");
  
  if (matrix_size==4)
     {
       MatrixtypeSquare(matrix_size);
       //first step is we swap the first and the last value of the matrix
       //so far we know that the last value is equal to arraysize minus one
       int last_mat = matrix_size-1;
       
       ArrayCopy(output_mat,Matrix);
       
       // first diagonal
       output_mat[0] = Matrix[last_mat]; //swap first array with last one
       output_mat[last_mat] = Matrix[0]; //swap the last array with the first one
       double first_diagonal = output_mat[0]*output_mat[last_mat];
       
       // second diagonal  //adiing negative signs >>>
       output_mat[1] = - Matrix[1]; 
       output_mat[2] = - Matrix[2];
       double second_diagonal = output_mat[1]*output_mat[2]; 
       
       if (m_debug)
        {
          Print("Diagonal already Swapped Matrix");
          MatrixPrint(output_mat,2,2);
        }
        
       //formula for inverse is 1/det(xTx) * (xtx)-1
       //determinant equals the product of the first diagonal minus the product of the second diagonal
       
       double det =  first_diagonal-second_diagonal;
       
       if (m_debug)
       Print("determinant =",det);
       
       
       for (int i=0; i<matrix_size; i++)
          { output_mat[i] = output_mat[i]*(1/det); DBL_MAX_MIN(output_mat[i]); }
       
     }
 }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
int CSimpleMatLinearRegression::MatrixtypeSquare(int sizearr)
 { 
//function for checking if the matrix is a square matrix or not
   
    int squarematrices[9] = {4,9,16,25,36,49,64,81,100}; //the squares of 2...10
    //int divident=0;
    int type=0;
    
     for (int i=0; i<9; i++)
       {
          if (sizearr % squarematrices[i] == 0)
             {
                //divident = sizearr/squarematrices[i];
                type = (int)sqrt(sizearr);
                printf("This is a %dx%d Matrix",type,type);
                break;
             }
             
          if (i==9 && sizearr % squarematrices[i] !=0 ) //if after 10 iterations the matrix size couldn't be found on the list then its not a square one
           Print("This is not a Square Matrix");
       }
    return (type);
 }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CSimpleMatLinearRegression::MatrixPrint(double &Matrix[],int rows,int cols,int digits=0)
 {
   Print("[ ");
   int start = 0;
   for (int i=0; i<cols; i++)
     {
       ArrayPrint(Matrix,digits,NULL,start,rows/cols);
       start += rows;     
     }
   Print("]");
 }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CSimpleMatLinearRegression::MatrixMultiply(double &A[],double &B[],double &output_arr[],int row1,int col1,int row2,int col2)
 {
//---   
   double MultPl_Mat[]; //where the multiplications will be stored
   
   if (col1 != row2)
        Alert("Matrix Multiplication Error, \n The number of columns in the first matrix is not equal to the number of rows in second matrix");
 
   else 
    { 
        ArrayResize(MultPl_Mat,row1*col2);
        
        int mat1_index, mat2_index;
        
        if (col1==1)  //Multiplication for 1D Array
         {
            for (int i=0; i<row1; i++)
              for(int k=0; k<row1; k++)
                 {
                   int index = k + (i*row1);
                   MultPl_Mat[index] = A[i] * B[k];
                 }
           //Print("Matrix Multiplication output");
           //ArrayPrint(MultPl_Mat);
         }
        else 
         {
         //if the matrix has more than 2 dimensionals
         for (int i=0; i<row1; i++)
          for (int j=0; j<col2; j++)
            { 
               int index = j + (i*col2);
               MultPl_Mat[index] = 0;
               
               for (int k=0; k<col1; k++)
                 {
                     mat1_index = k + (i*row2);   //k + (i*row2)
                     mat2_index = j + (k*col2);   //j + (k*col2)
                     
                     //Print("index out ",index," index a ",mat1_index," index b ",mat2_index);
                     
                       MultPl_Mat[index] += A[mat1_index] * B[mat2_index];
                       DBL_MAX_MIN(MultPl_Mat[index]);
                 }
               //Print(index," ",MultPl_Mat[index]);
             }
           ArrayCopy(output_arr,MultPl_Mat);
           ArrayFree(MultPl_Mat);
       }
    }
 }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CSimpleMatLinearRegression::MatrixDetectType(double &Matrix[],int rows,int &__r__,int &__c__)
{
  int size = ArraySize(Matrix);
     __c__ = size/rows;
     __r__ = size/__c__;
     
   //if (m_debug) 
   // printf("Matrix Type \n %dx%d Before Transpose/Original \n %dx%d After Transposed/Array Format",__r__,__c__,__c__,__r__);
}
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CSimpleMatLinearRegression::MatrixUnTranspose(double &Matrix[],int torows, int tocolumns)
 {
    int rows, columns;
    
    double Temp_Mat[]; //temporary array
  
      rows = torows;
      columns = tocolumns;
       
//--- UnTransposing Array Starting

          ArrayResize(Temp_Mat,ArraySize(Matrix));
          
          int index=0; int start_incr = 0;
          
           for (int C=0; C<columns; C++)
              {
                 start_incr= C; //the columns are the ones resposible for shaping the new array
                 
                  for (int R=0; R<rows; R++, index++) 
                     {
                       
                       Temp_Mat[index] = Matrix[start_incr];                       
                       
                       //if (m_debug)
                       //Print("Old Array Access key = ",index," New Array Access Key = ",start_incr);
                       
                       start_incr += columns;
                     }
                     
              }
       
       ArrayCopy(Matrix,Temp_Mat);
       ArrayFree(Temp_Mat);      
 }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+

