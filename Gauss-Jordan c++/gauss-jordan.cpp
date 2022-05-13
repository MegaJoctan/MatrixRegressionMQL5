#include <stdio.h>
#include <conio.h>
#include <math.h>

#define SIZE 10

int main()
{
		 float a[SIZE][SIZE], x[SIZE], ratio;
		 int i,j,k,n;
		 
		//  clrscr();

		 /* Inputs */
		 /* 1. Reading order of matrix */ 

		 printf("Enter order of matrix: ");
		 scanf("%d", &n);

		 /* 2. Reading Matrix */
		 printf("Enter coefficients of Matrix:\n");
		 for(i=1;i<=n;i++)
		 {
			  for(j=1;j<=n;j++)
			  {
				   printf("a[%d][%d] = ",i,j);
				   scanf("%f", &a[i][j]);
			  }
		 }
		 /* Augmenting Identity Matrix of Order n */
		 for(i=1;i<=n;i++)
		 {
			  for(j=1;j<=n;j++)
			  {
				   if(i==j)
				   {
				    	a[i][j+n] = 1;
				   }
				   else
				   {
				    	a[i][j+n] = 0;
				   }
			  }
		 }

/*
		printf("\n before gauss jordan elimination");
		for (i=1; i<=n; i++)
		   for (k=1; k<=n*2; k++)
		     {
				 printf("\n  a[%d][%d] = %.4f ",i,k, a[i][k]);
			 }
*/

		 /* Applying Gauss Jordan Elimination */
		 for(i=1;i<=n;i++)
		 {
			  if(a[i][i] == 0.0)
			  {
				   printf("Mathematical Error!");
				   exit(0);
			  }
			  for(j=1;j<=n;j++)
			  {    
				   if(i!=j)
				   {
					    ratio = a[j][i]/a[i][i];
					    for(k=1;k<=2*n;k++)
					    {
							// printf("\n a[%d][%d] val =%.2f a[%d][%d] val =%.2f ",j,k, a[j][k], i,k ,a[i][k]);
							
							// printf("\n jk val =%.4f, ratio = %.4f , ik val =%.4f ",a[j][k],ratio,a[i][k]);
					      	a[j][k] = a[j][k] - ratio*a[i][k];
							
						 	// printf(" \n a[%d][%d] = %0.5f ",j,k, a[j][k]);
					    }
				   }
			  }
		 }

/*
		printf("\n after gauss jordan elimination");
		for (i=1; i<=n; i++)
		   for (k=1; k<=n*2; k++)
		     {
				 printf("\n  a[%d][%d] = %.4f ",i,k, a[i][k]);
			 }
*/

		 //Row Operation to Make Principal Diagonal to 1 
		 for(i=1;i<=n;i++)
		 {
			  for(j=n+1;j<=2*n;j++)
			  {
			    //  printf("\n a[%d][%d] = %.2f a[%d][%d] = %.2f ",i,j, a[i][j], i,i, a[i][i]);
			   	 a[i][j] = a[i][j]/a[i][i];
				//  printf("\n mathematical operation = %.4f",a[i][j]);
			  }
		 }


		 //Displaying Inverse Matrix
		 printf("\nInverse Matrix is:\n");
		 for(i=1;i<=n;i++)
		 {
			  for(j=n+1;j<=2*n;j++)
			  {
			   	printf("%0.3f\t",a[i][j]);
			  }
			  printf("\n");
		 }


		 getch();
		 return(0);
}