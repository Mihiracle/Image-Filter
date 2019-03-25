#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include "bmplib.h"

using namespace std;

//============================Add function prototypes here======================
void dummy(unsigned char [SIZE][SIZE][RGB], unsigned char [SIZE][SIZE][RGB]); 
void convolve(unsigned char [SIZE][SIZE][RGB], unsigned char [SIZE][SIZE][RGB], int, double [11][11]); 
void sobel(unsigned char [SIZE][SIZE][RGB], unsigned char [SIZE][SIZE][RGB]); 
void gaussian(double [11][11], int, double); 
void gaussian_filter(unsigned char [SIZE][SIZE][RGB], unsigned char [SIZE][SIZE][RGB], int, double); 
void unsharp(unsigned char [SIZE][SIZE][RGB], unsigned char [SIZE][SIZE][RGB], int, double, double); 
float gaussian_equation(int , int , double); 


//============================Do not change code in main()======================

#ifndef AUTOTEST

int main(int argc, char* argv[])
{
   //First check argc
  if(argc < 3)
    {
      //we need at least ./filter <input file> <filter name> to continue
      cout << "usage: ./filter <input file> <filter name> <filter parameters>";
      cout << " <output file name>" << endl;
      return -1;
    }
   //then check to see if we can open the input file
   unsigned char input[SIZE][SIZE][RGB];
   unsigned char output[SIZE][SIZE][RGB];
   char* outfile;
   int N;
   double sigma, alpha;
   //double kernel[11][11];

   // read file contents into input array
   int status = readRGBBMP(argv[1], input); 
   if(status != 0)
   {
      cout << "unable to open " << argv[1] << " for input." << endl;
      return -1;
   }
   //Input file is good, now look at next argument
   if( strcmp("sobel", argv[2]) == 0)
   {
     sobel(output, input);
     outfile = argv[3];
   }
   else if( strcmp("blur", argv[2]) == 0)
   {
     if(argc < 6)
       {
	 cout << "not enough arguments for blur." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     outfile = argv[5];
     gaussian_filter(output, input, N, sigma);
   }
   else if( strcmp("unsharp", argv[2]) == 0)
   {
     if(argc < 7)
       {
	 cout << "not enough arguments for unsharp." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     alpha = atof(argv[5]);
     outfile = argv[6];
     unsharp(output, input, N, sigma, alpha);

   }
   else if( strcmp("dummy", argv[2]) == 0)
   {
     //do dummy
     dummy(output, input);
     outfile = argv[3];
   }
   else
   {
      cout << "unknown filter type." << endl;
      return -1;
   }

   if(writeRGBBMP(outfile, output) != 0)
   {
      cout << "error writing file " << argv[3] << endl;
   }

}   

#endif 

//=========================End Do not change code in main()=====================


// Creates an identity kernel (dummy kernel) that will simply
// copy input to output image and applies it via convolve()
void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
   double k[11][11];
   for (int i = 0; i < 3; i++)
   {
      for(int j = 0; j < 3; j++)
      {
         k[i][j] = 0;
      }
   }
   k[1][1] = 1;
   convolve(out, in, 3, k); 
}


// Convolves an input image with an NxN kernel to produce the output kernel
void convolve(unsigned char out[][SIZE][RGB],unsigned char in[][SIZE][RGB],int N, double kernel[][11])
{
 
  int halfK = (int) N / 2; 

  int padded[SIZE+11][SIZE+11][RGB];  // Use for input image with appropriate 
                                       // padding
  int temp[SIZE][SIZE][RGB];          // Use for the unclamped output pixel 
                                       // values then copy from temp to out, 
                                       // applying clamping 
  

   //first set all of padded to 0 (black)
  for (int i =0; i < SIZE+N; i++) { 
    for (int j =0; j < SIZE+N; j++) { 
      for (int k = 0; k < RGB; k++) { 
        padded[i][j][k] = 0; 
      }
    }
  }


   //now copy input into padding to appropriate locations
  for (int i = 0; i < SIZE; i++) { 
    for (int j =0; j < SIZE; j++) { 
      for (int k=0; k<RGB;k++) { 
        padded[i+halfK][j+halfK][k] = in[i][j][k];
      }
    }
  }
  
 //initialize temp pixels to 0 (black)
  for (int i =0; i < SIZE; i++) { 
    for (int j =0; j < SIZE; j++) { 
      for (int k = 0; k < RGB; k++) { 
        temp[i][j][k] = 0; 
      }
    }
  }

  //now perform convolve (using convolution equation on each pixel of the 
  // actual image) placing the results in temp (i.e. unclamped result)
  for(int y= halfK;y<SIZE+halfK ;y++) 
    for(int x= halfK;x<SIZE+halfK;x++) 
      for(int k=0;k<RGB;k++)
        for(int i= -halfK ; i<= halfK; i++)
            for(int j= -halfK; j<= halfK; j++)
                  temp[y-halfK][x-halfK][k] += padded[y+i+halfK-1][x+j+halfK-1][k]*kernel[halfK+i][halfK+j];
    
   //now clamp and copy to output
  for (int i = 0; i < SIZE; i++) { 
    for (int j = 0; j < SIZE; j++) { 
      for (int k = 0; k < RGB; k++) {
        if (temp[i][j][k] > 255) { 
          temp[i][j][k] = 255; 
        } 
        if (temp[i][j][k] < 0) { 
          temp[i][j][k] = 0; 
        }
      }
    }
  } 

  for (int i = 0; i < SIZE; i++) { 
    for (int j = 0; j < SIZE; j++) { 
      for (int k = 0; k < RGB; k++) {
        out[i][j][k] = (unsigned char) temp[i][j][k]; 
      }
    }
  } 

}

void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
   double k[11][11];
   double s_h1[3][3] = { {-1, 0, 1}, 
                         {-2, 0, 2}, 
                         {-1, 0, 1} };
   double s_h2[3][3] = { {1, 0, -1}, 
                         {2, 0, -2}, 
                         {1, 0, -1} };
   
   unsigned char h1_soble[SIZE][SIZE][RGB]; //hold intemediate images
   unsigned char h2_soble[SIZE][SIZE][RGB]; 

   for (int i = 0; i < 11; i++)
   {
      for(int j=0; j < 11; j++)
      {
         k[i][j] = 0;
      }
   }


   // Copy in 1st 3x3 horizontal sobel kernel (i.e. copy s_h1 into k)
  for (int i =0; i < 3; i++) 
    for (int j =0; j < 3; j++) 
      k[i][j] = s_h1[i][j]; 


   // Call convolve to apply horizontal sobel kernel with result in h1_soble
   convolve(h1_soble, in, 3, k); 


   // Copy in 2nd 3x3 horizontal sobel kernel
  for (int i =0; i < 3; i++) 
    for (int j =0; j < 3; j++) 
      k[i][j] = s_h2[i][j]; 


   // Call convolve to apply horizontal sobel kernel with result in h2_soble
   convolve(h2_soble, in, 3, k); 


   // Add the two results (applying clamping) to produce the final output 

  int temp[SIZE][SIZE][RGB]; 

  for (int i=0;i<SIZE;i++) { 
    for (int j=0;j<SIZE;j++) { 
      for(int q=0; q < RGB; q++) { 
        temp[i][j][q] = h1_soble[i][j][q] + h2_soble[i][j][q]; 
      }
    }
  }

  for (int i = 0; i < SIZE; i++) { 
    for (int j = 0; j < SIZE; j++) { 
      for (int q = 0; q < RGB; q++) {
        if (temp[i][j][q] > 255) { 
          temp[i][j][q] = 255; 
        } 
        if (temp[i][j][q] < 0) { 
          temp[i][j][q] = 0; 
        }
      }
    }
  } 

  for (int i = 0; i < SIZE; i++) { 
    for (int j = 0; j < SIZE; j++) { 
      for (int q = 0; q < RGB; q++) {
        out[i][j][q] = (unsigned char) temp[i][j][q]; 
      }
    }
  } 

  
}



void gaussian(double kernel[][11], int N, double sigma) { 
  for (int i = 0; i < 11; i++) { 
    for (int j = 0; j < 11; j++) { 
      kernel[i][j] = 0; 
    }
  } 
  float val = 0.0, totalSum = 0; 
  int gi = 0, gj = 0; 
  for (int i =0; i < N; i++) { 
    gi = i - (N/2); 
    for (int j =0; j<N;j++) { 
    gj = j - (N/2); 
    val = gaussian_equation(gi,gj,sigma);
    totalSum += val; 
    kernel[i][j] = val; //raw gaussian matrix 
    }
  }
  for (int i = 0; i <N;i++) { 
    for (int j=0;j<N;j++) { 
      kernel[i][j] = kernel[i][j] / totalSum; // now normalized 
    }
  }

  cout << "Gaussian Kernel: " << endl; 
  for (int i =0; i < N; i++) { 
    for (int j = 0; j < N; j++) { 
      cout << kernel[i][j] << " | ";
    }
    cout << endl; 
  }
}

float gaussian_equation(int x, int y, double sigma) { 
  float xVal = 0, yVal = 0, finalVal = 0; 
  xVal = (float) (pow(x,2.0)/ (float) (2.0*pow(sigma,2.0))); 
  yVal = (float) (pow(y,2.0)/ (float) (2.0*pow(sigma,2.0)));
  finalVal = exp((-(xVal+yVal)));
  return finalVal;
}

void gaussian_filter(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma) { 
  double kernel[11][11]; 
  gaussian(kernel,N,sigma); 
  convolve(out,in,N,kernel); 
}

void unsharp(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], int N, double sigma, double alpha) { 
  unsigned char blurredImage[SIZE][SIZE][RGB]; 
  int detailedImage[SIZE][SIZE][RGB]; 
  double temp[SIZE][SIZE][RGB]; 

  gaussian_filter(blurredImage,in,N,sigma);

  for (int i = 0; i < SIZE; i++) { 
    for (int j = 0; j < SIZE; j++) { 
      for (int k = 0; k < RGB; k++) { 
        detailedImage[i][j][k] = (int) in[i][j][k] - (int) blurredImage[i][j][k]; 
      }
    }
  }
  
  for (int i = 0; i < SIZE; i++) { 
    for (int j = 0; j < SIZE; j++) { 
      for (int k = 0; k < RGB; k++) { 
        temp[i][j][k] = (int) in[i][j][k] + (int) (alpha*detailedImage[i][j][k]); 
      }
    }
  }

  for (int i = 0; i < SIZE; i++) { 
    for (int j = 0; j < SIZE; j++) { 
      for (int k = 0; k < RGB; k++) {
        if (temp[i][j][k] > 255) { 
          temp[i][j][k] = 255; 
        } 
        if (temp[i][j][k] < 0) { 
          temp[i][j][k] = 0; 
        }
      }
    }
  } 

  for (int i = 0; i < SIZE; i++) { 
    for (int j = 0; j < SIZE; j++) { 
      for (int k = 0; k < RGB; k++) { 
        out[i][j][k] = (unsigned char) temp[i][j][k]; 
      }
    }
  }


}


