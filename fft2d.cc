// Distributed two-dimensional Discrete FFT transform
// Namrata Nadagouda

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include <chrono>
#include <vector>
#include <thread>

#include "Complex.h"
#include "InputImage.h"

constexpr unsigned int NUMTHREADS = 4;

using namespace std;
namespace sc = std::chrono;

int per_thread;
int extra;

int width;
int height;
Complex *data = new Complex[width * height];
double den;


//undergrad students can assume NUMTHREADS will evenly divide the number of rows in tested images
//graduate students should assume NUMTHREADS will not always evenly divide the number of rows in tested images.
// I will test with a different image than the one given

//dir ---> true:DFT, false:IFT
void Transform1D(Complex* h, int w, bool dir, Complex* H)
{
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.
  
  
  for (int i = 0; i < w; i++) {  // For each output element
    Complex *val = new Complex();
    
    for (int j = 0; j < w; j++) {  // For each input element
      double angle = 2 * M_PI * i * j / w;
      if (dir)
        angle = (-1)*angle;
      Complex factor(cos(angle), sin(angle));
      *val = *val + h[j] * factor;
    }
    H[i] = *val;
  }

}

void thread_func1(int threadID) 
{ 
  if (!per_thread){
    if(width>threadID){
      Complex *vec1 = new Complex[width];
      for (int j = 0; j < width; j++) 
        vec1[j] = data[threadID*width + j]; 

      Complex *vec2 = new Complex[width];
      Transform1D(vec1, width, true, vec2);
    
      for (int j = 0; j < width; j++) 
      data[threadID*width + j] = vec2[j];
    }
    return;
  }  

  int start = per_thread * threadID;

  for (int i = 0; i < per_thread; i++) {
    
    Complex *vec1 = new Complex[width];
    for (int j = 0; j < width; j++) 
      vec1[j] = data[(i+start)*width + j]; 

    Complex *vec2 = new Complex[width];
    Transform1D(vec1, width, true, vec2);
    
    for (int j = 0; j < width; j++) 
      data[(i+start)*width + j] = vec2[j];

  }

  int start2 = per_thread * NUMTHREADS + threadID;
  if (extra>threadID){
    Complex *vec1 = new Complex[width];
    for (int j = 0; j < width; j++) 
      vec1[j] = data[(threadID + start2)*width + j]; 

    Complex *vec2 = new Complex[width];
    Transform1D(vec1, width, true, vec2);
    
    for (int j = 0; j < width; j++) 
      data[(threadID + start2)*width + j] = vec2[j];
  }

 }

void thread_func2(int threadID) 
{ 
  if (!per_thread){
    if(width>threadID){
      Complex *vec1 = new Complex[width];
      for (int j = 0; j < width; j++) 
        vec1[j] = data[j*width + threadID]; 

      Complex *vec2 = new Complex[width];
      Transform1D(vec1, width, true, vec2);
    
      for (int j = 0; j < width; j++) 
      data[j*width + threadID] = vec2[j];
    }
    return;
  }  

  int start = per_thread * threadID;

  for (int i = 0; i < per_thread; i++) {
    
    Complex *vec1 = new Complex[width];
    for (int j = 0; j < width; j++) 
      vec1[j] = data[start + j*width + i]; 

    Complex *vec2 = new Complex[width];
    Transform1D(vec1, width, true, vec2);
    
    for (int j = 0; j < width; j++) 
      data[start + j*width + i] = vec2[j];

  }

  int start2 = per_thread * NUMTHREADS + threadID;
  if (extra>threadID){
    Complex *vec1 = new Complex[width];
    for (int j = 0; j < width; j++) 
      vec1[j] = data[start2 + j*width + threadID]; 

    Complex *vec2 = new Complex[width];
    Transform1D(vec1, width, true, vec2);
    
    for (int j = 0; j < width; j++) 
      data[start2 + j*width + threadID] = vec2[j];
  }

 }

void thread_func3(int threadID) 
{ 
  Complex mul(den, 0);

  if (!per_thread){
    if(width>threadID){
      Complex *vec1 = new Complex[width];
      for (int j = 0; j < width; j++) 
        vec1[j] = data[threadID*width + j]; 

      Complex *vec2 = new Complex[width];
      Transform1D(vec1, width, false, vec2);
    
      for (int j = 0; j < width; j++) 
      data[threadID*width + j] = mul*vec2[j];
    }
    return;
  }  

  int start = per_thread * threadID;
  
  for (int i = 0; i < per_thread; i++) {
    
    Complex *vec1 = new Complex[width];
    for (int j = 0; j < width; j++) 
      vec1[j] = data[(i+start)*width + j]; 

    Complex *vec2 = new Complex[width];
    Transform1D(vec1, width, false, vec2);
    
    for (int j = 0; j < width; j++) 
      data[(i+start)*width + j] = mul*vec2[j];

  }

  int start2 = per_thread * NUMTHREADS + threadID;
  if (extra>threadID){
    Complex *vec1 = new Complex[width];
    for (int j = 0; j < width; j++) 
      vec1[j] = data[(threadID + start2)*width + j]; 

    Complex *vec2 = new Complex[width];
    Transform1D(vec1, width, false, vec2);
    
    for (int j = 0; j < width; j++) 
      data[(threadID + start2)*width + j] = mul*vec2[j];
  }

 }

void thread_func4(int threadID) 
{ 
  Complex mul(den, 0);

  if (!per_thread){
    if(width>threadID){
      Complex *vec1 = new Complex[width];
      for (int j = 0; j < width; j++) 
        vec1[j] = data[j*width + threadID]; 

      Complex *vec2 = new Complex[width];
      Transform1D(vec1, width, false, vec2);
    
      for (int j = 0; j < width; j++) 
      data[j*width + threadID] = mul*vec2[j];
    }
    return;
  }  

  int start = per_thread * threadID;

  for (int i = 0; i < per_thread; i++) {
    
    Complex *vec1 = new Complex[width];
    for (int j = 0; j < width; j++) 
      vec1[j] = data[start + j*width + i]; 

    Complex *vec2 = new Complex[width];
    Transform1D(vec1, width, false, vec2);
    
    for (int j = 0; j < width; j++) 
      data[start + j*width + i] = mul*vec2[j];

  }

  int start2 = per_thread * NUMTHREADS + threadID;
  if (extra>threadID){
    Complex *vec1 = new Complex[width];
    for (int j = 0; j < width; j++) 
      vec1[j] = data[start2 + j*width + threadID]; 

    Complex *vec2 = new Complex[width];
    Transform1D(vec1, width, false, vec2);
    
    for (int j = 0; j < width; j++) 
      data[start2 + j*width + threadID] = mul*vec2[j];
  }

 }
 

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.
  // 2) Create a vector of complex objects of size width * height to hold
  //    values calculated
  // 3) Do the individual 1D transforms on the rows assigned to each thread
  // 4) Force each thread to wait until all threads have completed their row calculations
  //    prior to starting column calculations
  // 5) Perform column calculations
  // 6) Wait for all column calculations to complete
  // 7) Use SaveImageData() to output the final results

  InputImage image(inputFN);  // Create the helper object for reading the image
  // Step (1) in the comments is the line above.
  // Your code here, steps 2-7

  width = image.GetWidth();
  height = image.GetHeight();
  den = (double)1/width;

  data = image.GetImageData();

  per_thread = width/NUMTHREADS;
  extra = width % NUMTHREADS;

  thread *t = new thread[NUMTHREADS];

  //row transformation
  for (int i = 0; i <NUMTHREADS; i++) 
    t[i] = thread(thread_func1,i);
  
  for (int i = 0; i <NUMTHREADS; i++) 
    t[i].join();

  //column transformation
  for (int i = 0; i <NUMTHREADS; i++) 
    t[i] = thread(thread_func2, i);
  
  for (int i = 0; i <NUMTHREADS; i++) 
    t[i].join();

  image.SaveImageData("MyAfter2D.txt", data, width, height);

  //row inverse transformation
  for (int i = 0; i <NUMTHREADS; i++) 
    t[i] = thread(thread_func3, i);
  
  for (int i = 0; i <NUMTHREADS; i++) 
    t[i].join();
  
  //column inverse transformation
  for (int i = 0; i <NUMTHREADS; i++) 
    t[i] = thread(thread_func4, i);
  
  for (int i = 0; i <NUMTHREADS; i++) 
    t[i].join();


  image.SaveImageDataReal("MyAfterInverse.txt", data, width, height);

}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line

  auto start = sc::high_resolution_clock::now();

  Transform2D(fn.c_str()); // Perform the transform.

  auto end = sc::high_resolution_clock::now();

  cout << "Elapsed Time: " << sc::duration_cast<sc::milliseconds>(end - start).count() << "ms" << std::endl;

}  
  

  
