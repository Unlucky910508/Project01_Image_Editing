///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"

#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <vector>
#include <unordered_map>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char* d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
        data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    int index;
    double tmp;
    for (int i = 0; i < width * height; ++i) {
        index = i * 4;
        tmp = 0.299 * data[index] + 0.587 * data[index + 1] + 0.114 * data[index + 2];
        data[index] = tmp;
        data[index + 1] = tmp;
        data[index + 2] = tmp;
    }
    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    // red: 256 -> 8 , divide 32
    // green: 256 -> 8, divide 32
    // blue 256 -> 4, divide 64

    int tmp;
    int index;

    for (int i = 0; i < width * height; ++i) {
        index = i * 4;
        tmp = data[index] >> 5;
        data[index] = tmp << 5;
        tmp = data[index + 1] >> 5;
        data[index + 1] = tmp << 5;
        tmp = data[index + 2] >> 6;
        data[index + 2] = tmp << 6;
    }
    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity() {
    // find use most frequency color
    // 1 00000 00000 00000
    //   red   green blue
    std::unordered_map<int, int> Counter;
    std::vector<std::pair<int, int>> majorColor;

    int findColor;
    for (int i = 0; i < width * height; ++i) {
        findColor = 0;
        for (int j = 0; j < 3; ++j) {
            findColor += data[i * 4 + j] >> 3; // divid 8
            findColor <<= 5;
        }
        findColor >>= 5;
        //findColor--;
        ++Counter[findColor];
    }
    // write into vecotr
    for (const std::pair<int, int>& i : Counter) {
        std::pair<int, int> tmp; // <color, number>

        tmp.first = i.first;
        tmp.second = i.second;
        majorColor.push_back(tmp);
    }

    qsort(majorColor.data(), majorColor.size(), sizeof(std::pair<int, int>), helper::Pair_color_compare);
    if (majorColor.size() > 256)
        majorColor.resize(256);

    // Print
    //for (int i = 0; i < majorColor.size(); i++) { std::cout << i << ": " << majorColor[i].first << " " << majorColor[i].second << "\n"; }

    {// calculate most close color;
        int r_dis;
        int g_dis;
        int b_dis;

        int best_dis;
        int bestRGB[3];

        int index;
        int tmp_dis;
        int colorNumber = majorColor.size();

        for (int i = 0; i < width * height; ++i) {
            best_dis = 2147483647; // fuck
            for (int j = 0; j < colorNumber; ++j) {
                index = i * 4;
                int tmp = majorColor[j].first;
                int r = tmp & 0b111110000000000;
                int g = tmp & 0b000001111100000;
                int b = tmp & 0b000000000011111;
                
                r >>= 7; // fuck
                g >>= 2; // fuck
                b <<= 3; // fuck

                r_dis = data[index] - r;
                g_dis = data[index + 1] - g;
                b_dis = data[index + 2] - b;

                tmp_dis = r_dis * r_dis + g_dis * g_dis + b_dis * b_dis;

                if (tmp_dis <= best_dis) {
                    best_dis = tmp_dis;
                    bestRGB[0] = r;
                    bestRGB[1] = g;
                    bestRGB[2] = b;
                }
            }
            data[index] = bestRGB[0];
            data[index + 1] = bestRGB[1];
            data[index + 2] = bestRGB[2];
        }

        Counter.clear();
        majorColor.clear();

        return true;
    }// Quant_Populosity
}


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    // GrayScale
    To_Grayscale();

    {// main
        int index;
        const double threshold = 0.5;
        for (int i = 0; i < width * height; ++i) {
            index = i * 4;
            if ((double)data[index] / 256.0 < threshold) {
                data[index] = 0;
                data[index + 1] = 0;
                data[index + 2] = 0;
            }
            else {
                data[index] = 255;
                data[index + 1] = 255;
                data[index + 2] = 255;
            }
        }
    }// main

    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    // GrayScale
    To_Grayscale();

    srand(time(NULL));

    {// main
        int index;
        double threshold;

        for (int i = 0; i < width * height; ++i) {
            index = i * 4;
            double tmp = rand() % 40001;
            tmp /= 100000.0;
            tmp -= 0.2;
            threshold = 0.5 + tmp;
                if ((double)data[index] / 256.0 < threshold) {
                    data[index] = 0;
                    data[index + 1] = 0;
                    data[index + 2] = 0;
                }
                else {
                    data[index] = 255;
                    data[index + 1] =255;
                    data[index + 2] = 255;
                }
        }// main
    }
    return true;

}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    // GrayScale
    To_Grayscale();

    const double threshold = 0.5;
    double** err;

    {// initialize array;
        err = new double* [height];
        for (int i = 0; i < height; ++i) {
            err[i] = new double[width];
        }

        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                err[i][j] = 0;
            }
        }
    }// initialize array;

    {// FS
        bool dir = true;
        int index;

        double oldVal;
        double newVal;
        double transferErr;

        for (int i = 0; i < height; ++i) {
            if (dir) { // -->>
                // for
                for (int j = 0; j < width; ++j) {
                    index = (i * width + j) * 4;

                    //oldVal = data[index];
                    newVal = data[index] + err[i][j];

                    newVal = newVal > 255 ? 255.0 : newVal;
                    newVal = newVal < 0.0 ? 0.0 : newVal;

                    if (newVal / 256.0 < 0.5) {
                        data[index] = 0;
                        data[index + 1] = 0;
                        data[index + 2] = 0;
                        transferErr = newVal;
                    }
                    else {
                        data[index] = 255.0;
                        data[index + 1] = 255.0;
                        data[index + 2] = 255.0;
                        transferErr = newVal - 255.0;
                    }

                    { // TransferErr n
                        if (i == (height - 1) && j == (width - 1)) {
                            break;
                        }
                        else if (i == (height - 1)) {
                            err[i][j + 1] += transferErr * 7.0 / 16.0;
                        }
                        else if (j == (width - 1)) {
                            err[i + 1][j - 1] += transferErr * 3.0 / 16.0;
                            err[i + 1][j] += transferErr * 5.0 / 16.0;
                        }
                        else {
                            err[i][j + 1] += transferErr * 7.0 / 16.0;
                            err[i + 1][j - 1] += transferErr * 3.0 / 16.0;
                            err[i + 1][j] += transferErr * 5.0 / 16.0;
                            err[i + 1][j + 1] += transferErr * 1.0 / 16.0;
                        }

                    } // TransferErr
                } // for
                dir = false;
            }// -->>
            else { // <<--
                for (int j = width - 1; j >= 0; --j) {
                    index = (i * width + j) * 4;

                    newVal = data[index] + err[i][j];

                    newVal = newVal > 255 ? 255.0 : newVal;
                    newVal = newVal < 0.0 ? 0.0 : newVal;

                    if (newVal / 256.0 < 0.5) {
                        data[index] = 0;
                        data[index + 1] = 0;
                        data[index + 2] = 0;
                        transferErr = newVal;
                    }
                    else {
                        data[index] = 255.0;
                        data[index + 1] = 255.0;
                        data[index + 2] = 255.0;
                        transferErr = newVal - 255.0;
                    }

                    { // TransferErr
                        if (i == (height - 1) && j == (width - 1)) {
                            break;
                        }
                        else if (i == (height - 1)) {
                            err[i][j - 1] += transferErr * 7.0 / 16.0;
                        }
                        else if (j == 0) {
                            err[i + 1][j - 1] += transferErr * 3.0 / 16.0;
                            err[i + 1][j] += transferErr * 5.0 / 16.0;
                        }
                        else {
                            err[i][j - 1] += transferErr * 7.0 / 16.0;
                            err[i + 1][j - 1] += transferErr * 3.0 / 16.0;
                            err[i + 1][j] += transferErr * 5.0 / 16.0;
                            err[i + 1][j + 1] += transferErr * 1.0 / 16.0;
                        }

                    } // TransferErr 
                }// for
                dir = true;
            } // <<--
        }
    }
    // free memory
    for (int i = height-1; i <= 0; --i) {
        delete[] err[i];
    }
    delete[] err;
    return false;
} // Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    int index;
    this->To_Grayscale();
    unsigned int sum = 0;
    double threshold;

    {// find the threshold
        double avgIntensity = 0.0;
        int countColor[256]; // Store the distributed of color
        int flag; // the percentile rank of threshold

        // initialize array
        for (int i = 0; i < 256; ++i) {
            countColor[i] = 0;
        }

        // Count average intensity
        for (int i = 0; i < width * height; ++i) {
            index = i * 4;
            avgIntensity += (double)data[index] / 256;
            ++countColor[data[index]];
        }

        avgIntensity = avgIntensity / (width * height);
        std::cout << "avgIntensity: " << avgIntensity << "\n";
        flag = avgIntensity * width * height;

        // find where the flag one
        for (int i = 255; i >= 0; --i) {
            flag -= countColor[i];
            if (flag < 0) {
                threshold = (double)i / 256.0;
                break;
            }
        }
    }// find the threshold

    {// Count by threshold
        for (int i = 0; i < width * height; ++i) {
            index = i * 4;
            if ((double)data[index] / (double)256 < threshold) {
                data[index] = 0;
                data[index + 1] = 0;
                data[index + 2] = 0;
            }
            else {
                data[index] = 255;
                data[index + 1] = 255;
                data[index + 2] = 255;
            }
        }
    }// Count by threshold
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    this->To_Grayscale();

    double threshold;
    int index;

    double mask[4][4] = {
        0.7059, 0.3529, 0.5882, 0.2353,
        0.0588, 0.9412, 0.8235, 0.4118,
        0.4706, 0.7647, 0.8824, 0.1176,
        0.1765, 0.5249, 0.2491, 0.6471,
    };

    {// dither by cluster
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                threshold = mask[i % 4][j % 4];
                //std::cout << "threshold: " << threshold << "\n";
                index = (i * width + j) * 4;
                //std::cout << "index" << index << "\n";
                if ((double)data[index] / 256.0 < threshold) {
                    data[index] = 0;
                    data[index + 1] = 0;
                    data[index + 2] = 0;
                }
                else {
                    data[index] = 255;
                    data[index + 1] = 255;
                    data[index + 2] = 255;
                }
            }
        }
    }// dither by cluster
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    double*** err;

    {// initialize array;
        err = new double** [height + 1];
        for (int i = 0; i <= height; ++i) {
            err[i] = new double*[width + 1];
        }

        for (int i = 0; i <= height; ++i) {
            for (int j = 0; j <= width; ++j) {
                err[i][j] = new double[3];
            }
        }

        for (int i = 0; i <= height; ++i) {
            for (int j = 0; j <= width; ++j) {
                for (int k = 0; k < 3; ++k) {
                    err[i][j][k] = 0;
                }
            }
        }
    }// initialize array

    {// color
        bool dir = true;
        int index;

        double newVal[3];
        double transferErr[3];

        for (int i = 0; i < height; ++i) {
            if (dir) { // -->>
                for (int j = 0; j < width; ++j) {
                    index = (i * width + j) * 4;

                    // count new value by add all error
                    for (int k = 0; k < 3; ++k) {
                        newVal[k] = data[index + k] + err[i][j][k];
                        // in order not to excess 255 or 0
                        newVal[k] = newVal[k] > 255.0 ? 255.0 : newVal[k];
                        newVal[k] = newVal[k] < 0.0 ? 0.0 : newVal[k];
                    }

                    { // assign new data
                        data[index + 0] = (unsigned char)newVal[0] & 0b11100000;
                        data[index + 1] = (unsigned char)newVal[1] & 0b11100000;
                        data[index + 2] = (unsigned char)newVal[2] & 0b11000000;;
                    }

                    { // count transfer error
                        transferErr[0] = (unsigned char)newVal[0] & 0b11111;
                        transferErr[1] = (unsigned char)newVal[1] & 0b11111;
                        transferErr[2] = (unsigned char)newVal[2] & 0b111111;
                    }
                    //std::cout << transferErr[0] << "\n";
                    { // transfer the error
                        for (int k = 0; k < 3; ++k) {
                            if (j != 0) {
                                err[i + 1][j - 1][k] += transferErr[k] * 3.0 / 16.0;
                            }
                            err[i][j + 1][k] += transferErr[k] * 7.0 / 16.0;
                            err[i + 1][j][k] += transferErr[k] * 5.0 / 16.0;
                            err[i + 1][j + 1][k] += transferErr[k] * 1.0 / 16.0;
                        }
                    }
                } // for

                { // turn around
                    dir = false;
                } // turn around

            }// -->>
            else { // <<--
                for (int j = width - 1; j >= 0; --j) {
                    index = (i * width + j) * 4;

                    // count new value by add all error
                    for (int k = 0; k < 3; ++k) {
                        newVal[k] = data[index + k] + err[i][j][k];
                        // in order not to excess 255 or 0
                        newVal[k] = newVal[k] > 255.0 ? 255.0 : newVal[k];
                        newVal[k] = newVal[k] < 0.0 ? 0.0 : newVal[k];
                    }

                    { // diter into 8 color
                        data[0] = (unsigned char)newVal[0] & 0b11100000;
                        data[1] = (unsigned char)newVal[1] & 0b11100000;
                        data[2] = (unsigned char)newVal[2] & 0b11000000;
                    }

                    { // count transfer error
                        transferErr[0] = (unsigned char)newVal[0] & 0b11111;
                        transferErr[1] = (unsigned char)newVal[1] & 0b11111;
                        transferErr[2] = (unsigned char)newVal[2] & 0b111111;
                    }

                    { // transfer the error
                        for (int k = 0; k < 3; ++k) {
                            if (j != 0) {
                                err[i][j - 1][k] += transferErr[k] * 7.0 / 16.0;
                                err[i + 1][j - 1][k] += transferErr[k] * 1.0 / 16.0;
                            }
                            err[i + 1][j][k] += transferErr[k] * 5.0 / 16.0;
                            err[i + 1][j + 1][k] += transferErr[k] * 3.0 / 16.0;
                        }
                    }
                }

                { // turn around
                    dir = true;
                } // turn around

            } // <<--
        }
    }// color

    for (int i = height - 1; i >= 0; --i) {
        for (int j = width - 1; j >= 0; --j) {
            delete[] err[i][j];
        }
    }
    for (int i = height - 1; i >= 0; --i) {
        delete[] err[i];
    }
    delete[] err;

    return true;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    double mask[25];
    for (int i = 0; i < 25; ++i) {
            mask[i] = 0.04;
    }

    Comp_Imp(0.0, mask, 5, 5);

    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    double mask[25]{
        1, 2, 3, 2, 1,
        2, 4, 6, 4, 2,
        3, 6, 9, 6, 3, 
        2, 4, 6, 4, 2,
        1, 2, 3, 2, 1,
    };

    Comp_Imp(1.0/81.0, mask, 5, 5);

    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    double mask[25]{
       1, 4 , 6 , 4 , 1,
       4, 16, 24, 16, 4,
       6, 24, 36, 24, 6,
       4, 16, 24, 16, 4,
       1, 4 , 6 , 4 , 1,
    };

    Comp_Imp(1.0 / 256.0, mask, 5, 5);

    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N(unsigned int N) {
    if (N == 1) {
        return true;
    }
    else {

        double* arr1 = new double[N];
        double* arr2 = new double[N * N];

        double newVal;
        double store;
        double sum = 0.0;

        { //initialize
            for (int i = 0; i < N; ++i) {
                arr1[i] = 0;
            }
            for (int i = 0; i < N * N; ++i) {
                arr2[i] = 0;
            }
        }
        arr1[0] = 1;
        arr1[1] = 1;

        { // finish arr1
            for (int i = 2; i < N; ++i) {
                store = arr1[0];
                for (int j = 1; j < i; ++j) {
                    newVal = store + arr1[j];
                    store = arr1[j];
                    arr1[j] = newVal;
                }
                arr1[i] = 1;
            }
        }

        { // finish arr2
            for (int i = 0; i < N; i++) {
                arr2[i] = arr1[i];
                arr2[i * N] = arr1[i];
            }
            for (int i = 1; i < N; ++i) {
                for (int j = 1; j < N; ++j) {
                    arr2[i * N + j] = arr2[i * N] * arr1[j];
                }
            }
        }

        { // sum them
            for (int i = 0; i < N * N; ++i) {
                sum += arr2[i];
            }
        }

        this->Comp_Imp(1.0 / sum, arr2, N, N);

        for (int i = 0; i < N; ++i) {
            std::cout << arr2[i] << " ";
        }
        std::cout << "\n";

        delete[] arr1;
        delete[] arr2;

        return true;
    }
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    TargaImage tmp(width, height, data);
    this->Filter_Enhance();
    Difference(&tmp);
    /*
    double mask[25] = { 
        
        -1, -4 , -6 , -4 , -1,
        -4, -16, -24, -16, -4,
        -6, -24, 220, -24, -6,
        -4, -16, -24, -16, -4,
        -1, -4 , -6 , -4 , -1
    };

    double Mask[25] = {
         -1, -2, -3, -2, -1,
        -2, -4, -6, -4, -2,
        -3, -6, 72, -6, -3,
        -2, -4, -6, -4, -2,
        -1, -2, -3, -2, -1
    };

    //Comp_Imp(1.0 / 81.0, Mask, 5, 5);
    Comp_Imp(1.0 / 256.0, mask, 5, 5);
    Comp_Imp(1.0 / 256.0, mask, 5, 5);
    */
    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    double mask[25] = { /*-1, -2, -3, -2, -1,
        -2, -4, -6, -4, -2,
        -3, -6, 153, -6, -3,
        -2, -4, -6, -4, -2,
        -1, -2, -3, -2, -1*/
        -1, -4 , -6 , -4 , -1,
        -4, -16, -24, -16, -4,
        -6, -24, 476, -24, -6,
        -4, -16, -24, -16, -4,
        -1, -4 , -6 , -4 , -1
    };

    Comp_Imp(1.0 / 256.0, mask, 5, 5);

    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}

///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size() {
    double mask[9]{
        1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0,
        1.0 / 8.0 , 1.0 / 4.0, 1.0 / 8.0 ,
        1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0,
    };

    // TargaImage
    int h = height / 2;
    int w = width / 2;
    unsigned char* d = new unsigned char[h * w * 4];

    int runK;
    int runM;
    int index;

    double color[3];

    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            index = (i * w + j) * 4;

            { // initialize color
                color[0] = 0.0;
                color[1] = 0.0;
                color[2] = 0.0;
                d[index + 3] = data[(2 * i * w + 2 * j) * 4 + 3];
            }

            for (int k = 0; k < 3; ++k) {
                runK = 2 * i - 1 + k;
                //runK = runK > height ? height - 1 : runK;
                //runK = runK < 0 ? 0 : runK;
                runK = runK >= height ? 2 * height - runK - 1 : runK;
                runK = runK < 0 ? -runK : runK;

                for (int m = 0; m < 3; ++m) {
                    runM = 2 * j - 1 + m;
                    //runM = runM > height ? height - 1 : runM;
                    //runM = runM < 0 ? 0 : runM;
                    runM = runM >= width ? 2 * width - runM - 1 : runM;
                    runM = runM < 0 ? -runM : runM;

                    for (int l = 0; l < 3; ++l) {
                        color[l] += mask[k * 3 + m] * data[(runK * width + runM) * 4 + l];
                    }
                } // for m
            } // for k
            
            for (int k = 0; k < 3; ++k) {
                color[k] = color[k] > 255 ? 255.0 : color[k];
                color[k] = color[k] < 0 ? 0.0 : color[k];
                d[index + k] = (unsigned char)color[k];
            }

        } // for j
    } // for i

    delete[] this->data;

    this->height = h;
    this->width = w;
    this->data = d;

    return true;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size(){
    /*
    int h = height * 2;
    int w = width * 2;
    unsigned char* d = new unsigned char[h * w * 4];

    int tmpI = 0;
    int tmpJ = 0;
    int tmpIndex;
     
    for (int i = 0; i < h; i += 2) {
        tmpJ = 0;
        for (int j = 0; j < w; j += 2) {
            tmpIndex = (tmpI * width + tmpJ) * 4;
            for (int k = 0; k < 3; ++k) {
                d[(i * w + j) * 4 + k] = data[tmpIndex + k];
                d[(i * w + j + 1) * 4 + k] = data[tmpIndex + k];
                d[((i + 1) * w + j) * 4 + k] = data[tmpIndex + k];
                d[((i + 1) * w + j + 1) * 4 + k] = data[tmpIndex + k];
            }
            tmpJ++;
        }
        tmpI++;
    }

    this->height = h;
    this->width = w;
    this->data = d;
    */
    long double mask1[9]{
        1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0,
        1.0 / 8.0 , 1.0 / 4.0, 1.0 / 8.0 ,
        1.0 / 16.0, 1.0 / 8.0, 1.0 / 16.0,
    };
    long double mask2[16]{
        1.0 / 64.0, 3.0 / 64.0, 3.0 / 64.0, 1.0 / 64.0,
        3.0 / 64.0, 9.0 / 64.0, 9.0 / 64.0, 3.0 / 64.0,
        3.0 / 64.0, 9.0 / 64.0, 9.0 / 64.0, 3.0 / 64.0,
        1.0 / 64.0, 3.0 / 64.0, 3.0 / 64.0, 1.0 / 64.0,
    };
    long double mask3[12]{
        1.0 / 32.0, 2.0 / 32.0, 1.0 / 32.0,
        3.0 / 32.0, 6.0 / 32.0, 3.0 / 32.0,
        3.0 / 32.0, 6.0 / 32.0, 3.0 / 32.0,
        1.0 / 32.0, 2.0 / 32.0, 1.0 / 32.0,
    };
    long double mask4[12]{
        1.0 / 32.0, 3.0 / 32.0, 3.0 / 32.0, 1.0 / 32.0,
        2.0 / 32.0, 6.0 / 32.0, 6.0 / 32.0, 2.0 / 32.0,
        1.0 / 32.0, 3.0 / 32.0, 3.0 / 32.0, 1.0 / 32.0,
    };

    // TargaImage
    int h = height * 2;
    int w = width * 2;
    unsigned char* d = new unsigned char[h * w * 4];

    int runK;
    int runM;
    int index;
    int tmpI;
    int tmpJ;

    long double color[3];

    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            index = (i * w + j) * 4;
            tmpI = i / 2;
            tmpJ = j / 2;

            { // initialize color
                color[0] = 0.0;
                color[1] = 0.0;
                color[2] = 0.0;
                d[index + 3] = data[(tmpI * width + tmpJ) * 4 + 3];
            }
               
            // i even, j even
            if ((i % 2 == 0) && (j % 2 == 0)) {
                for (int k = 0; k < 3; k++) {
                    runK = tmpI - 1 + k;
                    runK = runK >= height ? height - 1 : runK;
                    runK = runK <= 0 ? 0 : runK;

                    for (int m = 0; m < 3; ++m) {
                        runM = tmpJ - 1 + m;
                        runM = runM >= width ? width - 1 : runM;
                        runM = runM <= 0 ? 0 : runM;

                        for (int l = 0; l < 3; ++l) {
                            color[l] += mask1[k * 3 + m] * data[(runK * width + runM) * 4 + l];
                        }
                    } // for m
                } // for k
            }
            // i even, j odd
            else if (i % 2 == 0) {
                for (int k = 0; k < 3; k++) {
                    runK = tmpI - 1 + k;
                    runK = runK >= height ? height - 1 : runK;
                    runK = runK <= 0 ? 0 : runK;

                    for (int m = 0; m < 4; ++m) {
                        runM = tmpJ - 1 + m;
                        runM = runM >= width ? width - 1 : runM;
                        runM = runM <= 0 ? 0 : runM;

                        for (int l = 0; l < 3; ++l) {
                            color[l] += mask4[k * 4 + m] * data[(runK * width + runM) * 4 + l];
                        }
                    } // for m
                } // for k
            }
            // i odd, j even
            else if (j % 2 == 0) {
                for (int k = 0; k < 4; ++k) {
                    runK = tmpI - 1 + k;
                    runK = runK >= height ? height - 1 : runK;
                    runK = runK <= 0 ? 0 : runK;

                    for (int m = 0; m < 3; ++m) {
                        runM = tmpJ - 1 + m;
                        runM = runM >= width ? width - 1 : runM;
                        runM = runM <= 0 ? 0 : runM;

                        for (int l = 0; l < 3; ++l) {
                            color[l] += mask3[k * 3 + m] * data[(runK * width + runM) * 4 + l];
                        }
                    } // form
                } // for k
            }
            // i odd, j odd
            else {
                for (int k = 0; k < 4; ++k) {
                    runK = tmpI - 1 + k;
                    runK = runK >= height ? height - 1 : runK;
                    runK = runK <= 0 ? 0 : runK;

                    for (int m = 0; m < 4; ++m) {
                        runM = tmpJ - 1 + m;
                        runM = runM >= width ? width - 1 : runM;
                        runM = runM <= 0 ? 0 : runM;

                        for (int l = 0; l < 3; ++l) {
                            color[l] += mask2[k * 4 + m] * data[(runK * width + runM) * 4 + l];
                        }
                    } // for m
                } // for k
            }

            for (int k = 0; k < 3; ++k) {
                color[k] = color[k] >= 255 ? 255.0 : color[k];
                color[k] = color[k] <= 0 ? 0.0 : color[k];
                d[index + k] = (unsigned char)color[k];
            }

        } // for j
    } // for i

    delete[] this->data;

    this->height = h;
    this->width = w;
    this->data = d;
    
    return true;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}

void TargaImage::Comp_Imp(double modify, double mask[], int h, int w) {
    long double* tmpMask = new long double[(long long)h * (long long)w * sizeof(long double)];
    long double countColor[3];

    int boundH = h / 2;
    int boundW = w / 2;
    int runH;
    int runW;

    { // copy mask value
        for (int i = 0; i < h * w; ++i) {
            tmpMask[i] = mask[i];
        }
        if (modify) {
            for (int i = 0; i < h * w; ++i) {
                tmpMask[i] *= modify;
            }
        }
    }

    { // count
        int index;
        int tmpIndex;

        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                index = (i * width + j) * 4;
                countColor[0] = 0.0;
                countColor[1] = 0.0;
                countColor[2] = 0.0;

                { // processing mask
                    for (int k = 0; k < h; ++k) {
                        runH = i - boundH + k;
                        runH = runH < 0 ? -runH : runH;
                        runH = runH >= height ? 2 * height - runH - 1 : runH;

                        for (int l = 0; l < w; ++l) {
                            runW = j - boundW + l;
                            runW = runW < 0 ? -runW : runW;
                            runW = runW >= width ? 2 * width - runW - 1 : runW;

                            // multi the mask
                            tmpIndex = (runH * width + runW) * 4;
                            for (int m = 0; m < 3; ++m) {
                                countColor[m] += tmpMask[k * w + l] * this->data[tmpIndex + m];
                            }
                        } // for l
                    } // for k

                    for (int k = 0; k < 3; ++k) {
                        countColor[k] = countColor[k] > 255.0 ? 255.0 : countColor[k];
                        countColor[k] = countColor[k] < 0.0 ? 0.0 : countColor[k];
                    }
                } // processing mask

                for (int k = 0; k < 3; ++k) {
                    this->data[index + k] = countColor[k];
                }

            } // for j
        } // for i
    } // count

    delete[] tmpMask;
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
    unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
    radius(iradius), x(ix), y(iy), r(ir), g(ig), b(ib), a(ia)
{

}

int helper::Pair_color_compare(const void* a, const void* b) {
    std::pair<int, int>* c = (std::pair<int, int>*) a;
    std::pair<int, int>* d = (std::pair<int, int>*) b;
    if (c->second > d->second)
        return -1;
    else if (c->second < d->second)
        return 1;
    else
        return 0;
}