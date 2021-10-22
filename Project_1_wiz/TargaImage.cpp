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
#include <vector>
#include <algorithm>

using namespace std;

// constants
const int           RED = 0;                // red channel
const int           GREEN = 1;                // green channel
const int           BLUE = 2;                // blue channel
const unsigned char BACKGROUND[3] = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1; i <= s; i++)
        res = (n - i + 1) * res / i;

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
    unsigned char* rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (!data)
        return NULL;

    // Divide out the alpha
    for (i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        int out_offset = i * width * 3;

        for (j = 0; j < width; j++)
        {
            RGBA_To_RGB(data + (in_offset + j * 4), rgb + (out_offset + j * 3));
        }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char* filename)
{
    TargaImage* out_image = Reverse_Rows();

    if (!out_image)
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
TargaImage* TargaImage::Load_Image(char* filename)
{
    unsigned char* temp_data;
    TargaImage* temp_image;
    TargaImage* result;
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


int TargaImage::index(int x, int y, int width) {
    int index = x * 4 + width * 4 * y;
    return index;
}


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    int pixel_amount = width * height;  // 1 pixel has r,g,b,a
    float pixel_gray = 0;
    for (int i = 0; i < pixel_amount; i++) {
        pixel_gray = float(data[i * 4] * 0.3 + data[i * 4 + 1] * 0.59 + data[i * 4 + 2] * 0.11);
        data[i * 4] = pixel_gray;
        data[i * 4 + 1] = pixel_gray;
        data[i * 4 + 2] = pixel_gray;
    }
    //ClearToBlack();
    return false;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    int pixel_amount = width * height;  // 1 pixel has r,g,b,a
    for (int i = 0; i < pixel_amount; i++) {
        data[i * 4] = ((data[i * 4] + 1) >> 5) << 5;            //r
        data[i * 4 + 1] = ((data[i * 4 + 1] + 1) >> 5) << 5;     //g
        data[i * 4 + 2] = ((data[i * 4 + 2] + 1) >> 6) << 6;     //b

    }
    //ClearToBlack();
    return false;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool cmp2(pair<int, int>a, pair<int, int>b)
{
    return a.second < b.second;
}
bool TargaImage::Quant_Populosity()
{
    int pixel_amount = width * height;  // 1 pixel has r,g,b,a
    int reg = 0;
    unsigned char color_counter[32][32][32] = { 0 }; //r,g,b,cnt
    vector<pair<long int, int>> color_cnt;


    int reg_r, reg_g, reg_b;
    for (int i = 0; i < pixel_amount; i++) {
        reg_r = int((data[i * 4] + 1) / 8);
        reg_g = int((data[i * 4 + 1] + 1) / 8);
        reg_b = int((data[i * 4 + 2] + 1) / 8);
        color_counter[reg_r][reg_g][reg_b] += 1;
    }
    //convert to 8 bits per color , from 0 to 32 and count color

    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) {
            for (int k = 0; k < 32; k++) {
                if (color_counter[i][j][k] > 0) {
                    reg = i * 10000 + j * 100 + k;
                    color_cnt.push_back({ reg,color_counter[i][j][k] });
                }
            }
        }
    }

    sort(color_cnt.rbegin(), color_cnt.rend(), cmp2);
    int popular_color[256][3];
    for (int i = 0; i < 256; i++) {
        popular_color[i][0] = color_cnt[i].first / 10000 * 8;
        popular_color[i][1] = color_cnt[i].first / 100 % 100 * 8;
        popular_color[i][2] = color_cnt[i].first % 100 * 8;
        //cout << popular_color[i][0]<<" "<< popular_color[i][1]<<" "<< popular_color[i][2] << endl;
    }


    int min, min_index, r2, g2, b2;
    for (int i = 0; i < pixel_amount; i++) {

        r2 = pow((data[i * 4] - popular_color[0][0]), 2);
        g2 = pow((data[i * 4 + 1] - popular_color[0][1]), 2);
        b2 = pow((data[i * 4 + 2] - popular_color[0][2]), 2);

        min_index = 0;
        min = r2 + g2 + b2;

        for (int j = 1; j < 256; j++) {

            r2 = pow((data[i * 4] - popular_color[j][0]), 2);
            g2 = pow((data[i * 4 + 1] - popular_color[j][1]), 2);
            b2 = pow((data[i * 4 + 2] - popular_color[j][2]), 2);

            reg = r2 + g2 + b2;
            if (reg < min) {
                min = reg;
                min_index = j;
            }
        }
        //cout << min_index << " " << min << endl;
        data[i * 4] = popular_color[min_index][0];
        data[i * 4 + 1] = popular_color[min_index][1];
        data[i * 4 + 2] = popular_color[min_index][2];
    }


    //ClearToBlack();
    return false;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    To_Grayscale();
    int pixel_amount = width * height;  // 1 pixel has r,g,b,a
    for (int i = 0; i < pixel_amount; i++) {
        data[i * 4] = (data[i * 4] + 1) / 128 * 255;
        data[i * 4 + 1] = (data[i * 4 + 1] + 1) / 128 * 255;
        data[i * 4 + 2] = (data[i * 4 + 2] + 1) / 128 * 255;
    }
    //ClearToBlack();
    return false;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    int pixel_amount = width * height;
    float pixel_gray = 0;
    float max, min;
    max = 0.2;
    min = -0.2;
    srand(time(NULL));

    for (int i = 0; i < pixel_amount; i++) {
        pixel_gray = (data[i * 4] * 0.3 + data[i * 4 + 1] * 0.59 + data[i * 4 + 2] * 0.11) / 256;
        pixel_gray = pixel_gray + ((max - min) * rand() / (RAND_MAX + 1.0) + min);
        if (pixel_gray > 0.5) {
            data[i * 4] = 255;
            data[i * 4 + 1] = 255;
            data[i * 4 + 2] = 255;
        }
        else {
            data[i * 4] = 0;
            data[i * 4 + 1] = 0;
            data[i * 4 + 2] = 0;
        }
    }

    //ClearToBlack();
    return false;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Dither_FS()
{
    To_Grayscale();

    for (int i = 0; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            int pos = i * width * 4 + j * 4;
            int o_r = data[pos],
                o_g = data[pos + 1],
                o_b = data[pos + 2];

            int factor = 1;

            int n_r = (factor * o_r / 128) * (255 / factor),
                n_g = (factor * o_g / 128) * (255 / factor),
                n_b = (factor * o_b / 128) * (255 / factor);

            data[pos] = n_r;
            data[pos + 1] = n_g;
            data[pos + 2] = n_b;

            float err_r = o_r - n_r,
                err_g = o_g - n_g,
                err_b = o_b - n_b;

            //right
            data[index(j + 1, i, width)] += err_r * 7.0 / 16;
            data[index(j + 1, i, width) + 1] += err_g * 7.0 / 16;
            data[index(j + 1, i, width) + 2] += err_b * 7.0 / 16;
            //left bottom
            data[index(j - 1, i + 1, width)] += err_r * 3.0 / 16;
            data[index(j - 1, i + 1, width) + 1] += err_g * 3.0 / 16;
            data[index(j - 1, i + 1, width) + 2] += err_b * 3.0 / 16;
            //bottom
            data[index(j, i + 1, width)] += err_r * 5.0 / 16;
            data[index(j, i + 1, width) + 1] += err_g * 5.0 / 16;
            data[index(j, i + 1, width) + 2] += err_b * 5.0 / 16;
            //right bottom
            data[index(j + 1, i + 1, width)] += err_r * 1.0 / 16;
            data[index(j + 1, i + 1, width) + 1] += err_g * 1.0 / 16;
            data[index(j + 1, i + 1, width) + 2] += err_b * 1.0 / 16;



        }
    }
    //ClearToBlack();
    return false;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    int pixel_amount = width * height;  // 1 pixel has r,g,b,a
    int threshold;
    float pixel_gray = 0;
    float avg_gray = 0;
    vector<int> gray_data;

    for (int i = 0; i < pixel_amount; i++) {
        pixel_gray = data[i * 4] * 0.3 + data[i * 4 + 1] * 0.59 + data[i * 4 + 2] * 0.11;
        data[i * 4] = pixel_gray;
        avg_gray += pixel_gray;
        gray_data.push_back(pixel_gray);
    }

    sort(gray_data.begin(), gray_data.begin() + pixel_amount);
    avg_gray /= pixel_amount * 256;
    threshold = (1 - avg_gray) * pixel_amount;
    //cout << pixel_amount << endl;
    cout << gray_data[threshold] << endl;
    //cout << "avg:" << avg_gray << endl;


    for (int i = 0; i < pixel_amount; i++) {
        if (data[i * 4] > gray_data[threshold]) {
            data[i * 4] = 255;
            data[i * 4 + 1] = 255;
            data[i * 4 + 2] = 255;
        }
        else {
            data[i * 4] = 0;
            data[i * 4 + 1] = 0;
            data[i * 4 + 2] = 0;
        }
    }
    //ClearToBlack();
    return false;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    int pixel_amount = width * height;  // 1 pixel has r,g,b,a
    float pixel_gray, threshold = 0;

    float mask[4][4] = {
        { 0.7059, 0.3529, 0.5882, 0.2353 },
        { 0.0588, 0.9412, 0.8235, 0.4118 },
        { 0.4706, 0.7647, 0.8824, 0.1176 },
        { 0.1765, 0.5294, 0.2941, 0.6471 }
    };


    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            pixel_gray = (data[(i * width + j) * 4] * 0.3 + data[(i * width + j) * 4 + 1] * 0.59 + data[(i * width + j) * 4 + 2] * 0.11) / 256;
            threshold = mask[i % 4][j % 4];
            if (pixel_gray > threshold) {
                data[(i * width + j) * 4] = 255;
                data[(i * width + j) * 4 + 1] = 255;
                data[(i * width + j) * 4 + 2] = 255;
            }
            else {
                data[(i * width + j) * 4] = 0;
                data[(i * width + j) * 4 + 1] = 0;
                data[(i * width + j) * 4 + 2] = 0;
            }
        }
    }



    //ClearToBlack();
    return false;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
int closest(int color, int value) {
    int red[8] = { 0, 36, 73, 109, 146, 182, 219, 255 },
        green[8] = { 0, 36, 73, 109, 146, 182, 219, 255 },
        blue[4] = { 0, 85, 170, 255 };

    int reg = 999;
    int close_value = 999;

    if (color == 0) {
        for (int i = 0; i < 8; i++) {
            if (abs(value - red[i]) < reg) {
                reg = abs(value - red[i]);
                close_value = red[i];
            }
        }
    }
    else if (color == 1) {
        for (int i = 0; i < 8; i++) {
            if (abs(value - green[i]) < reg) {
                reg = abs(value - green[i]);
                close_value = green[i];
            }
        }
    }
    else {
        for (int i = 0; i < 4; i++) {
            if (abs(value - blue[i]) < reg) {
                reg = abs(value - blue[i]);
                close_value = blue[i];
            }
        }
    }
    //cout << close_value << " ";
    return close_value;
}

bool TargaImage::Dither_Color() {

    for (int i = 0; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            int pos = i * width * 4 + j * 4;
            int o_r = data[pos],
                o_g = data[pos + 1],
                o_b = data[pos + 2];

            int factor = 1;

            int n_r = closest(0, o_r),
                n_g = closest(1, o_g),
                n_b = closest(2, o_b);

            data[pos] = n_r;
            data[pos + 1] = n_g;
            data[pos + 2] = n_b;

            float err_r = o_r - n_r,
                err_g = o_g - n_g,
                err_b = o_b - n_b;

            //right
            data[index(j + 1, i, width)] += err_r * 7.0 / 16;
            data[index(j + 1, i, width) + 1] += err_g * 7.0 / 16;
            data[index(j + 1, i, width) + 2] += err_b * 7.0 / 16;
            //left bottom
            data[index(j - 1, i + 1, width)] += err_r * 3.0 / 16;
            data[index(j - 1, i + 1, width) + 1] += err_g * 3.0 / 16;
            data[index(j - 1, i + 1, width) + 2] += err_b * 3.0 / 16;
            //bottom
            data[index(j, i + 1, width)] += err_r * 5.0 / 16;
            data[index(j, i + 1, width) + 1] += err_g * 5.0 / 16;
            data[index(j, i + 1, width) + 2] += err_b * 5.0 / 16;
            //right bottom
            data[index(j + 1, i + 1, width)] += err_r * 1.0 / 16;
            data[index(j + 1, i + 1, width) + 1] += err_g * 1.0 / 16;
            data[index(j + 1, i + 1, width) + 2] += err_b * 1.0 / 16;



        }
        //ClearToBlack();

        for (int i = 0; i < height - 1; i++) {
            for (int j = 1; j < width - 1; j++) {
                if (data[i * width * 4 + j * 4 + 0] > 255)data[i * width * 4 + j * 4 + 0] = 255;
                else if (data[i * width * 4 + j * 4 + 0] < 0)data[i * width * 4 + j * 4 + 0] = 0;
            }
        }
    }
    return false;
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
        cout << "Comp_Over: Images not the same size\n";
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

    for (int i = 0; i < width * height * 4; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i + 1] = abs(rgb1[1] - rgb2[1]);
        data[i + 2] = abs(rgb1[2] - rgb2[2]);
        data[i + 3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box() {
    int s;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int c = 0; c < 3; c++) {
                s = 0;
                for (int k = -2; k < 3; k++) {
                    if (i + k<0 || i + k>height - 1)continue;
                    for (int l = -2; l < 3; l++) {
                        if (j + l<0 || j + l>width - 1)continue;
                        s += data[(i + k) * width * 4 + (j + l) * 4 + c];
                    }
                }
                data[i * width * 4 + j * 4 + c] = s / 25;
            }
        }
    }
    //ClearToBlack();
    return false;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    int arr[5][5] = { { 1,2,3,2,1}, {2,4,6,4,2},{3,6,9,6,3},{2,4,6,4,2},{1,2,3,2,1} };
    int s;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int c = 0; c < 3; c++) {
                s = 0;
                for (int k = -2; k < 3; k++) {
                    if (i + k<0 || i + k>height - 1)continue;
                    for (int l = -2; l < 3; l++) {
                        if (j + l<0 || j + l>width - 1)continue;
                        s += data[(i + k) * width * 4 + (j + l) * 4 + c] * arr[k + 2][l + 2];
                    }
                }
                data[i * width * 4 + j * 4 + c] = s / 81;
            }
        }
    }
    //ClearToBlack();
    return false;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    int arr[5][5] = { { 1,4,6,4,1}, {4,16,24,16,4},{6,24,36,24,6},{4,16,24,16,4},{ 1,4,6,4,1} };
    int s;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int c = 0; c < 3; c++) {
                s = 0;
                for (int k = -2; k < 3; k++) {
                    if (i + k<0 || i + k>height - 1)continue;
                    for (int l = -2; l < 3; l++) {
                        if (j + l<0 || j + l>width - 1)continue;
                        s += data[(i + k) * width * 4 + (j + l) * 4 + c] * arr[k + 2][l + 2];
                    }
                }
                data[i * width * 4 + j * 4 + c] = s / 256;
            }
        }
    }
    //ClearToBlack();
    return false;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////


bool TargaImage::Filter_Gaussian_N(unsigned int N)
{
    Filter_Gaussian_N_(data, N);
    ////build array
    //int** arr;
    //arr = new int* [N];
    //for (int i = 0; i < N; i++) {
    //    arr[i] = new int[N];
    //}

    //for (int i = 0; i < ceil(N/2)+1; i++) {
    //    arr[0][i] = (i == 0 ? 1 : (i + 1) * 2);
    //    arr[i][0] = (i == 0 ? 1 : (i + 1) * 2);
    //}

    //for (int i = ceil(N / 2); i < N; i++) {
    //    arr[0][i] = arr[0][N - 1 - i];
    //    arr[i][0] = arr[N - 1 - i][0];
    //}
    //for (int i = 1; i < N; i++) {
    //    for (int j = 1; j < N; j++) {
    //        arr[i][j] = arr[i][0] * arr[0][j];
    //    }
    //}


    ////print array & sum_weight
    //int weight = 0;
    //for (int i = 0; i < N; i++) {
    //    for (int j = 0; j < N; j++) {
    //        cout << arr[i][j] << " ";
    //        weight += arr[i][j];
    //    }
    //    cout << endl;
    //}

    ////start filter
    //int s;
    //for (int i = 0; i < height; i++) {
    //    for (int j = 0; j < width; j++) {
    //        for (int c = 0; c < 3; c++) {
    //            s = 0;
    //            for (int k = -ceil(N / 2); k < ceil(N / 2)+1; k++) {
    //                if (i + k<0 || i + k>height - 1)continue;
    //                for (int l = -ceil(N / 2); l < ceil(N / 2)+1; l++) {
    //                    if (j + l<0 || j + l>width - 1)continue;
    //                    int n_k, n_l;
    //                    n_k = k + ceil(N / 2);
    //                    n_l = l + ceil(N / 2);
    //                    s += data[(i + k) * width * 4 + (j + l) * 4 + c] * arr[n_k][n_l];
    //                }
    //            }
    //            data[i * width * 4 + j * 4 + c] = s / weight;
    //        }
    //    }
    //}



    ////delete array
    //for (int i = 0; i < N; i++)
    //{
    //    delete[] arr[i];
    //}
    //delete[] arr;


    //ClearToBlack();
    return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    unsigned char *rgb = To_RGB();
    int arr[5][5] = { { 1,4,6,4,1}, {4,16,24,16,4},{6,24,-220,24,6},{4,16,24,16,4},{ 1,4,6,4,1} };
    double s;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int c = 0; c < 3; c++) {
                s = 0;
                for (int k = -2; k < 3; k++) {
                    if (i + k<0 || i + k>height - 1)continue;
                    for (int l = -2; l < 3; l++) {
                        if (j + l<0 || j + l>width - 1)continue;
                        s += rgb[(i + k) * width * 3 + (j + l) * 3 + c] * -(arr[k + 2][l + 2]);
                    }
                }
                s =  s / 256;
                if (s < 0)s = 0;
                if (s > 255)s = 255;
                data[i * width * 4 + j * 4 + c] = s ;
            }
        }
    }
    //ClearToBlack();
    return false;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    unsigned char* rgb = To_RGB();
    int arr[5][5] = { { 1,4,6,4,1}, {4,16,24,16,4},{6,24,36,24,6},{4,16,24,16,4},{ 1,4,6,4,1} };
    double s;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int c = 0; c < 3; c++) {
                s = 0;
                for (int k = -2; k < 3; k++) {
                    if (i + k<0 || i + k>height - 1)continue;
                    for (int l = -2; l < 3; l++) {
                        if (j + l<0 || j + l>width - 1)continue;
                        s += rgb[(i + k) * width * 3 + (j + l) * 3 + c] * arr[k + 2][l + 2];
                    }
                }
                s = rgb[i * width * 3 + j * 3 + c] + (rgb[i * width * 3 + j * 3 + c] - s / 256);
                if (s < 0)s = 0;
                if (s > 255)s = 255;
                data[i * width * 4 + j * 4 + c] = s;
            }
        }
    }
    return false;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void  TargaImage::Filter_Gaussian_N_(unsigned char* img,int N) {

    //store img 
    unsigned char* od = new unsigned char[width * height*4];
    for (int i = 0; i < height * width * 4; i++)od[i] = img[i];

    //build array
    int** arr;
    arr = new int* [N];
    for (int i = 0; i < N; i++) {
        arr[i] = new int[N];
    }

    for (int i = 0; i < ceil(N / 2) + 1; i++) {
        arr[0][i] = (i == 0 ? 1 : (i + 1) * 2);
        arr[i][0] = (i == 0 ? 1 : (i + 1) * 2);
    }

    for (int i = ceil(N / 2); i < N; i++) {
        arr[0][i] = arr[0][N - 1 - i];
        arr[i][0] = arr[N - 1 - i][0];
    }
    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            arr[i][j] = arr[i][0] * arr[0][j];
        }
    }


    //print array & sum_weight
    int weight = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            //cout << arr[i][j] << " ";
            weight += arr[i][j];
        }
        //cout << endl;
    }

    //start filter
    int s;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int c = 0; c < 3; c++) {
                s = 0;
                for (int k = -ceil(N / 2); k < ceil(N / 2) + 1; k++) {
                    if (i + k<0 || i + k>height - 1)continue;
                    for (int l = -ceil(N / 2); l < ceil(N / 2) + 1; l++) {
                        if (j + l<0 || j + l>width - 1)continue;
                        int n_k, n_l;
                        n_k = k + ceil(N / 2);
                        n_l = l + ceil(N / 2);
                        s += od[(i + k) * width * 4 + (j + l) * 4 + c] * arr[n_k][n_l];
                    }
                }
                img[i * width * 4 + j * 4 + c] = s / weight;
            }
        }
    }



    //delete array
    for (int i = 0; i < N; i++)
    {
        delete[] arr[i];
    }
    delete[] arr;

    delete[] od;

}
bool TargaImage::NPR_Paint()
{
    int radius[3] = { 1,3,7 };
    unsigned char* rgb = To_RGB();
    unsigned char* referenceImage = new unsigned char[width * height * 4];

    //Stroke s(300, 200, 200, 255, 255, 255, 255);
    //Paint_Stroke(s);

    //ref = data
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < 4; k++) {
                referenceImage[i * width * 4 + j * 4 + k] = data[i * width * 4 + j * 4 + k];
            }
        }
    }

    //let data be canvas
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < 3; k++) {
                data[i * width * 4 + j * 4 + k] = 255;
            }
        }
    }


    for (int r = 2; r >= 0; r--) {
        //reset color image

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                for (int k = 0; k < 3; k++) {
                    referenceImage[i * width * 4 + j * 4 + k] = rgb[i * width * 3 + j * 3 + k];
                }
            }
        }
        Filter_Gaussian_N_(referenceImage, 2*radius[r]+1);
        NPR_Paint_Layer(data, referenceImage, radius[r]);
    }

    //ClearToBlack();
    return false;
}


bool TargaImage::NPR_Paint_Layer(unsigned char* canvas, unsigned char* referenceImage, unsigned int r) {
    vector<Stroke> S;
    //canvas has channel 3
    //referenceImage has channel 4
    vector<double> diff(height * width);

    double reg,reg2;
    for (int i = 0; i < height; i+=r) {
        for (int j = 0; j < width; j+=r) {
            reg = 0;
            for (int k = 0; k < 3; k++) {

                reg2 =data[i * width * 4 + j * 4 + k] - referenceImage[i * width * 4 + j * 4 + k];
                //cout << reg2 << "\t";
                reg += pow((double)reg2, 2);
            }
            reg = sqrt(reg);
            diff.at(i*width+j) = reg;
            //cout << diff[i * width + j] << endl;
        }
    }//calculate diff

    //Stroke(unsigned int radius, unsigned int x, unsigned int y,
    //    unsigned char r, unsigned char g, unsigned char b, unsigned char a);
    double range = floor((double)r/2);
    //cout << range << endl;
    double sum, areaError, max;
    int  max_h, max_w, cnt, index;

    int T = 25;
    int fal = 0;
    int suc = 0;
    for (int i = 0; i < height; i+=r) {
        for (int j = 0; j < width; j+=r) {
            //sum region error and pick max point 
            
            //initialize
            sum = 0;
            max = -1;
            max_h = 0;
            max_w = 0;
            cnt = 0;
            for (int k = -range; k <= range; k++){
                if (i + k<0 || i + k>height - 1)continue;
                for (int l = -range; l <= range; l++) {
                    if (j + l<0 || j + l>width - 1)continue;
                    index = (i + k) * width + (j + l);

                   // cout << "diff:" << diff.at(index) << "\t";
                    sum += diff.at(index);
                    //cout << diff[index] << " ";
                    cnt += 1;
                    if (diff.at(index)> max) {
                        max = diff.at(index);;
                        max_h = i + k;
                        max_w = j + l;
                    }
                }
            }
            //cout<<sum<<" "<<cnt << endl;

            areaError = sum ;
            
            if (areaError >T) {
                S.push_back(Stroke(r, max_w, max_h,
                    referenceImage[max_h * width * 4 + max_w * 4 + 0],
                    referenceImage[max_h * width * 4 + max_w * 4 + 1],
                    referenceImage[max_h * width * 4 + max_w * 4 + 2],
                    referenceImage[max_h * width * 4 + max_w * 4 + 3]));
                suc += 1;
            }
            else {
                fal += 1;
            }
        }
    }
    //for(int i=0;i<height*width;i++)cout << diff.at(i) << " " ;
    cout <<">T: "<< suc << "\telse: " << fal << "\n";
    random_shuffle(S.begin(), S.end());
    for (int i = 0; i < S.size(); i++) {
        Paint_Stroke(S[i]);
    }
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{

    Resize(0.5);
    //ClearToBlack();
    return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{

    Resize(2);
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////

float Bilinear_Interpolation(int q11, int q12, int q21, int q22, int x1, int x2, int y1, int y2, float x, float y) {
    float fr1 = (x2 - x) / (x2 - x1) * q11 + (x - x1) / (x2 - x1) * q21;    //interpolation to x
    float fr2 = (x2 - x) / (x2 - x1) * q12 + (x - x1) / (x2 - x1) * q22;    //interpolation to x
    float fp = (y2 - y) / (y2 - y1) * fr1 + (y - y1) / (y2 - y1) * fr2;     //interpolation to y

    return fp;
}
bool TargaImage::Resize(float scale)
{
    int x, y, x1, y1, x2, y2, q11, q12, q21, q22;

    int old_height = height,
        old_width = width;

    unsigned char *old_data = To_RGB();

    //now is new size
    height *= scale;
    width *= scale;
    //data = new unsigned char[height * width * 4];

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float x = (float)old_height / height * i,
                y = (float)old_width / width * j;
            int x1 = floor(x),
                x2 = (x1 >= old_height) ? x1 : (x1 + 1),
                y1 = floor(y),
                y2 = (y1 >= old_width) ? y1 : (y1 + 1);
            int oq11 = x1 * old_width * 3 + y1 * 3,
                oq12 = x1 * old_width * 3 + y2 * 3,
                oq21 = x2 * old_width * 3 + y1 * 3,
                oq22 = x2 * old_width * 3 + y2 * 3;

            data[i * width * 4 + j * 4 + 0] = Bilinear_Interpolation(old_data[oq11 + 0], old_data[oq12 + 0], old_data[oq21 + 0], old_data[oq22 + 0], x1, x2, y1, y2, x, y);
            data[i * width * 4 + j * 4 + 1] = Bilinear_Interpolation(old_data[oq11 + 1], old_data[oq12 + 1], old_data[oq21 + 1], old_data[oq22 + 1], x1, x2, y1, y2, x, y);
            data[i * width * 4 + j * 4 + 2] = Bilinear_Interpolation(old_data[oq11 + 2], old_data[oq12 + 2], old_data[oq21 + 2], old_data[oq22 + 2], x1, x2, y1, y2, x, y);
            
        }
    }
    delete[] old_data;
    //ClearToBlack();
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
    unsigned char* rgb = To_RGB();
    float theta = (angleDegrees * c_pi / 180);

    int x, y;

    double vcos = cos(theta);
    double vsin = sin(theta);

    int x_0 = (width - 1) / 2;
    int y_0 = (height - 1) / 2;

    int nx2, ny2;

    for (int ny = 0; ny < height; ny++) {
        for (int nx = 0; nx < width; nx++) {

            nx2 = nx - x_0;
            ny2 = ny - y_0;

            x = (int)(nx2 * vcos + ny2 * vsin  + x_0);
            y = (int)(-nx2 * vsin + ny2 * vcos  + y_0);

            if (y >= 0 && y < height && x >= 0 && x < width) {
                data[ny * 4 * width + nx * 4 + 0] = rgb[y * 3 * width + x * 3 + 0];
                data[ny * 4 * width + nx * 4 + 1] = rgb[y * 3 * width + x * 3 + 1];
                data[ny * 4 * width + nx * 4 + 2] = rgb[y * 3 * width + x * 3 + 2];
            }
            else {
                data[ny * 4 * width + nx * 4 + 0] = data[ny * 4 * width + nx * 4 + 1] = data[ny * 4 * width + nx * 4 + 2] = 0;
            }
        }
    }


    delete[] rgb;
    //ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char* rgba, unsigned char* rgb)
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

        for (i = 0; i < 3; i++)
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
}// RGBA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char* dest = new unsigned char[width * height * 4];
    TargaImage* result;
    int 	        i, j;

    if (!data)
        return NULL;

    for (i = 0; i < height; i++)
    {
        int in_offset = (height - i - 1) * width * 4;
        int out_offset = i * width * 4;

        for (j = 0; j < width; j++)
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
                }
                else if (dist_squared == radius_squared + 1) {
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

