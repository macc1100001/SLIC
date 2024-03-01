#ifndef __COUNTNEIGHBOURS__
#define __COUNTNEIGHBOURS__

#include <opencv2/opencv.hpp>

using namespace cv;

//Count the number of adjacent pixels in 4- and 8- neighbourhood that are different. Those cells outside the image, count as diferent, so in 4 neighbourhood, corner pixels start counting on 2, lateral pixels on 1 and all the other on 0. in 8-neighbourhood, lateral and corner pixels start counting on 3 and all the rest on 0.

template <typename X>
int countNeighbours_4(Mat &M, int y, int x)
{
    X val, *apul, cont;
    int sz  = M.step/M.elemSize();

    assert (M.cols > 2 && M.rows > 2);

    val = M.at<X>(y,x);

    if (!y) 
    {
        if (!x)
        {
            cont = 2;
            apul = M.ptr<X>(0) + 1;
            if (*apul != val) cont++;
            apul+= sz -1;
            if (*apul != val) cont++;
            return cont;
        }
        if (x == sz-1)
        {
            cont = 2;
            apul = M.ptr<X>(0) + x - 1;
            if (*apul != val) cont++;
            apul+= sz + 1;
            if (*apul != val) cont++;
            return cont;
        }
        cont = 1;
        apul = M.ptr<X>(0) + x - 1;
        if (*apul != val) cont++;
        apul += 2;
        if (*apul != val) cont++;
        apul += sz - 1;
        if (*apul != val) cont++;
        return cont;
    }
    if (y == M.rows - 1)
    {
        if (!x)
        {
            cont = 2;
            apul = M.ptr<X>(y-1);
            if (*apul != val) cont++;
            apul += sz + 1;
            if (*apul != val) cont++;
            return cont;
        }
        if ( x == sz - 1)
        {
            cont = 2;
            apul = M.ptr<X>(y - 1) + x;
            if (*apul != val) cont++;
            apul += sz - 1;
            if (*apul != val) cont++;
            return cont;
    
        }
        cont = 1;
        apul = M.ptr<X>(y-1) + x;
        if (*apul != val) cont++;
        apul += sz - 1;
        if (*apul != val) cont++;
        apul += 2;
        if (*apul != val) cont++;
        return cont; 
   }
   if (!x)
   {
        cont = 1;
        apul = M.ptr<X>(y-1);
        if (*apul != val) cont++;
        apul += sz + 1;
        if (*apul != val) cont++;
        apul += sz - 1;
        if (*apul != val) cont++;
        return cont;
   }
   if (x==sz-1)
   {
        cont = 1;
        apul = M.ptr<X>(y-1) + x;
        if (*apul != val) cont++;
        apul += sz - 1;
        if (*apul != val) cont++;
        apul += sz + 1;
        if (*apul != val) cont++;
        return cont;
   }

   cont = 0;
   apul = M.ptr<X>(y) + x - sz;
   if (*apul != val) cont++;
   apul += sz - 1;
   if (*apul != val) cont++;
   apul += 2;
   if (*apul != val) cont++;
   apul += sz - 1;
   if (*apul != val) cont++;

   return cont;
}

 
template <typename X>
int countNeighbours_8(Mat &M, int y, int x)
{
    X val, *apul, cont = 3;
    int sz  = M.step/M.elemSize();

    assert (M.cols > 2 && M.rows > 2);

    val = M.at<X>(y,x);

    if (!y) 
    {
        if (!x)
        {
            apul = M.ptr<X>(0) + 1;
            if (*apul != val) cont++;
            apul+= sz -1;
            if (*apul != val) cont++;
            apul++;
            if (*apul != val) cont++;
            return cont;
        }
        if (x == sz-1)
        {
            cont = 3;
            apul = M.ptr<X>(0) + x - 1;
            if (*apul != val) cont++;
            apul+= sz;
            if (*apul != val) cont++;
            apul++;
            if (*apul != val) cont++;
            return cont;
        }
        apul = M.ptr<X>(0) + x - 1;
        if (*apul != val) cont++;
        apul+=2;
        if (*apul != val) cont++;
        apul += sz - 2;
        if (*apul != val) cont++;
        apul++;
        if (*apul != val) cont++;
        apul++;
        if (*apul != val) cont++;
        return cont;
    }
    if (y == M.rows-1)
    {
        if (!x)
        {
            apul = M.ptr<X>(y-1);
            if (*apul != val) cont++;
            apul++;
            if (*apul != val) cont++;
            apul += sz;
            if (*apul != val) cont++;
            return cont;
        }
        if ( x == sz-1)
        {
            apul = M.ptr<X>(y-1) + x - 1;
            if (*apul != val) cont++;
            apul++;
            if (*apul != val) cont++;
            apul += sz - 1;
            if (*apul != val) cont++;
            return cont;
    
        }
        apul = M.ptr<X>(y-1) + x - 1;
        if (*apul != val) cont++;
        apul++;
        if (*apul != val) cont++;
        apul++;
        if (*apul != val) cont++;
        apul += sz - 2;
        if (*apul != val) cont++;
        apul+=2;
        if (*apul != val) cont++;
        return cont; 
   }
   if (!x)
   {
        apul = M.ptr<X>(y-1);
        if (*apul != val) cont++;
        apul++;
        if (*apul != val) cont++;
        apul += sz;
        if (*apul != val) cont++;
        apul += sz - 1;
        if (*apul != val) cont++;
        apul++;
        if (*apul != val) cont++;
        return cont;
   }
   if (x == sz - 1)
   {
        apul = M.ptr<X>(y-1) + x - 1;
        if (*apul != val) cont++;
        apul++;
        if (*apul != val) cont++;
        apul += sz - 1;
        if (*apul != val) cont++;
        apul += sz;
        if (*apul != val) cont++;
        apul++;
        if (*apul != val) cont++;
        return cont;
   }

   cont = 0;
   apul = M.ptr<X>(y) + x - sz - 1;
   if (*apul != val) cont++;
   apul++;
   if (*apul != val) cont++;
   apul++;
   if (*apul != val) cont++;
   apul += sz - 2;
   if (*apul != val) cont++;
   apul += 2;
   if (*apul != val) cont++;
   apul += sz - 2;
   if (*apul != val) cont++;
   apul++;
   if (*apul != val) cont++;
   apul++;
   if (*apul != val) cont++;

   return cont;
}
#endif
