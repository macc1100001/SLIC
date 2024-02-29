#include "opencv2/opencv.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <fsort.h>
#include <iostream>
#include <countNeighbours.h>

using namespace std;
using namespace cv;


int main(int argc, char **argv)
{
    Mat M;

    M = imread("testimage.pgm", 0);

    for (int i=0;i<M.rows;++i)
        for (int j=0;j<M.cols;++j)
        {
            if (M.at<unsigned char>(i,j)!=0)
            {
            
            cout << (int)M.at<unsigned char>(i,j) << " | 8-vecinos & 4-vecinos@(" << i << "," << j << ") = [" << countNeighbours_8<unsigned char>(M,i,j) << ", " << countNeighbours_4<unsigned char>(M,i,j) << "]" << endl;
            }
        }
    

    return 0;
}
