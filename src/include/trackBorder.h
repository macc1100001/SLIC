#ifndef __TRACKBORDER__
#define __TRACKBORDER__

#include <opencv2/opencv.hpp>

using namespace std;

//Estructura de pila utilizado para almacenar logar que falta visitar en el
//algoritmo de etiquetación. La pila no es dinámica, sino estática por razones
//de velocidad.
struct stack
{
    struct coor 
    {
        int x, y;
    };

    coor *S;
    int sz, idx;
    stack (int t)
    {
        sz = t;
        idx = 0;
        if (sz > 0)
            S = new coor[sz];
        else
            S = NULL;
    }
    ~stack()
    {
       if (S)
        delete[] S;
    }
    bool push (int x, int y)
    {
        if (idx < sz)
        {
            S[idx].x = x;
            S[idx].y = y;
            idx++;
            return true;
        }
        return false;
    }
    bool pop(int &x, int &y)
    {
        if (idx > 0)
        {
            idx--;
            x = S[idx].x;
            y = S[idx].y;
            return true;
        }
        return false;
    }
    bool empty()
    {
        return !idx;
    }
};


template <typename X>
bool isBorder(Mat &I, Mat&O, int y, int x)
{
   X val, *apul;
   uchar *apuo;
   int sz;

   sz = I.step/I.elemSize();

   
   val = I.at<X>(y,x);
   apul = I.ptr<X>(y) + x - sz - 1;
   apuo = O.ptr<uchar>(y) + x - sz - 1;

    if (!y || !x || y >= I.rows - 1 || x >= I.cols - 1)
        return true;


   if (*apul != val && *apuo == 0) return true;
   apul++;
   apuo++;
   if (*apul != val && *apuo == 0) return true;
   apul++;
   apuo++;
   if (*apul != val && *apuo == 0) return true;
   apul += sz - 2;
   apuo += sz - 2;
   if (*apul != val && *apuo == 0) return true;
   apul += 2;
   apuo += 2;
   if (*apul != val && *apuo == 0) return true;
   apul += sz - 2;
   apuo += sz - 2;
   if (*apul != val && *apuo == 0) return true;
   apul++;
   apuo++;
   if (*apul != val && *apuo == 0) return true;
   apul++;
   apuo++;
   if (*apul != val && *apuo == 0) return true;

   return false;
}



//Algoritmo no recursivo de etiquetación.
template <typename X>
void trackBorder(Mat &I, Mat &O, X val)
{
    int i, j, k, l, x, y, xl, yk;

    stack S(I.rows*I.cols);

    x = y = 0;
    for (i=0;i<I.rows;++i)
    {
        for (j=0;j<I.cols;++j)
        {

            if (I.at<X>(i, j) == val &&  isBorder<X>(I, O, i, j) && O.at<uchar>(i,j) == 0)
            {
                O.at<uchar>(i, j) = 1;
                S.push(j, i);
                do
                {
                    S.pop(x, y);

                    for (k=-1;k<2;++k)
                    {
                        yk = y + k;
                        if (yk>=0 && yk < I.rows)
                            for (l=-1;l<2;++l)
                            {
                                if (!k && !l)
                                    continue;
                                xl = x + l;
                                if (xl>=0 && xl < I.cols)
                                    if (I.at<X>(yk, xl) == val &&  isBorder<X>(I, O, yk, xl) && O.at<uchar>(yk, xl) == 0 )
                                    {
                                        O.at<uchar>(yk, xl) = 1;
                                        S.push(xl, yk);
                                    }
                            }
                    }
                } while (!S.empty());
            }
        }
    }
}

#endif
