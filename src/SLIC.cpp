#include <cstdlib> // EXIT RETURN MESSAGES
#include <unistd.h> // getopt
//#include <stdio.h>	// maybe not needed if cout is used, change later

#include <opencv2/opencv.hpp>
//#include <opencv2/highgui/highgui.hpp> // not needed, included by Mosaic.h

#include <fsort.h>
#include <iostream>
#include <vector>
#include <cstdio>
#include <Mosaic.h>
#include <sys/timeb.h>
#include <vector>
#include <countNeighbours.h>
#include <trackBorder.h>

/*
	TODO:
	valgrind shows multiple errors about memory	
	FIX IT!
*/

using namespace std;
using namespace cv;

typedef enum {None=0, Simple, Mean} borderType;

bool verboseFlag = false;

/**
\fn void printMaxMin(Mat &M)
\brief This function print a matrix channels minimum and maximum values.
\param Mat &M The matris where the values are obtained.
**/
void printMaxMin(Mat &M)
{
    vector <Mat> Channels;

    split(M, Channels);
    for (unsigned int i=0;i<Channels.size();++i)
    {
        double max, min;

        minMaxIdx(Channels[i], &min, &max);
        cout << "Los valores minimos y maximos del canal " << i << " son [" << min << ", " << max << "] respectivamente\n" << endl; 
    }
}


/**
\struct Clusters
\brief A public class which stores the clusters states 
*/

struct Clusters 
{
    int S; ///< The grid separation space.

    int K;  ///< The number of superpixels (clases)

    float mS2;///< Keeps the square ratio of the weigth used to relative importance between the distance in the color space and the distance in the image. 

    Size szI, szC; ///< The size of the image to be processed.

    Mat colorCoors; ///< The (K,1) matrix where the cluster information is stored. Each element has 5 values needed to store the mean color (3), position(2) values that characterize each cluster.

    Mat clusterInfo; ///< The (K,1), that stores the cluster geometric information; each element has 2 members, area and perimeter.
    
    Mat distances; ///< Matrix that stores the minimum distance from an image pixel to a cluster.
    Mat labels; ///< Matrix that stores a the cluster number asociated with an image pixel.
    Mat __acum, __cont;

    /*!
    \brief Class constructor.
    \param Mat &I The matrix that contains the image to be processed.
    \param int s The grid separation space
    \param m Weight that control the contribution of the spatial and color distances in distance used.
    */
    Clusters (Mat &I, int s, float m)
    {
        S = s;
        mS2 = m*m/(S*S);
        szI = Size(I.cols, I.rows);
        szC = Size(szI.width/S, szI.height/S);
        K = szC.width * szC.height;

        colorCoors = Mat_<Vec<float, 5> >(Size(K, 1));
        clusterInfo= Mat_<Vec<int, 2> >(Size(K, 1));
        distances = INFINITY * Mat::ones(szI, CV_32FC1);
        labels = -Mat::ones(szI, CV_32SC1);
        __acum = Mat_<Vec<float, 5> >(Size(K, 1));
        __cont = Mat::zeros(Size(K,1), CV_32SC1);
        initPixels(I);
    }
    /*!
    \fn void initPixels(Mat &I)
    \brief This function initialize the super-pixels (i.e. clusters stored in colorCoors.
    \param Mat &I The matrix that contains the image to be processed.

    */
    void initPixels(Mat &I)
    {
        int i, j, k, l;
        int nC, nR, x, y;
        Vec<float, 5> *apuC;

        nC = S/2 + (szI.width % S)/2;
        nR = S/2 + (szI.height % S)/2;
        apuC = colorCoors.ptr<Vec<float, 5> >(0);
        for (l = 0, i = nR; i < szI.height; i += S, ++l)
        {
            for (k = 0, j = nC ;j < szI.width && k < K; j += S, k++, apuC++)
            {
                x = j;
                y = i;
                getCoorMinGradient(I, x, y);
                getMeanColor(I,x, y, (*apuC)[0], (*apuC)[1], (*apuC)[2]);
                (*apuC)[3] = (float)y; 
                (*apuC)[4] = (float)x; 
            }
        }
    }

    /*!
    \fn void getMeanColor (Mat &I, int x, int y, float &l, float &a, float &b)
    \brief Find the mean color vector in the S^2 region around coors x,y.
    \param Mat &I The matrix that contains the image to be processed.
    \param int x, y The region center coordinates
    \param float &lm, &a, &b the returned CIE-lab color vector components.
    */
    void getMeanColor (Mat &I, int x, int y, float &l, float &a, float &b)
    {
        Vec3f acum, *apui;
        int i, j, N, colLeft, colRight, rowUp, rowBottom;

        setRegionLimits(x, y, S/2, szI, rowUp, rowBottom, colLeft, colRight);

        for (i=rowUp; i<rowBottom;++i)
        {
            apui = I.ptr<Vec3f>(i);
            acum *= 0.0;
            for (j = colLeft; j < colRight; ++j, apui++)
                acum += *apui;
        }
        N = (rowBottom-rowUp)*(colRight-colLeft);
        if (N)
            acum *= 1.0 / N;
       l = acum[0];
       a = acum[1];
       b = acum[2];
    }
    
    /*!
    \fn void getCoorMinGradient(Mat &I, int &x, int &y)
    \brief find the minimum gradient value in the 3x3 neighbourhoor around coordinate x,y of image I.
    \param Mat &I The matrix that contains the image to be processed.
    \param int &x, &y the 3x3 region coordinat center. The minimum gradient coordinate on the way out.
    */
    void getCoorMinGradient(Mat &I, int &x, int &y)
    {
        int i, j, xmin, ymin, sz1;
        float g, gmin, G[9], *apug;
        Vec3f *apu, *end, dx2, dy2; 

        if (x <= 0 || x + 1 >= I.cols  || y <= 0 || y + 1 >= I.rows)
            return;

        sz1  = I.step/I.elemSize();
        apu = I.ptr<Vec3f>(y -1) + x - 1;

        dx2 = *(apu-1)-*(apu+1);
        dy2 = *(apu-sz1)-*(apu+sz1);
        gmin = dx2.dot(dx2) + dy2.dot(dy2);
        xmin = x - 1;
        ymin = y - 1;
        apug = G;
        for (i = y-1;i < y + 2; ++i)
        {
            apu = I.ptr<Vec3f>(i) + x - 1;
            end = apu + 3;
            for (j=x-1;apu < end;++apu, ++j, ++apug)
            {
                 dx2 = *(apu-1)-*(apu+1);
                 dy2 = *(apu-sz1)-*(apu+sz1);
                 *apug = g = dx2.dot(dx2) + dy2.dot(dy2);
                 if (g < gmin)
                 {
                    gmin = g;
                    ymin = i;
                    xmin = j;
                }
            }
        }
        x = xmin;
        y = ymin;
    }

    //ms2 should be m^2/S^2;
    /*!
    \fn float Distance(Vec<float, 5> &C, Vec3f &lab, int y, int x)
    \brief Computes the weighted squared distance between a cluster and CIE-lab color vector and a coordinate. Namely \f$D^2=d_c^2+m^2\frac{d_s}{S}^2\f$, where \f$d_c\f$ is the color-space distance, \f$ds\f$ is the image spatial distance, \f$S\f$ is the grid step size\$, and \f$m\f$ that determines the contribution between the color and spatial distances.
    \param Vec<float, 5> &C A reference to a coordinate.
    \param Vec3f &lab A CIElab vector.
    \param int x, y a coordinate.
    */
    float Distance(const Vec<float, 5> &C, const Vec3f &lab, int y, int x)
    {
        float tmp, acum, acum2;

        tmp = lab[0]-C[0];
        acum = tmp*tmp;
        tmp = lab[1]-C[1];
        acum += tmp*tmp;
        tmp = lab[2]-C[2];
        acum += tmp*tmp;

        tmp = y-C[3];
        acum2 = tmp*tmp;
        tmp = x -C[4];
        acum2 += tmp*tmp;
        acum2 *=mS2;
        return acum + acum2;
    }

    void segment(Mat &I, float threshold = 0.0001)
    {
        int  j, k, l, cont, S2, colLeft, colRight, rowUp, rowBottom;
        float *apud, D, E;
        int *apul;
        Vec<float, 5> *apuC;
        Vec<int, 2> *apuI;
        Vec3f *apui;
        
        S2 = 2 * S;
        cont = 0;
        do
        {
            apuC = colorCoors.ptr<Vec<float, 5> >(0);
            for (j = 0; j < K; ++j, apuC++)
            {
                setRegionLimits((int)(*apuC)[4], (int)(*apuC)[3], S2, szI, rowUp, rowBottom, colLeft, colRight);
                for (k = rowUp; k < rowBottom; ++k)
                {
                    apud = distances.ptr<float>(k) + colLeft;
                    apul = labels.ptr<int>(k) + colLeft;
                    apui = I.ptr<Vec3f>(k) + colLeft;
                    for (l = colLeft; l < colRight; ++l, apud++, apul++, apui++)
                    {
                        D = Distance (*apuC, *apui, k, l);
                        if (D < *apud)
                        {
                            *apud = D;
                            *apul = j;
                        }
                    }
                }
            }

            //Compute new Cluster Centers.
            E = computeCentroids(I);
            if(verboseFlag)
            	cout << "Iteration: " << cont << " ResidualError: " << E << endl;
            cont++;
        } while (E > threshold && cont < 100);
        reConnect();

        //Compute geometric information.
        memset ((void *)clusterInfo.ptr<Vec<int, 2> >(0), 0, clusterInfo.step);
        for (j=0;j<I.rows;++j)
        {
            apul = labels.ptr<int>(j);
            for (k=0;k<I.cols; ++k, ++apul)
            {
                apuI = clusterInfo.ptr<Vec<int, 2> >(0) + *apul;
                (*apuI)[0]++;
                (*apuI)[1] += countNeighbours_4<int>(labels, j, k);
            }
        }
    }

    void printInfoClusters()
    {
        int i;
        Vec<float, 5> *apuC;
        Vec<int, 2> *apuI;

        for (i=0;i<K;++i)
        {
            apuC = colorCoors.ptr<Vec<float, 5> >(0)+ i;
            apuI = clusterInfo.ptr<Vec<int, 2> >(0)+ i;
            cout << "Cluster " << i << ": ["
                 << (*apuC)[3] << ", "
                 << (*apuC)[4] << "; "
                 << (*apuI)[0] << ", "
                 << (*apuI)[1] << "]" << endl;
        }
        //cout << endl;
    }

    float computeCentroids(Mat &I)
    {
        int i, j;
        Vec<float, 5> *apuA, *apuB, tmp;
        Vec3f *apui;
        float E, *apuv, *apu;
        int *apuC, *apul;

        memset ((void *)__acum.ptr<Vec<float, 5> >(0), 0, __acum.step);
        memset ((void *)__cont.ptr<int>(0), 0, __cont.step);

        apuA = __acum.ptr<Vec<float, 5> >(0);
        apuC = __cont.ptr<int>(0);
        for (i=0;i<I.rows;++i)
        {
            apul = labels.ptr<int>(i);
            apui = I.ptr<Vec3f>(i);
            for (j = 0; j < I.cols; ++j, apul++, apui++)
            {
                apuA = __acum.ptr<Vec<float, 5> >(0)+*apul;
                apu = (float *)apui;
                apuv = (float *)apuA;
                *apuv += *apu;
                apuv++; apu++;
                *apuv += *apu;
                apuv++; apu++;
                *apuv += *apu;
                apuv++;
                *apuv += i;
                apuv++;
                *apuv += j;
                (*(apuC+*apul))++;
            }
        }
        apuA = __acum.ptr<Vec<float, 5> > (0);
        apuB = colorCoors.ptr<Vec<float, 5> > (0);
        apuC = __cont.ptr<int> (0);
        E = 0.0;
        for (i=0;i<K;++i, apuA++, apuB++, ++apuC)
        {
            if (*apuC)
            {
                *apuA /= *apuC;
                tmp = *apuA - *apuB;
                E += sqrt(tmp.dot(tmp));
            }
            else
                *apuA *= 0;
        }
        E /= K;
       __acum.copyTo(colorCoors); 
        return E;
    }

    void drawImage(Mat &I)
    {
	    int i, j;
        Vec <float, 5> *apuC;
        int *apul;
        float *apuf, *apuc;
        Vec3f *apui;
        
        I = Mat::zeros(szI, CV_32FC3);
        apuC = colorCoors.ptr<Vec<float, 5> >(0);
        for (i = 0; i < szI.height; ++i)
        {
            apul = labels.ptr<int>(i);
            apui = I.ptr<Vec3f>(i);
            for (j = 0; j < szI.width; ++j, apui++, apul++)
            {
                apuf = (float *)apui;
                apuc = (float *)(apuC+*apul); 
                *apuf = *apuc;
                apuf++; apuc++;
                *apuf = *apuc;
                apuf++; apuc++;
                *apuf = *apuc;
                apuf++; apuc++;
            }
        }
    }

    void drawBorders (Mat &I)
    {
        int i, j,  *apus;
        float *apuf;
        Vec<int, 2> *apuI;
        Vec3f *apui;
        Mat srt, idx, Brds;
        unsigned char *apub;

        //Compute an index on the clusters sorted by perimeter.
        srt = Mat_<int>(Size(K,1));
        Brds = Mat::zeros(szI, CV_8UC1);
        apuI = clusterInfo.ptr<Vec<int, 2> >(0);
        apus = srt.ptr<int>(0);
        for (i=0;i<K;++i, ++apus, ++apuI)
            *apus = (*apuI)[1];
        sortIdx (srt, idx, SORT_EVERY_ROW | SORT_DESCENDING);

        for (i=0;i<K;++i)
            trackBorder(labels, Brds, idx.at<int>(0,i));
        for (i=0;i<I.rows;++i)
        {
            apui = I.ptr<Vec3f>(i);
            apub = Brds.ptr<unsigned char>(i);
            for (j=0;j<I.cols;++j, ++apub, ++apui)
            {
                if (*apub)
                {
                    apuf = (float *)apui;
                    *apuf = 255.; apuf++;
                    *apuf = 255.; apuf++;
                    *apuf = 255.;
                }
            }
        }
    }

    void reConnect()
    {
        int i, j, k, sz;
        int *apul, *apun;
        int cont, ngh;
        int lbl[8], d, D, idx, orphIdx;

        sz  = labels.step/labels.elemSize();
        orphIdx = 0;
        for (i = 1; i < szI.height-1; ++i)
        {
            apul = labels.ptr<int>(i);
            for (j = 1; j < szI.width-1; ++j, apul++)
            {
                memset(lbl, 0, 8*sizeof(int));
                ngh = cont = 0;
                apun = apul - sz - 1;
                if (*apun == *apul)
                    ngh++;
                else
                    lbl[cont++] = *apun;
                apun++;
                if (*apun == *apul)
                    ngh++;
                else
                    lbl[cont++] = *apun;
                apun++;
                if (*apun == *apul)
                    ngh++;
                else
                    lbl[cont++] = *apun;
                if (*apul == *(apul-1))
                    ngh++;
                else
                    lbl[cont++] = *apun;
                if (*apul == *(apul+1))
                    ngh++;
                else
                    lbl[cont++] = *apun;
                apun = apul + sz - 1;
                if (*apun == *apul)
                    ngh++;
                else
                    lbl[cont++] = *apun;
                apun++;
                if (*apun == *apul)
                    ngh++;
                else
                    lbl[cont++] = *apun;
                apun++;
                if (*apun == *apul)
                    ngh++;
                else
                    lbl[cont++] = *apun;

                if (!ngh) //This is an orhpaned pixel.
                {
                    
                    switch (cont)
                    {
                        case 2: Sort2(lbl);
                                break;
                        case 3: Sort3(lbl);
                                break;
                        case 4: Sort4(lbl);
                                break;
                        case 5: Sort5(lbl);
                                break;
                        case 6: Sort6(lbl);
                                break;
                        case 7: Sort7(lbl);
                                break;
                    };
                    D = *(__cont.ptr<int>(0)+lbl[0]);
                    idx = 0;
                    for (k=1;k<cont;++k)
                        if (lbl[k] != lbl[k-1])
                        {
                            d = *(__cont.ptr<int>(0)+lbl[k]);
                            if (D > d)
                            {
                                D = d;
                                idx = k;
                            }
                         }
                    *apul = lbl[idx];
                    
                            
                    orphIdx++;
                    
                }
            }
        }
        if(verboseFlag)
        	cout << orphIdx  << " orphaned pixel were found" << endl;
    }

    void setRegionLimits(int x, int y, int S, Size sz,  int &rowUp, int &rowBottom, int &colLeft, int &colRight)
    {
        
        colLeft = x - S  >= 0 ? x - S : 0; 
        rowUp = y - S  >= 0 ? y - S : 0; 
        colRight = x + S <= sz.width ? x + S : sz.width;
        rowBottom = y + S <= sz.height ? y + S : sz.height;
    }
};

void uso(char* programa){
	fprintf(stderr, "Uso: %s [ -s VALOR ] [ -m VALOR ] [ -S / -M ] Archivo\n\
	-s VALOR: Separacion del grid\n\
	-m VALOR: Valor de ponderacion entre la distancia en el espacio de color y la distancia en la imagen\n\
	-S: Se dibuja el borde usando un borde simple.\n\
	-M: Se dibuja el borde usando el promedio.\n\
	-v Verbose\n\
	-h Muestra esta lista de ayuda.\n\
	Valores por defecto: s = 10, m = 8\n", programa);
}

int main(int argc, char **argv)
{
    Mat frame, fFrame, gFrame, labFrame, qframe, slicFrame;
    Mat Out;
    double iFact = 1. / 255.;
    int S = 10;
    float m = 8;
    Mosaic M;
    borderType border = None;
    
    
    int c;
    while((c = getopt(argc, argv, "s:m:vSMh")) != -1){
    	switch(c){
    		case 's':
    			S = atoi(optarg);
    			break;
    		case 'm':
    			m = strtof(optarg, NULL);
    			break;
			case 'S':
				border = Simple;
				break;
			case 'M':
				border = Mean;
				break;
    		case 'v':
    			verboseFlag = true;
    			break;
			case 'h':
    		default:    	
				uso(argv[0]);
				exit(EXIT_FAILURE);
    	}
    }
    if(!argv[optind]){
    	fprintf(stderr, "¡Se esperaba un nombre de archivo!\n");
    	uso(argv[0]);
    	exit(EXIT_FAILURE);
    }
    
    frame = imread(argv[optind], 1);
    frame.convertTo (fFrame, CV_32FC3);
    //Es necesario normalizar la image BGR al intervalo [0,1] antes de convertir a espacio CIE Lab; en este caso iFact = 1./255
    fFrame *= iFact;

    cvtColor (fFrame, labFrame, COLOR_BGR2Lab); // TODO: valgrind shows here
	if(verboseFlag)
	    printMaxMin(labFrame);
    
    Clusters superPixels(labFrame, S, m);

    superPixels.segment(labFrame);
    superPixels.drawImage(slicFrame); // TODO: valgrind shows error here
    switch(border)
    {
        case Mean:
            cvtColor (slicFrame, Out, COLOR_Lab2BGR);
            Out *= 255.0;
            superPixels.drawBorders(Out);
            Out.convertTo (Out, CV_8UC3);
            break;
        case Simple:
            fFrame *= 255.0;
            superPixels.drawBorders(fFrame);
            fFrame.convertTo (Out, CV_8UC3);
            break;
        case None:
            cvtColor (slicFrame, Out, COLOR_Lab2BGR);
            Out *= 255.0;
            Out.convertTo (Out, CV_8UC3);
            break;
    }
	if(verboseFlag)
    	superPixels.printInfoClusters();


   if (frame.cols > frame.rows)
        M.init(Size(frame.cols, frame.rows), 2, 1, 8, 8, CV_8UC3);
    else
        M.init(Size(frame.cols, frame.rows), 1, 2, 8, 8, CV_8UC3);

    M.init(Size(frame.cols, frame.rows), 1, 1, 8, 8, CV_8UC3);

    //Abrimos las ventanas para mostrar los resultados.
    namedWindow( "Mosaico", 1 );

    //Despliega la imagen capturada en una ventana y conviertela a 
    //una respresentacion de dobles.
    M.setFigure(Out, 0, 0);
    imwrite ("SLIC.png", Out);

    M.show("Mosaico");

    waitKeyEx( 0 );
          
    frame.release(); fFrame.release(); gFrame.release(); labFrame.release(); qframe.release();
    slicFrame.release();
    Out.release();
    

    exit(EXIT_SUCCESS);
}
