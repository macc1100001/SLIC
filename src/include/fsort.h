#ifndef __FSORT__
#define __FSORT__

/*
Sorting functions.

Algoritmos to sort small fixed-size arrays.

See:
http://en.wikipedia.org/wiki/Sorting_network#Constructing_sorting_networks
http://pages.ripco.net/~jgamble/nw.html
*/


#define Swap(x,y) {x^=y;y=x^y;x=y^x;}

#define iSwap(d, idx, a, b) if (d[idx[a]] > d[idx[b]]) Swap(idx[a], idx[b]);
#define sSwap(d, a, b) if (d[a] > d[b]) Swap(d[a], d[b]);


template <typename T>
void iSort2 (unsigned int *idx, T *d)
{
	 iSwap (d, idx, 0, 1);
}

template <typename T>
void iSort3 (unsigned int *idx, T *d)
{
	 iSwap (d, idx, 1, 2);
	 iSwap (d, idx, 0, 2);
	 iSwap (d, idx, 0, 1);
}

template <typename T>
void iSort4 (unsigned int *idx, T *d)
{
	 iSwap (d, idx, 0, 1);
	 iSwap (d, idx, 2, 3);
	 iSwap (d, idx, 0, 2);
	 iSwap (d, idx, 1, 3);
	 iSwap (d, idx, 1, 2);
}

template <typename T>
void iSort5 (unsigned int *idx, T *d)
{
	 iSwap (d, idx, 0, 1);
	 iSwap (d, idx, 3, 4);
	 iSwap (d, idx, 2, 4);
	 iSwap (d, idx, 2, 3);
	 iSwap (d, idx, 0, 3);
	 iSwap (d, idx, 0, 2);
	 iSwap (d, idx, 1, 4);
	 iSwap (d, idx, 1, 3);
	 iSwap (d, idx, 1, 2);
}

template <typename T>
void iSort6 (unsigned int *idx, T *d)
{
	 iSwap (d, idx, 1, 2);
	 iSwap (d, idx, 0, 2);
	 iSwap (d, idx, 0, 1);
	 iSwap (d, idx, 4, 5);
	 iSwap (d, idx, 3, 5);
	 iSwap (d, idx, 3, 4);
	 iSwap (d, idx, 0, 3);
	 iSwap (d, idx, 1, 4);
	 iSwap (d, idx, 2, 5);
	 iSwap (d, idx, 2, 4);
	 iSwap (d, idx, 1, 3);
	 iSwap (d, idx, 2, 3);
}

template <typename T>
void iSort7 (unsigned int *idx, T *d)
{
	 iSwap (d, idx, 1, 2);
	 iSwap (d, idx, 0, 2);
	 iSwap (d, idx, 0, 1);
	 iSwap (d, idx, 3, 4);
	 iSwap (d, idx, 5, 6);
	 iSwap (d, idx, 3, 5);
	 iSwap (d, idx, 4, 6);
	 iSwap (d, idx, 4, 5);
	 iSwap (d, idx, 0, 4);
	 iSwap (d, idx, 0, 3);
	 iSwap (d, idx, 1, 5);
	 iSwap (d, idx, 2, 6);
	 iSwap (d, idx, 2, 5);
	 iSwap (d, idx, 1, 3);
	 iSwap (d, idx, 2, 4);
	 iSwap (d, idx, 2, 3);
}

template <typename T>
void Sort2 (T *d)
{
	 sSwap (d, 0, 1);
}

template <typename T>
void Sort3 (T *d)
{
	 sSwap (d, 1, 2);
	 sSwap (d, 0, 2);
	 sSwap (d, 0, 1);
}

template <typename T>
void Sort4 (T *d)
{
	 sSwap (d, 0, 1);
	 sSwap (d, 2, 3);
	 sSwap (d, 0, 2);
	 sSwap (d, 1, 3);
	 sSwap (d, 1, 2);
}

template <typename T>
void Sort5 (T *d)
{
	 sSwap (d, 0, 1);
	 sSwap (d, 3, 4);
	 sSwap (d, 2, 4);
	 sSwap (d, 2, 3);
	 sSwap (d, 0, 3);
	 sSwap (d, 0, 2);
	 sSwap (d, 1, 4);
	 sSwap (d, 1, 3);
	 sSwap (d, 1, 2);
}

template <typename T>
void Sort6 (T *d)
{
	 sSwap (d, 1, 2);
	 sSwap (d, 0, 2);
	 sSwap (d, 0, 1);
	 sSwap (d, 4, 5);
	 sSwap (d, 3, 5);
	 sSwap (d, 3, 4);
	 sSwap (d, 0, 3);
	 sSwap (d, 1, 4);
	 sSwap (d, 2, 5);
	 sSwap (d, 2, 4);
	 sSwap (d, 1, 3);
	 sSwap (d, 2, 3);
}

template <typename T>
void Sort7 (T *d)
{
	 sSwap (d, 1, 2);
	 sSwap (d, 0, 2);
	 sSwap (d, 0, 1);
	 sSwap (d, 3, 4);
	 sSwap (d, 5, 6);
	 sSwap (d, 3, 5);
	 sSwap (d, 4, 6);
	 sSwap (d, 4, 5);
	 sSwap (d, 0, 4);
	 sSwap (d, 0, 3);
	 sSwap (d, 1, 5);
	 sSwap (d, 2, 6);
	 sSwap (d, 2, 5);
	 sSwap (d, 1, 3);
	 sSwap (d, 2, 4);
	 sSwap (d, 2, 3);
}


#endif
