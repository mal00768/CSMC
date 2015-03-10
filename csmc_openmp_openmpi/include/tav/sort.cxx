/** \file sort.cxx 
\brief Sorting, and sorted array operations
\ingroup gp_nr
 */
/*
August 08 2005 iswap via pointer in QuickSort
10/19/2005 Swap is from basic
11/03/2005 MinElement and MaxElement on array
06/20/2007 add find routine, which finds an exact match only
03/13/09 Modification to use reference and arbitrary template
06/25/09 Select, array cannot be constant. 
 */
/**
\defgroup gp_nr Numerical Recipes
*/
#ifndef __SORT_CXX__
#define __SORT_CXX__
#include <debug.h>
#include "basic.cxx" /*dependence: Swap function*/
namespace tav{
  /*!return location of smallest element*/
template <class clsT> 
unsigned long MinElement(unsigned long n, const clsT &arr)
{
  unsigned long ind=0;
  //typeof(arr[0]) x=arr[0];
    double x=arr[0];
  for (int i=0;i<n;i++) if (x>arr[i]) {x=arr[i];ind=i;}
  return ind;
}

  /* !return location of largest element*/
template <class clsT> 
unsigned long MaxElement(unsigned long n, const clsT &arr)
{
  unsigned long ind=0;
  //typeof(arr[0]) x=arr[0];
    double x=arr[0];
  for (int i=0;i<n;i++) if (x<arr[i]) {x=arr[i];ind=i;}
  return ind;
}



  //#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
/*! This function returns the k-th smallest element from an input arr of leng n. The input array is rearranged to have this value in 
arr[k], all smaller values are in arr[0] to arr[k-1]; all arger values are in arr[k+1..n-1]. Smaller and larger values are in arbitrary 
order.
Speed 2^sqrt(log_2(N))
*/
template <class clsT> 
clsT Select(unsigned long k, unsigned long n, clsT *arr)
{
	unsigned long i,ir,j,l,mid;
	//typeof(arr[0])  a,temp;
    double  a,temp;

	l=0;
	ir=n-1;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
			  Swap(arr[l],arr[ir]);
			}
			return arr[k-1];
		} else {
			mid=(l+ir-1) >> 1;
			Swap(arr[mid],arr[l+1]);
			if (arr[l] > arr[ir]) {
			  Swap(arr[l],arr[ir]);
			}
			if (arr[l+1] > arr[ir]) {
			  Swap(arr[l+1],arr[ir]);
			}
			if (arr[l] > arr[l+1]) {
			  Swap(arr[l],arr[l+1]);
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				Swap(arr[i],arr[j]);
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			if (j >= k-1) ir=j-1;
			if (j < k) l=i;
		}
	}
}
  //#undef SWAP

//#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 100
//#define ISWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
  //swap is ineffective but does not work in pointer, tested on integer
template <class clsT, class clsS> 
void QuickSort(unsigned long n, clsT &arr, clsS &brr)
{
	unsigned long i,ir=n,j,k,l=1;
	int *istack,jstack=0;
	//typeof(arr[0]) a,temp;
	//typeof(brr[0]) b,itemp;
    double a,temp;
	double b,itemp;
	//NSTACK>2*ln(N)/ln(2);
	istack=new int [NSTACK];
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j-1];
				b=brr[j-1];
				for (i=j-1;i>=l;i--) {
					if (arr[i-1] <= a) break;
					arr[i]=arr[i-1];
					brr[i]=brr[i-1];
				}
				arr[i]=a;
				brr[i]=b;
			}
			if (!jstack) {
			      
				delete []istack;
				return;
			}
			ir=istack[jstack-1];
			l=istack[jstack-2];
			jstack -= 2;
		} else {
			k=(l+ir) >> 1;
			Swap(arr[k-1],arr[l]);
			Swap(brr[k-1],brr[l]);
			if (arr[l-1] > arr[ir-1]) {
			  Swap(arr[l-1],arr[ir-1]);
			  Swap(brr[l-1],brr[ir-1]);
			}
			if (arr[l+1-1] > arr[ir-1]) {
			  Swap(arr[l],arr[ir-1]);
			  Swap(brr[l],brr[ir-1]);
			}
			if (arr[l-1] > arr[l]) {
			  Swap(arr[l-1],arr[l]);
			  Swap(brr[l-1],brr[l]);
			}
			i=l+1;
			j=ir;
			a=arr[l];
			b=brr[l];
			i--; //shift to speed
			j--; //shift to speed
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				Swap(arr[i],arr[j]);
				Swap(brr[i],brr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			brr[l]=brr[j];
			brr[j]=b;
			i++;
			j++;
			jstack += 2;
			if (jstack > NSTACK) FatalError("NSTACK too small in sort2.");
			if (ir-i+1 >= j-l) {
				istack[jstack-1]=ir;
				istack[jstack-2]=i;
				ir=j-1;
			} else {
				istack[jstack-1]=j-1;
				istack[jstack-2]=l;
				l=i;
			}
		}
	}
}
template <class clsT> 
void QuickSort(unsigned long n, clsT &arr)
{
	unsigned long i,ir=n,j,k,l=1;
	int *istack,jstack=0;
	//typeof(arr[0]) a,temp;
     double a,temp;

	istack=new int [NSTACK];
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j-1];
					for (i=j-1;i>=l;i--) {
					if (arr[i-1] <= a) break;
					arr[i]=arr[i-1];
					}
				arr[i]=a;
				}
			if (!jstack) {
			      
				delete []istack;
				return;
			}
			ir=istack[jstack-1];
			l=istack[jstack-2];
			jstack -= 2;
		} else {
			k=(l+ir) >> 1;
			Swap(arr[k-1],arr[l]);
				if (arr[l-1] > arr[ir-1]) {
				  Swap(arr[l-1],arr[ir-1]);
				}
			if (arr[l+1-1] > arr[ir-1]) {
			  Swap(arr[l],arr[ir-1]);
				}
			if (arr[l-1] > arr[l]) {
			  Swap(arr[l-1],arr[l]);
				}
			i=l+1;
			j=ir;
			a=arr[l];
				i--; //shift to speed
			j--; //shift to speed
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				Swap(arr[i],arr[j]);
				}
			arr[l]=arr[j];
			arr[j]=a;
				i++;
			j++;
			jstack += 2;
			if (jstack > NSTACK) FatalError("NSTACK too small in sort2.");
			if (ir-i+1 >= j-l) {
				istack[jstack-1]=ir;
				istack[jstack-2]=i;
				ir=j-1;
			} else {
				istack[jstack-1]=j-1;
				istack[jstack-2]=l;
				l=i;
			}
		}
	}
}
  //#undef ISWAP
#undef M
#undef NSTACK
//#undef SWAP


template <class clsT> 
void HeapSort(unsigned long n, clsT &ra)
{
	unsigned long i,ir,j,l;
	//typeof(ra[0])  rra;
    double rra;
    
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l-1];
		} else {
			rra=ra[ir-1];
			ra[ir-1]=ra[0];
			if (--ir == 1) {
				ra[0]=rra;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j-1] < ra[j]) j++;
			if (rra < ra[j-1]) {
				ra[i-1]=ra[j-1];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i-1]=rra;
	}
}

template <class clsT, class clsS> 
void HeapSort(unsigned long n, clsT & ra, clsS &rb)
{
	unsigned long i,ir,j,l;
	//typeof(ra[0])  rra;
   //typeof(rb[0]) rrb;
    double  rra;
    double rrb;
    
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l-1];
			rrb=rb[l-1];
		} else {
			rra=ra[ir-1];
			rrb=rb[ir-1];
			ra[ir-1]=ra[0];
			rb[ir-1]=rb[0];
			if (--ir == 1) {
				ra[0]=rra;
				rb[0]=rrb;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j-1] < ra[j]) j++;
			if (rra < ra[j-1]) {
				ra[i-1]=ra[j-1];
				rb[i-1]=rb[j-1];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i-1]=rra;
		rb[i-1]=rrb;
	}
}

/*! \brief Locate an element in a sorted list.

Given sorted array (either way ascending or descending) x[0]...x[n-1] 
and value x we search it
and return j so that x is between x[j] and x[j+1]; 
out of bound j=-1 or j=n-1
operator >= (cnumber, cnumber) is required
in cases with equal
\par Return values 
- ascending (123..) list return \f$j\,\, {\rm that}\,\, x\in [x[j],x[j+1])\quad   x[n]\equiv+\infty \f$
- descending (321..) list return \f$ x\in (x[j],x[j+1]]\quad  x[-1]\equiv-\infty \f$
 */
template <typename type_int, typename type_cnumber_array, typename cnumber> 
type_int Locate(const type_int n, const type_cnumber_array &xx, const cnumber x)

{
  //	if (x == xx[0]) return 0;
  //	if (x == xx[n-1]) return n-1; //this is still ok used to be n-2
  type_int ju,jm,jl;
	bool ascnd;
   
	//int n=xx.size();
	jl=-1;
	ju=n;
	ascnd=(xx[n-1] >= xx[0]);
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	return jl;
}
    /*! \brief Find an element in a sorted list.
     
     Given sorted array (either way ascending or descending) x[0]...x[n-1] 
     and value x we search it
    return only of exact match, else return n (out of bounds)
     */
    template <typename type_int, typename type_cnumber_array, typename cnumber> 
    type_int Find(const type_int n, const type_cnumber_array &xx, const cnumber x)
    
    {
 
        //	if (x == xx[0]) return 0;
        //	if (x == xx[n-1]) return n-1; //this is still ok used to be n-2
        type_int ju,jm,jl;
        bool ascnd;
        //int n=xx.size();

        ascnd=(xx[n-1] >= xx[0]); //true if ascending order
        jl=0;
        ju=n-1;
        if (ascnd) {if ((x< xx[jl])||(x>xx[ju]))  return n; } 
        else {if ((x> xx[jl])||(x<xx[ju])) return n; if (x==xx[jl]) return jl; }
        
        if (x==xx[ju]) return ju; 
        
        while (ju-jl > 1) {
            jm=(ju+jl) >> 1;
            if (x >= xx[jm] == ascnd)
                jl=jm;
            else
                ju=jm;
        }
        if (x==xx[jl]) return jl; else return n;
    }
/*
Given a monotonic array xx[0] to xx[n-1] and value x return a value jlo so that x is between x[jlo] x[jlo+1]
if jlo=-1 or jlo=n-1 the result  out of range, except x=xx[n-1], see my fix. Input jlo is used as initial guess. 
\todo there are possible note return of -1 may harm unsigned discuss later
*/
template <typename type_int, typename type_cnumber_array, typename cnumber> 
type_int Hunt(const type_int n, const type_cnumber_array &xx, const cnumber x, type_int jlo=0)
{
	type_int jm,jhi,inc;
	//note internal jlo is changed and returned. 
	bool ascnd;
	ascnd=(xx[n-1] >= xx[0]);
	if (jlo < 0 || jlo > n-1) {
		jlo=-1;
		jhi=n;
	} else {
		inc=1;
		if (x >= xx[jlo] == ascnd) {
			if (jlo == n-1) return jlo;
			jhi=jlo+1;
			while (x >= xx[jhi] == ascnd) {
				jlo=jhi;
				inc += inc;
				jhi=jlo+inc;
				if (jhi > n-1) {
					jhi=n;
					break;
				}
			}
		} else {
			if (jlo == 0) {
				jlo=-1;
				return jlo;
			}
			jhi=jlo--;
			while (x < xx[jlo] == ascnd) {
				jhi=jlo;
				inc <<= 1;
				if (inc >= jhi) {
					jlo=-1;
					break;
				}
				else jlo=jhi-inc;
			}
		}
	}
	while (jhi-jlo != 1) {
		jm=(jhi+jlo) >> 1;
		if (x >= xx[jm] == ascnd)
			jlo=jm;
		else
			jhi=jm;
	}
	if (x == xx[n-1]) jlo=n-1;//previously n-2
	if (x == xx[0]) jlo=0;
	return jlo;
}


}
#endif
