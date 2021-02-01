#define _GNU_SOURCE   // To usine getline

// std libraries to include
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>

/* Macro to assert a condition, print error message if condition is true and exit */
#define fatalif(pred,s,...) \
 (void)((pred)? \
    fprintf(stderr, \
      (s)[sizeof s-2]==':'?s" %s\n":s"%.0s\n",##__VA_ARGS__,strerror(errno)), \
    exit(1), 0: 0)


/* Macros used to create vector data structures */
// Pointer to header of vector
#define v_p(v)  ((int*)(v)-2)

// Number of elements in vector if allocated
#define v_n(v)  (v_p(v)[0])

// Max number of allocated eleemnts in vector
#define v_m(v)  (v_p(v)[1])

// Number of elements in vectory
#define v_count(v) ((v)?v_n(v):0)

// Size of memory allocated for elements and header
#define v_sz(v,n) ((n)*sizeof*(v)+sizeof(int[2]))

// Increase the size of the vector
#define v_grow(v,k) \
  ((v) = !(v)?            (int *)calloc(1,v_sz(v,k))+2 : \
         v_m(v)<v_n(v)+k? v_m(v) = 2 * v_m(v) + k + 1,   \
                          (void *)((int *)realloc(v_p(v),v_sz(v,v_m(v)))+2) : \
                          (v))

// Append a new element to the vector
#define v_adjoin(v,e) (v_grow(v,1), (v)[v_n(v)++] = e)

// Pop the last element of the vector
#define v_pop(v) (--v_n(v))

// Get the last element of the vector
#define v_last(v) ((v)[v_n(v)-1])

// Reset vector to 0 elements
#define v_reset(v) (void)((v)? v_n(v)=0: 0)

// Free memory allocated by vector
#define v_free(v) ((v)?free(v_p(v)):0)
