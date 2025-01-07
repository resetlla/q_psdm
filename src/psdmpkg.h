#ifndef PSDMPKG_H_
#define PSDMPKG_H_

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>

#include "fft.h"
#include "segy.h"

extern int xargc;
extern char **xargv;
extern int inited;

#define alloc1char(n1) malloc1c(n1)
#define alloc2char(n1, n2) malloc2c(n1, n2)
#define alloc1float(n1) malloc1f(n1)
#define alloc1int(n1) malloc1i(n1)
#define alloc2float(n1, n2) malloc2f(n1, n2)
#define alloc2int(n1, n2) malloc2i(n1, n2)
#define getparstring(key, val) getparstr(key, val)
#define free1char(p) free1c(p)
#define free2char(p) free2c(p)
#define free1int(p) free1i(p)
#define free2int(p) free2i(p)
#define free3int(p) free3i(p)
#define free1float(p) free1f(p)
#define free2float(p) free2f(p)
#define free3float(p) free3f(p)
#define free4float(p) free4f(p)
#define free1double(p) free1d(p)
#define free2double(p) free2d(p)
#define free1complex(p) free1x(p)
#define free2complex(p) free2x(p)

char *malloc1c(size_t n1);
char **malloc2c(size_t n1, size_t n2);

int *malloc1i(size_t n1);
int **malloc2i(size_t n1, size_t n2);
int ***malloc3i(size_t n1, size_t n2, size_t n3);

float *malloc1f(size_t n1);
float **malloc2f(size_t n1, size_t n2);
float ***malloc3f(size_t n1, size_t n2, size_t n3);
float ****malloc4f(size_t n1, size_t n2, size_t n3, size_t n4);

double *malloc1d(size_t n1);
double **malloc2d(size_t n1, size_t n2);
#define parseline(str, n) parseint(str, n)
int *parseint(const char *str, int *n);

void free1c(char *p);
void free2c(char **p);

void free1i(int *p);
void free2i(int **p);
void free3i(int ***p);

void free1f(float *p);
void free2f(float **p);
void free3f(float ***p);
void free4f(float ****p);

void free1d(double *p);
void free2d(double **p);

/* Parser the Command-Line Argument */
void initargs_inner(int argc, char *argv[], const char *sdoc[]);

/* Parser the Command-Line Argument */
#define initargs(argc, argv, sdoc)                                         \
  do {                                                                     \
    struct stat st;                                                        \
    int fd = fileno(stdin);                                                \
    fstat(fd, &st);                                                        \
    FILE *pp;                                                              \
    if (argc == 1 && S_ISCHR(st.st_mode)) {                                \
      pp = popen("more 1>&2", "w");                                        \
      for (int i = 0; sdoc[i] != NULL; ++i) {                              \
        fprintf(pp, "%s", sdoc[i]);                                        \
      }                                                                    \
      pclose(pp);                                                          \
      fprintf(stderr, "==================\n");                             \
      fprintf(stderr, "-built:      %s %s\n", __DATE__, __TIME__);         \
      fprintf(stderr, "-file:       %s\n", __FILE__);                      \
      fprintf(stderr, "-g++:        %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, \
              __GNUC_PATCHLEVEL__);                                        \
      fprintf(stderr, "-glibc:      %d.%d\n", __GLIBC__, __GLIBC_MINOR__); \
      fprintf(stderr, "-byte-order: %s\n",                                 \
              __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__ ? "little-endian"  \
                                                        : "big-endian");   \
      exit(0);                                                             \
    }                                                                      \
    initargs_inner(argc, argv, sdoc);                                      \
  } while (0)

#define egetparstring(name, value)                            \
  do {                                                        \
    if (!getparstring(name, value)) {                         \
      fprintf(stderr, "Error: %s must be specified\n", name); \
      exit(0);                                                \
    }                                                         \
  } while (0)

#define egetparint(name, value)                               \
  do {                                                        \
    if (!getparint(name, value)) {                            \
      fprintf(stderr, "Error: %s must be specified\n", name); \
      exit(0);                                                \
    }                                                         \
  } while (0)

#define egetparfloat(name, value)                             \
  do {                                                        \
    if (!getparfloat(name, value)) {                          \
      fprintf(stderr, "Error: %s must be specified\n", name); \
      exit(0);                                                \
    }                                                         \
  } while (0)

#define egetpardouble(name, value)                            \
  do {                                                        \
    if (!getpardouble(name, value)) {                         \
      fprintf(stderr, "Error: %s must be specified\n", name); \
      exit(0);                                                \
    }                                                         \
  } while (0)

size_t fsize(const char *path);

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64)
#include <fcntl.h>
#include <io.h>
#include <windows.h>

#define isatty(h) _isatty(h)

#ifndef access
#define access(f, m) _access((f), (m))
#endif

#undef popen
#define popen _popen
#undef pclose
#define pclose _pclose
#else
/* Make sure isatty() has a prototype.*/
#include <unistd.h>
/*extern int isatty(int);*/
/* popen and pclose are not C89 functions and so are sometimes omitted from
** the <stdio.h> header */
/*extern int fileno(FILE*);
extern FILE *popen(const char*,const char*);
extern int pclose(FILE*);*/
#endif

int getparint(const char *key, int *val);
int getparfloat(const char *key, float *val);
int getpardouble(const char *key, double *val);
int getparstr(const char *key, char **val); /* char *s */
char *strrev1(char *p);

#define exist(_p) fexist(_p)
int fexist(const char *path);

int fetch_and_taper(const float *data, int ns, int it1, int it2, float *data1);

int readparfile(const char *filename, const char *fmt, ...);

void time_offset_mute_santity(float *tomute, int ntomute);
char *strrtrim(char *src, int c);

void parsepath(char *path, char *projpath, char *projname, int *iproc,
               char *task);

float *parsetuple(const char *str, int *ntuple, int *nelem);

typedef enum {
  POLY6THCUT = 8,
  TANER4THCUT = 6,
  TANER2TH = 4,
  DSR = 2,
  OTHERS = 2
} TRAVELTIMEMETHOD;
typedef enum { AMPLITUDE = 1, AMPLITUDE_PHASE = 2 } CORRECTION;

typedef struct {
  int sx;
  int sy;
} coor_t;

#ifdef __cplusplus
extern "C" {
#endif

float getscale(short scalco);
void get_vel_coor(const char *velfilepath, int line, coor_t *vel_coor,
                  short *vel_scalel, int fci, int nxi, int fcv, int nxv,
                  int flv, int nyv);

/* Time-offset Mute Tuple */
float time_offset_mute(float *tomute, int ntomute, float offset);

#ifdef __cplusplus
}
#endif

#endif  // PSDMPKG_H_
