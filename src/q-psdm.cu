#include <cuda.h>
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include "fft.h"
#include "psdmpkg.h"
#include "queue.h"
#include "segy.h"

#include <sstream>
#include <string>
#include <vector>

#define QMIG3D
#define FULL_DERIVATIVE

// #define  WARN_ON

#ifdef WARN_ON
#define WARN(fmt, ...)                   \
  do {                                   \
    fprintf(stderr, fmt, ##__VA_ARGS__); \
  } while (0)
#else
#define WARN(fmt, ...)
#endif

const char *sdoc[] = {
#ifdef DMIG3D_WITHOUT_GATHER
    "Prestack Depth Migration: Dip Migration With Multithreading             "
    "\n",
#else
    "Prestack Depth Migration: Generate Dip-Angle Gather With Multithreading "
    "\n",
#endif
    "                          Using grid travel time table                  "
    "\n",
    " dgthr3dz\033[1;31m path=/..  parfilepath=.. nodename=**_[0|1..]         "
    "\n",
    "         ttpath=/..                                      \033[0m        "
    "\n",
    "        *[offspath=/..] *[noffpath=/..] *[velfilepath=/..]              "
    "\n",
    "        *[infopath=/..] *[dipdirpath=/..] *[outputpath=/..]             "
    "\n",
    "        *[logfilepath=/..] *[traveltime=2|4|6] *[ninv=20]               "
    "\n",
    "        *[f1=1] *[f2=3] *[f3=130.0] *[f4=140](Hz)                       "
    "\n",
    "        *fxs= dxs= nxs= fys= dys= nys=                                  "
    "\n",
    "        *dxt= nxt= dyt= nyt= dzt= nzt=                                  "
    "\n",
#ifdef DMIG3D_WITHOUT_GATHER
    "        *[opcrp=1] *[cdps=#cdp1-#cdp2] *[aspect=N]                      "
    "\n",
#else
    "        *[opcrp=0] *[cdps=#cdp1-#cdp2] *[aspect=A|X|Y|N]                "
    "\n",
    "        *[dipx=50] *[dipy=50] *[ddipx=1] *[ddipy=1](degree)             "
    "\n",
#endif
    "        *dyv=#header                                                  "
    "\n\n",
    "        curl=NULL                                                       "
    "\n",
    "    traveltime: 2-DSR 4-TANER2TH 6-TANER4THCUT                          "
    "\n",
    "          ninv: Time Ratio Used to Inverse the Interval Velocity        "
    "\n",
    "    f1,2,3,4  : Filter Parameter Applied to the Input Seismic Data      "
    "\n",
#ifdef DMIG3D_WITHOUT_GATHER
    "                                                                      "
    "\n\n",
#else
    "    opcrp     : 0, not output CRP                                       "
    "\n",
    "                1, output CRP                                           "
    "\n",
    "    cdps      : specify the CDPs output in CRP                          "
    "\n",
    "    aspect    : All,X+Y;  X; Y; N  for Dip-Gather                       "
    "\n",
    "    dipx,dipy : Maximum dip for X,Y Direction                           "
    "\n",
    "    ddipx,y   : Interval                                              "
    "\n\n",
#endif
    "\033[1;34m  SU Format     \033[0m                                     "
    "\n\n",
    " * indicate the path is derivate from the {path}                      "
    "\n\n",
    NULL};

segy tr, tr1;
FILE *fplog;
int is3D = 1;

typedef struct {
  int64_t idx;
  uint16_t nys;
  uint16_t nxs;
  uint16_t dys;
  uint16_t dxs;
} ttt_desc_t;

struct CPUParamsNoLine {
  // 从main中读取的参数
  float dyi;
  float dxi;
  int nxi;
  int nzi;
  int ntl;
  double scaled;
  int fci;
  float dzi;
  int ndipx;
  int ndipy;
  float ddipx;
  float ddipy;
  float dt;
  float fc;
  // 内部参数
  float eps;
  int *iaztab;
  float fxv1;
  float fyv1;
  float exv1;
  float eyv1;
  int nfft;
  int nffti;
  float *ww;
  float dt1;
  int ndipx2;
  int ndipy2;
  float lnG;
  float coef1;
  float coef2;
  float coef3;
  float coef4;
  float taperzone;
  float dzii;
};

struct CPUParamsWithLine {
  double fxi;
  double cyi;
  float dw;
  int ****aper;
  int **aznx1;
  int **aznx2;
  int **azny1;
  int **azny2;
  int ioffs;
  int itb1;
  int itb2;
  int nsb;
  float tdstart1;
  float tdstart2;
  int nf1;
  float tdmid;
  int nww;
  int nxg;
  int nzg;
  int nxi1;
  int nzi1;
  int ttt_cdp1;
  int ttt_cdp2;
  int qqq_cdp1;
  int qqq_cdp2;
  int ncoef1;
  int ixstart;
  int ixend;
  int nxb;
};

struct GPUParamsNoLine {
  int igpu;
  float *re;
  float *im;
  float *red;
  float *imd;
  int *itibeg;
  int *itibegd;
};

struct GPUParamsWithLine {
  int *nf3d;
  float *Qd;
  float *datav;
  fftwf_complex *wdatav;
  fftwf_plan planv;
  ttt_desc_t *qqq_descd;
  float *qqqd;
  fftwf_plan plan;
  fftwf_plan plani;
  fftwf_complex *wdata;
  float *data;
  float *datad;
  float *dipx1d;
  float *dipx2d;
  float *dipy1d;
  float *dipy2d;
  ttt_desc_t *ttt_descd;
  float *tttd;
  float *imgd;
  float *gthrxd;
  float *gthryd;
};

struct PROD_PARAMS {
  safe_queue<segy *> *datapool;
  int ioffl1;
  int ioffl2;
  char *offspath;
  int ioffs;
  int ntl;
  int ncons;
  char *nodename;
  uint64_t *t_read;
  off_t *s_read;
};

struct CONS_PARAMS {
  safe_queue<segy *> *datapool;
  CPUParamsNoLine *cpuParamsNoLine;
  CPUParamsWithLine *cpuParamsWithLine;
  GPUParamsNoLine *gpuParamsNoLine;
  GPUParamsWithLine *gpuParamsWithLine;
};

void *producer(void *arg);
void *consumer(void *arg);

void qpsdm(int ngpu, int ndata, int ncons, int poolsize, char *projpath,
           char *projtask, float dxi, int nxi, float dyi, int nyi, int nzi,
           float dzi, float dzii, float zmin, int fcv, double fxv, int flv,
           double fyv, int nxv, int nyv, int ntv, float dtv, int naxi, int nayi,
           float daxi, float dayi, int line1, int line2, int nyd, int fld,
           int ntl, float dt, float f1, float fc, float f3, float f4,
           int *noffs, int extxline, int *iline, int fci, int iyi1,
           float **offs, int *cdps, int ncdps, int nsmooth, double scaled,
           float **tmute, char *infopath, int aspect, int opcrp,
           char *outputpath, char *offspath, char *path, char *velfilepath,
           char *dipdirpath, char *f3dirpath, char *ttpath, char *nodename,
           float threshold, float contract, float taperzone, int nblock);

void psdm_kernel(segy tr,
                 // CPUParamsNoLine
                 float dyi, float dxi, int nxi, int nzi, int ntl, double scaled,
                 int fci, float dzi, int ndipx, int ndipy, float ddipx,
                 float ddipy, float dt, float fc, float eps, int *iaztab,
                 float fxv1, float fyv1, float exv1, float eyv1, int nfft,
                 int nffti, float *ww, float dt1, int ndipx2, int ndipy2,
                 float taperzone, float lnG, float coef1, float coef2,
                 float coef3, float coef4, float dzii,
                 // CPUParamsWithLine
                 int ****aper, int **aznx1, int **aznx2, int **azny1,
                 int **azny2, double fxi, double cyi, float dw, int ioffs,
#ifdef QMIG3D
                 int itb1, int itb2, int nsb, float tdstart1, float tdstart2,
                 int nf1, float tdmid,
#endif
                 int nww, int nxg, int nzg, int nxi1, int nzi1, int ttt_cdp1,
                 int ttt_cdp2, int qqq_cdp1, int qqq_cdp2, int ncoef1,
                 int ixstart, int ixend, int nxb,
// GPUParamsNoLine
#ifdef QMIG3D
                 float *re, float *im, float *red, float *imd,
#endif
                 int igpu, int *itibeg, int *itibegd,
// GPUParamsWithLine
#ifdef QMIG3D
                 int *nf3d, float *Qd, float *datav, fftwf_complex *wdatav,
                 fftwf_plan planv, ttt_desc_t *qqq_descd, float *qqqd,
#endif
                 float *data, float *datad, float *dipx1d, float *dipx2d,
                 float *dipy1d, float *dipy2d, ttt_desc_t *ttt_descd,
                 float *tttd, float *imgd);

__global__ void image_depth_gpu(
    float *img, int izimin, int nzi, float *data, int nt, float dt, int ixstart,
    int ixend, int nxb, float dzii,
#ifdef QMIG3D
    const float *__restrict__ red, const float *__restrict__ imd, int nf1,
    float fc, float dw, int nww, int *nf3, ttt_desc_t *qqq_desc, float *qqq,
    float tstart1, float tstart2, float tmid,
#endif
    float taperzone, float lnG, float coef1, float coef2, float coef3,
    float coef4, int ncoef1, int nxg, int nzg, int nxi1, int nzi1, int ttt_cdp1,
    int ttt_cdp2, float sx, float sy, float gx, float gy, float dxi, float dzi,
    int *itibegd, int ixi1, int ixi2, int ndipx, int ndipx2, float ddipx,
    int ndipy, int ndipy2, float ddipy, ttt_desc_t *ttt_desc, float *ttt,
    float *dipx1, float *dipx2, float *dipy1, float *dipy2);

int get_ttt_all(const char *ttpath, const char *torq, int line, int line1,
                int line2, int *nxg, int *nzg, int *nxi1, int *nzi1,
                int *ttt_cdp1, int *ttt_cdp2,
                std::vector<GPUParamsWithLine *> gpuParamsWithLineVec,
                int ncons);
void reset_ttt_all(ttt_desc_t *desc_out, float *data_out);

int main(int argc, char *argv[]) {
  /* 命令行参数 */
  char *path, *ttpath, *parfilepath, *nodename, *offspath, *noffpath,
      *velfilepath, *dipdirpath, *f3dirpath, *outputpath, *logfilepath,
      *infopath, *zcdps, *zaspect;
  int ninv;
  float f1, f2, f3, f4, fc;
  float threshold, contract, taperzone;

  /* 参数文件参数:数据空间 */
  int ntl, line1, line2; /* Time Sample Number,First Data Line,Last Data Line */
  float dt;              /* Time Sample Interval                              */

  /* 参数文件参数:成像空间 */
  int cdp1, cdp2; /* Image Range:cdp,lines,time */
  float dzi, zmin, zmax, dzii;
  int *iline, nyi;
  char *zline;

  /* 参数文件参数:切除和孔径 */
  char *zmute;
  int xline;

  /* 速度空间变量 */
  double fxv, fyv;
  float dxv, dyv, dtv, dzv;
  int flv, fcv, nxv, nyv, ntv;

  /* 成像空间变量 */
  int fci, nxi, nzi;
  /* AUX */
  float dyi, dxi;

  /* 数据空间变量 */
  int fld, nyd;

  /* 偏移距 */
  int *noffs;
  float offmin, offmax;
  float **offs = NULL;

  /* 切除信息 */
  float **tmute = NULL;

  /* 检查点信息变量 */
  int iyi1, ilb;

  /* 常规索引 */
  int idata, ioffs;

  /* MISC */
  char str[512];
  double scalev, scaled;
  int ret, nsmooth;
  int igpu, ngpu;
  FILE *fp;
  int i, line;
  time_t tic, toc;

  /* 命令行拼写参数 */
  int iproc;
  char projpath[1024];
  char projname[512];
  char projtask[512];

  /* 切除曲线 */
  float *tomute;
  int ntomute, npair;

  cudaError_t err;

  /* 倾角道集相关 */
  int opcrp, aspect;
  float dipx, dipy, ddipx, ddipy;
  int *cdps, ncdps;
  int ndipx, ndipy;

  /*	分块 */
  int nblock;

  /* 获取命令行参数:强制 */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  initargs(argc, argv, sdoc);

  if (!getparstring("path", &path)) {
    fprintf(stderr, "Error: path must be specified\n");
    exit(1);
  }
  /* 去掉path最后的一个或多个/(如果存在的话) */
  strrtrim(path, '/');

  if (!getparstring("parfilepath", &parfilepath)) {
    fprintf(stderr, "Error: parfilepath must be specified\n");
    exit(1);
  }

  if (!getparstring("nodename", &nodename)) {
    fprintf(stderr, "Error: nodename must be specified\n");
    exit(1);
  }

  if (!getparstring("ttpath", &ttpath)) {
    fprintf(stderr, "Error: ttpath must be specified\n");
    exit(1);
  }

  /*	多线程参数 */
  int ndata, ncons, poolsize;

  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  /* 解析path到各个子项 */
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  parsepath(path, projpath, projname, &iproc, projtask);
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  /* 再次获取命令行参数:弱获取 */
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (!getparstring("logfilepath", &logfilepath)) {
    logfilepath = alloc1char(1024);
    if (1 != readparfile(parfilepath, "logfilepath=%s", logfilepath)) {
      sprintf(logfilepath, "/tmp/%s_%s.log", nodename, projname);
    }
  }

  /* 日志系统开启 */
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  fplog = stderr;  // fopen(logfilepath,"w");
  if (fplog == NULL) {
    fprintf(stderr, "Warn: log file [%s] open error\n", logfilepath);
    fprintf(stderr, "      redirect it to stderr\n");
    fplog = stderr;
  }
  setbuf(fplog, NULL);

  time(&tic);
  fprintf(fplog,
          "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          "++\n");
  fprintf(fplog,
          "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          "++\n");
  fprintf(fplog,
          "+ Log File For Prestack Time Migration, Generate Dip Gather With "
          "Multithreading @ %s \n",
          nodename);
  fprintf(fplog,
          "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          "++\n");
  fprintf(fplog,
          "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
          "++\n");
  fprintf(fplog, "Time: %s\n", ctime(&tic));
  fflush(fplog);
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  if (!getparstring("offspath", &offspath)) {
    offspath = alloc1char(1024);
    if (1 != readparfile(parfilepath, "offspath=%s", offspath)) {
      sprintf(offspath, "%s/data%d", projpath, iproc);
    }
  }

  if (!getparstring("velfilepath", &velfilepath)) {
    velfilepath = alloc1char(1024);
    if (1 != readparfile(parfilepath, "velfilepath=%s", velfilepath)) {
      sprintf(velfilepath, "%s/vel/vel_loadin", projpath);
    }
  }

  if (!getparstring("infopath", &infopath)) {
    infopath = alloc1char(1024);
    if (1 != readparfile(parfilepath, "infopath=%s", infopath)) {
      strcpy(infopath, path);
    }
  }

  if (!getparstring("outputpath", &outputpath)) {
    outputpath = alloc1char(1024);
    if (1 != readparfile(parfilepath, "outputpath=%s", outputpath)) {
      sprintf(outputpath, "%s/result", path);
    }
  }

  if (!getparint("ninv", &ninv)) {
    if (1 != readparfile(parfilepath, "ninv=%d", &ninv)) {
      ninv = 20;
    }
  }

  if (!getparfloat("f1", &f1)) {
    if (1 != readparfile(parfilepath, "f1=%f", &f1)) {
      f1 = 1;
    }
  }

  if (!getparfloat("fc", &fc)) {
    if (1 != readparfile(parfilepath, "fc=%f", &fc)) {
      fc = 40;
    }
  }

  if (!getparfloat("contract", &contract)) {
    if (1 != readparfile(parfilepath, "contract=%f", &contract)) {
      contract = 3.0;
    }
  }

  if (!getparfloat("taperzone", &taperzone)) {
    if (1 != readparfile(parfilepath, "taperzone=%f", &taperzone)) {
      taperzone = 8.0;
    }
  }

  if (!getparfloat("threshold", &threshold)) {
    if (1 != readparfile(parfilepath, "threshold=%f", &threshold)) {
      threshold = 600.0;
    }
  }

  if (!getparfloat("f2", &f2)) {
    if (1 != readparfile(parfilepath, "f2=%f", &f2)) {
      f2 = 4;
    }
  }
  if (!getparfloat("f3", &f3)) {
    if (1 != readparfile(parfilepath, "f3=%f", &f3)) {
      f3 = 130;
    }
  }

  if (!getparfloat("f4", &f4)) {
    if (1 != readparfile(parfilepath, "f4=%f", &f4)) {
      f4 = 140;
    }
  }

  if (!getparint("ntl", &ntl)) {
    if (1 != readparfile(parfilepath, "ntl=%d", &ntl)) {
      fprintf(stderr, "Error: ntl must be specified\n");
      fprintf(fplog, "Error: ntl must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  if (!getparfloat("dt", &dt)) {
    if (1 != readparfile(parfilepath, "dt=%f", &dt)) {
      fprintf(stderr, "Error: dt must be specified\n");
      fprintf(fplog, "Error: dt must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  if (!getparint("line1", &line1)) {
    if (1 != readparfile(parfilepath, "line1=%d", &line1)) {
      fprintf(stderr, "Error: line1 must be specified\n");
      fprintf(fplog, "Error: line1 must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  if (!getparint("line2", &line2)) {
    if (1 != readparfile(parfilepath, "line2=%d", &line2)) {
      fprintf(stderr, "Error: line2 must be specified\n");
      fprintf(fplog, "Error: line2 must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  if (!getparfloat("zmin", &zmin)) {
    if (1 != readparfile(parfilepath, "zmin=%f", &zmin)) {
      fprintf(stderr, "Error: zmin must be specified\n");
      fprintf(fplog, "Error: zmin must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  if (!getparfloat("zmax", &zmax)) {
    if (1 != readparfile(parfilepath, "zmax=%f", &zmax)) {
      fprintf(stderr, "Error: zmax must be specified\n");
      fprintf(fplog, "Error: zmax must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  if (!getparfloat("dzi", &dzi)) {
    if (1 != readparfile(parfilepath, "dzi=%f", &dzi)) {
      fprintf(stderr, "Error: dzi must be specified\n");
      fprintf(fplog, "Error: dzi must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  if (!getparint("cdp1", &cdp1)) {
    if (1 != readparfile(parfilepath, "cdp1=%d", &cdp1)) {
      fprintf(stderr, "Error: cdp1 must be specified\n");
      fprintf(fplog, "Error: cdp1 must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  if (!getparint("cdp2", &cdp2)) {
    if (1 != readparfile(parfilepath, "cdp2=%d", &cdp2)) {
      fprintf(stderr, "Error: cdp2 must be specified\n");
      fprintf(fplog, "Error: cdp2 must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  if (!getparstring("line", &zline)) {
    zline = alloc1char(2048);
    if (1 != readparfile(parfilepath, "line=%s", zline)) {
      fprintf(stderr, "Error: line must be specified\n");
      fprintf(fplog, "Error: line must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  if (!getparstring("mute", &zmute)) {
    zmute = alloc1char(2048);
    if (1 != readparfile(parfilepath, "mute=%s", zmute)) {
      strcpy(zmute, "(1,9999),(100,9999)");
    }
  }

  if (!getparint("extxline", &xline)) {
    if (1 != readparfile(parfilepath, "extxline=%d", &xline)) {
      fprintf(stderr, "Error: extxline must be specified\n");
      fprintf(fplog, "Error: extxline must be specified\n");
      fflush(fplog);
      exit(1);
    }
  }

  /* Dip Gather Related */
  if (!getparint("opcrp", &opcrp)) {
    if (1 != readparfile(parfilepath, "opcrp=%d", &opcrp)) {
#ifdef DMIG3D_WITHOUT_GATHER
      opcrp = 1;
#else
      opcrp = 0;
#endif
    }
  }

  if (!getparstring("cdps", &zcdps)) {
    zcdps = alloc1char(1024);
    if (1 != readparfile(parfilepath, "cdps=%s", zcdps)) {
      sprintf(zcdps, "%d-%d", cdp1, cdp2);
    }
  }

  if (!getparstring("aspect", &zaspect)) {
    zaspect = alloc1char(16);
    if (1 != readparfile(parfilepath, "aspect=%s", zaspect)) {
#ifdef DMIG3D_WITHOUT_GATHER
      strcpy(zaspect, "None");
#else
      strcpy(zaspect, "All");
#endif
    }
  }

  if (!getparfloat("dipx", &dipx)) {
    if (1 != readparfile(parfilepath, "dipx=%f", &dipx)) {
      dipx = 50.0;
    }
  }

  if (!getparfloat("dipy", &dipy)) {
    if (1 != readparfile(parfilepath, "dipy=%f", &dipy)) {
      dipy = 50.0;
    }
  }

  if (!getparfloat("ddipx", &ddipx)) {
    if (1 != readparfile(parfilepath, "ddipx=%f", &ddipx)) {
      ddipx = 1.0;
    }
  }

  if (!getparfloat("ddipy", &ddipy)) {
    if (1 != readparfile(parfilepath, "ddipy=%f", &ddipy)) {
      ddipy = 1.0;
    }
  }

  if (!getparint("nblock", &nblock)) {
    if (1 != readparfile(parfilepath, "nblock=%d", &nblock)) {
      nblock = 1;
    }
  }

  /* 部分合理性检测 */
  if (f3 > 0.5 / dt - 10.0) /* freq(nyq)-10 hz */
  {
    f3 = 0.5 / dt - 10.0;
  }

  if (f4 > 0.5 / dt - 5.0) /* freq(nyq)-10 hz */
  {
    f4 = 0.5 / dt - 5.0;
  }

  if (f4 - f3 <= 5.0) {
    f3 = f4 - 5.0;
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /* 设置运行GPU参数定义 */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  ret = sscanf(nodename, "%*[^_]_%d", &igpu);
  if (ret == 1) {
    err = cudaGetDeviceCount(&ngpu);
    if (err != cudaSuccess) {
      fprintf(fplog, "===========================================\n");
      sprintf(str, "lsb_release -a     >>%s", logfilepath);
      if (system(str) != 0) {
        fprintf(stderr, "Command Executed Fail:%s\n", str);
        fprintf(fplog, "Command Executed Fail:%s\n", str);
        fflush(fplog);
      }

      sprintf(str, "lspci | grep 'VGA' >>%s", logfilepath);
      if (system(str) != 0) {
        fprintf(stderr, "Command Executed Fail:%s\n", str);
        fprintf(fplog, "Command Executed Fail:%s\n", str);
        fflush(fplog);
      }

      sprintf(str, "nvidia-smi         >>%s", logfilepath);
      if (system(str) != 0) {
        fprintf(stderr, "Command Executed Fail:%s\n", str);
        fprintf(fplog, "Command Executed Fail:%s\n", str);
        fflush(fplog);
      }

      sprintf(str, "free -g            >>%s", logfilepath);
      if (system(str) != 0) {
        fprintf(stderr, "Command Executed Fail:%s\n", str);
        fprintf(fplog, "Command Executed Fail:%s\n", str);
        fflush(fplog);
      }

      fprintf(fplog, "===========================================\n");

      fprintf(fplog, "Error: cudaGetDevice Error %d\n", ngpu);
      fprintf(fplog, "       %s\n", cudaGetErrorString(err));

      fflush(fplog);
    }
    if (igpu >= 0 && igpu < ngpu) {
      WARN("\n\nThis Machine has %d-GPU Card\n", ngpu);
    } else {
      igpu = igpu % ngpu;
      fprintf(stderr,
              "Warn: This Machine [%s] has no such GPU using id=%d instead\n",
              nodename, igpu);
      fprintf(fplog,
              "Warn: This Machine [%s] has no such GPU using id=%d instead\n",
              nodename, igpu);
    }
  } else {
    fprintf(stderr,
            "Error: Nodename [%s] Format Error or This Machine has no GPUs\n",
            nodename);
    fprintf(fplog,
            "Error: Nodename [%s] Format Error or This Machine has no GPUs\n",
            nodename);
    fflush(fplog);
    exit(1);
  }
  fflush(fplog);

  ngpu = 1;

  if (!getparint("ndata", &ndata)) {
    if (1 != readparfile(parfilepath, "ndata=%d", &ndata)) {
      ndata = ngpu;
    }
  }

  if (!getparint("ncons", &ncons)) {
    if (1 != readparfile(parfilepath, "ncons=%d", &ncons)) {
      ncons = ngpu;
    }
  }

  if (!getparint("poolsize", &poolsize)) {
    if (1 != readparfile(parfilepath, "poolsize=%d", &poolsize)) {
      poolsize = 10000;
    }
  }

  if (ndata > ngpu) ndata = ngpu;

  fprintf(stderr, "ndata=%d ncons=%d poolsize=%d\n", ndata, ncons, poolsize);
  fflush(fplog);

  /* 合理性检测以及变量内部化 */
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  cdps = parseline(zcdps, &ncdps);

  if (zaspect[0] == 'A' || zaspect[0] == 'a') {
    aspect = 1;
  } else if (zaspect[0] == 'X' || zaspect[0] == 'x') {
    aspect = 2;
  } else if (zaspect[0] == 'Y' || zaspect[0] == 'y') {
    aspect = 3;
  } else if (zaspect[0] == 'N' || zaspect[0] == 'n') {
    aspect = 4;
  } else {
    aspect = 1;
  }

  ndipx = ((int)(dipx / ddipx + 0.5));
  ndipy = ((int)(dipy / ddipy + 0.5));

  nsmooth = 2;

  /* 获取成像空间参数 */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  fci = cdp1;
  nxi = cdp2 - cdp1 + 1;
  iline = parseline(zline, &nyi);

  nzi = (int)(zmax / dzi + 0.5) + 1;

  /* 获取速度空间参数 */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  if (!exist(velfilepath)) {
    fprintf(stderr, "Error: On Node [%s], Velocity File [%s] Does Not Exist\n",
            nodename, velfilepath);
    fprintf(fplog, "Error: On Node [%s], Velocity File [%s] Does Not Exist\n",
            nodename, velfilepath);
    fflush(fplog);
    exit(1);
  }
  fp = fopen(velfilepath, "rb");

  fgettr(fp, &tr);
  fseek(fp, -(240L + tr.ns * sizeof(float)), SEEK_END);
  fgettr(fp, &tr1);
  fclose(fp);

  scalev = getscale(tr.scalco);

  fxv = tr.gx * scalev;
  fyv = tr.gy * scalev;
  dxv = (tr1.gx - tr.gx) * scalev / (tr1.cdp - tr.cdp);
  if (tr.fldr == tr1.fldr) {
    if (!getparfloat("dyv", &dyv)) {
      if (1 != readparfile(parfilepath, "dyv=%f", &dyv)) {
        fprintf(stderr,
                "since first fldr and last fldr in the velocity file is the "
                "same:%d please specify dyv!\n",
                tr.fldr);
        fprintf(fplog,
                "since first fldr and last fldr in the velocity file is the "
                "same:%d please specify dyv!\n",
                tr.fldr);
        fflush(fplog);
        exit(1);
      }
    }
  } else {
    dyv = (tr1.gy - tr.gy) * scalev / (tr1.fldr - tr.fldr);
    /* Make Sure the Commandline is always Vaild */
    if (!getparfloat("dyv", &dyv)) {
      readparfile(parfilepath, "dyv=%f", &dyv);
    }
  }
  nxv = tr1.cdp - tr.cdp + 1;
  nyv = tr1.fldr - tr.fldr + 1;
  fcv = tr.cdp;
  flv = tr.fldr;
  ntv = tr.ns;
  dtv = tr.dt * 5.0E-7;
  dzv = 8;
  dzii = dzv * round(50 / dzv);
  fprintf(fplog, "dzv=%f dzii=%f\n", dzv, dzii);

  /* 获取数据空间参数 */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  fld = line1;
  nyd = line2 - line1 + 1;

  /* 获取数据空间scalco->scaled */
  /* 合理性确认                 */
  line = iline[0];
  scaled = 1.0;
  for (i = line - 1000; i < line + 1000; ++i) {
    sprintf(str, "%s/off%d/line%d", offspath, 0, i);
    if (exist(str)) {
      fp = fopen(str, "rb");
      fgettr(fp, &tr);
      scaled = getscale(tr.scalco);

      if (tr.ns != ntl) {
        fprintf(stderr, "Error: In Data Space (tr.ns=%d)!=(ntl=%d)\n", tr.ns,
                ntl);
        fprintf(fplog, "Error: In Data Space (tr.ns=%d)!=(ntl=%d)\n", tr.ns,
                ntl);
        fflush(fplog);
        exit(1);
      }
      if (tr.dt != (int)(dt * 1.0E6)) {
        fprintf(stderr, "Error: In Data Space (tr.dt=%dms)!=(dt=%g*1.0E6)\n",
                tr.dt, dt);
        fprintf(fplog, "Error: In Data Space (tr.dt=%dms)!=(dt=%g*1.0E6)\n",
                tr.dt, dt);
        fflush(fplog);
        exit(1);
      }
      fclose(fp);
      break;
    }
  }

  /* 内部参数 */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  dxi = dxv;
  dyi = dyv;
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  char nodename0[128];
  noffpath = alloc1char(1024);
  sscanf(nodename, "%[^_]s", nodename0);
  noffs = alloc1int(ndata);
  memset(noffs, 0, ndata * sizeof(int));

  for (idata = 0; idata < ndata; ++idata) {
    sprintf(offspath, "%s/data%d", projpath, idata);
    fprintf(stderr, "offspath=%s\n", offspath);
    sprintf(nodename, "%s_%d", nodename0, idata);
    fprintf(stderr, "nodename=%s\n", nodename);
    sprintf(noffpath, "%s/nofffile.rcd", offspath);

    if (!exist(noffpath)) {
      fprintf(stderr, "Error: On %s,[%s] Does Not Exist\n", nodename, noffpath);
      fprintf(fplog, "Error: On %s,[%s] Does Not Exist\n", nodename, noffpath);
      fflush(fplog);
      exit(1);
    }

    fp = fopen(noffpath, "r");
    if (fscanf(fp, "nofffile=%d\n", &noffs[idata]) != 1) {
      fprintf(stderr, "Read noffs from %s fail! in LINE:%d\n", noffpath,
              __LINE__);
      fprintf(fplog, "Read noffs from %s fail! in LINE:%d\n", noffpath,
              __LINE__);
      fflush(fplog);
      exit(1);
    }
    fclose(fp);
  }

  offs = alloc2float(noffs[0] * 2, ndata);
  memset(offs[0], 0, noffs[0] * 2 * ndata * sizeof(float));
  for (idata = 0; idata < ndata; ++idata) {
    sprintf(offspath, "%s/data%d", projpath, idata);
    sprintf(nodename, "%s_%d", nodename0, idata);
    for (ioffs = 0; ioffs < noffs[idata]; ++ioffs) {
      sprintf(str, "%s/off%d", offspath, ioffs);
      if (!exist(str)) {
        fprintf(stderr, "Error: On %s,[%s] Does Not Exist\n", nodename, str);
        fprintf(fplog, "Error: On %s,[%s] Does Not Exist\n", nodename, str);
        fflush(fplog);
        exit(1);
      }
      sprintf(str, "%s/off%d/ell_off%d", offspath, ioffs, ioffs);
      if (!exist(str)) {
        fprintf(stderr, "Error: On %s,[%s] Does Not Exist\n", nodename, str);
        fprintf(fplog, "Error: On %s,[%s] Does Not Exist\n", nodename, str);
        fflush(fplog);
        exit(1);
      }

      fp = fopen(str, "r");
      if (fscanf(fp, "%f,%f\n", &offmin, &offmax) != 2) {
        fprintf(stderr, "Read offmin and offmax from %s fail! in LINE:%d\n",
                str, __LINE__);
        fprintf(fplog, "Read offmin and offmax from %s fail! in LINE:%d\n", str,
                __LINE__);
        fflush(fplog);
        exit(1);
      }
      fclose(fp);

      offs[idata][ioffs] = offmin;
      offs[idata][ioffs + noffs[idata]] = offmax;
    }
  }

  /* 读取CIG(CRP)切除参数,并插值出偏移距组所对应的切除时间点 */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  tomute = parsetuple(zmute, &ntomute, &npair);

  fprintf(fplog, "Mute Tuple:(DEPTH:m,OFFSET:m)\n");
  for (i = 0; i < ntomute; ++i) {
    fprintf(fplog, "%-5.2f %-5.0f\n", tomute[i * 2], tomute[i * 2 + 1]);
  }
  fprintf(fplog, "\n");

  tmute = alloc2float(noffs[0], ndata);
  memset(tmute[0], 0, noffs[0] * ndata * sizeof(float));

  time_offset_mute_santity(tomute, ntomute);

  for (idata = 0; idata < ndata; ++idata) {
    for (ioffs = 0; ioffs < noffs[idata]; ++ioffs) {
      tmute[idata][ioffs] = time_offset_mute(
          tomute, ntomute, offs[idata][ioffs]); /* offset-time pair */
      if (tmute[idata][ioffs] < 0) tmute[idata][ioffs] = 0;
      if (tmute[idata][ioffs] < zmin) tmute[idata][ioffs] = zmin;
    }
  }
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  sprintf(nodename, "%s_0", nodename0);
  fprintf(stderr, "%s - Start\n", nodename);
  fprintf(fplog, "%s - Start\n", nodename);
  fflush(fplog);

  time(&tic);

  qpsdm(ngpu, ndata, ncons, poolsize, projpath, projtask, dxi, nxi, dyi, nyi,
        nzi, dzi, dzii, zmin, fcv, fxv, flv, fyv, nxv, nyv, ntv, dtv, ndipx,
        ndipy, ddipx, ddipy, line1, line2, nyd, fld, ntl, dt, f1, fc, f3, f4,
        noffs, xline, iline, fci, iyi1, offs, cdps, ncdps, nsmooth, scaled,
        tmute, infopath, aspect, opcrp, outputpath, offspath, path, velfilepath,
        dipdirpath, f3dirpath, ttpath, nodename, threshold, contract, taperzone,
        nblock);

  time(&toc);

  fprintf(stderr, "%s,Time Elapse:%.0fs\n", nodename, difftime(toc, tic));
  fprintf(stderr, "%s:All Migration Work Done\n", nodename);

  fprintf(fplog, "%s,Time Elapse:%.0fs\n", nodename, difftime(toc, tic));
  fprintf(fplog, "%s:All Migration Work Done\n", nodename);
  fprintf(fplog, "\nTime: %s", ctime(&toc));
  fflush(fplog);

  // sprintf(str, "%s/break_point.txt", offspath);
  // unlink(str);

  // fprintf(fplog, "unlink BreakPoint:%s\n", str);
  // fflush(fplog);

  return 0;
}

void qpsdm(int ngpu, int ndata, int ncons, int poolsize, char *projpath,
           char *projtask, float dxi, int nxi, float dyi, int nyi, int nzi,
           float dzi, float dzii, float zmin, int fcv, double fxv, int flv,
           double fyv, int nxv, int nyv, int ntv, float dtv, int ndipx,
           int ndipy, float ddipx, float ddipy, int line1, int line2, int nyd,
           int fld, int ntl, float dt, float f1, float fc, float ff3, float f4,
           int *noffs, int extxline, int *iline, int fci, int iyi1,
           float **offs, int *cdps, int ncdps, int nsmooth, double scaled,
           float **tmute, char *infopath, int aspect, int opcrp,
           char *outputpath, char *offspath, char *path, char *velfilepath,
           char *dipdirpath, char *f3dirpath, char *ttpath, char *nodename,
           float threshold, float contract, float taperzone, int nblock) {
  (void)ngpu;
  (void)(nyd);
  (void)(fld);
  (void)(aspect);
  (void)(path);
  (void)(velfilepath);
  (void)(f3dirpath);
  (void)zmin;
  (void)ntv;
  (void)dtv;
  (void)f1;
  (void)nsmooth;
  (void)dipdirpath;

  char str[1024];

  float c1, c2;
  int i, i1;

  int ixi, iyi, iw;
  int ioffs;
  float odt1, df;
  double fxi, cyi;
  int n3, n4;

  int nfft, nffti;
  float dw, w;

  float *ww;
  int ioffl1, ioffl2;
  int line;
  FILE *fp;

  int cdp1;
  int izi;

  float **dipx1, **dipx2, **dipy1, **dipy2;
#ifdef QMIG3D
  float **f3;
  int **nf3;
#endif

  cudaError_t err;

  char nodename0[128];

  time_t toc, tic;

  float **img = NULL, **tmpimg = NULL;

#ifdef QMIG3D
  float **Q;
#endif

  float dt1;

  int izmin;

  float tanx1 = tanf(60.0 / 180.0 * M_PI);
  float tanx2 = tanf(75.0 / 180.0 * M_PI);
  float tany1 = tanf(10.0 / 180.0 * M_PI);
  float tany2 = tanf(15.0 / 180.0 * M_PI);

  int ndipx2, ndipy2;

  int nxb1, ixstart, ixend, nxb;
  nxb1 = int(nxi / nblock + 0.5);

  /*aper*/
  int ****aper;
  int **aznx1, **azny1, **aznx2, **azny2;

  float eps;
  int naztab, *iaztab, i2, iaz;
  float c3;

  /* travle time talbe store with grid */
  /* in each grid the z is contigous   */
  int nxg, nzg, nxi1, nzi1;

  // int izibeg;
  // izibeg = zmin / dzi;
  fprintf(stderr, "nxi=%d, nzi=%d\n", nxi, nzi);

  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  eps = 0.001f;
  naztab = (int)(2.0 / eps + 1.5);
  iaztab = alloc1int(naztab);

  c2 = (1.0 - cosf(10.0 * 0.0174532925199)) * M_PI * 0.5; /* 方位角 */
  c2 = 1.0 - cosf(c2 * 0.5);
  i2 = (int)(c2 / eps);
  for (i = 0; i < i2; ++i) {
    iaztab[i] = 0;
  }
  for (iaz = 1; iaz < 18; ++iaz) {
    c2 = iaz * 10.0;
    c1 = c2 - 10.0;
    c3 = c2 + 10.0;

    c1 = (1.0 - cosf(c1 * 0.0174532925199)) * M_PI * 0.5; /* 方位角 */
    c2 = (1.0 - cosf(c2 * 0.0174532925199)) * M_PI * 0.5;
    c3 = (1.0 - cosf(c3 * 0.0174532925199)) * M_PI * 0.5;
    c1 = 1.0 - cosf(c2 - (c2 - c1) * 0.5);
    c2 = 1.0 - cosf(c2 + (c3 - c2) * 0.5);

    i1 = (int)(c1 / eps);
    i2 = (int)(c2 / eps);

    for (i = i1; i < i2; ++i) iaztab[i] = iaz;
  }
  for (i = i2; i < naztab; ++i) iaztab[i] = 0;

  /* 有效的孔径起算时间 */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  ndipx2 = ndipx * 2 + 1;
  ndipy2 = ndipy * 2 + 1;

  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  sscanf(nodename, "%[^_]s", nodename0);
  cdp1 = fci;

  /* 申请数据空间以及重采样空间 */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  nfft = get_opt_n_fftw(ntl);
  nffti = get_opt_n_fftw(ntl * 2); /* 一般超过2倍的重采样率 */

  df = 1.0 / (nfft * dt);
  dw = 2.0 * M_PI * df;
  odt1 = nffti * df;
  dt1 = 1.0 / odt1;

  n3 = (int)(ff3 / df + 0.5);
  n4 = (int)(f4 / df + 0.5);

  ww = alloc1float(nffti / 2 + 1);

#if defined(HALF_DERIVATIVE) && defined(FULL_DERIVATIVE)
#error "Macro: HALF_DERIVATIVE & FULL_DERIVATIVE is both defined"
#elif !defined(HALF_DERIVATIVE) && !defined(FULL_DERIVATIVE)
#error "Macro: HALF_DERIVATIVE & FULL_DERIVATIVE is neither defined"
#endif

#ifdef HALF_DERIVATIVE
  for (iw = 0, w = 0.0; iw < n3; ++iw, w += dw) {
    ww[iw] = sqrtf(w);
  }
  c1 = n4 - n3;
  w = sqrtf(w);
  for (iw = n3; iw <= n4; ++iw) {
    c2 = cosf((iw - n3) / c1 * M_PI * 0.5);
    c2 *= c2;
    ww[iw] = w * c2;
  }
#elif defined(FULL_DERIVATIVE)
  for (iw = 0, w = 0.0; iw < n3; ++iw, w += dw) {
    ww[iw] = w;
  }
  c1 = n4 - n3;
  for (iw = n3; iw <= n4; ++iw) {
    c2 = cosf((iw - n3) / c1 * M_PI * 0.5);
    c2 *= c2;
    ww[iw] = w * c2;
  }
#endif
  for (iw = n4; iw < nffti / 2 + 1; ++iw) {
    ww[iw] = 0.0;
  }

  int nww = 110.0 / df;

  fprintf(stderr, "nww=%d\n", nww);
  fflush(fplog);

  /* 申请倾角场空间 */
  dipx1 = alloc2float(nzi, nxi);
  dipx2 = alloc2float(nzi, nxi);
  dipy1 = alloc2float(nzi, nxi);
  dipy2 = alloc2float(nzi, nxi);

#ifdef QMIG3D
  /* 申请nf3高通频空间 */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  f3 = alloc2float(nzi, nxi);
  nf3 = alloc2int(nzi, nxi);

  fprintf(fplog, "Create Space For F3\n");
  Q = alloc2float(nzi, nxi);
#endif

  /* 常量 */
  float lnGc, coef1c, coef2c, coef3c, coef4c;
  int ncoef1c;

  lnGc = logf(threshold);
  coef1c = threshold * (1.0 - lnGc - 2.5 * lnGc * lnGc);
  coef2c = threshold * (1.0 + 5.0 * lnGc);
  coef3c = -2.5 * threshold;
  coef4c = 1.1 * threshold;

  ncoef1c = (int)(10.0 / df);
  /* 倾角衰减模式系数拷贝到GPU */
  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  contract *= 0.01745329251994;
  taperzone *= 0.01745329251994;

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*++++++++++++++++++++++++    初始化GPU参数     ++++++++++++++++++++++*/
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

  // 存储参数
  std::vector<GPUParamsNoLine *> gpuParamsNoLineVec;
  std::vector<GPUParamsWithLine *> gpuParamsWithLineVec;
  fprintf(fplog, "Initialize Line-Independent CPU/GPU Parameters Start ...\n");
  /*	初始化每个线程的GPU参数 */
  for (int igpu = 0; igpu < ncons; ++igpu) {
    cudaSetDevice(igpu);

    GPUParamsNoLine *gpuParamsNoLine = new GPUParamsNoLine;
    GPUParamsWithLine *gpuParamsWithLine = new GPUParamsWithLine;

#ifdef QMIG3D
    float *re, *im;
    float *red, *imd;
    float *Qd;
    float *datav;
    fftwf_complex *wdatav;
    int *nf3d;
#endif

    int *itibeg;
    int *itibegd;

    float *dipx1d, *dipx2d, *dipy1d, *dipy2d;

    float *data = NULL;
    float *datad;

#ifdef QMIG3D
    re = alloc1float(ntl * 4);
    im = alloc1float(ntl * 4);
    cudaMalloc(&red, sizeof(float) * ntl * 4);
    cudaMalloc(&imd, sizeof(float) * ntl * 4);
    gpuParamsNoLine->re = re;
    gpuParamsNoLine->im = im;
    gpuParamsNoLine->red = red;
    gpuParamsNoLine->imd = imd;

    err = cudaMalloc(&nf3d, sizeof(int) * nzi * nxi);
    if (err != cudaSuccess) {
      fprintf(stderr, "Error: CUDA  Error in FILE:%s FUNCTION:%s LINE:%d\n",
              __FILE__, __FUNCTION__, __LINE__);
      fprintf(stderr, "       %s\n", cudaGetErrorString(err));
      exit(1);
    }
    cudaMalloc(&Qd, sizeof(float) * nxi * nzi);
    gpuParamsWithLine->nf3d = nf3d;
    gpuParamsWithLine->Qd = Qd;

    datav = fftwf_alloc_real(ntl * 4);
    wdatav = fftwf_alloc_complex(ntl * 4);
    gpuParamsWithLine->datav = datav;
    gpuParamsWithLine->wdatav = wdatav;
#endif

    itibeg = alloc1int(nxi);
    err = cudaMalloc(&itibegd, sizeof(int) * nxi);
    if (err != cudaSuccess) {
      fprintf(fplog, "Error: CUDA  Error in FILE:%s FUNCTION:%s LINE:%d\n",
              __FILE__, __FUNCTION__, __LINE__);
      fprintf(fplog, "       %s\n", cudaGetErrorString(err));
      exit(1);
    }
    gpuParamsNoLine->itibeg = itibeg;
    gpuParamsNoLine->itibegd = itibegd;

    cudaMalloc(&dipx1d, sizeof(float) * nxi * nzi);
    cudaMalloc(&dipx2d, sizeof(float) * nxi * nzi);
    cudaMalloc(&dipy1d, sizeof(float) * nxi * nzi);
    cudaMalloc(&dipy2d, sizeof(float) * nxi * nzi);
    gpuParamsWithLine->dipx1d = dipx1d;
    gpuParamsWithLine->dipx2d = dipx2d;
    gpuParamsWithLine->dipy1d = dipy1d;
    gpuParamsWithLine->dipy2d = dipy2d;

    data = alloc1float(nffti);
    cudaMalloc(&datad, sizeof(float) * nffti);
    gpuParamsWithLine->data = data;
    gpuParamsWithLine->datad = datad;

    gpuParamsNoLine->igpu = igpu;
    gpuParamsNoLineVec.push_back(gpuParamsNoLine);
    gpuParamsWithLineVec.push_back(gpuParamsWithLine);
  }

  /* 层速度模型空间 X方向两边各扩展1000m Y方向两边各扩展500m*/
  float fxv1, fyv1, exv1, eyv1;
  fxv1 = fxv - 2000;
  exv1 = fxv + nxv * dxi + 2000;
  fyv1 = fyv - 1000;
  eyv1 = fyv + nyv * dyi + 1000;

  int ttt_cdp1, ttt_cdp2;
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*+++++++++++++++++++++++++++    主循环     ++++++++++++++++++++++++*/
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  fprintf(fplog,
          "Initialize Line-dependent CPU/GPU Parameters And Start Migration "
          "With Multithreading ...\n");
  CPUParamsNoLine *cpuParamsNoLine = new CPUParamsNoLine;
  cpuParamsNoLine->dyi = dyi;
  cpuParamsNoLine->dxi = dxi;
  cpuParamsNoLine->nxi = nxi;
  cpuParamsNoLine->nzi = nzi;
  cpuParamsNoLine->ntl = ntl;
  cpuParamsNoLine->scaled = scaled;
  cpuParamsNoLine->fci = fci;
  cpuParamsNoLine->dzi = dzi;
  cpuParamsNoLine->ndipx = ndipx;
  cpuParamsNoLine->ndipy = ndipy;
  cpuParamsNoLine->ddipx = ddipx;
  cpuParamsNoLine->ddipy = ddipy;
  cpuParamsNoLine->dt = dt;
  cpuParamsNoLine->fc = fc;
  cpuParamsNoLine->eps = eps;
  cpuParamsNoLine->iaztab = iaztab;
  cpuParamsNoLine->fxv1 = fxv1;
  cpuParamsNoLine->fyv1 = fyv1;
  cpuParamsNoLine->exv1 = exv1;
  cpuParamsNoLine->eyv1 = eyv1;
  cpuParamsNoLine->nfft = nfft;
  cpuParamsNoLine->nffti = nffti;
  cpuParamsNoLine->ww = ww;
  cpuParamsNoLine->dt1 = dt1;
  cpuParamsNoLine->ndipx2 = ndipx2;
  cpuParamsNoLine->ndipy2 = ndipy2;
  cpuParamsNoLine->lnG = lnGc;
  cpuParamsNoLine->coef1 = coef1c;
  cpuParamsNoLine->coef2 = coef2c;
  cpuParamsNoLine->coef3 = coef3c;
  cpuParamsNoLine->coef4 = coef4c;
  cpuParamsNoLine->taperzone = taperzone;
  cpuParamsNoLine->dzii = dzii;
  /* 循环成像线 */
  for (iyi = iyi1; iyi < nyi; ++iyi) {
    time(&tic);

    line = iline[iyi];
    fxi = fxv + (fci - fcv) * dxi;
    cyi = fyv + (line - flv) * dyi;

    fprintf(fplog, "Start Image Line:%d\n", line);

    /* 确定一个偏移距文件夹下,对该成像线有贡献的数据域线范围 */
    ioffl1 = line - extxline;
    ioffl2 = line + extxline;
    if (ioffl1 < line1) ioffl1 = line1;
    if (ioffl2 > line2) ioffl2 = line2;

    /* store the original coor */
    coor_t *vel_coor = (coor_t *)malloc(sizeof(coor_t) * nxi);
    short vel_scalel;
    get_vel_coor(velfilepath, line, &vel_coor[0], &vel_scalel, fci, nxi, fcv,
                 nxv, flv, nyv);

    /* get travel time table */
    get_ttt_all(ttpath, "t", line, line1, line2, &nxg, &nzg, &nxi1, &nzi1,
                &ttt_cdp1, &ttt_cdp2, gpuParamsWithLineVec, ncons);
#ifdef QMIG3D
    int qqq_cdp1, qqq_cdp2;
    get_ttt_all(ttpath, "q", line, line1, line2, &nxg, &nzg, &nxi1, &nzi1,
                &qqq_cdp1, &qqq_cdp2, gpuParamsWithLineVec, ncons);
    if (ttt_cdp1 != qqq_cdp1 || ttt_cdp2 != qqq_cdp2 || ttt_cdp1 >= ttt_cdp2) {
      fprintf(
          fplog,
          "travel time table cdp range error, or not confirm with qfield\n");
      exit(1);
    }
#endif

    /* 将倾角转存为其tan值 */
    for (ixi = 0; ixi < nxi; ++ixi) {
      for (izi = 0; izi < nzi; ++izi) {
        dipx1[ixi][izi] =
            tanf(-10 * 0.0174532925199 + contract); /* *PI/180.0 */
        dipx2[ixi][izi] = tanf(10 * 0.0174532925199 - contract);
        dipy1[ixi][izi] = tanf(-10 * 0.0174532925199 + contract);
        dipy2[ixi][izi] = tanf(10 * 0.0174532925199 - contract);
      }
    }
#ifdef QMIG3D
    /* 读取高通频f3,插值 */
    for (ixi = 0; ixi < nxi; ++ixi) {
      for (izi = 0; izi < nzi; ++izi) {
        if (f3[ixi][izi] > ff3) {
          f3[ixi][izi] = 90;
        }

        nf3[ixi][izi] = (int)(f3[ixi][izi] / df + 0.5);
      }
    }
    fprintf(fplog, "Interpolation of F3 for line:%d\n", line);
#endif

#ifdef QMIG3D
    for (ixi = 0; ixi < nxi; ++ixi) {
      for (izi = 0; izi < nzi; ++izi) {
        Q[ixi][izi] = dw / (2.0 * (100.0 + izi * 100.0 / nzi) /*tr.data[izi]*/);
      }
    }
#endif

    for (int igpu = 0; igpu < ncons; ++igpu) {
      cudaSetDevice(igpu);

      GPUParamsWithLine *gpuParamsWithLine = gpuParamsWithLineVec[igpu];

      cudaMemcpy(gpuParamsWithLine->dipx1d, dipx1[0], sizeof(float) * nzi * nxi,
                 cudaMemcpyHostToDevice);
      cudaMemcpy(gpuParamsWithLine->dipx2d, dipx2[0], sizeof(float) * nzi * nxi,
                 cudaMemcpyHostToDevice);
      cudaMemcpy(gpuParamsWithLine->dipy1d, dipy1[0], sizeof(float) * nzi * nxi,
                 cudaMemcpyHostToDevice);
      cudaMemcpy(gpuParamsWithLine->dipy2d, dipy2[0], sizeof(float) * nzi * nxi,
                 cudaMemcpyHostToDevice);

#ifdef QMIG3D
      cudaMemcpy(gpuParamsWithLine->nf3d, nf3[0], sizeof(int) * nxi * nzi,
                 cudaMemcpyHostToDevice);
      cudaMemcpy(gpuParamsWithLine->Qd, Q[0], sizeof(float) * nxi * nzi,
                 cudaMemcpyHostToDevice);
#endif
    }

    CPUParamsWithLine *cpuParamsWithLine = new CPUParamsWithLine;
    for (int iblock = 0; iblock < nblock; ++iblock) {
      ixstart = iblock * nxb1;
      ixend = (iblock + 1) * nxb1;
      if (iblock == nblock - 1) ixend = nxi;
      nxb = ixend - ixstart;
      cpuParamsWithLine->ixstart = ixstart;
      cpuParamsWithLine->ixend = ixend;
      cpuParamsWithLine->nxb = nxb;

      img = alloc2float(nzi, nxb);
      tmpimg = alloc2float(nzi, nxb);

      for (int igpu = 0; igpu < ncons; ++igpu) {
        cudaSetDevice(igpu);

        float *imgd = NULL;

        cudaMalloc(&imgd, sizeof(float) * nxb * nzi);
        gpuParamsWithLineVec[igpu]->imgd = imgd;
      }

      for (int idata = 0; idata < ndata; ++idata) {
        fprintf(fplog, "Create Apeture Azimuth Index Table\n");

        aznx1 = alloc2int(18, noffs[idata]);
        aznx2 = alloc2int(18, noffs[idata]);
        azny1 = alloc2int(18, noffs[idata]);
        azny2 = alloc2int(18, noffs[idata]);

        aper = NULL;
        fprintf(fplog, "Aperture For Line %d\n", line);
        fprintf(fplog, "Start Offset Group Loop\n");
        fflush(fplog);

        cpuParamsWithLine->aper = aper;
        cpuParamsWithLine->aznx1 = aznx1;
        cpuParamsWithLine->aznx2 = aznx2;
        cpuParamsWithLine->azny1 = azny1;
        cpuParamsWithLine->azny2 = azny2;

        sprintf(outputpath, "%s/%d/%s/result", projpath, idata, projtask);
        sprintf(nodename, "%s_%d", nodename0, idata);
        sprintf(offspath, "%s/data%d", projpath, idata);
        sprintf(infopath, "%s/%d/%s", projpath, idata, projtask);
        /* 循环偏移距组文件夹 */
        for (ioffs = 0; ioffs < noffs[idata]; ++ioffs) {
          fprintf(fplog, "Input Offset Group: %d of %d in Block%d \n",
                  ioffs + 1, noffs[idata], iblock);
          fflush(fplog);

          sprintf(str, "%s/%s.off%d.L%d.su", outputpath, nodename, ioffs, line);

          if (iblock == 0) {
            if (exist(str)) {
              fprintf(fplog, "File:%s Exist\n", str);
              if (0 == unlink(str))
                fprintf(fplog, "File:%s been unlinked\n", str);
              else
                fprintf(fplog, "File:%s unlink failure\n", str);
            }
          }

          /* 设置切除参数 */
          izmin = (int)(tmute[idata][ioffs] / dzi);

          if (izmin > nzi) {
            fprintf(
                fplog,
                "Warn: Offset=%f is greater than the offset mute parameter\n",
                tmute[idata][ioffs]);
            fprintf(fplog, "      No image is generate\n");
            fflush(fplog);
            continue;
          }

#ifdef QMIG3D
          /* according to offset mute -> split to two segment */
          float t00, t11;
          int it00pre = 0.48 / dt;
          int it00suf = 0.36 / dt;
          float tdstart1, tdstart2, tdmid;
          int itb1, itb2, nsb, nf1, nf4;
          int iti;

#if 0
          t00 = ts01[ioffs][0] - it00pre * dt;
          t00 = t00 > 0.0 ? t00 : 0.0;
          t11 = ts01[ioffs][1] + it00suf * dt;
          t11 = t11 < (ntl - 1) * dt ? t11 : (ntl - 1) * dt;
#endif
          t00 = 0;
          t11 = (ntl - 1) * dt;

          nsb = (t11 - t00 - it00pre * dt - it00suf * dt) / 2.0 / dt;

          itb1 = t00 / dt;
          tdstart1 = itb1 * dt;

          tdmid = t00 + it00pre * dt + nsb * dt;

          nsb = nsb + it00pre + it00suf;
          nsb = get_opt_n_fftw(nsb);

          itb2 = t11 / dt - nsb + 1;
          tdstart2 = itb2 * dt;

          fprintf(fplog, "offset =%f\n", offs[idata][ioffs]);
          fprintf(fplog, "t00 =%f\n", t00);
          fprintf(fplog, "t11 =%f\n", t11);
          fprintf(fplog, "tdstart1 =%f\n", tdstart1);
          fprintf(fplog, "tdmid    =%f\n", tdmid);
          fprintf(fplog, "tdstart2 =%f\n", tdstart2);
          fprintf(fplog, "nsb =%d\n", nsb);

          df = 1.0 / nsb / dt;
          dw = 2.0 * M_PI * df;
          nf1 = (int)(f1 / df);
          nf4 = (int)(f4 / df);
          nww = nf4 - nf1 + 1;

          ncoef1c = (int)(10.0 / df);

          for (ixi = 0; ixi < nxi; ++ixi) {
            for (iti = 0; iti < nzi; ++iti) {
              nf3[ixi][iti] = (int)(f3[ixi][iti] / df + 0.5);
            }
          }
#endif

          /* 初始化offset剖面 */
          if (opcrp) {
            memset(img[0], 0, sizeof(float) * nxb * nzi);
            memset(tmpimg[0], 0, sizeof(float) * nxb * nzi);
          }

          for (int igpu = 0; igpu < ncons; ++igpu) {
            cudaSetDevice(igpu);

            GPUParamsWithLine *gpuParamsWithLine = gpuParamsWithLineVec[igpu];
            if (opcrp) {
              cudaMemset(gpuParamsWithLine->imgd, 0, sizeof(float) * nxb * nzi);
            }

#ifdef QMIG3D
            gpuParamsWithLine->planv =
                fftwf_plan_dft_r2c_1d(nsb, gpuParamsWithLine->datav,
                                      gpuParamsWithLine->wdatav, FFTW_MEASURE);
            cudaMemcpy(gpuParamsWithLine->nf3d, nf3[0], sizeof(int) * nxi * nzi,
                       cudaMemcpyHostToDevice);
#endif
          }
          cpuParamsWithLine->fxi = fxi;
          cpuParamsWithLine->cyi = cyi;
          cpuParamsWithLine->dw = dw;
          cpuParamsWithLine->ioffs = ioffs;
#ifdef QMIG3D
          cpuParamsWithLine->itb1 = itb1;
          cpuParamsWithLine->itb2 = itb2;
          cpuParamsWithLine->nsb = nsb;
          cpuParamsWithLine->tdstart1 = tdstart1;
          cpuParamsWithLine->tdstart2 = tdstart2;
          cpuParamsWithLine->nf1 = nf1;
          cpuParamsWithLine->tdmid = tdmid;
          cpuParamsWithLine->qqq_cdp1 = qqq_cdp1;
          cpuParamsWithLine->qqq_cdp2 = qqq_cdp2;
#endif
          cpuParamsWithLine->nww = nww;
          cpuParamsWithLine->nxg = nxg;
          cpuParamsWithLine->nzg = nzg;
          cpuParamsWithLine->nxi1 = nxi1;
          cpuParamsWithLine->nzi1 = nzi1;
          cpuParamsWithLine->ttt_cdp1 = ttt_cdp1;
          cpuParamsWithLine->ttt_cdp2 = ttt_cdp2;
          cpuParamsWithLine->ncoef1 = ncoef1c;

          safe_queue<segy *> datapool(poolsize);
          pthread_t prod_thread;

          PROD_PARAMS *prod_params = new PROD_PARAMS;
          prod_params->datapool = &datapool;
          prod_params->ioffl1 = ioffl1;
          prod_params->ioffl2 = ioffl2;
          prod_params->offspath = offspath;
          prod_params->ioffs = ioffs;
          prod_params->ntl = ntl;
          prod_params->ncons = ncons;
          prod_params->nodename = nodename;

          if (pthread_create(&prod_thread, NULL, producer, prod_params) != 0) {
            fprintf(stderr, "Error: Producer Thread Create failed\n");
            exit(1);
          }

          pthread_t cons_threads[ncons];
          std::vector<CONS_PARAMS *> consParamVec;
          for (int igpu = 0; igpu < ncons; ++igpu) {
            GPUParamsNoLine *gpuParamsNoLine = gpuParamsNoLineVec[igpu];
            GPUParamsWithLine *gpuParamsWithLine = gpuParamsWithLineVec[igpu];

            /* consumer	*/
            CONS_PARAMS *cons_params = new CONS_PARAMS;
            cons_params->datapool = &datapool;
            cons_params->cpuParamsNoLine = cpuParamsNoLine;
            cons_params->cpuParamsWithLine = cpuParamsWithLine;
            cons_params->gpuParamsNoLine = gpuParamsNoLine;
            cons_params->gpuParamsWithLine = gpuParamsWithLine;

            if (pthread_create(&cons_threads[igpu], NULL, consumer,
                               cons_params) != 0) {
              fprintf(stderr, "Error: Consumer-%d Thread Create failed\n",
                      igpu);
              exit(1);
            }
            consParamVec.push_back(cons_params);
          }

          pthread_join(prod_thread, NULL);

          for (i = 0; i < ncons; ++i) {
            pthread_join(cons_threads[i], NULL);
          }

          delete prod_params;

          if (opcrp) {
            fprintf(fplog, "Get Image(crp) from GPU For Line %d\n", line);

            for (int igpu = 0; igpu < ncons; ++igpu) {
              cudaSetDevice(igpu);

              GPUParamsWithLine *gpuParamsWithLine = gpuParamsWithLineVec[igpu];
              cudaMemcpy(tmpimg[0], gpuParamsWithLine->imgd,
                         sizeof(float) * nxb * nzi, cudaMemcpyDeviceToHost);

              for (i = 0; i < nxb * nzi; ++i) {
                img[0][i] += tmpimg[0][i];
              }
              memset(tmpimg[0], 0, sizeof(float) * nxb * nzi);
              delete consParamVec[igpu];
            }
            consParamVec.clear();

#ifdef QMIG3D
            for (ixi = 0; ixi < nxb; ++ixi) {
              for (iti = 0; iti < nzi; ++iti) {
                img[ixi][iti] /= nsb;
              }
            }
#endif
            fprintf(fplog, "Output:%s in Block%d \n", str, iblock);

            /* append binary */
            fp = fopen(str, "ab");

            if (fp != NULL) {
              memset(&tr, 0, sizeof(tr));
              tr.ns = nzi;
              tr.dt = (int)(dzi * 1.0E3);
              tr.fldr = line;
              tr.offset = (int)(offs[idata][ioffs]);
              tr.trid = 1;
              tr.counit = 1;
              /*  tr.delrt=(int)(itmin*dti*2.0E3); */
              tr.nvs = 1;
              i1 = 0;
              for (i = 0; i < ncdps; ++i) {
                ixi = cdps[i] - cdp1;
                if (ixi >= ixstart && ixi < ixend) {
                  ++i1;
                  tr.tracl = i1;
                  tr.cdp = fci + ixi;
                  /* output original coor in vel */
                  tr.sx = vel_coor[ixi].sx;
                  tr.sy = vel_coor[ixi].sy;
                  tr.gx = vel_coor[ixi].sx;
                  tr.gy = vel_coor[ixi].sy;
                  tr.scalco = vel_scalel;

                  memcpy((tr.data), img[ixi - ixstart], sizeof(float) * nzi);
                  fputtr(fp, &tr);
                }
              }
              fclose(fp);

              fprintf(fplog, "Write Image(crp) To Disk For Line %d:%s\n", line,
                      str);
              fflush(fplog);
            } else {
              fprintf(fplog, "Error: On [%s], File [%s] Open Error\n", nodename,
                      str);
              fprintf(fplog, "       %s\n", strerror(errno));

              exit(1);
            }
          }
        } /* 循环offs 结束 */
      }
    }

    delete cpuParamsWithLine;

    /* 输出局部检查点信息 */
    sprintf(str, "%s/break_point.txt", offspath);
    fp = fopen(str, "w");
    fprintf(fp, "%d\n", line);
    fclose(fp);

    sprintf(str, "%s/%s_%d", infopath, nodename, line);
    fp = fopen(str, "w");
    if (fp != NULL) {
      time(&toc);
      fprintf(fplog, "%s\n", ctime(&toc));
      fclose(fp);
    }

    for (int igpu = 0; igpu < ncons; ++igpu) {
      cudaSetDevice(igpu);

      GPUParamsWithLine *gpuParamsWithLine = gpuParamsWithLineVec[igpu];
      reset_ttt_all(gpuParamsWithLine->ttt_descd, gpuParamsWithLine->tttd);
#ifdef QMIG3D
      reset_ttt_all(gpuParamsWithLine->qqq_descd, gpuParamsWithLine->qqqd);
#endif
    }

    time(&toc);
    fprintf(stderr, "Time Elapse of Line-%d:%.0fs\n", line, difftime(toc, tic));
  } /* 循环成像线 结束 */

  sprintf(infopath, "%s/%d/%s", projpath, 0, projtask);
  sprintf(nodename, "%s", nodename0);

  /* 输出任务完成信息(共享盘) */
  sprintf(str, "%s/%s", infopath, nodename);
  fp = fopen(str, "w");
  if (fp != NULL) {
    time(&toc);
    fprintf(fplog, "%s\n", ctime(&toc));
    fclose(fp);
  }
  /* 释放内存和显存 */
  /* ... */
}

void *producer(void *arg) {
  PROD_PARAMS *prod_params = (PROD_PARAMS *)arg;
  safe_queue<segy *> *datapool = prod_params->datapool;
  int ioffl1 = prod_params->ioffl1;
  int ioffl2 = prod_params->ioffl2;
  char *offspath = prod_params->offspath;
  int ioffs = prod_params->ioffs;
  int ntl = prod_params->ntl;
  int ncons = prod_params->ncons;
  char *nodename = prod_params->nodename;

  char str[1024];
  struct stat st;
  FILE *fp;
  int ntr;
  int itr;
  char nodename0[128];
  sscanf(nodename, "%[^_]s", nodename0);

  unsigned char *buf;
  int bufsize = (100 << 20); /* 100MB */

  buf = (unsigned char *)malloc(sizeof(unsigned char) * bufsize);

  for (int ioffl = ioffl1; ioffl <= ioffl2; ++ioffl) {
    sprintf(str, "%s/off%d/line%d", offspath, ioffs, ioffl);
    if (!exist(str)) continue; /* 如果文件不存在进行下一个 */

    /* 读取某一条数据域的有贡献线 进行成像 */
    stat(str, &st);
    if (st.st_size > bufsize) {
      bufsize = st.st_size;
      free(buf);
      buf = (unsigned char *)malloc(bufsize);
      WARN("Buffer Reallocate Triggered nodename=%s ioffs=%d ioffl=%d\n",
           nodename, ioffs, ioffl);
    }

    fp = fopen(str, "rb");
    if (fread(buf, 1L, st.st_size, fp) != (size_t)st.st_size) {
      fprintf(stderr, "Read %s fail:%s in LINE:%d\n", str, strerror(errno),
              __LINE__);
      fprintf(fplog, "Read %s fail:%s in LINE:%d\n", str, strerror(errno),
              __LINE__);
      fflush(fplog);
      exit(1);
    }
    fclose(fp);

    ntr = st.st_size / (240L + sizeof(float) * ntl);
    for (itr = 0; itr < ntr; ++itr) {
      segy *trace = new segy;
      swaptr(&buf[itr * (240L + sizeof(float) * ntl)], trace);
      datapool->push(trace);
    }
  }

  free(buf);
  segy *poison = NULL;
  for (int igpu = 0; igpu < ncons; ++igpu) {
    datapool->push(poison);
  }

  return NULL;
}

void *consumer(void *arg) {
  CONS_PARAMS *cons_params = (CONS_PARAMS *)arg;
  safe_queue<segy *> *datapool = cons_params->datapool;

  CPUParamsNoLine *cpuParamsNoLine = cons_params->cpuParamsNoLine;
  float dyi = cpuParamsNoLine->dyi;
  float dxi = cpuParamsNoLine->dxi;
  int nxi = cpuParamsNoLine->nxi;
  int nzi = cpuParamsNoLine->nzi;
  int ntl = cpuParamsNoLine->ntl;
  double scaled = cpuParamsNoLine->scaled;
  int fci = cpuParamsNoLine->fci;
  float dzi = cpuParamsNoLine->dzi;
  int ndipx = cpuParamsNoLine->ndipx;
  int ndipy = cpuParamsNoLine->ndipy;
  float ddipx = cpuParamsNoLine->ddipx;
  float ddipy = cpuParamsNoLine->ddipy;
  float dt = cpuParamsNoLine->dt;
  float fc = cpuParamsNoLine->fc;
  float eps = cpuParamsNoLine->eps;
  int *iaztab = cpuParamsNoLine->iaztab;
  float fxv1 = cpuParamsNoLine->fxv1;
  float fyv1 = cpuParamsNoLine->fyv1;
  float exv1 = cpuParamsNoLine->exv1;
  float eyv1 = cpuParamsNoLine->eyv1;
  int nfft = cpuParamsNoLine->nfft;
  int nffti = cpuParamsNoLine->nffti;
  float *ww = cpuParamsNoLine->ww;
  float dt1 = cpuParamsNoLine->dt1;
  int ndipx2 = cpuParamsNoLine->ndipx2;
  int ndipy2 = cpuParamsNoLine->ndipy2;
  float lnG = cpuParamsNoLine->lnG;
  float coef1 = cpuParamsNoLine->coef1;
  float coef2 = cpuParamsNoLine->coef2;
  float coef3 = cpuParamsNoLine->coef3;
  float coef4 = cpuParamsNoLine->coef4;
  float taperzone = cpuParamsNoLine->taperzone;
  float dzii = cpuParamsNoLine->dzii;

  CPUParamsWithLine *cpuParamsWithLine = cons_params->cpuParamsWithLine;
  int ****aper = cpuParamsWithLine->aper;
  int **aznx1 = cpuParamsWithLine->aznx1;
  int **aznx2 = cpuParamsWithLine->aznx2;
  int **azny1 = cpuParamsWithLine->azny1;
  int **azny2 = cpuParamsWithLine->azny2;
  double fxi = cpuParamsWithLine->fxi;
  double cyi = cpuParamsWithLine->cyi;
  float dw = cpuParamsWithLine->dw;
  int ioffs = cpuParamsWithLine->ioffs;
#ifdef QMIG3D
  int itb1 = cpuParamsWithLine->itb1;
  int itb2 = cpuParamsWithLine->itb2;
  int nsb = cpuParamsWithLine->nsb;
  float tdstart1 = cpuParamsWithLine->tdstart1;
  float tdstart2 = cpuParamsWithLine->tdstart2;
  int nf1 = cpuParamsWithLine->nf1;
  float tdmid = cpuParamsWithLine->tdmid;
#endif
  int nww = cpuParamsWithLine->nww;
  int nxg = cpuParamsWithLine->nxg;
  int nzg = cpuParamsWithLine->nzg;
  int nxi1 = cpuParamsWithLine->nxi1;
  int nzi1 = cpuParamsWithLine->nzi1;
  int ttt_cdp1 = cpuParamsWithLine->ttt_cdp1;
  int ttt_cdp2 = cpuParamsWithLine->ttt_cdp2;
  int qqq_cdp1 = cpuParamsWithLine->qqq_cdp1;
  int qqq_cdp2 = cpuParamsWithLine->qqq_cdp2;
  int ncoef1 = cpuParamsWithLine->ncoef1;
  int ixstart = cpuParamsWithLine->ixstart;
  int ixend = cpuParamsWithLine->ixend;
  int nxb = cpuParamsWithLine->nxb;

  GPUParamsNoLine *gpuParamsNoLine = cons_params->gpuParamsNoLine;
#ifdef QMIG3D
  float *re = gpuParamsNoLine->re;
  float *im = gpuParamsNoLine->im;
  float *red = gpuParamsNoLine->red;
  float *imd = gpuParamsNoLine->imd;
#endif
  int igpu = gpuParamsNoLine->igpu;
  int *itibeg = gpuParamsNoLine->itibeg;
  int *itibegd = gpuParamsNoLine->itibegd;

  GPUParamsWithLine *gpuParamsWithLine = cons_params->gpuParamsWithLine;
#ifdef QMIG3D
  int *nf3d = gpuParamsWithLine->nf3d;
  float *Qd = gpuParamsWithLine->Qd;
  float *datav = gpuParamsWithLine->datav;
  fftwf_complex *wdatav = gpuParamsWithLine->wdatav;
  fftwf_plan planv = gpuParamsWithLine->planv;
  ttt_desc_t *qqq_descd = gpuParamsWithLine->qqq_descd;
  float *qqqd = gpuParamsWithLine->qqqd;
#endif
  float *data = gpuParamsWithLine->data;
  float *datad = gpuParamsWithLine->datad;
  float *dipx1d = gpuParamsWithLine->dipx1d;
  float *dipx2d = gpuParamsWithLine->dipx2d;
  float *dipy1d = gpuParamsWithLine->dipy1d;
  float *dipy2d = gpuParamsWithLine->dipy2d;
  ttt_desc_t *ttt_descd = gpuParamsWithLine->ttt_descd;
  float *tttd = gpuParamsWithLine->tttd;
  float *imgd = gpuParamsWithLine->imgd;

  while (true) {
    segy *trace;
    datapool->pop(trace);
    if (trace == NULL) {
      break;
    }
    psdm_kernel(
        *trace, dyi, dxi, nxi, nzi, ntl, scaled, fci, dzi, ndipx, ndipy, ddipx,
        ddipy, dt, fc, eps, iaztab, fxv1, fyv1, exv1, eyv1, nfft, nffti, ww,
        dt1, ndipx2, ndipy2, taperzone, lnG, coef1, coef2, coef3, coef4, dzii,
        aper, aznx1, aznx2, azny1, azny2, fxi, cyi, dw, ioffs,
#ifdef QMIG3D
        itb1, itb2, nsb, tdstart1, tdstart2, nf1, tdmid,
#endif
        nww, nxg, nzg, nxi1, nzi1, ttt_cdp1, ttt_cdp2, qqq_cdp1, qqq_cdp2,
        ncoef1, ixstart, ixend, nxb,
#ifdef QMIG3D
        re, im, red, imd,
#endif
        igpu, itibeg, itibegd,
#ifdef QMIG3D
        nf3d, Qd, datav, wdatav, planv, qqq_descd, qqqd,
#endif
        data, datad, dipx1d, dipx2d, dipy1d, dipy2d, ttt_descd, tttd, imgd);
    delete trace;
  }

  fprintf(fplog, "GPU-%d Finish Migration. RETURN\n", igpu);

  return NULL;
}

void psdm_kernel(segy tr,
                 // CPUParamsNoLine
                 float dyi, float dxi, int nxi, int nzi, int ntl, double scaled,
                 int fci, float dzi, int ndipx, int ndipy, float ddipx,
                 float ddipy, float dt, float fc, float eps, int *iaztab,
                 float fxv1, float fyv1, float exv1, float eyv1, int nfft,
                 int nffti, float *ww, float dt1, int ndipx2, int ndipy2,
                 float taperzone, float lnG, float coef1, float coef2,
                 float coef3, float coef4, float dzii,
                 // CPUParamsWithLine
                 int ****aper, int **aznx1, int **aznx2, int **azny1,
                 int **azny2, double fxi, double cyi, float dw, int ioffs,
#ifdef QMIG3D
                 int itb1, int itb2, int nsb, float tdstart1, float tdstart2,
                 int nf1, float tdmid,
#endif
                 int nww, int nxg, int nzg, int nxi1, int nzi1, int ttt_cdp1,
                 int ttt_cdp2, int qqq_cdp1, int qqq_cdp2, int ncoef1,
                 int ixstart, int ixend, int nxb,
// GPUParamsNoLine
#ifdef QMIG3D
                 float *re, float *im, float *red, float *imd,
#endif
                 int igpu, int *itibeg, int *itibegd,
// GPUParamsWithLine
#ifdef QMIG3D
                 int *nf3d, float *Qd, float *datav, fftwf_complex *wdatav,
                 fftwf_plan planv, ttt_desc_t *qqq_descd, float *qqqd,
#endif
                 float *data, float *datad, float *dipx1d, float *dipx2d,
                 float *dipy1d, float *dipy2d, ttt_desc_t *ttt_descd,
                 float *tttd, float *imgd) {
  cudaSetDevice(igpu);
  (void)nxi;
  (void)dt;
  (void)fc;
  (void)aper;
  (void)nww;
  (void)qqq_cdp1;
  (void)qqq_cdp2;
  (void)Qd;
  (void)ww;
  (void)nfft;

  double sx, sy, gx, gy, ddx, ddy, offset;
  int cdpx, cdpy;
  float cosa, cosx, c1;
  int i, iaz, izimin, ixi, ixi1, ixi2, iti;
  int ntibeg, naznx, naznx1, naznx2, nazny1, nazny2;

  cudaError_t err;

  sx = tr.sx * scaled;
  sy = tr.sy * scaled;
  gx = tr.gx * scaled;
  gy = tr.gy * scaled;

  if (sx < fxv1 || sx > exv1) return;
  if (gx < fxv1 || gx > exv1) return;
  if (sy < fyv1 || sy > eyv1) return;
  if (gy < fyv1 || gy > eyv1) return;

  ddx = gx - sx;
  ddy = gy - sy;

  offset = sqrtf(ddx * ddx + ddy * ddy);

  cdpx = (int)(((sx + gx) * 0.5 - fxi) / dxi + 0.5);
  cdpy = (int)(((sy + gy) * 0.5 - cyi) / dyi + 0.5);

  /* 由坐标求取其对应的方位角 */
  if (offset > 1.0) {
    cosa = ddx / offset;
    cosx = fabs(cosa);

    if (ddx * ddy < 0.0) cosx = -cosx;
    c1 = 1.0 - cosx;
    i = (int)(c1 / eps);
    iaz = iaztab[i];
  } else {
    iaz = 0;
  }

  /* 获取孔径参数 */
  naznx1 = aznx1[ioffs][iaz];
  naznx2 = aznx2[ioffs][iaz];
  nazny1 = azny1[ioffs][iaz];
  nazny2 = azny2[ioffs][iaz];
  naznx = naznx1 + naznx2;

  if (cdpy > nazny1 || cdpy < -nazny2) return;

  izimin = nzi;
  for (ixi = 0, ntibeg = 0; ixi < naznx; ++ixi) {
    iti = 0;
    if (iti < nzi - 10 && cdpx - naznx1 + ixi >= ixstart &&
        cdpx - naznx1 + ixi < ixend) {
      if (iti < izimin) {
        izimin = iti;
      }
      itibeg[ntibeg] = iti;
      ixi2 = cdpx - naznx1 + ixi;
      ++ntibeg;
    }
  }

  ixi1 = ixi2 - ntibeg + 1;

  if (ntibeg == 0) return;

#ifdef QMIG3D
  int iw;
  float w;
  for (int it = 0; it < ntl; ++it) {
    if (it * dt < 0.1)
      tr.data[it] *= 0.1 * 0.1;
    else
      tr.data[it] *= (it * dt) * (it * dt);
  }
  fetch_and_taper(tr.data, ntl, 0, ntl - 1, data);
  fetch_and_taper(data, ntl, itb1, itb1 + nsb - 1, datav);
  fftwf_execute(planv); /* tdata->wdata */

  if (is3D) /* 3-Dimension: Full Derivate */
  {
    for (iw = 0, w = nf1 * dw; iw < nww; ++iw, w += dw) {
      re[iw] = +wdatav[iw + nf1][1] * w;
      im[iw] = -wdatav[iw + nf1][0] * w;
    }
  } else /* 2-Dimension: Half Derivate */
  {
    for (iw = 0, w = nf1 * dw; iw < nww; ++iw, w += dw) {
      float cc = sqrtf(w) / sqrtf(2.0);
      re[iw] = cc * (wdatav[iw + nf1][0] + wdatav[iw + nf1][1]);
      im[iw] = cc * (wdatav[iw + nf1][1] - wdatav[iw + nf1][0]);
    }
  }

  fetch_and_taper(data, ntl, itb2, itb2 + nsb - 1, datav);
  fftwf_execute(planv); /* tdata->wdata */

  if (is3D) /* 3-Dimension: derivate */
  {
    for (iw = 0, w = nf1 * dw; iw < nww; ++iw, w += dw) {
      re[iw + nww] = +wdatav[iw + nf1][1] * w;
      im[iw + nww] = -wdatav[iw + nf1][0] * w;
    }
  } else /* 2-Dimension: derivate */
  {
    for (iw = 0, w = nf1 * dw; iw < nww; ++iw, w += dw) {
      float cc = sqrtf(w) / sqrtf(2.0);
      re[iw + nww] = cc * (wdatav[iw + nf1][1] + wdatav[iw + nf1][1]);
    }
  }

  cudaMemcpy(red, re, sizeof(float) * 2 * nww, cudaMemcpyHostToDevice);
  cudaMemcpy(imd, im, sizeof(float) * 2 * nww, cudaMemcpyHostToDevice);
#endif
  sx -= fxi;
  gx -= fxi;
  sy -= cyi;
  gy -= cyi;

  cudaMemcpy(itibegd, itibeg, sizeof(int) * ntibeg, cudaMemcpyHostToDevice);

  dim3 block(32, 16);
  dim3 grid((int)(ceil((nzi - izimin) / 32.0)), (int)(ceil((ntibeg) / 16.0)));

  image_depth_gpu<<<grid, block, nww * 4 * sizeof(float)>>>(
      imgd, izimin, nzi, datad, nffti, dt1, ixstart, ixend, nxb, dzii, red, imd,
      nf1, fc, dw, nww, nf3d, qqq_descd, qqqd, tdstart1, tdstart2, tdmid,
      taperzone, lnG, coef1, coef2, coef3, coef4, ncoef1, nxg, nzg, nxi1, nzi1,
      ttt_cdp1 - fci, ttt_cdp2 - fci, sx, sy, gx, gy, dxi, dzi, itibegd, ixi1,
      ixi2, ndipx, ndipx2, ddipx, ndipy, ndipy2, ddipy, ttt_descd, tttd, dipx1d,
      dipx2d, dipy1d, dipy2d);
  if ((err = cudaGetLastError()) != cudaSuccess) {
    fprintf(fplog, "Error: CUDA  Error in FILE:%s FUNCTION:%s LINE:%d\n",
            __FILE__, __FUNCTION__, __LINE__);
    fprintf(fplog, "       %s\n", cudaGetErrorString(err));
    exit(1);
  }
}

__global__ void image_depth_gpu(
    float *img, int izimin, int nzi, float *data, int nt, float dt, int ixstart,
    int ixend, int nxb, float dzii,
#ifdef QMIG3D
    const float *__restrict__ red, const float *__restrict__ imd, int nf1,
    float fc, float dw, int nww, int *nf3, ttt_desc_t *qqq_desc, float *qqq,
    float tstart1, float tstart2, float tmid,
#endif
    float taperzone, float lnG, float coef1, float coef2, float coef3,
    float coef4, int ncoef1, int nxg, int nzg, int nxi1, int nzi1, int ttt_cdp1,
    int ttt_cdp2, float sx, float sy, float gx, float gy, float dxi, float dzi,
    int *itibegd, int ixi1, int ixi2, int ndipx, int ndipx2, float ddipx,
    int ndipy, int ndipy2, float ddipy, ttt_desc_t *ttt_desc, float *ttt,
    float *dipx1, float *dipx2, float *dipy1, float *dipy2) {
  int izi = blockIdx.x * blockDim.x + threadIdx.x + izimin;
  int ix = blockIdx.y * blockDim.y + threadIdx.y;

  float ts, tg, t;
  int it1, it2;

  float x, y, z, dsx, dgx;
  float tanx, tany, tanx1, tanx2, tany1, tany2, tanx1c, tanx2c, tany1c, tany2c;

  float w, tr;

  float c1, c2;
  float ds, dg;
  float v;
  float z2;

  int ixg, izg;
  float tst;

#ifdef QMIG3D
  extern __shared__ float shared[];
  float *re, *im;

  re = &shared[0];
  im = &shared[nww * 2];
  for (int i = threadIdx.y * blockDim.x + threadIdx.x; i < nww * 2;
       i += blockDim.x * blockDim.y) {
    re[i] = red[i];
    im[i] = imd[i];
  }

  __syncthreads();
#endif

  int izii, izii1, izg1;

  int ixi = ixi1 + ix;

  if (ix >= nxb) return;
  if (izi >= nzi || izi < itibegd[ix] || ixi > ixi2 || ixi < ttt_cdp1 ||
      ixi > ttt_cdp2)
    return;

  int ixzi = ixi * nzi + izi;

  z = izi * dzi;
  z2 = z * z;
  x = ixi * dxi;

  dsx = sx - x;
  dgx = gx - x;

  /* dip */
  ds = sqrtf(dsx * dsx + sy * sy + z2);
  dg = sqrtf(dgx * dgx + gy * gy + z2);
  x = dsx / ds + dgx / dg;
  y = sy / ds + gy / dg;
  z = z / ds + z / dg;

  tanx = x / z;
  tany = y / z;

  tanx1 = dipx1[ixzi];
  tanx1c = tanx1 - taperzone * (1.0f + tanx1 * tanx1);
  tanx2 = dipx2[ixzi];
  tanx2c = tanx2 + taperzone * (1.0f + tanx2 * tanx2);

  tany1 = dipy1[ixzi];
  tany1c = tany1 - taperzone * (1.0f + tany1 * tany1);
  tany2 = dipy2[ixzi];
  tany2c = tany2 + taperzone * (1.0f + tany2 * tany2);

  if (tanx <= tanx1c || tanx2c <= tanx || tany <= tany1c || tany2c <= tany)
    return;

  w = 1.0f;

  if (tanx < tanx1) w *= (tanx - tanx1c) / (tanx1 - tanx1c);
  if (tanx > tanx2) w *= (tanx2c - tanx) / (tanx2c - tanx2);
  if (tany < tany1) w *= (tany - tany1c) / (tany1 - tany1c);
  if (tany > tany2) w *= (tany2c - tany) / (tany2c - tany2);

  /* find the grid index */
  z = izi * dzi;
  ixg = (ixi - ttt_cdp1) / nxi1;
  izii = (int)z / dzii;
  izg = izii / nzi1;

  c2 = z / dzii - izii;
  c1 = 1 - c2;
  ts = 1.0;
  tg = 1.0;
  t = ts + tg;

  it1 = t / dt;
  it2 = it1 + 1;

  if (it1 > 0 && it2 < nt) {
#ifdef QMIG3D
    t = ts + tg;

    if (t > tmid) {
      re = &shared[nww];
      im = &shared[3 * nww];
      tst = t - tstart2;
    }
    v = 100;
    if (t > 0.1)
      v /= t * t;
    else
      v *= 100.0f;
#endif
    tr = (ds / dg * ds / dg * ds / dg * tg / ts);
    tr = w * (tr + 1.0 / tr);
    v *= tr;

    ixi -= ixstart;
    ixzi = ixi * nzi + izi;
    img[ixzi] += v;
  }
}

int get_ttt_all(const char *ttpath, const char *torq, int line, int line1,
                int line2, int *nxg, int *nzg, int *nxi1, int *nzi1,
                int *ttt_cdp1, int *ttt_cdp2,
                std::vector<GPUParamsWithLine *> gpuParamsWithLineVec,
                int ncons) {
  char index_file[4096], data_file[4096];
  char index_suff[10], data_suff[10];
  ttt_desc_t *ttt_desc;
  FILE *fp;
  int nxg0, nzg0, nxi0, nzi0;
  ssize_t ret;

  float *ttt = NULL;
  float *ttt1;
  size_t size = 0;

  int cdp1, cdp2;
  float dx, dz;

  // choose ttt or qqq
  if (strcmp(torq, "t") == 0) {
    sprintf(index_suff, "tidx");
    sprintf(data_suff, "ttt");
  } else {
    sprintf(index_suff, "qidx");
    sprintf(data_suff, "qqq");
  }

  // get index
  // =======================================================
  sprintf(index_file, "%s/Lxxx.%s", ttpath, index_suff);
  printf("load index file: %s\n", index_file);
  fp = fopen(index_file, "rb");

  if (fread(&cdp1, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "Error: Read cdp1 in %s\n", index_file);
  };
  if (fread(&cdp2, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "Error: Read cdp2 in %s\n", index_file);
  };
  if (fread(&dx, sizeof(float), 1, fp) != 1) {
    fprintf(stderr, "Error: Read dx in %s\n", index_file);
  };
  if (fread(&dz, sizeof(float), 1, fp) != 1) {
    fprintf(stderr, "Error: Read dz in %s\n", index_file);
  };
  if (fread(&nzg0, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "Error: Read nzg0 in %s\n", index_file);
  };
  if (fread(&nxg0, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "Error: Read nxg0 in %s\n", index_file);
  };
  if (fread(&nzi0, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "Error: Read nzi0 in %s\n", index_file);
  };
  if (fread(&nxi0, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "Error: Read nxi0 in %s\n", index_file);
  };

  *nxg = nxg0;
  *nzg = nzg0;
  *nxi1 = nxi0;
  *nzi1 = nzi0;

  *ttt_cdp1 = cdp1;
  *ttt_cdp2 = cdp2;

  printf("line = %d, nzg=%d, nxg=%d, nzi0=%d, nxi0=%d\n", line, nzg0, nxg0,
         nzi0, nxi0);
  printf("ttt_cdp1 = %d, ttt_cdp2 = %d\n", cdp1, cdp2);

  ttt_desc = (ttt_desc_t *)malloc(sizeof(ttt_desc_t) * nxg0 * nzg0);
  ret = fread(ttt_desc, sizeof(ttt_desc_t), nxg0 * nzg0, fp);
  fclose(fp);

  if (ret != (ssize_t)(nxg0 * nzg0)) {
    fprintf(stderr, " travle time table or qfield index size error %s\n",
            index_file);
    exit(1);
  }

  for (int igpu = 0; igpu < ncons; ++igpu) {
    cudaSetDevice(igpu);

    GPUParamsWithLine *gpuParamsWithLine = gpuParamsWithLineVec[igpu];

    if (strcmp(torq, "t") == 0) {
      ttt_desc_t *ttt_descd;
      cudaMalloc(&ttt_descd, sizeof(ttt_desc_t) * nxg0 * nzg0);
      cudaMemcpy(ttt_descd, ttt_desc, sizeof(ttt_desc_t) * nxg0 * nzg0,
                 cudaMemcpyHostToDevice);
      gpuParamsWithLine->ttt_descd = ttt_descd;
    } else {
      ttt_desc_t *qqq_descd;
      cudaMalloc(&qqq_descd, sizeof(ttt_desc_t) * nxg0 * nzg0);
      cudaMemcpy(qqq_descd, ttt_desc, sizeof(ttt_desc_t) * nxg0 * nzg0,
                 cudaMemcpyHostToDevice);
      gpuParamsWithLine->qqq_descd = qqq_descd;
    }
  }

  free(ttt_desc);

  // load travel time table bin/or qfield
  // ========================================================

  int line1find = -1, line2find = -1;
  // find line -> line1
  for (int iline = line; iline >= line1; --iline) {
    sprintf(data_file, "%s/L%d.%s", ttpath, iline, data_suff);

    if (fexist(data_file)) {
      line1find = iline;
      break;
    }
  }

  // find line -> line2
  for (int iline = line; iline <= line2; ++iline) {
    sprintf(data_file, "%s/L%d.%s", ttpath, iline, data_suff);

    if (fexist(data_file)) {
      line2find = iline;
      break;
    }
  }

  printf("%d -> %d -> %d\n", line1find, line, line2find);

  fprintf(stderr, "load travel time table success: line %d\n", line);

  for (int igpu = 0; igpu < ncons; ++igpu) {
    cudaSetDevice(igpu);

    GPUParamsWithLine *gpuParamsWithLine = gpuParamsWithLineVec[igpu];

    if (strcmp(torq, "t") == 0) {
      float *tttd;
      cudaMalloc(&tttd, size);
      cudaMemcpy(tttd, ttt, size, cudaMemcpyHostToDevice);
      gpuParamsWithLine->tttd = tttd;
    } else {
      float *qqqd;
      cudaMalloc(&qqqd, size);
      cudaMemcpy(qqqd, ttt, size, cudaMemcpyHostToDevice);
      gpuParamsWithLine->qqqd = qqqd;
    }
  }

  if (line1find == line && line2find == line) {
    // load travel time table line1find only
    sprintf(data_file, "%s/L%d.%s", ttpath, line, data_suff);
    printf("data_file = %s\n", data_file);

    size = fsize(data_file);
    ttt = (float *)malloc(size);

    fp = fopen(data_file, "rb");
    ret = fread(ttt, sizeof(char), size, fp);
    fclose(fp);

    if (ret != (ssize_t)size) {
      fprintf(stderr, " travle time table size error\n");
      exit(1);
    }
  } else if (line1find != -1 && line2find != -1) {
    // load travel time table line1find
    sprintf(data_file, "%s/L%d.%s", ttpath, line1find, data_suff);
    printf("data_file = %s\n", data_file);

    size = fsize(data_file);
    ttt = (float *)malloc(size);

    fp = fopen(data_file, "rb");
    ret = fread(ttt, sizeof(char), size, fp);
    fclose(fp);

    if (ret != (ssize_t)size) {
      fprintf(stderr, " travle time table size error\n");
      exit(1);
    }

    // load travel time table line2find
    sprintf(data_file, "%s/L%d.%s", ttpath, line2find, data_suff);
    printf("data_file = %s\n", data_file);

    if (size != fsize(data_file)) {
      printf("travel time table size not confirm\n");
    }
    ttt1 = (float *)malloc(size);

    fp = fopen(data_file, "rb");
    ret = fread(ttt1, sizeof(char), size, fp);
    fclose(fp);

    if (ret != (ssize_t)size) {
      fprintf(stderr, " travle time table size error\n");
      exit(1);
    }

    float c1, c2, t, t1;

    c2 = 1.0 * (line - line1find) / (line2find - line1find);
    c1 = 1.0 - c2;

    for (size_t i = 0; i < size / 4; ++i) {
      t = ttt[i];
      t1 = ttt1[i];

      t = t * c1 + t1 * c2;

      ttt[i] = t;
    }

    free(ttt1);
  } else {
    printf("no travel time table exists\n");
    exit(1);
  }

  fprintf(stderr, "load travel time table success: line %d\n", line);

  for (int igpu = 0; igpu < ncons; ++igpu) {
    cudaSetDevice(igpu);

    GPUParamsWithLine *gpuParamsWithLine = gpuParamsWithLineVec[igpu];

    if (strcmp(torq, "t") == 0) {
      cudaMalloc(&(gpuParamsWithLine->tttd), size);
      cudaMemcpy(gpuParamsWithLine->tttd, ttt, size, cudaMemcpyHostToDevice);
    } else {
      cudaMalloc(&(gpuParamsWithLine->qqqd), size);
      cudaMemcpy(gpuParamsWithLine->qqqd, ttt, size, cudaMemcpyHostToDevice);
    }
  }

  free(ttt);

  return 0;
  // ========================================================
}

void reset_ttt_all(ttt_desc_t *desc_out, float *data_out) {
  cudaFree(desc_out);
  cudaFree(data_out);
}
