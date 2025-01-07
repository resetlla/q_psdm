#include "segy.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>
#include "omp.h"
#include "fft.h"
// #include "psdmpkg.h"

#include <vector>

#define MSK_MANT_IBM 0XFFFFFF
#define MSK_EXP_IBM 0X7F000000
#define MSK_SIGN_IBM 0X80000000
#define MSK_NORM_IEEE 0X800000
#define MSK_NO_SIGN 0X7FFFFFFF
#define ALFA 0.1
#define EV0_MAX 4000  // Effective velocity maximum value
#define EV0_MIN 2000  // Effective velocity minimun value
#define ETA_MAX 0.2   // Effective anisotropic parameter maximum value
#define ETA_MIN 0.0   // Effective anisotropic parameter minimun value
#define EV0_APERTURE EV0_MAX - EV0_MIN
#define ETA_APERTURE ETA_MAX - ETA_MIN

#define signal(s) ((s < 0) ? (-1.) : (1.))

/* Internal Structure */
/* Only Accessed Internal */
typedef union {
  int i;
  unsigned char s[4];
} int_union;
typedef union {
  short i;
  unsigned char s[2];
} short_union;
typedef union {
  unsigned short i;
  unsigned char s[2];
} ushort_union;
typedef union {
  float f;
  int i;
  unsigned char s[4];
} float_union;

int keyindex_inited;

static struct {
  const char *key;
  const char *type;
  int offs;
} hdr[] = {
    {"tracl", "i", 0},    {"tracr", "i", 4},      {"fldr", "i", 8},
    {"tracf", "i", 12},   {"ep", "i", 16},        {"cdp", "i", 20},
    {"cdpt", "i", 24},    {"trid", "h", 28},      {"nvs", "h", 30},
    {"nhs", "h", 32},     {"duse", "h", 34},      {"offset", "i", 36},
    {"gelev", "i", 40},   {"selev", "i", 44},     {"sdepth", "i", 48},
    {"gdel", "i", 52},    {"sdel", "i", 56},      {"swdep", "i", 60},
    {"gwdep", "i", 64},   {"scalel", "h", 68},    {"scalco", "h", 70},
    {"sx", "i", 72},      {"sy", "i", 76},        {"gx", "i", 80},
    {"gy", "i", 84},      {"counit", "h", 88},    {"wevel", "h", 90},
    {"swevel", "h", 92},  {"sut", "h", 94},       {"gut", "h", 96},
    {"sstat", "h", 98},   {"gstat", "h", 100},    {"tstat", "h", 102},
    {"laga", "h", 104},   {"lagb", "h", 106},     {"delrt", "h", 108},
    {"muts", "h", 110},   {"mute", "h", 112},     {"ns", "u", 114},
    {"dt", "u", 116},     {"gain", "h", 118},     {"igc", "h", 120},
    {"igi", "h", 122},    {"corr", "h", 124},     {"sfs", "h", 126},
    {"sfe", "h", 128},    {"slen", "h", 130},     {"styp", "h", 132},
    {"stas", "h", 134},   {"stae", "h", 136},     {"tatyp", "h", 138},
    {"afilf", "h", 140},  {"afils", "h", 142},    {"nofilf", "h", 144},
    {"nofils", "h", 146}, {"lcf", "h", 148},      {"hcf", "h", 150},
    {"lcs", "h", 152},    {"hcs", "h", 154},      {"year", "h", 156},
    {"day", "h", 158},    {"hour", "h", 160},     {"minute", "h", 162},
    {"sec", "h", 164},    {"timbas", "h", 166},   {"trwf", "h", 168},
    {"grnors", "h", 170}, {"grnofr", "h", 172},   {"grnlof", "h", 174},
    {"gaps", "h", 176},   {"otrav", "h", 178},    {"cdpx", "i", 180},
    {"cdpy", "i", 184},   {"iline", "i", 188},    {"xline", "i", 192},
    {"sp", "i", 196},     {"unscale", "i", 200},  {"ntr", "i", 204},
    {"mark", "h", 208},   {"shortpad", "h", 210},
};

void fsgethd(FILE *fp, hd3600 *hd) {
  static int_union i0, i1;
  static short_union h0, h1;
  static unsigned char ebcdic2ascii[256] = {
      0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20,
      0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20,
      0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20,
      0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20,
      0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20,
      0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20,
      0X20, 0X20, 0X5B, 0X2E, 0X3C, 0X28, 0X2B, 0X21, 0X26, 0X20, 0X20, 0X20,
      0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X5D, 0X24, 0X2A, 0X29, 0X3B, 0X5E,
      0X2D, 0X2F, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0X7C, 0X2C,
      0X25, 0X5F, 0X3E, 0X3F, 0X20, 0X20, 0X20, 0X20, 0X20, 0X20, 0XEE, 0XA0,
      0XA1, 0X60, 0X3A, 0X23, 0X40, 0X27, 0X3D, 0X22, 0XE6, 0X61, 0X62, 0X63,
      0X64, 0X65, 0X66, 0X67, 0X68, 0X69, 0XA4, 0XA5, 0XE4, 0XA3, 0XE5, 0XA8,
      0XA9, 0X6A, 0X6B, 0X6C, 0X6D, 0X6E, 0X6F, 0X70, 0X71, 0X72, 0XAA, 0XAB,
      0XAC, 0XAD, 0XAE, 0XAF, 0XEF, 0X7E, 0X73, 0X74, 0X75, 0X76, 0X77, 0X78,
      0X79, 0X7A, 0XE0, 0XE1, 0XE2, 0XE3, 0XA6, 0XA2, 0XEC, 0XEB, 0XA7, 0XE8,
      0XED, 0XE9, 0XE7, 0XEA, 0X9E, 0X80, 0X81, 0X96, 0X84, 0X85, 0X94, 0X83,
      0X7B, 0X41, 0X42, 0X43, 0X44, 0X45, 0X46, 0X47, 0X48, 0X49, 0X95, 0X88,
      0X89, 0X8A, 0X8B, 0X8C, 0X7D, 0X4A, 0X4B, 0X4C, 0X4D, 0X4E, 0X4F, 0X50,
      0X51, 0X52, 0X8D, 0X8E, 0X8F, 0X9F, 0X90, 0X91, 0X5C, 0X20, 0X53, 0X54,
      0X55, 0X56, 0X57, 0X58, 0X59, 0X5A, 0X92, 0X93, 0X86, 0X82, 0X9C, 0X9B,
      0X30, 0X31, 0X32, 0X33, 0X34, 0X35, 0X36, 0X37, 0X38, 0X39, 0X87, 0X98,
      0X9D, 0X99, 0X97, 0X20};

  int i;

  int r = fread(hd, 1, 3600, fp);
  if (r != 3600) {
    fprintf(
        stderr,
        "read segy 3600 header fail, we try to get 3600 byte, but only got %d",
        r);
  }

  /* Convert EBCDIC to ASCII */
  for (i = 0; i < 3200; ++i) {
    hd->text[i] = ebcdic2ascii[hd->text[i]];
  }

  /* Swap byte order of Binary File Header */
#if defined(BHDR_ALL) || defined(BHDR_JOBID)
  i0.i = hd->jobid;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  hd->jobid = i1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_LINO)
  i0.i = hd->lino;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  hd->lino = i1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_RENO)
  i0.i = hd->reno;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  hd->reno = i1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_NTRPR)
  h0.i = hd->ntrpr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->ntrpr = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_NART)
  h0.i = hd->nart;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->nart = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HDT)
  h0.i = hd->hdt;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->hdt = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_DTO)
  h0.i = hd->dto;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->dto = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HNS)
  h0.i = hd->hns;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->hns = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_NSO)
  h0.i = hd->nso;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->nso = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_FORMAT)
  h0.i = hd->format;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->format = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_FOLD)
  h0.i = hd->fold;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->fold = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_TSORT)
  h0.i = hd->tsort;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->tsort = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_VSCODE)
  h0.i = hd->vscode;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->vscode = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSFS)
  h0.i = hd->hsfs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->hsfs = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSFE)
  h0.i = hd->hsfe;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->hsfe = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSLEN)
  h0.i = hd->hslen;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->hslen = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSTYP)
  h0.i = hd->hstyp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->hstyp = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_SCHN)
  h0.i = hd->schn;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->schn = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSTAS)
  h0.i = hd->hstas;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->hstas = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSTAE)
  h0.i = hd->hstae;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->hstae = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HTATYP)
  h0.i = hd->htatyp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->htatyp = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HCORR)
  h0.i = hd->hcorr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->hcorr = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_BGRCV)
  h0.i = hd->bgrcv;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->bgrcv = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_RCVM)
  h0.i = hd->rcvm;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->rcvm = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_MFEET)
  h0.i = hd->mfeet;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->mfeet = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_POLYT)
  h0.i = hd->polyt;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->polyt = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_VPOL)
  h0.i = hd->vpol;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd->vpol = h1.i;
#endif
}

void initkeyindex(keyindex *keyidx) {
  int i;

  for (i = 0; i < SEGY_TOTAL_KEY_NUM; ++i) {
    keyidx[i].start = hdr[i].offs + 1;
    keyidx[i].bytes = hdr[i + 1].offs - hdr[i].offs;
  }
  keyindex_inited = 1;
}

void setkeyindex(keyindex *keyidx, const char *key, int start, int bytes) {
  int i;

  for (i = 0; i < 80; ++i) {
    if (strncmp(key, hdr[i].key, strlen(key)) == 0) {
      keyidx[i].start = start;
      keyidx[i].bytes = bytes;
      break;
    }
  }
  if (i == 81) {
    fprintf(stderr, "Error: key=%s is not found in hdr\n", key);
    exit(0);
  }
  keyindex_inited = 1;
}

/*http://crack.seismo.unr.edu/ftp/pub/louie/proposals/scec/ibm2ieee.c*/
int fsgettr(FILE *fp, segy *tr, short int format, const keyindex *keyidx) {
  int_union i0, i1;
  short_union h0;
  ushort_union u0;
  static unsigned char buf[240];
  register int itr, ibuf, exponent, mantissa, sign, fpn;
  int ns, ret;
  int *p;

  if (!keyindex_inited) {
    fprintf(stderr,
            "Error: initkeyindex must be called before calling fsgettr\n");
    exit(0);
  }

  ret = fread(buf, 240, 1, fp);
  if (ret != 1) return 0;

#if defined(SEGY_ALL) || defined(SEGY_TRACL)
  ibuf = keyidx[0].start - 1;
  if (keyidx[0].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->tracl = i0.i;
  } else if (keyidx[0].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->tracl = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACR)
  ibuf = keyidx[1].start - 1;
  if (keyidx[1].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->tracr = i0.i;
  } else if (keyidx[1].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->tracr = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_FLDR)
  ibuf = keyidx[2].start - 1;
  if (keyidx[2].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->fldr = i0.i;
  } else if (keyidx[2].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->fldr = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACF)
  ibuf = keyidx[3].start - 1;
  if (keyidx[3].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->tracf = i0.i;
  } else if (keyidx[3].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->tracf = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_EP)
  ibuf = keyidx[4].start - 1;
  if (keyidx[4].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->ep = i0.i;
  } else if (keyidx[4].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->ep = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDP)
  ibuf = keyidx[5].start - 1;
  if (keyidx[5].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->cdp = i0.i;
  } else if (keyidx[5].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->cdp = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPT)
  ibuf = keyidx[6].start - 1;
  if (keyidx[6].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->cdpt = i0.i;
  } else if (keyidx[6].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->cdpt = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRID)
  ibuf = keyidx[7].start - 1;
  if (keyidx[7].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->trid = h0.i;
  } else if (keyidx[7].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->trid = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NVS)
  ibuf = keyidx[8].start - 1;
  if (keyidx[8].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->nvs = h0.i;
  } else if (keyidx[8].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->nvs = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NHS)
  ibuf = keyidx[9].start - 1;
  if (keyidx[9].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->nhs = h0.i;
  } else if (keyidx[9].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->nhs = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DUSE)
  ibuf = keyidx[10].start - 1;
  if (keyidx[10].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->duse = h0.i;
  } else if (keyidx[10].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->duse = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_OFFSET)
  ibuf = keyidx[11].start - 1;
  if (keyidx[11].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->offset = i0.i;
  } else if (keyidx[11].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->offset = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GELEV)
  ibuf = keyidx[12].start - 1;
  if (keyidx[12].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gelev = i0.i;
  } else if (keyidx[12].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gelev = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SELEV)
  ibuf = keyidx[13].start - 1;
  if (keyidx[13].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->selev = i0.i;
  } else if (keyidx[13].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->selev = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEPTH)
  ibuf = keyidx[14].start - 1;
  if (keyidx[14].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sdepth = i0.i;
  } else if (keyidx[14].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sdepth = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GDEL)
  ibuf = keyidx[15].start - 1;
  if (keyidx[15].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gdel = i0.i;
  } else if (keyidx[15].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gdel = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEL)
  ibuf = keyidx[16].start - 1;
  if (keyidx[16].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sdel = i0.i;
  } else if (keyidx[16].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sdel = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWDEP)
  ibuf = keyidx[17].start - 1;
  if (keyidx[17].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->swdep = i0.i;
  } else if (keyidx[17].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->swdep = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GWDEP)
  ibuf = keyidx[18].start - 1;
  if (keyidx[18].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gwdep = i0.i;
  } else if (keyidx[18].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gwdep = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALEL)
  ibuf = keyidx[19].start - 1;
  if (keyidx[19].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->scalel = h0.i;
  } else if (keyidx[19].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->scalel = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALCO)
  ibuf = keyidx[20].start - 1;
  if (keyidx[20].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->scalco = h0.i;
  } else if (keyidx[20].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->scalco = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SX)
  ibuf = keyidx[21].start - 1;
  if (keyidx[21].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sx = i0.i;
  } else if (keyidx[21].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sx = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SY)
  ibuf = keyidx[22].start - 1;
  if (keyidx[22].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sy = i0.i;
  } else if (keyidx[22].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sy = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GX)
  ibuf = keyidx[23].start - 1;
  if (keyidx[23].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gx = i0.i;
  } else if (keyidx[23].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gx = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GY)
  ibuf = keyidx[24].start - 1;
  if (keyidx[24].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gy = i0.i;
  } else if (keyidx[24].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gy = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_COUNIT)
  ibuf = keyidx[25].start - 1;
  if (keyidx[25].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->counit = h0.i;
  } else if (keyidx[25].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->counit = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_WEVEL)
  ibuf = keyidx[26].start - 1;
  if (keyidx[26].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->wevel = h0.i;
  } else if (keyidx[26].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->wevel = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWEVEL)
  ibuf = keyidx[27].start - 1;
  if (keyidx[27].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->swevel = h0.i;
  } else if (keyidx[27].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->swevel = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SUT)
  ibuf = keyidx[28].start - 1;
  if (keyidx[28].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sut = h0.i;
  } else if (keyidx[28].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sut = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GUT)
  ibuf = keyidx[29].start - 1;
  if (keyidx[29].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gut = h0.i;
  } else if (keyidx[29].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gut = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SSTAT)
  ibuf = keyidx[30].start - 1;
  if (keyidx[30].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sstat = h0.i;
  } else if (keyidx[30].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sstat = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GSTAT)
  ibuf = keyidx[31].start - 1;
  if (keyidx[31].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gstat = h0.i;
  } else if (keyidx[31].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gstat = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TSTAT)
  ibuf = keyidx[32].start - 1;
  if (keyidx[32].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->tstat = h0.i;
  } else if (keyidx[32].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->tstat = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGA)
  ibuf = keyidx[33].start - 1;
  if (keyidx[33].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->laga = h0.i;
  } else if (keyidx[33].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->laga = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGB)
  ibuf = keyidx[34].start - 1;
  if (keyidx[34].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->lagb = h0.i;
  } else if (keyidx[34].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->lagb = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DELRT)
  ibuf = keyidx[35].start - 1;
  if (keyidx[35].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->delrt = h0.i;
  } else if (keyidx[35].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->delrt = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTS)
  ibuf = keyidx[36].start - 1;
  if (keyidx[36].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->muts = h0.i;
  } else if (keyidx[36].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->muts = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTE)
  ibuf = keyidx[37].start - 1;
  if (keyidx[37].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->mute = h0.i;
  } else if (keyidx[37].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->mute = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NS)
  ibuf = keyidx[38].start - 1;
  if (keyidx[38].bytes == 2) {
    u0.s[1] = buf[ibuf];
    u0.s[0] = buf[ibuf + 1];
    tr->ns = u0.i;
  } else if (keyidx[38].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->ns = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DT)
  ibuf = keyidx[39].start - 1;
  if (keyidx[39].bytes == 2) {
    u0.s[1] = buf[ibuf];
    u0.s[0] = buf[ibuf + 1];
    tr->dt = u0.i;
  } else if (keyidx[39].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->dt = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAIN)
  ibuf = keyidx[40].start - 1;
  if (keyidx[40].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gain = h0.i;
  } else if (keyidx[40].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gain = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGC)
  ibuf = keyidx[41].start - 1;
  if (keyidx[41].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->igc = h0.i;
  } else if (keyidx[41].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->igc = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGI)
  ibuf = keyidx[42].start - 1;
  if (keyidx[42].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->igi = h0.i;
  } else if (keyidx[42].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->igi = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CORR)
  ibuf = keyidx[43].start - 1;
  if (keyidx[43].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->corr = h0.i;
  } else if (keyidx[43].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->corr = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFS)
  ibuf = keyidx[44].start - 1;
  if (keyidx[44].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sfs = h0.i;
  } else if (keyidx[44].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sfs = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFE)
  ibuf = keyidx[45].start - 1;
  if (keyidx[45].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sfe = h0.i;
  } else if (keyidx[45].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sfe = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SLEN)
  ibuf = keyidx[46].start - 1;
  if (keyidx[46].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->slen = h0.i;
  } else if (keyidx[46].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->slen = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_STYP)
  ibuf = keyidx[47].start - 1;
  if (keyidx[47].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->styp = h0.i;
  } else if (keyidx[47].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->styp = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAS)
  ibuf = keyidx[48].start - 1;
  if (keyidx[48].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->stas = h0.i;
  } else if (keyidx[48].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->stas = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAE)
  ibuf = keyidx[49].start - 1;
  if (keyidx[49].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->stae = h0.i;
  } else if (keyidx[49].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->stae = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TATYP)
  ibuf = keyidx[50].start - 1;
  if (keyidx[50].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->tatyp = h0.i;
  } else if (keyidx[50].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->tatyp = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILF)
  ibuf = keyidx[51].start - 1;
  if (keyidx[51].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->afilf = h0.i;
  } else if (keyidx[51].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->afilf = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILS)
  ibuf = keyidx[52].start - 1;
  if (keyidx[52].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->afils = h0.i;
  } else if (keyidx[52].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->afils = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILF)
  ibuf = keyidx[53].start - 1;
  if (keyidx[53].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->nofilf = h0.i;
  } else if (keyidx[53].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->nofilf = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILS)
  ibuf = keyidx[54].start - 1;
  if (keyidx[54].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->nofils = h0.i;
  } else if (keyidx[54].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->nofils = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCF)
  ibuf = keyidx[55].start - 1;
  if (keyidx[55].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->lcf = h0.i;
  } else if (keyidx[55].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->lcf = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCF)
  ibuf = keyidx[56].start - 1;
  if (keyidx[56].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->hcf = h0.i;
  } else if (keyidx[56].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->hcf = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCS)
  ibuf = keyidx[57].start - 1;
  if (keyidx[57].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->lcs = h0.i;
  } else if (keyidx[57].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->lcs = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCS)
  ibuf = keyidx[58].start - 1;
  if (keyidx[58].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->hcs = h0.i;
  } else if (keyidx[58].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->hcs = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_YEAR)
  ibuf = keyidx[59].start - 1;
  if (keyidx[59].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->year = h0.i;
  } else if (keyidx[59].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->year = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DAY)
  ibuf = keyidx[60].start - 1;
  if (keyidx[60].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->day = h0.i;
  } else if (keyidx[60].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->day = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_HOUR)
  ibuf = keyidx[61].start - 1;
  if (keyidx[61].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->hour = h0.i;
  } else if (keyidx[61].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->hour = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_MINUTE)
  ibuf = keyidx[62].start - 1;
  if (keyidx[62].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->minute = h0.i;
  } else if (keyidx[62].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->minute = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SEC)
  ibuf = keyidx[63].start - 1;
  if (keyidx[63].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sec = h0.i;
  } else if (keyidx[63].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sec = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TIMBAS)
  ibuf = keyidx[64].start - 1;
  if (keyidx[64].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->timbas = h0.i;
  } else if (keyidx[64].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->timbas = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRWF)
  ibuf = keyidx[65].start - 1;
  if (keyidx[65].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->trwf = h0.i;
  } else if (keyidx[65].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->trwf = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNORS)
  ibuf = keyidx[66].start - 1;
  if (keyidx[66].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->grnors = h0.i;
  } else if (keyidx[66].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->grnors = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNOFR)
  ibuf = keyidx[67].start - 1;
  if (keyidx[67].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->grnofr = h0.i;
  } else if (keyidx[67].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->grnofr = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNLOF)
  ibuf = keyidx[68].start - 1;
  if (keyidx[68].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->grnlof = h0.i;
  } else if (keyidx[68].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->grnlof = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAPS)
  ibuf = keyidx[69].start - 1;
  if (keyidx[69].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gaps = h0.i;
  } else if (keyidx[69].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gaps = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_OTRAV)
  ibuf = keyidx[70].start - 1;
  if (keyidx[70].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->otrav = h0.i;
  } else if (keyidx[70].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->otrav = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPX)
  ibuf = keyidx[71].start - 1;
  if (keyidx[71].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->cdpx = h0.i;
  } else if (keyidx[71].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->cdpx = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPY)
  ibuf = keyidx[72].start - 1;
  if (keyidx[72].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->cdpy = h0.i;
  } else if (keyidx[72].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->cdpy = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_ILINE)
  ibuf = keyidx[73].start - 1;
  if (keyidx[73].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->iline = h0.i;
  } else if (keyidx[73].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->iline = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_XLINE)
  ibuf = keyidx[74].start - 1;
  if (keyidx[74].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->xline = h0.i;
  } else if (keyidx[74].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->xline = i0.i;
  }
#endif

  ns = tr->ns;

  ret = fread((tr->data), sizeof(float), ns, fp);
  if (ret != ns) return 0;

  p = (int *)(tr->data);
  /* 4-byte IBM float-point */
  if (format == 1) {
    for (itr = 0; itr < ns; ++itr) {
      i0.i = p[itr];
      i1.s[0] = i0.s[3];
      i1.s[1] = i0.s[2];
      i1.s[2] = i0.s[1];
      i1.s[3] = i0.s[0];

      fpn = i1.i;

      if ((fpn & MSK_MANT_IBM) == 0) {
        fpn = 0;
      } else {
        mantissa = (fpn & MSK_MANT_IBM);
        exponent = ((((fpn & MSK_EXP_IBM) >> 24) - 64) << 2);
        sign = (fpn & MSK_SIGN_IBM);
        while ((mantissa & MSK_NORM_IEEE) == 0) {
          mantissa <<= 1; /* normalize */
          exponent--;
        }
        mantissa = mantissa & 0x7fffff; /* shift understood one out */
        exponent += 126;
        if ((exponent < 0) || (exponent > 255)) {
          fprintf(stderr, "IBM floating point exponent out of range \n");
          exponent = 0;
        }
        exponent <<= 23;
        fpn = sign | exponent | mantissa;
      }
      p[itr] = fpn;
    }
  } else if (format == 5) /* 4-byte IEEE float-point */
  {
    for (itr = 0; itr < ns; ++itr) {
      i0.i = p[itr];
      i1.s[0] = i0.s[3];
      i1.s[1] = i0.s[2];
      i1.s[2] = i0.s[1];
      i1.s[3] = i0.s[0];
      p[itr] = i1.i;
    }
  } else {
    fprintf(stderr, "Unsupported SEGY format=%d\n", format);
    return 0;
  }

  return ns;
}

/*http://crack.seismo.unr.edu/ftp/pub/louie/proposals/scec/ibm2ieee.c*/
int ssgettr(unsigned char *buf, segy *tr, short int format,
            const keyindex *keyidx) {
  int_union i0, i1;
  short_union h0;
  ushort_union u0;
  /*	static unsigned char buf[240];*/
  register int itr, ibuf, exponent, mantissa, sign, fpn;
  int ns;
  int *p;

  if (!keyindex_inited) {
    fprintf(stderr,
            "Error: initkeyindex must be called before calling fsgettr\n");
    exit(0);
  }

  /*	ret=fread(buf,240,1,fp);
          if (ret!=1) return 0;
  */

#if defined(SEGY_ALL) || defined(SEGY_TRACL)
  ibuf = keyidx[0].start - 1;
  if (keyidx[0].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->tracl = i0.i;
  } else if (keyidx[0].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->tracl = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACR)
  ibuf = keyidx[1].start - 1;
  if (keyidx[1].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->tracr = i0.i;
  } else if (keyidx[1].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->tracr = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_FLDR)
  ibuf = keyidx[2].start - 1;
  if (keyidx[2].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->fldr = i0.i;
  } else if (keyidx[2].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->fldr = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACF)
  ibuf = keyidx[3].start - 1;
  if (keyidx[3].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->tracf = i0.i;
  } else if (keyidx[3].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->tracf = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_EP)
  ibuf = keyidx[4].start - 1;
  if (keyidx[4].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->ep = i0.i;
  } else if (keyidx[4].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->ep = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDP)
  ibuf = keyidx[5].start - 1;
  if (keyidx[5].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->cdp = i0.i;
  } else if (keyidx[5].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->cdp = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPT)
  ibuf = keyidx[6].start - 1;
  if (keyidx[6].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->cdpt = i0.i;
  } else if (keyidx[6].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->cdpt = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRID)
  ibuf = keyidx[7].start - 1;
  if (keyidx[7].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->trid = h0.i;
  } else if (keyidx[7].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->trid = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NVS)
  ibuf = keyidx[8].start - 1;
  if (keyidx[8].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->nvs = h0.i;
  } else if (keyidx[8].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->nvs = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NHS)
  ibuf = keyidx[9].start - 1;
  if (keyidx[9].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->nhs = h0.i;
  } else if (keyidx[9].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->nhs = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DUSE)
  ibuf = keyidx[10].start - 1;
  if (keyidx[10].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->duse = h0.i;
  } else if (keyidx[10].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->duse = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_OFFSET)
  ibuf = keyidx[11].start - 1;
  if (keyidx[11].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->offset = i0.i;
  } else if (keyidx[11].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->offset = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GELEV)
  ibuf = keyidx[12].start - 1;
  if (keyidx[12].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gelev = i0.i;
  } else if (keyidx[12].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gelev = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SELEV)
  ibuf = keyidx[13].start - 1;
  if (keyidx[13].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->selev = i0.i;
  } else if (keyidx[13].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->selev = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEPTH)
  ibuf = keyidx[14].start - 1;
  if (keyidx[14].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sdepth = i0.i;
  } else if (keyidx[14].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sdepth = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GDEL)
  ibuf = keyidx[15].start - 1;
  if (keyidx[15].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gdel = i0.i;
  } else if (keyidx[15].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gdel = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEL)
  ibuf = keyidx[16].start - 1;
  if (keyidx[16].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sdel = i0.i;
  } else if (keyidx[16].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sdel = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWDEP)
  ibuf = keyidx[17].start - 1;
  if (keyidx[17].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->swdep = i0.i;
  } else if (keyidx[17].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->swdep = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GWDEP)
  ibuf = keyidx[18].start - 1;
  if (keyidx[18].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gwdep = i0.i;
  } else if (keyidx[18].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gwdep = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALEL)
  ibuf = keyidx[19].start - 1;
  if (keyidx[19].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->scalel = h0.i;
  } else if (keyidx[19].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->scalel = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALCO)
  ibuf = keyidx[20].start - 1;
  if (keyidx[20].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->scalco = h0.i;
  } else if (keyidx[20].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->scalco = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SX)
  ibuf = keyidx[21].start - 1;
  if (keyidx[21].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sx = i0.i;
  } else if (keyidx[21].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sx = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SY)
  ibuf = keyidx[22].start - 1;
  if (keyidx[22].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sy = i0.i;
  } else if (keyidx[22].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sy = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GX)
  ibuf = keyidx[23].start - 1;
  if (keyidx[23].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gx = i0.i;
  } else if (keyidx[23].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gx = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GY)
  ibuf = keyidx[24].start - 1;
  if (keyidx[24].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gy = i0.i;
  } else if (keyidx[24].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gy = h0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_COUNIT)
  ibuf = keyidx[25].start - 1;
  if (keyidx[25].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->counit = h0.i;
  } else if (keyidx[25].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->counit = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_WEVEL)
  ibuf = keyidx[26].start - 1;
  if (keyidx[26].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->wevel = h0.i;
  } else if (keyidx[26].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->wevel = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWEVEL)
  ibuf = keyidx[27].start - 1;
  if (keyidx[27].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->swevel = h0.i;
  } else if (keyidx[27].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->swevel = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SUT)
  ibuf = keyidx[28].start - 1;
  if (keyidx[28].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sut = h0.i;
  } else if (keyidx[28].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sut = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GUT)
  ibuf = keyidx[29].start - 1;
  if (keyidx[29].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gut = h0.i;
  } else if (keyidx[29].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gut = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SSTAT)
  ibuf = keyidx[30].start - 1;
  if (keyidx[30].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sstat = h0.i;
  } else if (keyidx[30].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sstat = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GSTAT)
  ibuf = keyidx[31].start - 1;
  if (keyidx[31].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gstat = h0.i;
  } else if (keyidx[31].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gstat = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TSTAT)
  ibuf = keyidx[32].start - 1;
  if (keyidx[32].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->tstat = h0.i;
  } else if (keyidx[32].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->tstat = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGA)
  ibuf = keyidx[33].start - 1;
  if (keyidx[33].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->laga = h0.i;
  } else if (keyidx[33].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->laga = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGB)
  ibuf = keyidx[34].start - 1;
  if (keyidx[34].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->lagb = h0.i;
  } else if (keyidx[34].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->lagb = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DELRT)
  ibuf = keyidx[35].start - 1;
  if (keyidx[35].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->delrt = h0.i;
  } else if (keyidx[35].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->delrt = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTS)
  ibuf = keyidx[36].start - 1;
  if (keyidx[36].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->muts = h0.i;
  } else if (keyidx[36].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->muts = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTE)
  ibuf = keyidx[37].start - 1;
  if (keyidx[37].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->mute = h0.i;
  } else if (keyidx[37].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->mute = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NS)
  ibuf = keyidx[38].start - 1;
  if (keyidx[38].bytes == 2) {
    u0.s[1] = buf[ibuf];
    u0.s[0] = buf[ibuf + 1];
    tr->ns = u0.i;
  } else if (keyidx[38].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->ns = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DT)
  ibuf = keyidx[39].start - 1;
  if (keyidx[39].bytes == 2) {
    u0.s[1] = buf[ibuf];
    u0.s[0] = buf[ibuf + 1];
    tr->dt = u0.i;
  } else if (keyidx[39].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->dt = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAIN)
  ibuf = keyidx[40].start - 1;
  if (keyidx[40].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gain = h0.i;
  } else if (keyidx[40].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gain = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGC)
  ibuf = keyidx[41].start - 1;
  if (keyidx[41].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->igc = h0.i;
  } else if (keyidx[41].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->igc = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGI)
  ibuf = keyidx[42].start - 1;
  if (keyidx[42].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->igi = h0.i;
  } else if (keyidx[42].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->igi = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CORR)
  ibuf = keyidx[43].start - 1;
  if (keyidx[43].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->corr = h0.i;
  } else if (keyidx[43].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->corr = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFS)
  ibuf = keyidx[44].start - 1;
  if (keyidx[44].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sfs = h0.i;
  } else if (keyidx[44].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sfs = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFE)
  ibuf = keyidx[45].start - 1;
  if (keyidx[45].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sfe = h0.i;
  } else if (keyidx[45].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sfe = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SLEN)
  ibuf = keyidx[46].start - 1;
  if (keyidx[46].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->slen = h0.i;
  } else if (keyidx[46].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->slen = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_STYP)
  ibuf = keyidx[47].start - 1;
  if (keyidx[47].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->styp = h0.i;
  } else if (keyidx[47].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->styp = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAS)
  ibuf = keyidx[48].start - 1;
  if (keyidx[48].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->stas = h0.i;
  } else if (keyidx[48].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->stas = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAE)
  ibuf = keyidx[49].start - 1;
  if (keyidx[49].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->stae = h0.i;
  } else if (keyidx[49].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->stae = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TATYP)
  ibuf = keyidx[50].start - 1;
  if (keyidx[50].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->tatyp = h0.i;
  } else if (keyidx[50].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->tatyp = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILF)
  ibuf = keyidx[51].start - 1;
  if (keyidx[51].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->afilf = h0.i;
  } else if (keyidx[51].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->afilf = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILS)
  ibuf = keyidx[52].start - 1;
  if (keyidx[52].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->afils = h0.i;
  } else if (keyidx[52].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->afils = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILF)
  ibuf = keyidx[53].start - 1;
  if (keyidx[53].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->nofilf = h0.i;
  } else if (keyidx[53].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->nofilf = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILS)
  ibuf = keyidx[54].start - 1;
  if (keyidx[54].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->nofils = h0.i;
  } else if (keyidx[54].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->nofils = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCF)
  ibuf = keyidx[55].start - 1;
  if (keyidx[55].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->lcf = h0.i;
  } else if (keyidx[55].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->lcf = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCF)
  ibuf = keyidx[56].start - 1;
  if (keyidx[56].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->hcf = h0.i;
  } else if (keyidx[56].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->hcf = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCS)
  ibuf = keyidx[57].start - 1;
  if (keyidx[57].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->lcs = h0.i;
  } else if (keyidx[57].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->lcs = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCS)
  ibuf = keyidx[58].start - 1;
  if (keyidx[58].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->hcs = h0.i;
  } else if (keyidx[58].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->hcs = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_YEAR)
  ibuf = keyidx[59].start - 1;
  if (keyidx[59].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->year = h0.i;
  } else if (keyidx[59].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->year = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DAY)
  ibuf = keyidx[60].start - 1;
  if (keyidx[60].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->day = h0.i;
  } else if (keyidx[60].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->day = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_HOUR)
  ibuf = keyidx[61].start - 1;
  if (keyidx[61].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->hour = h0.i;
  } else if (keyidx[61].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->hour = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_MINUTE)
  ibuf = keyidx[62].start - 1;
  if (keyidx[62].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->minute = h0.i;
  } else if (keyidx[62].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->minute = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SEC)
  ibuf = keyidx[63].start - 1;
  if (keyidx[63].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->sec = h0.i;
  } else if (keyidx[63].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->sec = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TIMBAS)
  ibuf = keyidx[64].start - 1;
  if (keyidx[64].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->timbas = h0.i;
  } else if (keyidx[64].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->timbas = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRWF)
  ibuf = keyidx[65].start - 1;
  if (keyidx[65].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->trwf = h0.i;
  } else if (keyidx[65].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->trwf = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNORS)
  ibuf = keyidx[66].start - 1;
  if (keyidx[66].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->grnors = h0.i;
  } else if (keyidx[66].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->grnors = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNOFR)
  ibuf = keyidx[67].start - 1;
  if (keyidx[67].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->grnofr = h0.i;
  } else if (keyidx[67].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->grnofr = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNLOF)
  ibuf = keyidx[68].start - 1;
  if (keyidx[68].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->grnlof = h0.i;
  } else if (keyidx[68].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->grnlof = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAPS)
  ibuf = keyidx[69].start - 1;
  if (keyidx[69].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->gaps = h0.i;
  } else if (keyidx[69].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->gaps = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_OTRAV)
  ibuf = keyidx[70].start - 1;
  if (keyidx[70].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->otrav = h0.i;
  } else if (keyidx[70].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->otrav = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPX)
  ibuf = keyidx[71].start - 1;
  if (keyidx[71].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->cdpx = h0.i;
  } else if (keyidx[71].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->cdpx = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPY)
  ibuf = keyidx[72].start - 1;
  if (keyidx[72].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->cdpy = h0.i;
  } else if (keyidx[72].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->cdpy = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_ILINE)
  ibuf = keyidx[73].start - 1;
  if (keyidx[73].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->iline = h0.i;
  } else if (keyidx[73].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->iline = i0.i;
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_XLINE)
  ibuf = keyidx[74].start - 1;
  if (keyidx[74].bytes == 2) {
    h0.s[1] = buf[ibuf];
    h0.s[0] = buf[ibuf + 1];
    tr->xline = h0.i;
  } else if (keyidx[74].bytes == 4) {
    i0.s[3] = buf[ibuf];
    i0.s[2] = buf[ibuf + 1];
    i0.s[1] = buf[ibuf + 2];
    i0.s[0] = buf[ibuf + 3];
    tr->xline = i0.i;
  }
#endif

  ns = tr->ns;
  /*
          ret=fread((tr->data),sizeof(float),ns,fp);
          if (ret!=ns) return 0;
  */
  memcpy((tr->data), buf + 240L, sizeof(float) * ns);

  p = (int *)(tr->data);
  /* 4-byte IBM float-point */
  if (format == 1) {
    for (itr = 0; itr < ns; ++itr) {
      i0.i = p[itr];
      i1.s[0] = i0.s[3];
      i1.s[1] = i0.s[2];
      i1.s[2] = i0.s[1];
      i1.s[3] = i0.s[0];

      fpn = i1.i;

      if ((fpn & MSK_MANT_IBM) == 0) {
        fpn = 0;
      } else {
        mantissa = (fpn & MSK_MANT_IBM);
        exponent = ((((fpn & MSK_EXP_IBM) >> 24) - 64) << 2);
        sign = (fpn & MSK_SIGN_IBM);
        while ((mantissa & MSK_NORM_IEEE) == 0) {
          mantissa <<= 1; /* normalize */
          exponent--;
        }
        mantissa = mantissa & 0x7fffff; /* shift understood one out */
        exponent += 126;
        if ((exponent < 0) || (exponent > 255)) {
          fprintf(stderr, "IBM floating point exponent out of range \n");
          exponent = 0;
        }
        exponent <<= 23;
        fpn = sign | exponent | mantissa;
      }
      p[itr] = fpn;
    }
  } else if (format == 5) /* 4-byte IEEE float-point */
  {
    for (itr = 0; itr < ns; ++itr) {
      i0.i = p[itr];
      i1.s[0] = i0.s[3];
      i1.s[1] = i0.s[2];
      i1.s[2] = i0.s[1];
      i1.s[3] = i0.s[0];
      p[itr] = i1.i;
    }
  } else {
    fprintf(stderr, "Unsupported SEGY format=%d\n", format);
    return 0;
  }

  return ns;
}

void inithd3600(hd3600 *hd, int ns, int dt_um, short int format) {
  const char *s =
      "\
C    Institute of Geology & Geophysics, Chinese Academy of Sciences             \
C         Theory and methods of seismic exploration                             \
C         Old Building(6th) Room 117                                            \
C         +86-010-82998160";

  memset(hd, 0, 3600);

  hd->jobid = 1;
  hd->lino = 1;
  hd->reno = 1;

  hd->hns = ns;
  hd->hdt = dt_um;

  if (format == 1 || format == 5) {
    hd->format = format;
  } else {
    fprintf(stderr, "Error: Unsupport File Format=%d in inithd3600\n", format);
    exit(0);
  }

  sethd3600text(hd, s);
}
void sethd3600text(hd3600 *hd, const char *text) {
  static unsigned char ascii2ebcdic[256] = {
#if 0
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,	/*         	*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,/*        	*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,	/*         	*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,	/*         	*/
	0x40,0x4F,0x7F,0x7B,0x5B,0x6C,0x50,0x7D,	/*  !"#$%&'	*/
	0x4D,0x5D,0x5C,0x4E,0x6B,0x60,0x4B,0x61,	/* ()*+,-./	*/
	0xF0,0xF1,0xF2,0xF3,0xF4,0xF5,0xF6,0xF7,	/* 01234567	*/
	0xF8,0xF9,0x7A,0x5E,0x4C,0x7E,0x6E,0x6F,	/* 89:;<=>?	*/
	0x7C,0xC1,0xC2,0xC3,0xC4,0xC5,0xC6,0xC7,	/* @ABCDEFG	*/
	0xC8,0xC9,0xD1,0xD2,0xD3,0xD4,0xD5,0xD6,	/* HIJKLMNO	*/
	0xD7,0xD8,0xD9,0xE2,0xE3,0xE4,0xE5,0xE6,	/* PQRSTUVW	*/
	0xE7,0xE8,0xE9,0x4A,0xE0,0x5A,0x5F,0x6D,	/* XYZ[\]^_	*/
	0x79,0x81,0x82,0x83,0x84,0x85,0x86,0x87,	/* `abcdefg	*/
	0x88,0x89,0x91,0x92,0x93,0x94,0x95,0x96,	/* hijklmno	*/
	0x97,0x98,0x99,0xA2,0xA3,0xA4,0xA5,0xA6,	/* pqrstuvw	*/
	0xA7,0xA8,0xA9,0xC0,0x6A,0xD0,0xA1,0x40,	/* xyz{|}~ 	*/
	0xB9,0xBA,0xED,0xBF,0xBC,0xBD,0xEC,0xFA,	/*              */
	0xCB,0xCC,0xCD,0xCE,0xCF,0xDA,0xDB,0xDC,	/*             	*/
	0xDE,0xDF,0xEA,0xEB,0xBE,0xCA,0xBB,0xFE,	/*		*/
	0xFB,0xFD,0x7d,0xEF,0xEE,0xFC,0xB8,0xDD,	/*		*/
	0x77,0x78,0xAF,0x8D,0x8A,0x8B,0xAE,0xB2,	/*		*/
	0x8F,0x90,0x9A,0x9B,0x9C,0x9D,0x9E,0x9F,	/*		*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,	/*	       	*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,	/*	       	*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,	/*	       	*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,	/*	       	*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,	/*	       	*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,	/*	       	*/
	0xAA,0xAB,0xAC,0xAD,0x8C,0x8E,0x80,0xB6,	/* ��������	*/
	0xB3,0xB5,0xB7,0xB1,0xB0,0xB4,0x76,0xA0,	/* ������?	*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40,	/*        	*/
	0x40,0x40,0x40,0x40,0x40,0x40,0x40,0x40}; 	/*        	*/
#endif
    0,
    1,
    2,
    3,
    55,
    45,
    46,
    47,
    22,
    5,
    37,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    60,
    61,
    50,
    38,
    24,
    25,
    63,
    39,
    28,
    29,
    30,
    31,
    64,
    79,
    127,
    123,
    91,
    108,
    80,
    125,
    77,
    93,
    92,
    78,
    107,
    96,
    75,
    97,
    240,
    241,
    242,
    243,
    244,
    245,
    246,
    247,
    248,
    249,
    122,
    94,
    76,
    126,
    110,
    111,
    124,
    193,
    194,
    195,
    196,
    197,
    198,
    199,
    200,
    201,
    209,
    210,
    211,
    212,
    213,
    214,
    215,
    216,
    217,
    226,
    227,
    228,
    229,
    230,
    231,
    232,
    233,
    74,
    224,
    90,
    95,
    109,
    121,
    129,
    130,
    131,
    132,
    133,
    134,
    135,
    136,
    137,
    145,
    146,
    147,
    148,
    149,
    150,
    151,
    152,
    153,
    162,
    163,
    164,
    165,
    166,
    167,
    168,
    169,
    192,
    106,
    208,
    161,
    7,
    32,
    33,
    34,
    35,
    36,
    21,
    6,
    23,
    40,
    41,
    42,
    43,
    44,
    9,
    10,
    27,
    48,
    49,
    26,
    51,
    52,
    53,
    54,
    8,
    56,
    57,
    58,
    59,
    4,
    20,
    62,
    225,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    81,
    82,
    83,
    84,
    85,
    86,
    87,
    88,
    89,
    98,
    99,
    100,
    101,
    102,
    103,
    104,
    105,
    112,
    113,
    114,
    115,
    116,
    117,
    118,
    119,
    120,
    128,
    138,
    139,
    140,
    141,
    142,
    143,
    144,
    154,
    155,
    156,
    157,
    158,
    159,
    160,
    170,
    171,
    172,
    173,
    174,
    175,
    176,
    177,
    178,
    179,
    180,
    181,
    182,
    183,
    184,
    185,
    186,
    187,
    188,
    189,
    190,
    191,
    202,
    203,
    204,
    205,
    206,
    207,
    218,
    219,
    220,
    221,
    222,
    223,
    234,
    235,
    236,
    237,
    238,
    239,
    250,
    251,
    252,
    253,
    254,
    255
  };

  int i, len;
  len = strlen(text);

  memset(hd->text, 0X40, 3200);

  len = len < 3600 ? len : 3600;
  for (i = 0; i < len; ++i) {
    hd->text[i] = ascii2ebcdic[(int)(text[i])];
  }
  for (i = 0; i < 3200; i += 80) {
    hd->text[i] = 0XC3;
  }
}

void fsputhd(FILE *fp, const hd3600 *hd) {
  static hd3600 hd_s;
  int_union i0, i1;
  short_union h0, h1;

  memcpy(&hd_s, hd, 3600);

#if defined(BHDR_ALL) || defined(BHDR_JOBID)
  i0.i = hd->jobid;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  hd_s.jobid = i1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_LINO)
  i0.i = hd->lino;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  hd_s.lino = i1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_RENO)
  i0.i = hd->reno;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  hd_s.reno = i1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_NTRPR)
  h0.i = hd->ntrpr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.ntrpr = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_NART)
  h0.i = hd->nart;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.nart = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HDT)
  h0.i = hd->hdt;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.hdt = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_DTO)
  h0.i = hd->dto;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.dto = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HNS)
  h0.i = hd->hns;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.hns = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_NSO)
  h0.i = hd->nso;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.nso = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_FORMAT)
  h0.i = hd->format;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.format = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_FOLD)
  h0.i = hd->fold;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.fold = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_TSORT)
  h0.i = hd->tsort;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.tsort = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_VSCODE)
  h0.i = hd->vscode;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.vscode = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSFS)
  h0.i = hd->hsfs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.hsfs = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSFE)
  h0.i = hd->hsfe;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.hsfe = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSLEN)
  h0.i = hd->hslen;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.hslen = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSTYP)
  h0.i = hd->hstyp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.hstyp = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_SCHN)
  h0.i = hd->schn;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.schn = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSTAS)
  h0.i = hd->hstas;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.hstas = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HSTAE)
  h0.i = hd->hstae;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.hstae = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HTATYP)
  h0.i = hd->htatyp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.htatyp = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_HCORR)
  h0.i = hd->hcorr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.hcorr = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_BGRCV)
  h0.i = hd->bgrcv;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.bgrcv = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_RCVM)
  h0.i = hd->rcvm;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.rcvm = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_MFEET)
  h0.i = hd->mfeet;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.mfeet = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_POLYT)
  h0.i = hd->polyt;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.polyt = h1.i;
#endif

#if defined(BHDR_ALL) || defined(BHDR_VPOL)
  h0.i = hd->vpol;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  hd_s.vpol = h1.i;
#endif

  fwrite(&hd_s, 3600, 1, fp);
}

void fsputtr(FILE *fp, const segy *tr, short int format,
             const keyindex *keyidx) {
  static int_union i0, i1;
  static short_union h0;
  static ushort_union u0;

  static segy tr_s;
  register int itr, ibuf, exponent, mantissa, fpn;
  int ns;
  int *p, *ps;
  unsigned char *buf = (unsigned char *)(&tr_s);

  if (!keyindex_inited) {
    fprintf(stderr,
            "Error: initkeyindex must be called before calling fsputtr\n");
    exit(0);
  }

  ns = tr->ns;

#if defined(SEGY_ALL) || defined(SEGY_TRACL)
  ibuf = keyidx[0].start - 1;
  if (keyidx[0].bytes == 4) {
    i0.i = tr->tracl;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[0].bytes == 2) {
    h0.i = tr->tracl;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACR)
  ibuf = keyidx[1].start - 1;
  if (keyidx[1].bytes == 4) {
    i0.i = tr->tracr;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[1].bytes == 2) {
    h0.i = tr->tracr;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_FLDR)
  ibuf = keyidx[2].start - 1;
  if (keyidx[2].bytes == 4) {
    i0.i = tr->fldr;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[2].bytes == 2) {
    h0.i = tr->fldr;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACF)
  ibuf = keyidx[3].start - 1;
  if (keyidx[3].bytes == 4) {
    i0.i = tr->tracf;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[3].bytes == 2) {
    h0.i = tr->tracf;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_EP)
  ibuf = keyidx[4].start - 1;
  if (keyidx[4].bytes == 4) {
    i0.i = tr->ep;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[4].bytes == 2) {
    h0.i = tr->ep;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDP)
  ibuf = keyidx[5].start - 1;
  if (keyidx[5].bytes == 4) {
    i0.i = tr->cdp;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[5].bytes == 2) {
    h0.i = tr->cdp;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPT)
  ibuf = keyidx[6].start - 1;
  if (keyidx[6].bytes == 4) {
    i0.i = tr->cdpt;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[6].bytes == 2) {
    h0.i = tr->cdpt;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRID)
  ibuf = keyidx[7].start - 1;
  if (keyidx[7].bytes == 2) {
    h0.i = tr->trid;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[7].bytes == 4) {
    i0.i = tr->trid;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NVS)
  ibuf = keyidx[8].start - 1;
  if (keyidx[8].bytes == 2) {
    h0.i = tr->nvs;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[8].bytes == 4) {
    i0.i = tr->nvs;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NHS)
  ibuf = keyidx[9].start - 1;
  if (keyidx[9].bytes == 2) {
    h0.i = tr->nhs;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[9].bytes == 4) {
    i0.i = tr->nhs;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DUSE)
  ibuf = keyidx[10].start - 1;
  if (keyidx[10].bytes == 2) {
    h0.i = tr->duse;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[10].bytes == 4) {
    i0.i = tr->duse;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_OFFSET)
  ibuf = keyidx[11].start - 1;
  if (keyidx[11].bytes == 4) {
    i0.i = tr->offset;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[11].bytes == 2) {
    h0.i = tr->offset;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GELEV)
  ibuf = keyidx[12].start - 1;
  if (keyidx[12].bytes == 4) {
    i0.i = tr->gelev;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[12].bytes == 2) {
    h0.i = tr->gelev;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SELEV)
  ibuf = keyidx[13].start - 1;
  if (keyidx[13].bytes == 4) {
    i0.i = tr->selev;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[13].bytes == 2) {
    h0.i = tr->selev;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEPTH)
  ibuf = keyidx[14].start - 1;
  if (keyidx[14].bytes == 4) {
    i0.i = tr->sdepth;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[14].bytes == 2) {
    h0.i = tr->sdepth;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GDEL)
  ibuf = keyidx[15].start - 1;
  if (keyidx[15].bytes == 4) {
    i0.i = tr->gdel;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[15].bytes == 2) {
    h0.i = tr->gdel;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEL)
  ibuf = keyidx[16].start - 1;
  if (keyidx[16].bytes == 4) {
    i0.i = tr->sdel;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[16].bytes == 2) {
    h0.i = tr->sdel;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWDEP)
  ibuf = keyidx[17].start - 1;
  if (keyidx[17].bytes == 4) {
    i0.i = tr->swdep;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[17].bytes == 2) {
    h0.i = tr->swdep;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GWDEP)
  ibuf = keyidx[18].start - 1;
  if (keyidx[18].bytes == 4) {
    i0.i = tr->gwdep;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[18].bytes == 2) {
    h0.i = tr->gwdep;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALEL)
  ibuf = keyidx[19].start - 1;
  if (keyidx[19].bytes == 2) {
    h0.i = tr->scalel;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[19].bytes == 4) {
    i0.i = tr->scalel;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALCO)
  ibuf = keyidx[20].start - 1;
  if (keyidx[20].bytes == 2) {
    h0.i = tr->scalco;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[20].bytes == 4) {
    i0.i = tr->scalco;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SX)
  ibuf = keyidx[21].start - 1;
  if (keyidx[21].bytes == 4) {
    i0.i = tr->sx;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[21].bytes == 2) {
    h0.i = tr->sx;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SY)
  ibuf = keyidx[22].start - 1;
  if (keyidx[22].bytes == 4) {
    i0.i = tr->sy;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[22].bytes == 2) {
    h0.i = tr->sy;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GX)
  ibuf = keyidx[23].start - 1;
  if (keyidx[23].bytes == 4) {
    i0.i = tr->gx;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[23].bytes == 2) {
    h0.i = tr->gx;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GY)
  ibuf = keyidx[24].start - 1;
  if (keyidx[24].bytes == 4) {
    i0.i = tr->gy;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  } else if (keyidx[24].bytes == 2) {
    h0.i = tr->gy;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_COUNIT)
  ibuf = keyidx[25].start - 1;
  if (keyidx[25].bytes == 2) {
    h0.i = tr->counit;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[25].bytes == 4) {
    i0.i = tr->counit;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_WEVEL)
  ibuf = keyidx[26].start - 1;
  if (keyidx[26].bytes == 2) {
    h0.i = tr->wevel;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[26].bytes == 4) {
    i0.i = tr->wevel;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWEVEL)
  ibuf = keyidx[27].start - 1;
  if (keyidx[27].bytes == 2) {
    h0.i = tr->swevel;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[27].bytes == 4) {
    i0.i = tr->swevel;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SUT)
  ibuf = keyidx[28].start - 1;
  if (keyidx[28].bytes == 2) {
    h0.i = tr->sut;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[28].bytes == 4) {
    i0.i = tr->sut;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GUT)
  ibuf = keyidx[29].start - 1;
  if (keyidx[29].bytes == 2) {
    h0.i = tr->gut;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[29].bytes == 4) {
    i0.i = tr->gut;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SSTAT)
  ibuf = keyidx[30].start - 1;
  if (keyidx[30].bytes == 2) {
    h0.i = tr->sstat;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[30].bytes == 4) {
    i0.i = tr->sstat;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GSTAT)
  ibuf = keyidx[31].start - 1;
  if (keyidx[31].bytes == 2) {
    h0.i = tr->gstat;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[31].bytes == 4) {
    i0.i = tr->gstat;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TSTAT)
  ibuf = keyidx[32].start - 1;
  if (keyidx[32].bytes == 2) {
    h0.i = tr->tstat;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[32].bytes == 4) {
    i0.i = tr->tstat;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGA)
  ibuf = keyidx[33].start - 1;
  if (keyidx[33].bytes == 2) {
    h0.i = tr->laga;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[33].bytes == 4) {
    i0.i = tr->laga;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGB)
  ibuf = keyidx[34].start - 1;
  if (keyidx[34].bytes == 2) {
    h0.i = tr->lagb;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[34].bytes == 4) {
    i0.i = tr->lagb;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DELRT)
  ibuf = keyidx[35].start - 1;
  if (keyidx[35].bytes == 2) {
    h0.i = tr->delrt;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[35].bytes == 4) {
    i0.i = tr->delrt;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTS)
  ibuf = keyidx[36].start - 1;
  if (keyidx[36].bytes == 2) {
    h0.i = tr->muts;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[36].bytes == 4) {
    i0.i = tr->muts;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTE)
  ibuf = keyidx[37].start - 1;
  if (keyidx[37].bytes == 2) {
    h0.i = tr->mute;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[37].bytes == 4) {
    i0.i = tr->mute;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NS)
  ibuf = keyidx[38].start - 1;
  if (keyidx[38].bytes == 2) {
    u0.i = tr->ns;
    buf[ibuf] = u0.s[1];
    buf[ibuf + 1] = u0.s[0];
  } else if (keyidx[38].bytes == 4) {
    i0.i = tr->ns;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DT)
  ibuf = keyidx[39].start - 1;
  if (keyidx[39].bytes == 2) {
    u0.i = tr->dt;
    buf[ibuf] = u0.s[1];
    buf[ibuf + 1] = u0.s[0];
  } else if (keyidx[39].bytes == 4) {
    i0.i = tr->dt;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAIN)
  ibuf = keyidx[40].start - 1;
  if (keyidx[40].bytes == 2) {
    h0.i = tr->gain;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[40].bytes == 4) {
    i0.i = tr->gain;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGC)
  ibuf = keyidx[41].start - 1;
  if (keyidx[41].bytes == 2) {
    h0.i = tr->igc;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[41].bytes == 4) {
    i0.i = tr->igc;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGI)
  ibuf = keyidx[42].start - 1;
  if (keyidx[42].bytes == 2) {
    h0.i = tr->igi;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[42].bytes == 4) {
    i0.i = tr->igi;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CORR)
  ibuf = keyidx[43].start - 1;
  if (keyidx[43].bytes == 2) {
    h0.i = tr->corr;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[43].bytes == 4) {
    i0.i = tr->corr;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFS)
  ibuf = keyidx[44].start - 1;
  if (keyidx[44].bytes == 2) {
    h0.i = tr->sfs;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[44].bytes == 4) {
    i0.i = tr->sfs;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFE)
  ibuf = keyidx[45].start - 1;
  if (keyidx[45].bytes == 2) {
    h0.i = tr->sfe;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[45].bytes == 4) {
    i0.i = tr->sfe;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SLEN)
  ibuf = keyidx[46].start - 1;
  if (keyidx[46].bytes == 2) {
    h0.i = tr->slen;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[46].bytes == 4) {
    i0.i = tr->slen;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_STYP)
  ibuf = keyidx[47].start - 1;
  if (keyidx[47].bytes == 2) {
    h0.i = tr->styp;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[47].bytes == 4) {
    i0.i = tr->styp;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAS)
  ibuf = keyidx[48].start - 1;
  if (keyidx[48].bytes == 2) {
    h0.i = tr->stas;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[48].bytes == 4) {
    i0.i = tr->stas;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAE)
  ibuf = keyidx[49].start - 1;
  if (keyidx[49].bytes == 2) {
    h0.i = tr->stae;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[49].bytes == 4) {
    i0.i = tr->stae;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TATYP)
  ibuf = keyidx[50].start - 1;
  if (keyidx[50].bytes == 2) {
    h0.i = tr->tatyp;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[50].bytes == 4) {
    i0.i = tr->tatyp;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILF)
  ibuf = keyidx[51].start - 1;
  if (keyidx[51].bytes == 2) {
    h0.i = tr->afilf;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[51].bytes == 4) {
    i0.i = tr->afilf;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILS)
  ibuf = keyidx[52].start - 1;
  if (keyidx[52].bytes == 2) {
    h0.i = tr->afils;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[52].bytes == 4) {
    i0.i = tr->afils;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILF)
  ibuf = keyidx[53].start - 1;
  if (keyidx[53].bytes == 2) {
    h0.i = tr->nofilf;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[53].bytes == 4) {
    i0.i = tr->nofilf;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILS)
  ibuf = keyidx[54].start - 1;
  if (keyidx[54].bytes == 2) {
    h0.i = tr->nofils;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[54].bytes == 4) {
    i0.i = tr->nofils;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCF)
  ibuf = keyidx[55].start - 1;
  if (keyidx[55].bytes == 2) {
    h0.i = tr->lcf;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[55].bytes == 4) {
    i0.i = tr->lcf;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCF)
  ibuf = keyidx[56].start - 1;
  if (keyidx[56].bytes == 2) {
    h0.i = tr->hcf;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[56].bytes == 4) {
    i0.i = tr->hcf;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCS)
  ibuf = keyidx[57].start - 1;
  if (keyidx[57].bytes == 2) {
    h0.i = tr->lcs;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[57].bytes == 4) {
    i0.i = tr->lcs;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCS)
  ibuf = keyidx[58].start - 1;
  if (keyidx[58].bytes == 2) {
    h0.i = tr->hcs;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[58].bytes == 4) {
    i0.i = tr->hcs;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_YEAR)
  ibuf = keyidx[59].start - 1;
  if (keyidx[59].bytes == 2) {
    h0.i = tr->year;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[59].bytes == 4) {
    i0.i = tr->year;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_DAY)
  ibuf = keyidx[60].start - 1;
  if (keyidx[60].bytes == 2) {
    h0.i = tr->day;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[60].bytes == 4) {
    i0.i = tr->day;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_HOUR)
  ibuf = keyidx[61].start - 1;
  if (keyidx[61].bytes == 2) {
    h0.i = tr->hour;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[61].bytes == 4) {
    i0.i = tr->hour;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_MINUTE)
  ibuf = keyidx[62].start - 1;
  if (keyidx[62].bytes == 2) {
    h0.i = tr->minute;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[62].bytes == 4) {
    i0.i = tr->minute;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_SEC)
  ibuf = keyidx[63].start - 1;
  if (keyidx[63].bytes == 2) {
    h0.i = tr->sec;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[63].bytes == 4) {
    i0.i = tr->sec;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TIMBAS)
  ibuf = keyidx[64].start - 1;
  if (keyidx[64].bytes == 2) {
    h0.i = tr->timbas;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[64].bytes == 4) {
    i0.i = tr->timbas;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRWF)
  ibuf = keyidx[65].start - 1;
  if (keyidx[65].bytes == 2) {
    h0.i = tr->trwf;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[65].bytes == 4) {
    i0.i = tr->trwf;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNORS)
  ibuf = keyidx[66].start - 1;
  if (keyidx[66].bytes == 2) {
    h0.i = tr->grnors;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[66].bytes == 4) {
    i0.i = tr->grnors;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNOFR)
  ibuf = keyidx[67].start - 1;
  if (keyidx[67].bytes == 2) {
    h0.i = tr->grnofr;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[67].bytes == 4) {
    i0.i = tr->grnofr;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNLOF)
  ibuf = keyidx[68].start - 1;
  if (keyidx[68].bytes == 2) {
    h0.i = tr->grnlof;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[68].bytes == 4) {
    i0.i = tr->grnlof;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAPS)
  ibuf = keyidx[69].start - 1;
  if (keyidx[69].bytes == 2) {
    h0.i = tr->gaps;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[69].bytes == 4) {
    i0.i = tr->gaps;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_OTRAV)
  ibuf = keyidx[70].start - 1;
  if (keyidx[70].bytes == 2) {
    h0.i = tr->otrav;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[70].bytes == 4) {
    i0.i = tr->otrav;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPX)
  ibuf = keyidx[71].start - 1;
  if (keyidx[71].bytes == 2) {
    h0.i = tr->cdpx;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[71].bytes == 4) {
    i0.i = tr->cdpx;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPY)
  ibuf = keyidx[72].start - 1;
  if (keyidx[72].bytes == 2) {
    h0.i = tr->cdpy;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[72].bytes == 4) {
    i0.i = tr->cdpy;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_ILINE)
  ibuf = keyidx[73].start - 1;
  if (keyidx[73].bytes == 2) {
    h0.i = tr->iline;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[73].bytes == 4) {
    i0.i = tr->iline;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

#if defined(SEGY_ALL) || defined(SEGY_XLINE)
  ibuf = keyidx[74].start - 1;
  if (keyidx[74].bytes == 2) {
    h0.i = tr->xline;
    buf[ibuf] = h0.s[1];
    buf[ibuf + 1] = h0.s[0];
  } else if (keyidx[74].bytes == 4) {
    i0.i = tr->xline;
    buf[ibuf] = i0.s[3];
    buf[ibuf + 1] = i0.s[2];
    buf[ibuf + 2] = i0.s[1];
    buf[ibuf + 3] = i0.s[0];
  }
#endif

  p = (int *)(tr->data);
  ps = (int *)(tr_s.data);

  /* 4-byte IBM float-point */
  if (format == 1) {
    for (itr = 0; itr < ns; ++itr) {
      fpn = p[itr];
      if (fpn) {
        mantissa = (0x007fffff & fpn) | 0x00800000;
        exponent = (int)((0x7f800000 & fpn) >> 23) - 126;
        while (exponent & 0x3) {
          ++exponent;
          mantissa >>= 1;
        }
        fpn = (0x80000000 & fpn) | (((exponent >> 2) + 64) << 24) | mantissa;
      }
      i0.i = fpn;
      i1.s[0] = i0.s[3];
      i1.s[1] = i0.s[2];
      i1.s[2] = i0.s[1];
      i1.s[3] = i0.s[0];
      ps[itr] = i1.i;
    }

  } else if (format == 5) /* 4-byte IEEE float-point */
  {
    for (itr = 0; itr < ns; ++itr) {
      i0.i = p[itr];
      i1.s[0] = i0.s[3];
      i1.s[1] = i0.s[2];
      i1.s[2] = i0.s[1];
      i1.s[3] = i0.s[0];
      ps[itr] = i1.i;
    }
  } else {
    fprintf(stderr, "Unsupported SEGY format=%d\n", format);
    return;
  }

  fwrite(&tr_s, 240 + sizeof(float) * ns, 1, fp);
}

/* SU format file read & write */
int fgettr(FILE *fp, segy *tr) {
  int_union i0, i1;
  short_union h0, h1;
  ushort_union u0, u1;
  float_union f0, f1;

  register int itr;
  int ns;
  int ret;

  ret = fread(tr, 240, 1, fp);
  if (ret != 1) return 0;

#if defined(SEGY_ALL) || defined(SEGY_TRACL)
  i0.i = tr->tracl;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->tracl = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACR)
  i0.i = tr->tracr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->tracr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_FLDR)
  i0.i = tr->fldr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->fldr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACF)
  i0.i = tr->tracf;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->tracf = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_EP)
  i0.i = tr->ep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->ep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDP)
  i0.i = tr->cdp;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->cdp = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPT)
  i0.i = tr->cdpt;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->cdpt = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRID)
  h0.i = tr->trid;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->trid = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NVS)
  h0.i = tr->nvs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->nvs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NHS)
  h0.i = tr->nhs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->nhs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DUSE)
  h0.i = tr->duse;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->duse = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_OFFSET)
  i0.i = tr->offset;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->offset = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GELEV)
  i0.i = tr->gelev;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->gelev = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SELEV)
  i0.i = tr->selev;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->selev = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEPTH)
  i0.i = tr->sdepth;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->sdepth = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GDEL)
  i0.i = tr->gdel;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->gdel = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEL)
  i0.i = tr->sdel;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->sdel = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWDEP)
  i0.i = tr->swdep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->swdep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GWDEP)
  i0.i = tr->gwdep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->gwdep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALEL)
  h0.i = tr->scalel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->scalel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALCO)
  h0.i = tr->scalco;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->scalco = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SX)
  i0.i = tr->sx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->sx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SY)
  i0.i = tr->sy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->sy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GX)
  i0.i = tr->gx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->gx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GY)
  i0.i = tr->gy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->gy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_COUNIT)
  h0.i = tr->counit;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->counit = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_WEVEL)
  h0.i = tr->wevel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->wevel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWEVEL)
  h0.i = tr->swevel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->swevel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SUT)
  h0.i = tr->sut;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->sut = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GUT)
  h0.i = tr->gut;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->gut = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SSTAT)
  h0.i = tr->sstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->sstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GSTAT)
  h0.i = tr->gstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->gstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TSTAT)
  h0.i = tr->tstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->tstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGA)
  h0.i = tr->laga;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->laga = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGB)
  h0.i = tr->lagb;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->lagb = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DELRT)
  h0.i = tr->delrt;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->delrt = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTS)
  h0.i = tr->muts;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->muts = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTE)
  h0.i = tr->mute;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->mute = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NS)
  u0.i = tr->ns;
  u1.s[0] = u0.s[1];
  u1.s[1] = u0.s[0];
  tr->ns = u1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DT)
  u0.i = tr->dt;
  u1.s[0] = u0.s[1];
  u1.s[1] = u0.s[0];
  tr->dt = u1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAIN)
  h0.i = tr->gain;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->gain = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGC)
  h0.i = tr->igc;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->igc = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGI)
  h0.i = tr->igi;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->igi = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CORR)
  h0.i = tr->corr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->corr = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFS)
  h0.i = tr->sfs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->sfs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFE)
  h0.i = tr->sfe;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->sfe = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SLEN)
  h0.i = tr->slen;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->slen = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STYP)
  h0.i = tr->styp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->styp = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAS)
  h0.i = tr->stas;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->stas = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAE)
  h0.i = tr->stae;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->stae = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TATYP)
  h0.i = tr->tatyp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->tatyp = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILF)
  h0.i = tr->afilf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->afilf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILS)
  h0.i = tr->afils;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->afils = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILF)
  h0.i = tr->nofilf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->nofilf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILS)
  h0.i = tr->nofils;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->nofils = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCF)
  h0.i = tr->lcf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->lcf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCF)
  h0.i = tr->hcf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->hcf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCS)
  h0.i = tr->lcs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->lcs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCS)
  h0.i = tr->hcs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->hcs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_YEAR)
  h0.i = tr->year;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->year = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DAY)
  h0.i = tr->day;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->day = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HOUR)
  h0.i = tr->hour;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->hour = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MINUTE)
  h0.i = tr->minute;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->minute = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SEC)
  h0.i = tr->sec;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->sec = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TIMBAS)
  h0.i = tr->timbas;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->timbas = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRWF)
  h0.i = tr->trwf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->trwf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNORS)
  h0.i = tr->grnors;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->grnors = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNOFR)
  h0.i = tr->grnofr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->grnofr = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNLOF)
  h0.i = tr->grnlof;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->grnlof = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAPS)
  h0.i = tr->gaps;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->gaps = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_OTRAV)
  h0.i = tr->otrav;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->otrav = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPX)
  i0.i = tr->cdpx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->cdpx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPY)
  i0.i = tr->cdpy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->cdpy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_ILINE)
  i0.i = tr->iline;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->iline = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_XLINE)
  i0.i = tr->xline;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->xline = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SP)
  i0.i = tr->sp;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->sp = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_UNSCALE)
  i0.i = tr->unscale;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->unscale = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NTR)
  i0.i = tr->ntr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->ntr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MARK)
  h0.i = tr->mark;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->mark = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SHORTPAD)
  h0.i = tr->shortpad;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->shortpad = h1.i;
#endif

  ns = tr->ns;
  ret = fread((tr->data), sizeof(float), ns, fp);
  if (ret != ns) return 0;

  for (itr = 0; itr < ns; ++itr) {
    memcpy(&f0, &((tr->data)[itr]), sizeof(float));
    f1.s[0] = f0.s[3];
    f1.s[1] = f0.s[2];
    f1.s[2] = f0.s[1];
    f1.s[3] = f0.s[0];
    tr->data[itr] = f1.f;
  }

  return ns;
}

void fputtr(FILE *fp, const segy *tr) {
  static int_union i0, i1;
  static short_union h0, h1;
  static ushort_union u0, u1;
  static float_union f0, f1;

  static segy tr_s;

  register int itr;
  static int ns;

  ns = tr->ns;

#if defined(SEGY_ALL) || defined(SEGY_TRACL)
  i0.i = tr->tracl;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.tracl = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACR)
  i0.i = tr->tracr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.tracr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_FLDR)
  i0.i = tr->fldr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.fldr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACF)
  i0.i = tr->tracf;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.tracf = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_EP)
  i0.i = tr->ep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.ep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDP)
  i0.i = tr->cdp;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.cdp = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPT)
  i0.i = tr->cdpt;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.cdpt = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRID)
  h0.i = tr->trid;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.trid = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NVS)
  h0.i = tr->nvs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.nvs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NHS)
  h0.i = tr->nhs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.nhs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DUSE)
  h0.i = tr->duse;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.duse = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_OFFSET)
  i0.i = tr->offset;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.offset = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GELEV)
  i0.i = tr->gelev;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.gelev = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SELEV)
  i0.i = tr->selev;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.selev = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEPTH)
  i0.i = tr->sdepth;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.sdepth = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GDEL)
  i0.i = tr->gdel;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.gdel = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEL)
  i0.i = tr->sdel;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.sdel = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWDEP)
  i0.i = tr->swdep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.swdep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GWDEP)
  i0.i = tr->gwdep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.gwdep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALEL)
  h0.i = tr->scalel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.scalel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALCO)
  h0.i = tr->scalco;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.scalco = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SX)
  i0.i = tr->sx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.sx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SY)
  i0.i = tr->sy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.sy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GX)
  i0.i = tr->gx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.gx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GY)
  i0.i = tr->gy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.gy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_COUNIT)
  h0.i = tr->counit;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.counit = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_WEVEL)
  h0.i = tr->wevel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.wevel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWEVEL)
  h0.i = tr->swevel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.swevel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SUT)
  h0.i = tr->sut;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.sut = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GUT)
  h0.i = tr->gut;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.gut = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SSTAT)
  h0.i = tr->sstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.sstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GSTAT)
  h0.i = tr->gstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.gstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TSTAT)
  h0.i = tr->tstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.tstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGA)
  h0.i = tr->laga;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.laga = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGB)
  h0.i = tr->lagb;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.lagb = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DELRT)
  h0.i = tr->delrt;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.delrt = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTS)
  h0.i = tr->muts;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.muts = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTE)
  h0.i = tr->mute;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.mute = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NS)
  u0.i = tr->ns;
  u1.s[0] = u0.s[1];
  u1.s[1] = u0.s[0];
  tr_s.ns = u1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DT)
  u0.i = tr->dt;
  u1.s[0] = u0.s[1];
  u1.s[1] = u0.s[0];
  tr_s.dt = u1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAIN)
  h0.i = tr->gain;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.gain = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGC)
  h0.i = tr->igc;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.igc = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGI)
  h0.i = tr->igi;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.igi = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CORR)
  h0.i = tr->corr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.corr = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFS)
  h0.i = tr->sfs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.sfs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFE)
  h0.i = tr->sfe;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.sfe = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SLEN)
  h0.i = tr->slen;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.slen = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STYP)
  h0.i = tr->styp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.styp = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAS)
  h0.i = tr->stas;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.stas = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAE)
  h0.i = tr->stae;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.stae = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TATYP)
  h0.i = tr->tatyp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.tatyp = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILF)
  h0.i = tr->afilf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.afilf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILS)
  h0.i = tr->afils;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.afils = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILF)
  h0.i = tr->nofilf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.nofilf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILS)
  h0.i = tr->nofils;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.nofils = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCF)
  h0.i = tr->lcf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.lcf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCF)
  h0.i = tr->hcf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.hcf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCS)
  h0.i = tr->lcs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.lcs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCS)
  h0.i = tr->hcs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.hcs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_YEAR)
  h0.i = tr->year;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.year = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DAY)
  h0.i = tr->day;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.day = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HOUR)
  h0.i = tr->hour;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.hour = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MINUTE)
  h0.i = tr->minute;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.minute = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SEC)
  h0.i = tr->sec;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.sec = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TIMBAS)
  h0.i = tr->timbas;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.timbas = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRWF)
  h0.i = tr->trwf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.trwf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNORS)
  h0.i = tr->grnors;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.grnors = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNOFR)
  h0.i = tr->grnofr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.grnofr = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNLOF)
  h0.i = tr->grnlof;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.grnlof = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAPS)
  h0.i = tr->gaps;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.gaps = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_OTRAV)
  h0.i = tr->otrav;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.otrav = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPX)
  i0.i = tr->cdpx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.cdpx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPY)
  i0.i = tr->cdpy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.cdpy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_ILINE)
  i0.i = tr->iline;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.iline = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_XLINE)
  i0.i = tr->xline;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.xline = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SP)
  i0.i = tr->sp;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.sp = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_UNSCALE)
  i0.i = tr->unscale;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.unscale = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NTR)
  i0.i = tr->ntr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.ntr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MARK)
  h0.i = tr->mark;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.mark = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SHORTPAD)
  h0.i = tr->shortpad;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.shortpad = h1.i;
#endif

  for (itr = 0; itr < ns; ++itr) {
    memcpy(&f0, &(tr->data[itr]), sizeof(float));
    f1.s[0] = f0.s[3];
    f1.s[1] = f0.s[2];
    f1.s[2] = f0.s[1];
    f1.s[3] = f0.s[0];
    tr_s.data[itr] = f1.f;
  }

  fwrite(&tr_s, 240 + sizeof(float) * ns, 1, fp);
}

void sputtr(unsigned char *buf, const segy *tr) {
  static int_union i0, i1;
  static short_union h0, h1;
  static ushort_union u0, u1;
  static float_union f0, f1;

  static segy tr_s;

  register int itr;
  static int ns;

  ns = tr->ns;

#if defined(SEGY_ALL) || defined(SEGY_TRACL)
  i0.i = tr->tracl;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.tracl = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACR)
  i0.i = tr->tracr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.tracr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_FLDR)
  i0.i = tr->fldr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.fldr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACF)
  i0.i = tr->tracf;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.tracf = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_EP)
  i0.i = tr->ep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.ep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDP)
  i0.i = tr->cdp;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.cdp = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPT)
  i0.i = tr->cdpt;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.cdpt = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRID)
  h0.i = tr->trid;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.trid = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NVS)
  h0.i = tr->nvs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.nvs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NHS)
  h0.i = tr->nhs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.nhs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DUSE)
  h0.i = tr->duse;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.duse = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_OFFSET)
  i0.i = tr->offset;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.offset = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GELEV)
  i0.i = tr->gelev;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.gelev = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SELEV)
  i0.i = tr->selev;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.selev = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEPTH)
  i0.i = tr->sdepth;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.sdepth = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GDEL)
  i0.i = tr->gdel;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.gdel = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEL)
  i0.i = tr->sdel;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.sdel = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWDEP)
  i0.i = tr->swdep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.swdep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GWDEP)
  i0.i = tr->gwdep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.gwdep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALEL)
  h0.i = tr->scalel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.scalel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALCO)
  h0.i = tr->scalco;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.scalco = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SX)
  i0.i = tr->sx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.sx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SY)
  i0.i = tr->sy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.sy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GX)
  i0.i = tr->gx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.gx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GY)
  i0.i = tr->gy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.gy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_COUNIT)
  h0.i = tr->counit;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.counit = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_WEVEL)
  h0.i = tr->wevel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.wevel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWEVEL)
  h0.i = tr->swevel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.swevel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SUT)
  h0.i = tr->sut;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.sut = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GUT)
  h0.i = tr->gut;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.gut = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SSTAT)
  h0.i = tr->sstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.sstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GSTAT)
  h0.i = tr->gstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.gstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TSTAT)
  h0.i = tr->tstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.tstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGA)
  h0.i = tr->laga;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.laga = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGB)
  h0.i = tr->lagb;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.lagb = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DELRT)
  h0.i = tr->delrt;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.delrt = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTS)
  h0.i = tr->muts;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.muts = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTE)
  h0.i = tr->mute;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.mute = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NS)
  u0.i = tr->ns;
  u1.s[0] = u0.s[1];
  u1.s[1] = u0.s[0];
  tr_s.ns = u1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DT)
  u0.i = tr->dt;
  u1.s[0] = u0.s[1];
  u1.s[1] = u0.s[0];
  tr_s.dt = u1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAIN)
  h0.i = tr->gain;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.gain = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGC)
  h0.i = tr->igc;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.igc = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGI)
  h0.i = tr->igi;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.igi = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CORR)
  h0.i = tr->corr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.corr = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFS)
  h0.i = tr->sfs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.sfs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFE)
  h0.i = tr->sfe;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.sfe = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SLEN)
  h0.i = tr->slen;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.slen = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STYP)
  h0.i = tr->styp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.styp = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAS)
  h0.i = tr->stas;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.stas = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAE)
  h0.i = tr->stae;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.stae = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TATYP)
  h0.i = tr->tatyp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.tatyp = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILF)
  h0.i = tr->afilf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.afilf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILS)
  h0.i = tr->afils;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.afils = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILF)
  h0.i = tr->nofilf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.nofilf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILS)
  h0.i = tr->nofils;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.nofils = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCF)
  h0.i = tr->lcf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.lcf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCF)
  h0.i = tr->hcf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.hcf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCS)
  h0.i = tr->lcs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.lcs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCS)
  h0.i = tr->hcs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.hcs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_YEAR)
  h0.i = tr->year;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.year = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DAY)
  h0.i = tr->day;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.day = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HOUR)
  h0.i = tr->hour;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.hour = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MINUTE)
  h0.i = tr->minute;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.minute = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SEC)
  h0.i = tr->sec;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.sec = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TIMBAS)
  h0.i = tr->timbas;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.timbas = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRWF)
  h0.i = tr->trwf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.trwf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNORS)
  h0.i = tr->grnors;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.grnors = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNOFR)
  h0.i = tr->grnofr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.grnofr = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNLOF)
  h0.i = tr->grnlof;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.grnlof = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAPS)
  h0.i = tr->gaps;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.gaps = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_OTRAV)
  h0.i = tr->otrav;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.otrav = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPX)
  i0.i = tr->cdpx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.cdpx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPY)
  i0.i = tr->cdpy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.cdpy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_ILINE)
  i0.i = tr->iline;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.iline = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_XLINE)
  i0.i = tr->xline;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.xline = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SP)
  i0.i = tr->sp;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.sp = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_UNSCALE)
  i0.i = tr->unscale;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.unscale = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NTR)
  i0.i = tr->ntr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr_s.ntr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MARK)
  h0.i = tr->mark;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.mark = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SHORTPAD)
  h0.i = tr->shortpad;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr_s.shortpad = h1.i;
#endif

  for (itr = 0; itr < ns; ++itr) {
    memcpy(&f0, &(tr->data[itr]), sizeof(float));
    f1.s[0] = f0.s[3];
    f1.s[1] = f0.s[2];
    f1.s[2] = f0.s[1];
    f1.s[3] = f0.s[0];
    tr_s.data[itr] = f1.f;
  }

  memcpy(buf, &tr_s, 240L + sizeof(float) * ns);
}

int sgettr(unsigned char *buf, segy *tr) {
  int_union i0, i1;
  short_union h0, h1;
  ushort_union u0, u1;

  register int itr;
  int ns;
  int *p, *ps;
  /*	int ret; */

  /*	ret=fread(tr,240,1,fp);
          if (ret!=1) return 0;
  */
  memcpy(tr, buf, 240);

#if defined(SEGY_ALL) || defined(SEGY_TRACL)
  i0.i = tr->tracl;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->tracl = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACR)
  i0.i = tr->tracr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->tracr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_FLDR)
  i0.i = tr->fldr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->fldr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACF)
  i0.i = tr->tracf;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->tracf = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_EP)
  i0.i = tr->ep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->ep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDP)
  i0.i = tr->cdp;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->cdp = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPT)
  i0.i = tr->cdpt;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->cdpt = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRID)
  h0.i = tr->trid;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->trid = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NVS)
  h0.i = tr->nvs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->nvs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NHS)
  h0.i = tr->nhs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->nhs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DUSE)
  h0.i = tr->duse;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->duse = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_OFFSET)
  i0.i = tr->offset;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->offset = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GELEV)
  i0.i = tr->gelev;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->gelev = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SELEV)
  i0.i = tr->selev;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->selev = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEPTH)
  i0.i = tr->sdepth;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->sdepth = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GDEL)
  i0.i = tr->gdel;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->gdel = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEL)
  i0.i = tr->sdel;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->sdel = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWDEP)
  i0.i = tr->swdep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->swdep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GWDEP)
  i0.i = tr->gwdep;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->gwdep = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALEL)
  h0.i = tr->scalel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->scalel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALCO)
  h0.i = tr->scalco;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->scalco = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SX)
  i0.i = tr->sx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->sx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SY)
  i0.i = tr->sy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->sy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GX)
  i0.i = tr->gx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->gx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GY)
  i0.i = tr->gy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->gy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_COUNIT)
  h0.i = tr->counit;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->counit = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_WEVEL)
  h0.i = tr->wevel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->wevel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWEVEL)
  h0.i = tr->swevel;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->swevel = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SUT)
  h0.i = tr->sut;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->sut = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GUT)
  h0.i = tr->gut;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->gut = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SSTAT)
  h0.i = tr->sstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->sstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GSTAT)
  h0.i = tr->gstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->gstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TSTAT)
  h0.i = tr->tstat;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->tstat = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGA)
  h0.i = tr->laga;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->laga = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGB)
  h0.i = tr->lagb;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->lagb = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DELRT)
  h0.i = tr->delrt;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->delrt = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTS)
  h0.i = tr->muts;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->muts = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTE)
  h0.i = tr->mute;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->mute = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NS)
  u0.i = tr->ns;
  u1.s[0] = u0.s[1];
  u1.s[1] = u0.s[0];
  tr->ns = u1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DT)
  u0.i = tr->dt;
  u1.s[0] = u0.s[1];
  u1.s[1] = u0.s[0];
  tr->dt = u1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAIN)
  h0.i = tr->gain;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->gain = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGC)
  h0.i = tr->igc;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->igc = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGI)
  h0.i = tr->igi;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->igi = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CORR)
  h0.i = tr->corr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->corr = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFS)
  h0.i = tr->sfs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->sfs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFE)
  h0.i = tr->sfe;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->sfe = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SLEN)
  h0.i = tr->slen;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->slen = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STYP)
  h0.i = tr->styp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->styp = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAS)
  h0.i = tr->stas;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->stas = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAE)
  h0.i = tr->stae;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->stae = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TATYP)
  h0.i = tr->tatyp;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->tatyp = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILF)
  h0.i = tr->afilf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->afilf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILS)
  h0.i = tr->afils;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->afils = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILF)
  h0.i = tr->nofilf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->nofilf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILS)
  h0.i = tr->nofils;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->nofils = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCF)
  h0.i = tr->lcf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->lcf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCF)
  h0.i = tr->hcf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->hcf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCS)
  h0.i = tr->lcs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->lcs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCS)
  h0.i = tr->hcs;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->hcs = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_YEAR)
  h0.i = tr->year;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->year = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_DAY)
  h0.i = tr->day;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->day = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_HOUR)
  h0.i = tr->hour;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->hour = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MINUTE)
  h0.i = tr->minute;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->minute = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SEC)
  h0.i = tr->sec;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->sec = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TIMBAS)
  h0.i = tr->timbas;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->timbas = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRWF)
  h0.i = tr->trwf;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->trwf = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNORS)
  h0.i = tr->grnors;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->grnors = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNOFR)
  h0.i = tr->grnofr;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->grnofr = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNLOF)
  h0.i = tr->grnlof;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->grnlof = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAPS)
  h0.i = tr->gaps;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->gaps = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_OTRAV)
  h0.i = tr->otrav;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->otrav = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPX)
  i0.i = tr->cdpx;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->cdpx = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPY)
  i0.i = tr->cdpy;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->cdpy = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_ILINE)
  i0.i = tr->iline;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->iline = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_XLINE)
  i0.i = tr->xline;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->xline = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SP)
  i0.i = tr->sp;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->sp = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_UNSCALE)
  i0.i = tr->unscale;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->unscale = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_NTR)
  i0.i = tr->ntr;
  i1.s[0] = i0.s[3];
  i1.s[1] = i0.s[2];
  i1.s[2] = i0.s[1];
  i1.s[3] = i0.s[0];
  tr->ntr = i1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_MARK)
  h0.i = tr->mark;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->mark = h1.i;
#endif

#if defined(SEGY_ALL) || defined(SEGY_SHORTPAD)
  h0.i = tr->shortpad;
  h1.s[0] = h0.s[1];
  h1.s[1] = h0.s[0];
  tr->shortpad = h1.i;
#endif

  ns = tr->ns;
  /*
          ret=fread((tr->data),sizeof(float),ns,fp);
          if (ret!=ns) return 0;
  */
  memcpy((tr->data), buf + 240, sizeof(float) * ns);

  p = (int *)(&buf[240]);
  ps = (int *)(tr->data);
  for (itr = 0; itr < ns; ++itr) {
    i0.i = p[itr];
    i1.s[0] = i0.s[3];
    i1.s[1] = i0.s[2];
    i1.s[2] = i0.s[1];
    i1.s[3] = i0.s[0];

    ps[itr] = i1.i;
  }

  return ns;
}

void swaptr_macro(unsigned char *buf, segy *tr) {
  float_union f0;
  register int itr;
  unsigned char *p;
  int ns;
  memcpy(tr, buf, 240);

#if defined(SEGY_ALL) || defined(SEGY_TRACL)
  tr->tracl = SWAP4(tr->tracl);
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACR)
  tr->tracr = SWAP4(tr->tracr);
#endif

#if defined(SEGY_ALL) || defined(SEGY_FLDR)
  tr->fldr = SWAP4(tr->fldr);
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRACF)
  tr->tracf = SWAP4(tr->tracf);
#endif

#if defined(SEGY_ALL) || defined(SEGY_EP)
  tr->ep = SWAP4(tr->ep);
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDP)
  tr->cdp = SWAP4(tr->cdp);
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPT)
  tr->cdpt = SWAP4(tr->cdpt);
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRID)
  tr->trid = SWAP2(tr->trid);
#endif

#if defined(SEGY_ALL) || defined(SEGY_NVS)
  tr->nvs = SWAP2(tr->nvs);
#endif

#if defined(SEGY_ALL) || defined(SEGY_NHS)
  tr->nhs = SWAP2(tr->nhs);
#endif

#if defined(SEGY_ALL) || defined(SEGY_DUSE)
  tr->duse = SWAP2(tr->duse);
#endif

#if defined(SEGY_ALL) || defined(SEGY_OFFSET)
  tr->offset = SWAP4(tr->offset);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GELEV)
  tr->gelev = SWAP4(tr->gelev);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SELEV)
  tr->selev = SWAP4(tr->selev);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEPTH)
  tr->sdepth = SWAP4(tr->sdepth);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GDEL)
  tr->gdel = SWAP4(tr->gdel);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SDEL)
  tr->sdel = SWAP4(tr->sdel);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWDEP)
  tr->swdep = SWAP4(tr->swdep);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GWDEP)
  tr->gwdep = SWAP4(tr->gwdep);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALEL)
  tr->scalel = SWAP2(tr->scalel);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SCALCO)
  tr->scalco = SWAP2(tr->scalco);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SX)
  tr->sx = SWAP4(tr->sx);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SY)
  tr->sy = SWAP4(tr->sy);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GX)
  tr->gx = SWAP4(tr->gx);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GY)
  tr->gy = SWAP4(tr->gy);
#endif

#if defined(SEGY_ALL) || defined(SEGY_COUNIT)
  tr->counit = SWAP2(tr->counit);
#endif

#if defined(SEGY_ALL) || defined(SEGY_WEVEL)
  tr->wevel = SWAP2(tr->wevel);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SWEVEL)
  tr->swevel = SWAP2(tr->swevel);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SUT)
  tr->sut = SWAP2(tr->sut);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GUT)
  tr->gut = SWAP2(tr->gut);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SSTAT)
  tr->sstat = SWAP2(tr->sstat);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GSTAT)
  tr->gstat = SWAP2(tr->gstat);
#endif

#if defined(SEGY_ALL) || defined(SEGY_TSTAT)
  tr->tstat = SWAP2(tr->tstat);
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGA)
  tr->laga = SWAP2(tr->laga);
#endif

#if defined(SEGY_ALL) || defined(SEGY_LAGB)
  tr->lagb = SWAP2(tr->lagb);
#endif

#if defined(SEGY_ALL) || defined(SEGY_DELRT)
  tr->delrt = SWAP2(tr->delrt);
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTS)
  tr->muts = SWAP2(tr->muts);
#endif

#if defined(SEGY_ALL) || defined(SEGY_MUTE)
  tr->mute = SWAP2(tr->mute);
#endif

#if defined(SEGY_ALL) || defined(SEGY_NS)
  tr->ns = SWAP2(tr->ns);
#endif

#if defined(SEGY_ALL) || defined(SEGY_DT)
  tr->dt = SWAP2(tr->dt);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAIN)
  tr->gain = SWAP2(tr->gain);
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGC)
  tr->igc = SWAP2(tr->igc);
#endif

#if defined(SEGY_ALL) || defined(SEGY_IGI)
  tr->igi = SWAP2(tr->igi);
#endif

#if defined(SEGY_ALL) || defined(SEGY_CORR)
  tr->corr = SWAP2(tr->corr);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFS)
  tr->sfs = SWAP2(tr->sfs);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SFE)
  tr->sfe = SWAP2(tr->sfe);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SLEN)
  tr->slen = SWAP2(tr->slen);
#endif

#if defined(SEGY_ALL) || defined(SEGY_STYP)
  tr->styp = SWAP2(tr->styp);
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAS)
  tr->stas = SWAP2(tr->stas);
#endif

#if defined(SEGY_ALL) || defined(SEGY_STAE)
  tr->stae = SWAP2(tr->stae);
#endif

#if defined(SEGY_ALL) || defined(SEGY_TATYP)
  tr->tatyp = SWAP2(tr->tatyp);
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILF)
  tr->afilf = SWAP2(tr->afilf);
#endif

#if defined(SEGY_ALL) || defined(SEGY_AFILS)
  tr->afils = SWAP2(tr->afils);
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILF)
  tr->nofilf = SWAP2(tr->nofilf);
#endif

#if defined(SEGY_ALL) || defined(SEGY_NOFILS)
  tr->nofils = SWAP2(tr->nofils);
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCF)
  tr->lcf = SWAP2(tr->lcf);
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCF)
  tr->hcf = SWAP2(tr->hcf);
#endif

#if defined(SEGY_ALL) || defined(SEGY_LCS)
  tr->lcs = SWAP2(tr->lcs);
#endif

#if defined(SEGY_ALL) || defined(SEGY_HCS)
  tr->hcs = SWAP2(tr->hcs);
#endif

#if defined(SEGY_ALL) || defined(SEGY_YEAR)
  tr->year = SWAP2(tr->year);
#endif

#if defined(SEGY_ALL) || defined(SEGY_DAY)
  tr->day = SWAP2(tr->day);
#endif

#if defined(SEGY_ALL) || defined(SEGY_HOUR)
  tr->hour = SWAP2(tr->hour);
#endif

#if defined(SEGY_ALL) || defined(SEGY_MINUTE)
  tr->minute = SWAP2(tr->minute);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SEC)
  tr->sec = SWAP2(tr->sec);
#endif

#if defined(SEGY_ALL) || defined(SEGY_TIMBAS)
  tr->timbas = SWAP2(tr->timbas);
#endif

#if defined(SEGY_ALL) || defined(SEGY_TRWF)
  tr->trwf = SWAP2(tr->trwf);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNORS)
  tr->grnors = SWAP2(tr->grnors);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNOFR)
  tr->grnofr = SWAP2(tr->grnofr);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GRNLOF)
  tr->grnlof = SWAP2(tr->grnlof);
#endif

#if defined(SEGY_ALL) || defined(SEGY_GAPS)
  tr->gaps = SWAP2(tr->gaps);
#endif

#if defined(SEGY_ALL) || defined(SEGY_OTRAV)
  tr->otrav = SWAP2(tr->otrav);
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPX)
  tr->cdpx = SWAP4(tr->cdpx);
#endif

#if defined(SEGY_ALL) || defined(SEGY_CDPY)
  tr->cdpy = SWAP4(tr->cdpy);
#endif

#if defined(SEGY_ALL) || defined(SEGY_ILINE)
  tr->iline = SWAP4(tr->iline);
#endif

#if defined(SEGY_ALL) || defined(SEGY_XLINE)
  tr->xline = SWAP4(tr->xline);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SP)
  tr->sp = SWAP4(tr->sp);
#endif

#if defined(SEGY_ALL) || defined(SEGY_UNSCALE)
  tr->unscale = SWAP4(tr->unscale);
#endif

#if defined(SEGY_ALL) || defined(SEGY_NTR)
  tr->ntr = SWAP4(tr->ntr);
#endif

#if defined(SEGY_ALL) || defined(SEGY_MARK)
  tr->mark = SWAP2(tr->mark);
#endif

#if defined(SEGY_ALL) || defined(SEGY_SHORTPAD)
  tr->shortpad = SWAP2(tr->shortpad);
#endif

  ns = tr->ns;
  p = buf + 240;
  for (itr = 0; itr < ns; ++itr) {
    f0.s[3] = *p++;
    f0.s[2] = *p++;
    f0.s[1] = *p++;
    f0.s[0] = *p++;
    tr->data[itr] = f0.f;
  }
}