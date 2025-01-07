#include "psdmpkg.h"

int  xargc;
char **xargv;
int  inited;

void initargs_inner(int argc, char *argv[], const char *sdoc[]) {
  int max_args = 128 * 1024;
  for (int i = 0; i < argc; ++i) {
    int len = strlen(argv[i]);
    max_args = max_args > len ? max_args : len;
  }

  xargc = argc;
  xargv = alloc2char(max_args, xargc + 1);

  for (int i = 0; i < xargc; ++i) {
    memcpy(xargv[i], argv[i], strlen(argv[i]));
    xargv[i][strlen(argv[i])] = 0;
  }
  inited = 1;
}

char *malloc1c(size_t n1) {
  char *p;
  p = (char *)malloc(sizeof(char) * n1);
  if (p == NULL) {
    fprintf(stderr, "malloc1c error when call malloc1c(%zu)\n", n1);
    exit(0);
  }
  return p;
}
int *malloc1i(size_t n1) {
  int *p;
  p = (int *)malloc(sizeof(int) * n1);
  if (p == NULL) {
    fprintf(stderr, "malloc1i error when call malloc1i(%zu) p==NULL\n", n1);
    exit(0);
  }
  return p;
}
int **malloc2i(size_t n1, size_t n2) {
  int **p;
  size_t i;

  p = (int **)malloc(sizeof(int *) * n2);
  if (p == NULL) {
    fprintf(stderr, "malloc2i error when call malloc2i(%zu,%zu) p==NULL\n", n1,
            n2);
    exit(0);
  }
  p[0] = (int *)malloc(sizeof(int) * n1 * n2);
  if (p[0] == NULL) {
    fprintf(stderr, "malloc2i error when call malloc2i(%zu,%zu) p[0]==NULL\n",
            n1, n2);
    exit(0);
  }
  for (i = 1; i < n2; ++i) {
    p[i] = &p[0][i * n1];
  }
  return p;
}
float *malloc1f(size_t n1) {
  float *p;
  p = (float *)malloc(sizeof(float) * n1);
  if (p == NULL) {
    fprintf(stderr, "malloc1f error when call malloc1f(%zu) p==NULL\n", n1);
    exit(0);
  }
  return p;
}
float **malloc2f(size_t n1, size_t n2) {
  float **p;
  size_t i;

  p = (float **)malloc(sizeof(float *) * n2);
  if (p == NULL) {
    fprintf(stderr, "malloc2f error when call malloc2f(%zu,%zi) p==NULL\n", n1,
            n2);
    exit(0);
  }
  p[0] = (float *)malloc(sizeof(float) * n1 * n2);
  if (p[0] == NULL) {
    fprintf(stderr, "malloc2f error when call malloc2f(%zu,%zu) p[0]==NULL\n",
            n1, n2);
    exit(0);
  }
  for (i = 1; i < n2; ++i) {
    p[i] = &p[0][i * n1];
  }
  return p;
}
char **malloc2c(size_t n1,size_t n2)
{
	char **p;
	size_t i;
	
	p=(char**)malloc(sizeof(char*)*n2);
	if (p==NULL) {fprintf(stderr,"malloc2c error when call malloc2c(%zu,%zu) p==NULL\n",n1,n2);exit(0);}
	p[0]=(char*)malloc(sizeof(char)*n1*n2);
	if (p[0]==NULL) {fprintf(stderr,"malloc2c error when call malloc2c(%zu,%zu) p[0]==NULL\n",n1,n2);exit(0);}
	for (i=1;i<n2;++i)
	{
		p[i]=&p[0][i*n1];
	}
	return p;
}
int *parseline(const char *str, int *n) {
  int nmark, imark;
  int *mloc = NULL, *mval = NULL;

  int len;
  int i, j, ii, npt;

  int *pt;
  int *p;

  len = strlen(str);

  /* Statistic the split symbols */
  nmark = 0;
  for (i = 0; i < len; ++i) {
    if (str[i] == ',' || str[i] == '-' || str[i] == ':') ++nmark;
  }

  /* Allocate Space For the symbols  */
  if (nmark > 0) {
    mloc = (int *)malloc(sizeof(int) * nmark);
    mval = (int *)malloc(sizeof(int) * nmark);
  }
  pt = (int *)malloc(sizeof(int) * (nmark + 1));

  /* Get the marker's position & value */
  imark = 0;
  for (i = 0; i < len; ++i) {
    if (str[i] == ',' || str[i] == '-' || str[i] == ':') {
      mloc[imark] = i;
      mval[imark] = str[i];
      ++imark;
    }
  }

  /* Get the item between the markers  */
  sscanf(str, "%d[^,-:]", &pt[0]);
  for (i = 0; i < nmark; ++i) {
    sscanf(&str[mloc[i] + 1], "%d[^,-:]", &pt[i + 1]);
  }

  npt = 1;
  for (i = 0; i < nmark; ++i) {
    if (mval[i] == ',') {
      ++npt;
    }
    if (mval[i] == '-') {
      npt += pt[i + 1] - pt[i];
    }
    if (mval[i] == ':') {
      npt += (pt[i + 2] - pt[i]) / pt[i + 1];
      ++i;
    }
  }

  *n = npt;
  p = (int *)malloc(sizeof(int) * npt);

  p[0] = pt[0];
  j = 1;
  for (i = 0; i < nmark; ++i) {
    if (mval[i] == ',') {
      p[j] = pt[i + 1];
      ++j;
    }
    if (mval[i] == '-') {
      for (ii = pt[i] + 1; ii <= pt[i + 1]; ++ii) {
        p[j] = ii;
        ++j;
      }
    }
    if (mval[i] == ':') {
      --j;
      if (pt[i + 1] > 0)
        for (ii = pt[i]; ii <= pt[i + 2]; ii += pt[i + 1]) {
          p[j] = ii;
          ++j;
        }

      if (pt[i + 1] < 0)
        for (ii = pt[i]; ii >= pt[i + 2]; ii += pt[i + 1]) {
          p[j] = ii;
          ++j;
        }
      ++i;
    }
  }

  /* Release the internal allocated memory  */
  if (nmark > 0) {
    free(mloc);
    free(mval);
  }
  free(pt);

  return p;
}

int fexist(const char *path) { return (!access(path, 06)); }

int fetch_and_taper(const float *data, int ns, int it1, int it2, float *data1) {
  int it, ito;
  int ns1;
  float c;

  /* clip-off one segment */
  for (it = it1; it <= it2; ++it) {
    ito = it - it1;
    if (it < 0 || it >= ns)
      data1[ito] = 0.0;
    else
      data1[ito] = data[it];
  }

  ns1 = it2 - it1 + 1;
  /* taper the two end */

  for (it = 0; it < 200 && it < ns1; ++it) {
    c = 0.5 - 0.5 * cosf(M_PI * it / 200.0);
    data1[it] *= c;
    data1[ns1 - 1 - it] *= c;
  }
  return 1;
}

int getparint(const char *key, int *val) {
  int i;

  if (!inited) {
    fprintf(stderr, "argument is not initialized \n");
    exit(0);
  }

  for (i = 0; i < xargc; ++i) {
    if (0 == strncmp(key, xargv[i], strlen(key)) &&
        xargv[i][strlen(key)] == '=') {
      break;
    }
  }
  if (i < xargc) {
    val[0] = strtol(&xargv[i][strlen(key) + 1], NULL, 10);
    return 1;
  }
  return 0;
}

int getparfloat(const char *key, float *val) {
  int i;

  if (!inited) {
    fprintf(stderr, "argument is not initialized \n");
    exit(0);
  }

  for (i = 0; i < xargc; ++i) {
    if (0 == strncmp(key, xargv[i], strlen(key)) &&
        xargv[i][strlen(key)] == '=') {
      break;
    }
  }
  if (i < xargc) {
    val[0] = strtod(&xargv[i][strlen(key) + 1], NULL);
    return 1;
  }
  return 0;
}

int getpardouble(const char *key, double *val) {
  int i;

  if (!inited) {
    fprintf(stderr, "argument is not initialized \n");
    exit(0);
  }

  for (i = 0; i < xargc; ++i) {
    if (0 == strncmp(key, xargv[i], strlen(key)) &&
        xargv[i][strlen(key)] == '=') {
      break;
    }
  }
  if (i < xargc) {
    val[0] = strtod(&xargv[i][strlen(key) + 1], NULL);
    return 1;
  }
  return 0;
}

int getparstr(const char *key, char **val) {
  int i;

  if (!inited) {
    fprintf(stderr, "argument is not initialized \n");
    exit(0);
  }

  for (i = 0; i < xargc; ++i) {
    if (0 == strncmp(key, xargv[i], strlen(key)) &&
        xargv[i][strlen(key)] == '=') {
      break;
    }
  }
  if (i < xargc) {
    val[0] = &xargv[i][strlen(key) + 1];
    return 1;
  }
  return 0;
}

int getparstr(const char *key, const char **val) {
  int i;

  if (!inited) {
    fprintf(stderr, "argument is not initialized \n");
    exit(0);
  }

  for (i = 0; i < xargc; ++i) {
    if (0 == strncmp(key, xargv[i], strlen(key)) &&
        xargv[i][strlen(key)] == '=') {
      break;
    }
  }
  if (i < xargc) {
    val[0] = &xargv[i][strlen(key) + 1];
    return 1;
  }
  return 0;
}

int readparfile(const char *filename, const char *fmt, ...) {
  va_list ap;
  int ncomma, nequal, npercent;
  int i, len, i1, i2, i12, ret, nret, j;

  int ii;
  float ff;
  char ss[1024];

  FILE *fp;

  char buf[1024];
  char pattern[256];
#define NCOMMAPOS 200
  int commapos[NCOMMAPOS + 4];
  int npattern[NCOMMAPOS + 4] = {0};
  char *s;
  char *p;

  int innn, n, nn, nnn, totalsize = 1024;
  char *fbuf, *fbuf1;
  char *pp;

  /* Static The Numbers of Comma,equal,percent */
  ncomma = 0;
  nequal = 0;
  npercent = 0;
  commapos[0] = -1;
  len = strlen(fmt);
  for (i = 0; i < len; ++i) {
    switch (fmt[i]) {
      case ',': {
        ncomma++;
        commapos[ncomma] = i;
        if (ncomma >= NCOMMAPOS) {
          fprintf(stderr,
                  "Error: In FUNCTION: %s fmt[%s] is too long to parse\n",
                  __FUNCTION__, fmt);
          fprintf(stderr,
                  "       Set Macro NCOMMAPOS Larger and recompile it \n");
          exit(0);
        }
        break;
      }
      case '=':
        nequal++;
        break;
      case '%':
        npercent++;
        break;
      default:
        break;
    }
  }
  commapos[ncomma + 1] = len;

  /* Check Its Sanity */
  if ((ncomma + 1) != nequal || nequal != npercent) {
    fprintf(stderr, "Error: In FUNCTION: %s fmt[%s] [, = %%] do not Match\n",
            __FUNCTION__, fmt);
    exit(0);
  }

  /* Open File & Read Record */
  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: Parameter File [%s] Open Error\n", filename);
    exit(0);
  }

  /* Seek File to Find the Recod According to the Pattern Specified by fmt */
  nret = 0;
  va_start(ap, fmt);

  fbuf = (char *)malloc(sizeof(char) * totalsize);
  fbuf1 = (char *)malloc(sizeof(char) * totalsize);
  fbuf[0] = 0;
  n = 0;
  for (i = 0; i < nequal; ++i) {
    i1 = commapos[i];
    i2 = commapos[i + 1];
    i12 = i2 - i1 - 1;
    strncpy(pattern, &fmt[i1 + 1], i12);
    pattern[i12] = 0; /*eg a=%d */

    rewind(fp);
    fbuf[0] = 0;
    fbuf1[0] = 0;
    n = 0;

    while (!feof(fp)) {
      p = fgets(buf, 1024, fp);
      if (p == NULL) {
        continue;
      }
      nn = strlen(buf);
      if (nn == 0) {
        n = 0;
        fbuf[0] = 0;
        continue;
      }
      n += nn;

      if (n >= totalsize) {
        fbuf = (char *)realloc(fbuf, n + 16);
        fbuf1 = (char *)realloc(fbuf1, n + 16);
        totalsize = n;
        if (fbuf == NULL || fbuf1 == NULL) {
          fprintf(stderr, "Error: realloc error in readparfile\n");
          exit(0);
        }
      }

      if (buf[nn - 1] == '\n') /* current line stop */
      {
        strcat(fbuf, buf);
        n = 0;
      } else /* concatation  */
      {
        strcat(fbuf, buf);
        continue;
      }
      nnn = strlen(fbuf);
      if (fbuf[nnn - 1] == '\n' && nnn >= 2 &&
          fbuf[nnn - 2] == '\\') /* line goes on  */
      {
        fbuf[nnn - 2] = 0;
        continue;
      }
      if (fbuf[nnn - 1] == '\n') /* remove '\n' */
      {
        fbuf[nnn - 1] = 0;
      }
      /* Remove The Left-Side WhiteSpace */
      j = 0;
      while (isspace(fbuf[j])) {
        j++;
      }
      strcpy(fbuf1, &fbuf[j]);

      /* Justice Comment */
      if (fbuf1[0] == '#') {
        fbuf[0] = 0;
        n = 0;
        continue;
      }
      nnn = strlen(fbuf1);
      /* Trim Last Comment */
      for (innn = nnn - 1; innn >= 0; --innn) {
        if (fbuf1[innn] == '#') fbuf1[innn] = 0;
      }
      nnn = strlen(fbuf1);
      for (innn = nnn - 1; innn >= 0; --innn) {
        if (isspace(fbuf1[innn]))
          fbuf1[innn] = 0;
        else
          break;
      }
      if (strlen(fbuf1) == 0) {
        fbuf[0] = 0;
        n = 0;
        continue;
      }
      switch (pattern[i12 - 1]) {
        case 'd': {
          ret = sscanf(fbuf1, pattern, &ii);
          if (ret == 1) {
            ++npattern[i];
            if (nret >= nequal) return -1;
            *va_arg(ap, int *) = ii;
            ++nret;
          }
          break;
        }
        case 'g':
        case 'f': {
          ret = sscanf(fbuf1, pattern, &ff);
          if (ret == 1) {
            ++npattern[i];
            if (nret >= nequal) return -1;
            *va_arg(ap, float *) = ff;
            ++nret;
          }

          break;
        }
        case 's': {
          ret = sscanf(fbuf1, pattern, ss);
          if (ret == 1) {
            ++npattern[i];
            if (nret >= nequal) return -1;
            s = va_arg(ap, char *);
            j = 0;
            while (ss[j] && j < 1023) {
              s[j] = ss[j];
              j++;
            }
            s[j] = 0;

            ++nret;
          }
          break;
        }
        case 'S': {
          nn = strlen(pattern);
          if (0 == strncmp(fbuf1, pattern, nn - 2)) {
            ret = 1;
            ++npattern[i];
            if (nret >= nequal) return -1;
            pp = (char *)malloc(sizeof(char) * totalsize);
            strcpy(pp, &fbuf1[nn - 2]);
            (*va_arg(ap, char **)) = pp;
            ++nret;
          }
        }
      }
      fbuf[0] = 0;
      n = 0;
      if (nret > nequal) return -1;
    }
  }
  va_end(ap);

  fclose(fp);

  free(fbuf);
  free(fbuf1);

  for (i = 0; i < nequal; ++i) {
    if (npattern[i] != 1) {
      return -1;
    }
  }

  return nret;
#undef NCOMMAPOS
}

void time_offset_mute_santity(float *tomute, int ntomute) {
  int i, n;
  n = ntomute;

  for (i = 0; i < n - 1; ++i) {
    if (tomute[i * 2] > tomute[(i + 1) * 2]) {
      tomute[(i + 1) * 2] = 2.0 * tomute[i * 2] - tomute[(i + 1) * 2];
    }

    if (tomute[i * 2 + 1] > tomute[(i + 1) * 2 + 1]) {
      tomute[(i + 1) * 2 + 1] =
          2.0 * tomute[i * 2 + 1] - tomute[(i + 1) * 2 + 1];
    }
  }

  for (i = 0; i < n - 1; ++i) {
    if (tomute[i * 2] >= tomute[(i + 1) * 2]) {
      tomute[(i + 1) * 2] = tomute[i * 2] + 0.01;
    }
    if (tomute[i * 2 + 1] >= tomute[(i + 1) * 2 + 1]) {
      tomute[(i + 1) * 2 + 1] = tomute[i * 2 + 1] + 10.0;
    }
  }
}

char *strrev1(char *p) {
  char *q = p;
  if (!p || !*p) return p;
  while (q && *q) ++q;
  for (--q; p < q; ++p, --q) *p = *p ^ *q, *q = *p ^ *q, *p = *p ^ *q;

  return p;
}

size_t fsize(const char *path) {
  struct stat st;
  stat(path, &st);
  return st.st_size;
}

char *strrtrim(char *src, int c) {
  int i, len;

  if (src)
    len = strlen(src);
  else
    return src;

  i = len - 1;
  while (i >= 0 && src[i] == c) {
    src[i] = 0;
    i--;
  }

  return src;
}

void parsepath(char *path, char *projpath, char *projname, int *iproc,
               char *task) {
  int i, j, len;
  char buf[16];

  /* Project Task */
  len = strlen(path);

  j = 0;
  for (i = len - 1; i >= 0; --i) {
    if (path[i] != '/') {
      task[j] = path[i];
      ++j;
    } else {
      break;
    }
  }
  task[j] = 0;
  strrev1(task);
  --i;

  j = 0;
  for (; i >= 0; --i) {
    if (path[i] != '/') {
      buf[j] = path[i];
      ++j;
    } else {
      break;
    }
  }
  strcpy(projpath, path);
  projpath[i] = 0;

  buf[j] = 0;
  strrev1(buf);
  --i;

  sscanf(buf, "%d", iproc);

  j = 0;
  for (; i >= 0; --i) {
    if (path[i] != '/') {
      projname[j] = path[i];
      ++j;
    } else {
      break;
    }
  }
  projname[j] = 0;
  strrev1(projname);
  --i;
}

float *parsetuple(const char *str, int *ntuple, int *nelem) {
  int nleft, nright;
  int nlen;
  int nel1, nelv, ntup;

  int i, j, k, ret, ileft, iright;
  float *p;
  char buf[128];

  nlen = strlen(str);
  nleft = 0;
  nright = 0;
  for (i = 0; i < nlen; ++i) {
    if (str[i] == '(') ++nleft;
    if (str[i] == ')') ++nright;
  }

  if (nleft != nright || nleft == 0 || nright == 0) {
    *ntuple = 0;
    *nelem = 0;

    fprintf(stderr,
            "Warn: String [%s]'s brackets do not match,or has no element\n",
            str);

    return NULL;
  }

  /* Search For the First Left Bracket */
  for (ileft = 0; ileft < nlen; ++ileft) {
    if (str[ileft] == '(') break;
  }

  /* Search For the First Right Bracket */
  for (iright = 0; iright < nlen; ++iright) {
    if (str[iright] == ')') break;
  }

  /* Check How Many Elements in the First Brackets */
  nel1 = 0;
  for (i = ileft + 1; i < iright; ++i) {
    if ((str[i] >= '0' && str[i] <= '9') || (str[i] == '.') ||
        (str[i] == 'E') || (str[i] == 'e')) {
      while ((str[i] >= '0' && str[i] <= '9') || (str[i] == '.') ||
             (str[i] == 'E') || (str[i] == 'e')) {
        ++i;
      }
      --i;
      ++nel1;
    } else {
      while (str[i] == ' ' || str[i] == '\t' || str[i] == ',' ||
             str[i] == ';') {
        ++i;
      }
      --i;
    }
  }
  ++ileft;
  ++iright;

  ntup = 1;
  nelv = 0;
  /* Check Each Tuple's Elements */
  while (ileft < nlen && iright < nlen) {
    for (; ileft < nlen; ++ileft) {
      if (str[ileft] == '(') break;
    }
    for (; iright < nlen; ++iright) {
      if (str[iright] == ')') break;
    }

    if (ileft >= nlen || iright >= nlen) {
      break;
    }

    nelv = 0;
    for (i = ileft + 1; i < iright; ++i) {
      if ((str[i] >= '0' && str[i] <= '9') || (str[i] == '.') ||
          (str[i] == 'E') || (str[i] == 'e')) {
        while ((str[i] >= '0' && str[i] <= '9') || str[i] == '.' ||
               (str[i] == 'E') || (str[i] == 'e')) {
          ++i;
        }
        --i;
        ++nelv;
      } else {
        while (str[i] == ' ' || str[i] == '\t' || str[i] == ',' ||
               str[i] == ';') {
          ++i;
        }
        --i;
      }
    }
    ++ileft;
    ++iright;

    if (nelv != nel1) {
      *ntuple = 0;
      *nelem = 0;
      fprintf(stderr, "Error: String [%s]'s Elements Count do not match\n",
              str);
      return NULL;
    }

    ++ntup;
  }

  /* Allocate Memory For the Values(Elements) */
  p = (float *)malloc(sizeof(float) * ntup * nelv);
  if (p == NULL) {
    fprintf(stderr, "Error: Memory Allocation Error in %s\n", __FUNCTION__);
    *ntuple = 0;
    *nelem = 0;
    return NULL;
  }

  ileft = 0;
  iright = 0;
  k = 0;
  /* Get Values in the Elements */
  while (ileft < nlen && iright < nlen) {
    for (; ileft < nlen; ++ileft) {
      if (str[ileft] == '(') break;
    }
    for (; iright < nlen; ++iright) {
      if (str[iright] == ')') break;
    }

    if (ileft >= nlen || iright >= nlen) {
      break;
    }

    for (i = ileft + 1; i < iright; ++i) {
      if ((str[i] >= '0' && str[i] <= '9') || (str[i] == '.') ||
          (str[i] == 'E') || (str[i] == 'e')) {
        j = 0;
        while ((str[i] >= '0' && str[i] <= '9') || str[i] == '.' ||
               (str[i] == 'E') || (str[i] == 'e')) {
          buf[j] = str[i];
          ++j;
          ++i;
        }
        --i;
        buf[j] = 0;

        ret = sscanf(buf, "%f", &p[k]);
        if (ret != 1) {
          fprintf(stderr, "Error: Elements Parse Error\n");
          *ntuple = 0;
          *nelem = 0;
          return NULL;
        }
        k++;
      } else {
        while (str[i] == ' ' || str[i] == '\t' || str[i] == ',' ||
               str[i] == ';') {
          ++i;
        }
        --i;
      }
    }
    ++ileft;
    ++iright;
  }

  if (k != ntup * nelv) {
    *ntuple = 0;
    *nelem = 0;
    fprintf(stderr, "Error: String [%s]'s Parse Error\n", str);
    return NULL;
  }

  *ntuple = ntup;
  *nelem = nelv;

  return p;
}

void free1c(char *p) {
  if (p == NULL) return;
  free(p);
}

void free2c(char **p) {
  if (p == NULL) return;
  if (p[0] == NULL) {
    free(p);
    return;
  }

  free(p[0]);
  free(p);
}

void free1i(int *p) {
  if (p == NULL) return;
  free(p);
}
void free2i(int **p) {
  if (p == NULL) return;
  if (p[0] == NULL) {
    free(p);
    return;
  }

  free(p[0]);
  free(p);
}
void free3i(int ***p) {
  if (p == NULL) return;
  if (p[0] == NULL) {
    free(p);
    return;
  }
  if (p[0][0] == NULL) {
    free(p[0]);
    free(p);
    return;
  }

  free(p[0][0]);
  free(p[0]);
  free(p);
}

void free1f(float *p) {
  if (p == NULL) return;
  free(p);
}
void free2f(float **p) {
  if (p == NULL) return;
  if (p[0] == NULL) {
    free(p);
    return;
  }
  free(p[0]);
  free(p);
}
void free3f(float ***p) {
  if (p == NULL) return;
  if (p[0] == NULL) {
    free(p);
    return;
  }
  if (p[0][0] == NULL) {
    free(p[0]);
    free(p);
    return;
  }
  free(p[0][0]);
  free(p[0]);
  free(p);
}

void free4f(float ****p) {
  if (p == NULL) return;
  if (p[0] == NULL) {
    free(p);
    return;
  }
  if (p[0][0] == NULL) {
    free(p[0]);
    free(p);
    return;
  }
  if (p[0][0][0] == NULL) {
    free(p[0][0]);
    free(p[0]);
    free(p);
    return;
  }
  free(p[0][0][0]);
  free(p[0][0]);
  free(p[0]);
  free(p);
}

void free1d(double *p) {
  if (p == NULL) return;
  free(p);
}
void free2d(double **p) {
  if (p == NULL) return;
  if (p[0] == NULL) {
    free(p);
    return;
  }
  free(p[0]);
  free(p);
}

/* get the actual coor from the velocity file */
/* it is stored in gx,gy, scalel              */
void get_vel_coor(const char *velfilepath, int line, coor_t *vel_coor,
                  short *vel_scalel, int fci, int nxi, int fcv, int nxv,
                  int flv, int nyv) {
  FILE *fp;
  segy tr;

  fp = fopen(velfilepath, "rb");
  if (NULL == fp) {
    fprintf(stderr, "Velocity File open Error\n");
    exit(1);
  }

  fgettr(fp, &tr);
  //*vel_scalel = tr.scalel;
  *vel_scalel = tr.scalco;

  if (line >= flv + nyv) line = flv + nyv - 1;

  fseek(fp, (sizeof(float) * tr.ns + 240L) * ((line - flv) * nxv + fci - fcv),
        SEEK_SET);

  for (int itr = 0; itr < nxi; ++itr) {
    fgettr(fp, &tr);

    vel_coor[itr].sx = tr.sx;
    vel_coor[itr].sy = tr.sy;
  }

  fclose(fp);
}
/**********************************/
/* 由道头中的scalco参数获得乘因子 */
/* sx_actual=tr.sx*scale          */
/**********************************/
float getscale(short scalco) {
  float scale;

  if (scalco == 0) {
    scale = 1.0;
  } else if (scalco > 0) {
    scale = scalco;
  } else if (scalco < 0) {
    scale = -1.0 / scalco;
  }

  return scale;
}

/* Calculate the time at given offset value,with   */
/* inner  point: linear interpolation              */
/* outter point: using the nearest boundary point  */
float time_offset_mute(float *tomute, int ntomute, float offset) {
  int i, n;
  float o1, o2, t1, t2, o, t;
  n = ntomute;

  o = offset;

  if (n >= 2) {
    if (o < tomute[1]) {
      /* 线性插值 */
      o1 = tomute[1];
      o2 = tomute[3];
      t1 = tomute[0];
      t2 = tomute[2];

      t = t1 + (t2 - t1) / (o2 - o1) * (o - o1);
    } else if (o >= tomute[(n - 1) * 2 + 1]) {
      o1 = tomute[(n - 2) * 2 + 1];
      o2 = tomute[(n - 1) * 2 + 1];
      t1 = tomute[(n - 2) * 2];
      t2 = tomute[(n - 1) * 2];

      t = t1 + (t2 - t1) / (o2 - o1) * (o - o1);
    } else {
      for (i = 0; i < n - 1; ++i) {
        if (tomute[i * 2 + 1] <= o && o < tomute[(i + 1) * 2 + 1]) {
          o1 = tomute[i * 2 + 1];
          o2 = tomute[(i + 1) * 2 + 1];
          t1 = tomute[i * 2];
          t2 = tomute[(i + 1) * 2];

          t = t1 + (t2 - t1) / (o2 - o1) * (o - o1);
        }
      }
    }
  } else {
    t = tomute[0];
  }

  return t;
}