#ifndef SUFFT_H_
#define SUFFT_H_

#define FFTW_MEASURE 1
#define FFTW_ESTIMATE (1U << 6)

typedef struct {
  float *input;
  float *output;
  int num;
  int flag;
  int sign;
  int mode;  // 201: c2r; 102: r2c
} fftwf_plan;

typedef float fftwf_complex[2];

int get_opt_n_fftw(int n);
fftwf_plan fftwf_plan_dft_r2c_1d(int n, void *input, void *output, int flag);
fftwf_plan fftwf_plan_dft_c2r_1d(int n, void *input, void *output, int flag);

float *fftwf_alloc_real(int n);
fftwf_complex *fftwf_alloc_complex(int n);
void fftwf_free(void *ptr);

void fftwf_execute(fftwf_plan plan);
void fftwf_destroy_plan(fftwf_plan plan);

#endif  // SUFFT_H_
