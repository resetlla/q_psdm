#ifndef SEGY_H_
#define SEGY_H_

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <string>
#include <vector>

#define SU_NFLTS 65535 /*   Arbitrary limit on data array size	*/

/* Define Segy Format IO Field in its Header */
#define BHDR_ALL

/* total pased keys */
#define SEGY_TOTAL_KEY_NUM 75

/*#define SEGY_ALL*/

#define SEGY_TRACL
#define SEGY_TRACF
#define SEGY_TRACR

#define SEGY_FLDR
#define SEGY_CDP
#define SEGY_CDPT
#define SEGY_OFFSET
#define SEGY_NHS
#define SEGY_NVS
#define SEGY_EP

#define SEGY_GELEV
#define SEGY_SELEV
#define SEGY_SDEPTH
#define SEGY_GDEL
#define SEGY_SDEL

#define SEGY_SCALEL
#define SEGY_SCALCO
#define SEGY_SX
#define SEGY_SY
#define SEGY_GX
#define SEGY_GY

#define SEGY_COUNIT
#define SEGY_DELRT

#define SEGY_TSTAT
#define SEGY_MUTS
#define SEGY_MUTE

#define SEGY_NS
#define SEGY_DT

#define SEGY_TRID
#define SEGY_MINUTE
#define SEGY_WEVEL

#define SEGY_CDPX
#define SEGY_CDPY
#define SEGY_ILINE
#define SEGY_XLINE
#define SEGY_TRWF

/* End Define Segy Format IO Field in its Header */

#define fvgettr(fp, tr) fgettr(fp, tr)
#define fvputtr(fp, tr) fputtr(fp, tr)
#define swaptr(buf, tr) sgettr(buf, tr)

#define SWAP2(n) ((((n)&0XFF) << 8) | (((n)&0XFF00) >> 8))
#define SWAP4(n)                                                      \
  ((((n)&0XFF) << 24) | (((n)&0XFF00) << 8) | (((n) >> 8) & 0XFF00) | \
   (((n)&0XFF000000) >> 24))

typedef struct segy segy;
typedef struct hd3600 hd3600;
typedef struct keyindex keyindex;

/* keyidx[75] */

extern int keyindex_inited;

/* Standard Segy Format File IO */

void initkeyindex(keyindex *keyidx);

/* setkeyindex(keyidx,"sx",73,4); from 1 not 0 */
void setkeyindex(keyindex *keyidx, const char *key, int start, int bytes);

void inithd3600(hd3600 *hd, int ns, int dt_um, short int format);
void sethd3600text(hd3600 *hd, const char *text);

void fsgethd(FILE *fp, hd3600 *hd);
void fsputhd(FILE *fp, const hd3600 *hd);

int fsgettr(FILE *fp, segy *tr, short int format, const keyindex *keyidx);
void fsputtr(FILE *fp, const segy *tr, short int format,
             const keyindex *keyidx);
int ssgettr(unsigned char *buf, segy *tr, short int format,
            const keyindex *keyidx);

/* SU Format File IO */
int fgettr(FILE *fp, segy *tr);
void fputtr(FILE *fp, const segy *tr);

int sgettr(unsigned char *buf, segy *tr);
void sputtr(unsigned char *buf, const segy *tr);
/* some test shows,the swaptr_macro is much slower than swaptr */
void swaptr_macro(unsigned char *buf, segy *tr);

struct SegyInfo {
  int ns;
  int ntr;
  int dt;  // micro-seconds
};

/* Function Pototype End */

/* TYPEDEFS */
struct segy { /* segy - trace identification header */

  int tracl; /* Trace sequence number within line
                --numbers continue to increase if the
                same line continues across multiple
                SEG Y files.
                byte# 1-4
              */

  int tracr; /* Trace sequence number within SEG Y file
                ---each file starts with trace sequence
                one
                byte# 5-8
              */

  int fldr; /* Original field record number
               byte# 9-12
            */

  int tracf; /* Trace number within original field record
                byte# 13-16
             */

  int ep; /* energy source point number
             ---Used when more than one record occurs
             at the same effective surface location.
             byte# 17-20
           */

  int cdp; /* Ensemble number (i.e. CDP, CMP, CRP,...)
              byte# 21-24
           */

  int cdpt; /* trace number within the ensemble
               ---each ensemble starts with trace number one.
               byte# 25-28
             */

  short trid; /* trace identification code:
              -1 = Other
               0 = Unknown
               1 = Seismic data
               2 = Dead
               3 = Dummy
               4 = Time break
               5 = Uphole
               6 = Sweep
               7 = Timing
               8 = Water break
               9 = Near-field gun signature
              10 = Far-field gun signature
              11 = Seismic pressure sensor
              12 = Multicomponent seismic sensor
                      - Vertical component
              13 = Multicomponent seismic sensor
                      - Cross-line component
              14 = Multicomponent seismic sensor
                      - in-line component
              15 = Rotated multicomponent seismic sensor
                      - Vertical component
              16 = Rotated multicomponent seismic sensor
                      - Transverse component
              17 = Rotated multicomponent seismic sensor
                      - Radial component
              18 = Vibrator reaction mass
              19 = Vibrator baseplate
              20 = Vibrator estimated ground force
              21 = Vibrator reference
              22 = Time-velocity pairs
              23 ... N = optional use
                      (maximum N = 32,767)

              Following are CWP id flags:

              109 = autocorrelation
              110 = Fourier transformed - no packing
                   xr[0],xi[0], ..., xr[N-1],xi[N-1]
              111 = Fourier transformed - unpacked Nyquist
                   xr[0],xi[0],...,xr[N/2],xi[N/2]
              112 = Fourier transformed - packed Nyquist
                   even N:
                   xr[0],xr[N/2],xr[1],xi[1], ...,
                      xr[N/2 -1],xi[N/2 -1]
                      (note the exceptional second entry)
                   odd N:
                   xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
                      xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
                      (note the exceptional second & last entries)
              113 = Complex signal in the time domain
                   xr[0],xi[0], ..., xr[N-1],xi[N-1]
              114 = Fourier transformed - amplitude/phase
                   a[0],p[0], ..., a[N-1],p[N-1]
              115 = Complex time signal - amplitude/phase
                   a[0],p[0], ..., a[N-1],p[N-1]
              116 = Real part of complex trace from 0 to Nyquist
              117 = Imag part of complex trace from 0 to Nyquist
              118 = Amplitude of complex trace from 0 to Nyquist
              119 = Phase of complex trace from 0 to Nyquist
              121 = Wavenumber time domain (k-t)
              122 = Wavenumber frequency (k-omega)
              123 = Envelope of the complex time trace
              124 = Phase of the complex time trace
              125 = Frequency of the complex time trace
              126 = log amplitude
              127 = real cepstral domain F(t_c)= invfft[log[fft(F(t)]]
              130 = Depth-Range (z-x) traces
              201 = Seismic data packed to bytes (by supack1)
              202 = Seismic data packed to 2 bytes (by supack2)
                 byte# 29-30
              */

  short nvs; /* Number of vertically summed traces yielding
                this trace. (1 is one trace,
                2 is two summed traces, etc.)
                byte# 31-32
              */

  short nhs; /* Number of horizontally summed traces yielding
                this trace. (1 is one trace
                2 is two summed traces, etc.)
                byte# 33-34
              */

  short duse; /* Data use:
                      1 = Production
                      2 = Test
                 byte# 35-36
               */

  int offset; /* Distance from the center of the source point
                 to the center of the receiver group
                 (negative if opposite to direction in which
                 the line was shot).
                 byte# 37-40
               */

  int gelev; /* Receiver group elevation from sea level
                (all elevations above the Vertical datum are
                positive and below are negative).
                byte# 41-44
              */

  /*	union
          {*/
  int selev;  /* Surface elevation at source.
         byte# 45-48
       */
              /*		int esx;
                              int egx;
                      };
            
                      union
                      {*/
  int sdepth; /* Source depth below surface (a positive number).
         byte# 49-52
       */
              /*		int esy;
                              int egy;
                      };*/

  int gdel; /* Datum elevation at receiver group.
               byte# 53-56
            */

  int sdel; /* Datum elevation at source.
               byte# 57-60
            */

  int swdep; /* Water depth at source.
                byte# 61-64
             */

  int gwdep; /* Water depth at receiver group.
                byte# 65-68
             */

  short scalel; /* Scalar to be applied to the previous 7 entries
                   to give the real value.
                   Scalar = 1, +10, +100, +1000, +10000.
                   If positive, scalar is used as a multiplier,
                   if negative, scalar is used as a divisor.
                   byte# 69-70
                 */

  short scalco; /* Scalar to be applied to the next 4 entries
                   to give the real value.
                   Scalar = 1, +10, +100, +1000, +10000.
                   If positive, scalar is used as a multiplier,
                   if negative, scalar is used as a divisor.
                   byte# 71-72
                 */

  int sx; /* Source coordinate - X
             byte# 73-76
          */

  int sy; /* Source coordinate - Y
             byte# 77-80
          */

  int gx; /* Group coordinate - X
             byte# 81-84
          */

  int gy; /* Group coordinate - Y
             byte# 85-88
          */

  short counit; /* Coordinate units: (for previous 4 entries and
                        for the 7 entries before scalel)
                   1 = Length (meters or feet)
                   2 = Seconds of arc
                   3 = Decimal degrees
                   4 = Degrees, minutes, seconds (DMS)

                In case 2, the X values are longitude and
                the Y values are latitude, a positive value designates
                the number of seconds east of Greenwich
                        or north of the equator

                In case 4, to encode +-DDDMMSS
                counit = +-DDD*10^4 + MM*10^2 + SS,
                with scalco = 1. To encode +-DDDMMSS.ss
                counit = +-DDD*10^6 + MM*10^4 + SS*10^2
                with scalco = -100.
                   byte# 89-90
                */

  short wevel; /* Weathering velocity.
                  byte# 91-92
               */

  short swevel; /* Subweathering velocity.
                   byte# 93-94
                */

  short sut; /* Uphole time at source in milliseconds.
                byte# 95-96
             */

  short gut; /* Uphole time at receiver group in milliseconds.
                byte# 97-98
             */

  short sstat; /* Source static correction in milliseconds.
                  byte# 99-100
               */

  short gstat; /* Group static correction  in milliseconds.
                  byte# 101-102
               */

  short tstat; /* Total static applied  in milliseconds.
                  (Zero if no static has been applied.)
                  byte# 103-104
               */

  short laga; /* Lag time A, time in ms between end of 240-
                 byte trace identification header and time
                 break, positive if time break occurs after
                 end of header, time break is defined as
                 the initiation pulse which maybe recorded
                 on an auxiliary trace or as otherwise
                 specified by the recording system
                 byte# 105-106
              */

  short lagb; /* lag time B, time in ms between the time break
                 and the initiation time of the energy source,
                 may be positive or negative
                 byte# 107-108
              */

  short delrt; /* delay recording time, time in ms between
                  initiation time of energy source and time
                  when recording of data samples begins
                  (for deep water work if recording does not
                  start at zero time)
                  byte# 109-110
               */

  short muts; /* mute time--start
                 byte# 111-112
              */

  short mute; /* mute time--end
                 byte# 113-114
              */

  unsigned short ns; /* number of samples in this trace
                byte# 115-116
             */

  unsigned short dt; /* sample interval; in micro-seconds
                byte# 117-118
             */

  short gain; /* gain type of field instruments code:
                      1 = fixed
                      2 = binary
                      3 = floating point
                      4 ---- N = optional use
                 byte# 119-120
              */

  short igc; /* instrument gain constant
                byte# 121-122
             */

  short igi; /* instrument early or initial gain
                byte# 123-124
             */

  short corr; /* correlated:
                      1 = no
                      2 = yes
                 byte# 125-126
              */

  short sfs; /* sweep frequency at start
                byte# 127-128
             */

  short sfe; /* sweep frequency at end
                byte# 129-130
             */

  short slen; /* sweep length in ms
                 byte# 131-132
              */

  short styp; /* sweep type code:
                      1 = linear
                      2 = cos-squared
                      3 = other
                 byte# 133-134
              */

  short stas; /* sweep trace length at start in ms
                 byte# 135-136
              */

  short stae; /* sweep trace length at end in ms
                 byte# 137-138
              */

  short tatyp; /* taper type: 1=linear, 2=cos^2, 3=other
                  byte# 139-140
               */

  short afilf; /* alias filter frequency if used
                  byte# 141-142
               */

  short afils; /* alias filter slope
                  byte# 143-144
               */

  short nofilf; /* notch filter frequency if used
                   byte# 145-146
                */

  short nofils; /* notch filter slope
                   byte# 147-148
                */

  short lcf; /* low cut frequency if used
                byte# 149-150
             */

  short hcf; /* high cut frequncy if used
                byte# 151-152
             */

  short lcs; /* low cut slope
                byte# 153-154
             */

  short hcs; /* high cut slope
                byte# 155-156
             */

  short year; /* year data recorded
                 byte# 157-158
              */

  short day; /* day of year
                byte# 159-160
             */

  short hour; /* hour of day (24 hour clock)
                 byte# 161-162
              */

  short minute; /* minute of hour
                   byte# 163-164
                */

  short sec; /* second of minute
                byte# 165-166
             */

  short timbas; /* time basis code:
                        1 = local
                        2 = GMT
                        3 = other
                   byte# 167-168
                */

  short trwf; /* trace weighting factor, defined as 1/2^N
                 volts for the least sigificant bit
                 byte# 169-170
              */

  short grnors; /* geophone group number of roll switch
                   position one
                   byte# 171-172
                */

  short grnofr; /* geophone group number of trace one within
                   original field record
                   byte# 173-174
                */

  short grnlof; /* geophone group number of last trace within
                   original field record
                   byte# 175-176
                */

  short gaps; /* gap size (total number of groups dropped)
                 byte# 177-178
              */

  short otrav; /* overtravel taper code:
                       1 = down (or behind)
                       2 = up (or ahead)
                  byte# 179-180
               */
               /* industrial used position */
  int cdpx;    /* sample spacing for non-seismic data
                  byte# 181-184
               */

  int cdpy; /* first sample location for non-seismic data
               byte# 185-188
            */

  int iline; /* sample spacing between traces
                byte# 189-192
             */

  int xline; /* first trace location
                byte# 193-196
             */
             /* industrial used position end */
  int sp;    /* negative of power used for dynamic
                        range compression
                        byte# 197-200
                     */

  int unscale; /* reciprocal of scaling factor to normalize
                  range
                  byte# 201-204
               */

  int ntr; /* number of traces
              byte# 205-208
           */

  short mark; /* mark selected traces
                 byte# 209-210
              */

  short shortpad; /* alignment padding
                         byte# 211-212
                      */

  short unass[14]; /* unassigned--NOTE: last entry causes
              a break in the word alignment, if we REALLY
              want to maintain 240 bytes, the following
              entry should be an odd number of short/UINT2
              OR do the insertion above the "mark" keyword
              entry
              byte# 213-240
           */

  float data[SU_NFLTS];
};

typedef struct { /* segy - trace identification header */

  int tracl; /* Trace sequence number within line
                --numbers continue to increase if the
                same line continues across multiple
                SEG Y files.
                byte# 1-4
              */

  int tracr; /* Trace sequence number within SEG Y file
                ---each file starts with trace sequence
                one
                byte# 5-8
              */

  int fldr; /* Original field record number
               byte# 9-12
            */

  int tracf; /* Trace number within original field record
                byte# 13-16
             */

  int ep; /* energy source point number
             ---Used when more than one record occurs
             at the same effective surface location.
             byte# 17-20
           */

  int cdp; /* Ensemble number (i.e. CDP, CMP, CRP,...)
              byte# 21-24
           */

  int cdpt; /* trace number within the ensemble
               ---each ensemble starts with trace number one.
               byte# 25-28
             */

  short trid; /* trace identification code:
              -1 = Other
               0 = Unknown
               1 = Seismic data
               2 = Dead
               3 = Dummy
               4 = Time break
               5 = Uphole
               6 = Sweep
               7 = Timing
               8 = Water break
               9 = Near-field gun signature
              10 = Far-field gun signature
              11 = Seismic pressure sensor
              12 = Multicomponent seismic sensor
                      - Vertical component
              13 = Multicomponent seismic sensor
                      - Cross-line component
              14 = Multicomponent seismic sensor
                      - in-line component
              15 = Rotated multicomponent seismic sensor
                      - Vertical component
              16 = Rotated multicomponent seismic sensor
                      - Transverse component
              17 = Rotated multicomponent seismic sensor
                      - Radial component
              18 = Vibrator reaction mass
              19 = Vibrator baseplate
              20 = Vibrator estimated ground force
              21 = Vibrator reference
              22 = Time-velocity pairs
              23 ... N = optional use
                      (maximum N = 32,767)

              Following are CWP id flags:

              109 = autocorrelation
              110 = Fourier transformed - no packing
                   xr[0],xi[0], ..., xr[N-1],xi[N-1]
              111 = Fourier transformed - unpacked Nyquist
                   xr[0],xi[0],...,xr[N/2],xi[N/2]
              112 = Fourier transformed - packed Nyquist
                   even N:
                   xr[0],xr[N/2],xr[1],xi[1], ...,
                      xr[N/2 -1],xi[N/2 -1]
                      (note the exceptional second entry)
                   odd N:
                   xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
                      xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
                      (note the exceptional second & last entries)
              113 = Complex signal in the time domain
                   xr[0],xi[0], ..., xr[N-1],xi[N-1]
              114 = Fourier transformed - amplitude/phase
                   a[0],p[0], ..., a[N-1],p[N-1]
              115 = Complex time signal - amplitude/phase
                   a[0],p[0], ..., a[N-1],p[N-1]
              116 = Real part of complex trace from 0 to Nyquist
              117 = Imag part of complex trace from 0 to Nyquist
              118 = Amplitude of complex trace from 0 to Nyquist
              119 = Phase of complex trace from 0 to Nyquist
              121 = Wavenumber time domain (k-t)
              122 = Wavenumber frequency (k-omega)
              123 = Envelope of the complex time trace
              124 = Phase of the complex time trace
              125 = Frequency of the complex time trace
              126 = log amplitude
              127 = real cepstral domain F(t_c)= invfft[log[fft(F(t)]]
              130 = Depth-Range (z-x) traces
              201 = Seismic data packed to bytes (by supack1)
              202 = Seismic data packed to 2 bytes (by supack2)
                 byte# 29-30
              */

  short nvs; /* Number of vertically summed traces yielding
                this trace. (1 is one trace,
                2 is two summed traces, etc.)
                byte# 31-32
              */

  short nhs; /* Number of horizontally summed traces yielding
                this trace. (1 is one trace
                2 is two summed traces, etc.)
                byte# 33-34
              */

  short duse; /* Data use:
                      1 = Production
                      2 = Test
                 byte# 35-36
               */

  int offset; /* Distance from the center of the source point
                 to the center of the receiver group
                 (negative if opposite to direction in which
                 the line was shot).
                 byte# 37-40
               */

  int gelev; /* Receiver group elevation from sea level
                (all elevations above the Vertical datum are
                positive and below are negative).
                byte# 41-44
              */
             /*union
             {*/
  int selev; /* Surface elevation at source.
        byte# 45-48
      */
             /*	int esx;
                     int egx;
             };*/

  /*union
  {*/
  int sdepth; /* Source depth below surface (a positive number).
         byte# 49-52
       */
  /*	int esy;
          int egy;
  };*/

  int gdel; /* Datum elevation at receiver group.
               byte# 53-56
            */

  int sdel; /* Datum elevation at source.
               byte# 57-60
            */

  int swdep; /* Water depth at source.
                byte# 61-64
             */

  int gwdep; /* Water depth at receiver group.
                byte# 65-68
             */

  short scalel; /* Scalar to be applied to the previous 7 entries
                   to give the real value.
                   Scalar = 1, +10, +100, +1000, +10000.
                   If positive, scalar is used as a multiplier,
                   if negative, scalar is used as a divisor.
                   byte# 69-70
                 */

  short scalco; /* Scalar to be applied to the next 4 entries
                   to give the real value.
                   Scalar = 1, +10, +100, +1000, +10000.
                   If positive, scalar is used as a multiplier,
                   if negative, scalar is used as a divisor.
                   byte# 71-72
                 */

  int sx; /* Source coordinate - X
             byte# 73-76
          */

  int sy; /* Source coordinate - Y
             byte# 77-80
          */

  int gx; /* Group coordinate - X
             byte# 81-84
          */

  int gy; /* Group coordinate - Y
             byte# 85-88
          */

  short counit; /* Coordinate units: (for previous 4 entries and
                        for the 7 entries before scalel)
                   1 = Length (meters or feet)
                   2 = Seconds of arc
                   3 = Decimal degrees
                   4 = Degrees, minutes, seconds (DMS)

                In case 2, the X values are longitude and
                the Y values are latitude, a positive value designates
                the number of seconds east of Greenwich
                        or north of the equator

                In case 4, to encode +-DDDMMSS
                counit = +-DDD*10^4 + MM*10^2 + SS,
                with scalco = 1. To encode +-DDDMMSS.ss
                counit = +-DDD*10^6 + MM*10^4 + SS*10^2
                with scalco = -100.
                   byte# 89-90
                */

  short wevel; /* Weathering velocity.
                  byte# 91-92
               */

  short swevel; /* Subweathering velocity.
                   byte# 93-94
                */

  short sut; /* Uphole time at source in milliseconds.
                byte# 95-96
             */

  short gut; /* Uphole time at receiver group in milliseconds.
                byte# 97-98
             */

  short sstat; /* Source static correction in milliseconds.
                  byte# 99-100
               */

  short gstat; /* Group static correction  in milliseconds.
                  byte# 101-102
               */

  short tstat; /* Total static applied  in milliseconds.
                  (Zero if no static has been applied.)
                  byte# 103-104
               */

  short laga; /* Lag time A, time in ms between end of 240-
                 byte trace identification header and time
                 break, positive if time break occurs after
                 end of header, time break is defined as
                 the initiation pulse which maybe recorded
                 on an auxiliary trace or as otherwise
                 specified by the recording system
                 byte# 105-106
              */

  short lagb; /* lag time B, time in ms between the time break
                 and the initiation time of the energy source,
                 may be positive or negative
                 byte# 107-108
              */

  short delrt; /* delay recording time, time in ms between
                  initiation time of energy source and time
                  when recording of data samples begins
                  (for deep water work if recording does not
                  start at zero time)
                  byte# 109-110
               */

  short muts; /* mute time--start
                 byte# 111-112
              */

  short mute; /* mute time--end
                 byte# 113-114
              */

  unsigned short ns; /* number of samples in this trace
                byte# 115-116
             */

  unsigned short dt; /* sample interval; in micro-seconds
                byte# 117-118
             */

  short gain; /* gain type of field instruments code:
                      1 = fixed
                      2 = binary
                      3 = floating point
                      4 ---- N = optional use
                 byte# 119-120
              */

  short igc; /* instrument gain constant
                byte# 121-122
             */

  short igi; /* instrument early or initial gain
                byte# 123-124
             */

  short corr; /* correlated:
                      1 = no
                      2 = yes
                 byte# 125-126
              */

  short sfs; /* sweep frequency at start
                byte# 127-128
             */

  short sfe; /* sweep frequency at end
                byte# 129-130
             */

  short slen; /* sweep length in ms
                 byte# 131-132
              */

  short styp; /* sweep type code:
                      1 = linear
                      2 = cos-squared
                      3 = other
                 byte# 133-134
              */

  short stas; /* sweep trace length at start in ms
                 byte# 135-136
              */

  short stae; /* sweep trace length at end in ms
                 byte# 137-138
              */

  short tatyp; /* taper type: 1=linear, 2=cos^2, 3=other
                  byte# 139-140
               */

  short afilf; /* alias filter frequency if used
                  byte# 141-142
               */

  short afils; /* alias filter slope
                  byte# 143-144
               */

  short nofilf; /* notch filter frequency if used
                   byte# 145-146
                */

  short nofils; /* notch filter slope
                   byte# 147-148
                */

  short lcf; /* low cut frequency if used
                byte# 149-150
             */

  short hcf; /* high cut frequncy if used
                byte# 151-152
             */

  short lcs; /* low cut slope
                byte# 153-154
             */

  short hcs; /* high cut slope
                byte# 155-156
             */

  short year; /* year data recorded
                 byte# 157-158
              */

  short day; /* day of year
                byte# 159-160
             */

  short hour; /* hour of day (24 hour clock)
                 byte# 161-162
              */

  short minute; /* minute of hour
                   byte# 163-164
                */

  short sec; /* second of minute
                byte# 165-166
             */

  short timbas; /* time basis code:
                        1 = local
                        2 = GMT
                        3 = other
                   byte# 167-168
                */

  short trwf; /* trace weighting factor, defined as 1/2^N
                 volts for the least sigificant bit
                 byte# 169-170
              */

  short grnors; /* geophone group number of roll switch
                   position one
                   byte# 171-172
                */

  short grnofr; /* geophone group number of trace one within
                   original field record
                   byte# 173-174
                */

  short grnlof; /* geophone group number of last trace within
                   original field record
                   byte# 175-176
                */

  short gaps; /* gap size (total number of groups dropped)
                 byte# 177-178
              */

  short otrav; /* overtravel taper code:
                       1 = down (or behind)
                       2 = up (or ahead)
                  byte# 179-180
               */

  /* begin Unocal SU segy.h differences */

  /* cwp local assignments */
  float d1; /* sample spacing for non-seismic data
               byte# 181-184
            */

  float f1; /* first sample location for non-seismic data
               byte# 185-188
            */

  float d2; /* sample spacing between traces
               byte# 189-192
            */

  float f2; /* first trace location
               byte# 193-196
            */

  float ungpow; /* negative of power used for dynamic
                   range compression
                   byte# 197-200
                */

  float unscale; /* reciprocal of scaling factor to normalize
                    range
                    byte# 201-204
                 */

  short mark; /* mark selected traces
                 byte# 205-206
              */

  /* SLTSU local assignments */
  short mutb; /* mute time at bottom (start time)
                 bottom mute ends at last sample
                 byte# 207-208
              */
  float dz;   /* depth sampling interval in (m or ft)
              if =0.0, input are time samples
                 byte# 209-212
              */

  float fz; /* depth of first sample in (m or ft)
               byte# 213-116
            */

  short n2; /* number of traces per cdp or per shot
               byte# 217-218
            */

  short shortpad; /* alignment padding
                     byte# 219-220
                  */

  int ntr; /* number of traces
              byte# 221-224
           */

  /* SLTSU local assignments end */

  short unass[8]; /* unassigned
                     byte# 225-240
                  */
} segyh;

struct hd3600 {
  /* EBCDIC Header */
  unsigned char text[3200];

  /* Binary File Header */
  int jobid;        /* job identification number */
  int lino;         /* line number (only one line per reel) */
  int reno;         /* reel number */
  short int ntrpr;  /* number of data traces per record */
  short int nart;   /* number of auxiliary traces per record */
  short int hdt;    /* sample interval (microsecs) for this reel */
  short int dto;    /* same for original field recording */
  short int hns;    /* number of samples per trace for this reel */
  short int nso;    /* same for original field recording */
  short int format; /* data sample format code:
                                                  1 = 4-byte IBM Floating-point
                                             *2 = fixed point (4 bytes)
                                             *3 = fixed point (2 bytes)
                                             *4 = fixed point w/gain code (4
                       bytes) 5 = 4-byte IEEE Floating-point
                                             * Not Support */
  short int fold;   /* CDP fold expected per CDP ensemble */
  short int tsort;  /* trace sorting code:
                                                  1 = as recorded (no sorting)
                                                  2 = CDP ensemble
                                                  3 = single fold continuous
                       profile  4 = horizontally stacked */
  short int vscode; /* vertical sum code:
                                                  1 = no sum
                                                  2 = two sum ...
                                                  N = N sum (N = 32,767) */
  short int hsfs;   /* sweep frequency at start */
  short int hsfe;   /* sweep frequency at end */
  short int hslen;  /* sweep length (ms) */
  short int hstyp;  /* sweep type code:
                                                  1 = linear
                                                  2 = parabolic
                                                  3 = exponential
                                                  4 = other */
  short int schn;   /* trace number of sweep channel */
  short int hstas;  /* sweep trace taper length at start if
                                     tapered (the taper starts at zero time
                                     and is effective for this length) */
  short int hstae;  /* sweep trace taper length at end (the ending
                                     taper starts at sweep length minus the taper
                                     length at end) */
  short int htatyp; /* sweep trace taper type code:
                                                  1 = linear
                                                  2 = cos-squared
                                                  3 = other */
  short int hcorr;  /* correlated data traces code:
                                                  1 = no
                                                  2 = yes */
  short int bgrcv;  /* binary gain recovered code:
                                                  1 = yes
                                                  2 = no */
  short int rcvm;   /* amplitude recovery method code:
                                                  1 = none
                                                  2 = spherical divergence
                                                  3 = AGC
                                                  4 = other */
  short int mfeet;  /* measurement system code:
                                                  1 = meters
                                                  2 = feet */
  short int polyt;  /* impulse signal polarity code:
                                                  1 = increase in pressure or
                       upward  geophone case movement gives  negative number on
                       tape  2 = increase in pressure or upward  geophone case
                       movement gives  positive number on tape */
  short int
      vpol;                  /* vibratory polarity code:
                                                           code    seismic signal lags pilot by
                                                           1       337.5 to  22.5 degrees
                                                           2        22.5 to  67.5 degrees
                                                           3        67.5 to 112.5 degrees
                                                           4       112.5 to 157.5 degrees
                                                           5       157.5 to 202.5 degrees
                                                           6       202.5 to 247.5 degrees
                                                           7       247.5 to 292.5 degrees
                                                           8       293.5 to 337.5 degrees */
  unsigned char hunass[340]; /* unassigned */
};

struct keyindex {
  int start;
  int bytes;
};

#endif  // SEGY_H_
