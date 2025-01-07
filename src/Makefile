NVCC ?= nvcc

all: q-psdm

segy.o: segy.cc
	g++ -O2 -std=c++98 -c $< -o $@

psdmpkg.o: psdmpkg.cc
	g++ -O2 -std=c++98 -c $< -o $@

fft.o: fft.cc
	g++ -O2 -std=c++98 -c $< -o $@

q-psdm: q-psdm.cu segy.o psdmpkg.o fft.o
	${NVCC} -o $@ $< segy.o psdmpkg.o fft.o -O2 -Xcompiler "-W -Wall" --use_fast_math -Xptxas -v ${GENCODE_FLAGS} -D_DEFAULT_SOURCE
	cp $@ ../bin/

clean:
	rm -rf q-psdm segy.o psdmpkg.o fft.o ../bin/q-psdm