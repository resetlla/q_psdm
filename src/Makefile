NVCC ?= nvcc

all: qmig3dd_multithreading

segy.o: segy.cc
	g++ -O2 -std=c++98 -c $< -o $@

psdmpkg.o: psdmpkg.cc
	g++ -O2 -std=c++98 -c $< -o $@

fft.o: fft.cc
	g++ -O2 -std=c++98 -c $< -o $@

qmig3dd_multithreading: qmig3dd_multithreading.cu segy.o psdmpkg.o fft.o
	${NVCC} -o $@ $< segy.o psdmpkg.o fft.o -O2 -Xcompiler "-W -Wall" --use_fast_math -Xptxas -v ${GENCODE_FLAGS} -D_DEFAULT_SOURCE
	cp $@ ../bin/

clean:
	rm -rf qmig3dd_multithreading segy.o psdmpkg.o fft.o ../bin/qmig3dd_multithreading