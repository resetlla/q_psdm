# 1. Project Overview, File Structure, System Requirements
## Overview
    Q-PSDM is a Q-compensated prestack depth migration with GPU/CPU multithreading approach.
## File Structure
    .
    ├── README.md
    ├── Makefile
    ├── src
    │   ├── q-psdm.cu        # Core CUDA implementation
    │   ├── fft.cc           # FFT-related implementation
    │   ├── psdmpkg.cc       # PSDM core package
    │   ├── segy.cc          # SEG-Y data reading/writing
    │   ├── fft.h
    │   ├── sufft.h
    │   ├── psdmpkg.h
    │   ├── segy.h
    │   ├── queue.h
    ├── example              # Example data and scripts
    └── bin                  # Compiled executables (generated)
## System Requirements
	- CUDA 7.5 or higher
	- GCC 7 or compatible C++ compiler
	- Linux

# 2. Compilation
	1.	Clone the repository: git clone git@github.com:resetlla/q_psdm.git
    2.  Compile the code: make 
        The executable will be generated at bin/q-psdm.
    3.  Clean up compiled files: make clean

# 3. Usage Instructions
    go to example directory: bash run.sh

# 4. Notes and Precautions
	- Ensure your GPU driver and CUDA version are compatible.
	- FFT and PSDM core functions depend on fft.cc and psdmpkg.cc; they must be compiled together with q-psdm.cu.
	- When using custom SEG-Y data, verify that the file format is correct.
	- For large datasets, ensure sufficient GPU memory is available.
