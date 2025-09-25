# 1. Project Overview, File Structure, System Requirements
## Overview
### Q-PSDM (**Q-compensated Prestack Depth Migration**) is a 3D seismic imaging tool designed to perform deabsorption prestack depth migration. This implementation is optimized with **CUDA** for GPU acceleration and supports **multithreading on CPU**, making it highly efficient for large-scale seismic datasets.  
## File Structure
### The repository is organized as follows:
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
- **src/** contains all CUDA and C++ source files.  
- **example/** includes example SEG-Y datasets and a simple script to validate installation.  
- **bin/** will contain the compiled executable (`q-psdm`) after building.  
- **Makefile** handles the compilation process and cleanup. 
## System Requirements
To successfully compile and run Q-PSDM, the following requirements should be met:

1. **CUDA Toolkit 7.5 or higher**  
   Required for GPU acceleration. Newer versions (CUDA 9.x, 10.x, 11.x) are generally supported but may require minor adjustments.  

2. **GCC 7 or compatible C++ compiler**  
   Ensures compatibility with the CUDA toolkit and modern C++ features used in the code.  

3. **Linux operating system**  
   The project is developed and tested on Linux. Other operating systems are not officially supported.  

4. **GPU hardware**  
   A CUDA-capable NVIDIA GPU with sufficient memory (at least 4 GB for small examples; 16 GB or more recommended for large datasets).  

# 2. Compilation
Follow the steps below to build the project:
- 1. Clone the repository: `git clone git@github.com:resetlla/q_psdm.git`
- 2. Compile the code: `make`
  - 1. The executable will be generated at bin/q-psdm.
- 3. Clean up compiled files: `make clean`

# 3. Usage Instructions
After compilation, you can test the program using the provided example:
- 1. `cd examples`
- 2. `bash run.sh`
# 4. Notes and Precautions
- 1. GPU and CUDA Compatibility:Ensure that your GPU driver and CUDA version are properly installed and compatible. Incorrect driver/toolkit versions may cause compilation or runtime errors.
- 2. Code Dependencies:Core migration functions depend on fft.cc and psdmpkg.cc. These files must always be compiled together with q-psdm.cu.
- 3. Data Format:Input seismic data must be in SEG-Y format. When using custom SEG-Y data, verify that the file headers and traces conform to the expected standard.
- 4. Memory Usage:For large-scale datasets, make sure your GPU has enough memory. If memory is insufficient, consider dividing the dataset into smaller chunks.
- 5. Performance Tuning:
     1. Adjust thread block and grid sizes in q-psdm.cu for optimal GPU utilization.
     2. Consider compiling with optimization flags (-O2 or -O3) for faster runtime.
- 6. Support and Extensions:This project is designed to be modular. Advanced users can extend the implementation by:
     1. Adding new migration algorithms.
     2. Integrating with other geophysical processing tools.
     3. Modifying the FFT routines for custom applications.
