## RifDock
## Compilation Guide for Ubuntu

### Prerequisites

This guide assumes you are on an **Ubuntu**-based system. Before you begin, ensure you have the following dependencies:

* **CMake** & **Ninja-build**: Core build system tools.
* **GCC/G++ version 5**: RifDock requires an older compiler version.
* **Boost Library version 1.65.0**: A specific version of the Boost C++ libraries.
* **Rosetta**: The Rosetta software suite must be available and compiled.

### Installation

The installation is a multi-step process that involves setting up the environment, compiling dependencies, and finally, building RifDock itself.

#### 1. Install Build Tools & Compilers

First, install the necessary build utilities and the specific GCC-5 compiler required by RifDock.

```bash
# Update package list and install cmake/ninja
sudo apt update
sudo apt install -y cmake ninja-build

# Add the Ubuntu 16.04 (Xenial) repository to access older packages
echo "deb http://archive.ubuntu.com/ubuntu/ xenial main" | sudo tee /etc/apt/sources.list.d/xenial.list

# Update package list again and install gcc-5
sudo apt-get update
sudo apt install -y gcc-5 g++-5
```

> [!NOTE]
> RifDock was developed with older compilers. These commands add a repository from an older Ubuntu release (Xenial) specifically to allow `apt` to find and install the required `gcc-5` and `g++-5` packages on a modern system.

#### 2. Install Boost Library (v1.65.0)

Next, download and compile the required version of the Boost library.

```bash
# Download and extract Boost
wget -O boost_1_65_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.65.0/boost_1_65_0.tar.gz/download
tar -xzvf boost_1_65_0.tar.gz
cd boost_1_65_0/

# Bootstrap and build the library. Adjust -j20 based on your CPU cores.
./bootstrap.sh
./b2 -j20

# Return to the parent directory
cd ..
```

#### 3. Compile Rosetta

You must have a compiled version of the `rosetta_scripts` binary. Navigate to your Rosetta source directory to build it.

```bash
# Navigate to Rosetta's source directory
cd /path/to/your/rosetta_src_bundle/main/source

# Build rosetta_scripts (choose one)
# For a standard release build:
./ninja_build.py cxx11_omp -t rosetta_scripts -remake

# Or for a debug build:
./ninja_build.py cxx11_omp_debug -t rosetta_scripts -remake
```

#### 4. Compile RifDock

Finally, with all dependencies in place, you can compile RifDock.

1. **Clone the repository:**

```bash
git clone https://github.com/rifdock/rifdock
cd rifdock
mkdir build && cd build
```

2. **Set the compiler and run CMake:** Before running `cmake`, you must tell the system to use the `gcc-5` compiler you installed. Replace the placeholder paths with the **absolute paths** to your Boost and Rosetta directories.

```bash
# Set gcc-5 as the compiler for this terminal session
export CC=gcc-5 CXX=g++-5

# Run CMake with paths to dependencies
cmake .. \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_PREFIX_PATH=/path/to/your/boost_1_65_0/stage/lib \
-DCMAKE_ROSETTA_PATH=/path/to/your/rosetta_src_bundle/main \
-DCXXFLAGS="-isystem /path/to/your/boost_1_65_0"
```

3. **Compile the binaries:**

```bash
# Adjust -j20 based on your CPU cores
make -j20 rif_dock_test rifgen
```

> [!NOTE]
> The final executables, `rif_dock_test` and `rifgen`, will be located in your `rifdock/build/` directory. Congratulations, your installation is complete!
