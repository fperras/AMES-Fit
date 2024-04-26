# AMES-Fit
Automatic Multiple Experiment Simulation and Fitting

**Windows**
To use the program in Windows, you can simply download the binary and paste it in the directory with the input and spectral files.
Using a Ccmd (command prompt) window, cd into this directory and run the program with the command:

AMES-Fit.exe input_filename

**Unix**
No binaries are provided for a UNIX-like environment (Mac or Linux) but these can be created by compiling the program with the command:
g++ AMES-Fit.cpp -o AMES-Fit -fopenmp -lm -O3 -march=native

And then run in the terminal (or schedular, such as SLURM) with the command:

./AMES-Fit input_filename

Copyright 2024. Iowa State University. This material was produced under U.S. Government contract DE-AC02-07CH11358 for the Ames National Laboratory, which is operated by Iowa State University for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software. NEITHER THE GOVERNMENT, AMES NATIONAL LABORATORY, NOR IOWA STATE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE. If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from Ames National Laboratory.

Additionally, this program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. Accordingly, this program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. (https://www.gnu.org/licenses/gpl-3.0.en.html#license-text)
