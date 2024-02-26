# AMES-Fit
Automatic Multiple Experiment Simulation and Fitting

**Windows**
To use the program in Windows, you can simply download the binary and paste it in the directory with the input and spectral files.
Using a Ccmd (command prompt) window, cd into this directory and run the program with the command:
AMES-Fit.exe <input filename>

**Unix**
No biraries are provided for a UNIX-like environment (Mac or Linux) but these can be created by compiling the program with the command:
g++ AMES-Fit.cpp -o AMES-Fit -fopenmp -lm -O3 -march=native

And then run in the terminal (or schedular, such as SLURM) with the command:
./AMES-Fit <input filename>
