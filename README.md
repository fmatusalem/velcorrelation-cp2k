
BY FILIPE MATUSALEM, DEC 2020     filipematus@gmail.com 
Program to compute normalized velocity correlation function from CP2K xyz velocity file

Compilation: g++ -o velcorrelation-cp2k.x velcorrelation-cp2k-v2.c

Usage: ./velcorrelation-cp2k.x cp2k-velocity-file.xyz

The maximum time lag can be entered as a second argument (Default is half number of MD steps on cp2k-velocity-file.xyz file)
