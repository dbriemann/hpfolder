hpfolder
========

A protein folding genetic algorithm. Works on a simplified HP 2D model.


USAGE: 

   ./hpfolder  [-e <pos. int>] [-y <neg. int>] -p <string> [-s <pos. int>]
               [-m <float 0..1>] [-c <float 0..1>] [-g] [--] [--version]
               [-h]


Where: 

   -e <pos. int>,  --maxeval <pos. int>
     Sets the maximum number of energy evaluations

   -y <neg. int>,  --minenergy <neg. int>
     Sets the minimum energy

   -p <string>,  --protein <string>
     (required)  Sets the protein

   -s <pos. int>,  --popsize <pos. int>
     Sets the population size

   -m <float 0..1>,  --mutateprob <float 0..1>
     Sets the mutation probability

   -c <float 0..1>,  --crossprob <float 0..1>
     Sets the crossover probability

   -g,  --graphics
     Enables OpenGL Window

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Folding Visualizer - HP 2d Model



Example Proteins
----------------
1. BWBWWBBWBWWBWBBWWBWB -9
2. WWWBBWWBBWWWWWBBBBBBBWWBBWWWWBBWWBWW -14
3. WWBWWBBWWBBWWWWWBBBBBBBBBBWWWWWWBBWWBBWWBWWBBBBB -23
4. BBBBBBBBBBBBWBWBWWBBWWBBWWBWWBBWWBBWWBWWBBWWBBWWBWBWBBBBBBBBBBBB -40
5. BBBBWWWWBBBBBBBBBBBBWWWWWWBBBBBBBBBBBBWWWBBBBBBBBBBBBWWWBBBBBBBBBBBBWWWBWWBBWWBBWWBWB -52

