Instructions for running Michael Wang's Simulations for MAE5730:

1) Go into MAE5730_Final_Project.m file

2) Run 1st block of code for Lagrange Method-derived pendulum simulation. Set animation variable to character 't' to run animation, Set it to 'f' to plot total energy of system over time. Edit variable n for desired number of links on the pendulum. Currently, the folder supports links of: 2, 3, 4 10, and 20. Run 2nd block of code to display data.

3) Run 3rd block of code for Newton-Euler Method-derived pendulum simulation. Set animation variable to character 't' to run animation, Set it to 'f' to plot total energy of system over time.

4) Run 4th block of code for DAE Method-derived pendulum simulation. Set animation variable to character 't' to run animation, Set it to 'f' to plot total energy of system over time.

5) Run 5th block of code for DAE Method-derived 4-bar-linkage simulation. Set animation variable to character 't' to run animation, Set it to 'f' to plot total energy of system over time.

6) Run 6th block to find periodic solutions to pendulum problem. Edit variable n for desired number of links on the pendulum. Currently, the folder supports links of: 2, 3, 4 10, and 20. Run 2nd block of code to display data.

Special:
 • createFunctionPendulumN_Lagrange(N): Function that is able to generate RHS file of N-linked Pendulum using Lagrange's Method. Run this command in command window. 