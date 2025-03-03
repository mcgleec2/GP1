serial.c:

- serial.c is a C program which is supposed to simulate electrical conductivity for a 2D wire.
- currently only set to 10 electrons and time step of 10.
- Need to use -lm at the end of gcc ...... -lm compiling.
- 


plot.py:

- This code was run using Colab Notebooks and not on the Command Prompy.
- Had to copy the output data from serial.c manually and create own txt file with the data and 
headers. Then upload that txt file to the same Google Drive folder that contains the python program.
- python code to plot the x vs y positions of electrons at time steps and also 
to plot the average velocty vs the time steps.
- Average velcoity vs time step is just a linearly decreasing line ???.
- 


parallel.c:

- Parallelised version of serial.c
- Processes handle own sections and pass electrons along with messages.
- Last process at end of boundary wraps bacj to the root for electrons to pass through again.


trouble_python_code.py
- code to run simulation on python code
- hard understanding of the final graphs
- due to lack of time to fix the errors and understand the plots we decided to focus on the c code


Note:
- Original repo was set to be private, therefore this new one had to be create to let only view acces.
