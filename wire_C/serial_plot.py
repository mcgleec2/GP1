# *** Note ran on Colab.
# Data file manually extracted and made into txt file
# then uploaded int drive folder.

# Average velocity vs time step is linearly decreasing, needs fixed -> Something in C program data?


# ************************************************************************************

## Text:
# Mount Google Drive. Create working directory and image directory.

## Code:
# Import Google Drive from Google Colab library and import operating system
from google.colab import drive
import os

# Mount Google Drive
drive.mount('/content/drive')

# Create directory of working files
work_dir = "/content/drive/My Drive/Colab Notebooks"

# Create directory of image data
image_dir = os.path.join(work_dir, 'images')



## Text:
# List the files in the working directory

## Code:
!pwd
%cd '/content/drive/My Drive/Colab Notebooks'
!ls



## Text:
# Program  to read & plot the data

## Code:
# Import relevant libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Function to read the data file and return it as a Pandas DataFrame
def read_data(filename):
    return pd.read_csv(filename, delimiter=', ') # The deimiter used is a comma followed by a space (', ').


# Function to plot x vs y positions of electrons
def plot_positions(data):

    # Set figure dimensions
    plt.figure(figsize=(10, 6))

    # Use a scatter plot to plot x vs y data.
    plt.scatter(data['X'], data['Y'], s=1) # Set size of the points (s=1).

    # Set the x and y labels.
    plt.xlabel("X Position (m)")
    plt.ylabel("Y Position (m)")

    # Set title.
    plt.title("Electron Positions")

    # Add grid
    plt.grid()

    # Display plot.
    plt.show()


# Function to calculate velocities, then plot the average velocitiy vs time step
def plot_avg_velocity_vs_time(data):

    # Calculate velocity of each electron at each time step (Pythagoreus: sqrt(VX^2 + VY^2)
    data['Velocity'] = np.sqrt(data['VX']**2 + data['VY']**2)

    # Group the data by 'Time' and compute the average velocity for each time step.
    # Group data by 'Time', then find the average velocity at each time step.
    avg_velocity = data.groupby('Time')['Velocity'].mean()

    # Set figure dimensions.
    plt.figure(figsize=(8, 5))

    # Plot average velocity vs time step.
    plt.plot(avg_velocity.index, avg_velocity.values, marker='o', linestyle='-') # marker='o' for circle points & linestyle='-' for dashed connection

    # Set the x and y labels.
    plt.xlabel("Time Step")
    plt.ylabel("Average Velocity (m/s)")

    # Set title.
    plt.title("Average Velocity vs Time Step")

    # Add grid
    plt.grid()

    # Display plot.
    plt.show()


# Function to run the program: Read data & create plots.
def main():
    # Set the filename the data file.
    filename = "electron_data1.txt"

    # Read the data using read_data function.
    data = read_data(filename)

    # Plot the electron positions using plot_positions function.
    plot_positions(data)

    # Plot the average velocity vs time using plot_avg_velocity_vs_time function.
    plot_avg_velocity_vs_time(data)


# Call the main() function
main()
