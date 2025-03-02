// *** Notes:
// Running the code:
// mcgleec2@cheetah:~$ gcc GP1/wire_C/serial.c -o bin/GP1/wire_C/serial -lm   //-lm is needed at end
// mcgleec2@cheetah:~$ bin/gp1/work/serial
// ***


// Serial C Program to simulate electrical conductivity in a 2D wire

// Headers
#include <stdio.h>  // Standard input/output functions -> file & printf()
#include <stdlib.h> // Standard library -> random number generator [rand()]
#include <math.h>   // Math library
#include <time.h>   // Library -> seed random number generator = different random values every run

// Define constants
// Physical parameters of electrons -> Charge, Mass
// Grid partitioning -> Grid Size, Cell Size
// Simulation constraints -> Number of Electrons, Number of Time Steps, Time Step Size
#define NUM_ELECTRONS 10    // Number of electrons in the simulation
#define TIME_STEPS 10       // Number of time steps for the simulation
#define DT 1e-15            // Time step size in seconds
#define L 1e-6              // Length of the simulation area (1 um)
#define GRID_SIZE 10        // Number of divisions per dimension in the spatial partition grid
#define CELL_SIZE (L / GRID_SIZE) // Size of each grid cell
#define E_CHARGE 1.602e-19  // Elementary charge of an electron (Coulombs [C])
#define E_MASS 9.109e-31    // Mass of an electron (kg)
#define EPSILON_0 8.854e-12 // Vacuum permittivity (Farads per meter [F/m])

// Electron Structure
// Each electron has an x & y position component -> 2D space.
// Each electron has a vx & vy velocity component -> 2D space.
typedef struct {
    double x, y;   // Position coordinates (meters)
    double vx, vy; // Velocity components (meters/second)
} Electron;

// Data Structure for Partitioning
// First, store all electrons.
// Second, store pointers to electrons in each grid cell.
// Lastly, store the number of electrons in each grid cell.
// Checking nearby cells and not all electrons -> Computational: O(N^2) —> O(N).
Electron electrons[NUM_ELECTRONS]; // Array to store all electrons
Electron *grid[GRID_SIZE][GRID_SIZE][NUM_ELECTRONS]; // 3D array for grid-based spatial partitioning
int grid_count[GRID_SIZE][GRID_SIZE]; // Array to store the count of electrons in each grid cell

// Function to initialise electron positions and velocities randomly
// Randomly place electrons into the simulation.
// Then assigns random positions and velocities.
void initialize_electrons() {
    srand(time(NULL)); // Seed the random number generator with the current time
    for (int i = 0; i < NUM_ELECTRONS; i++) {
        electrons[i].x = ((double)rand() / RAND_MAX) * L;  // Assign random x position within system bounds
        electrons[i].y = ((double)rand() / RAND_MAX) * L;  // Assign random y position within system bounds
        electrons[i].vx = ((double)rand() / RAND_MAX - 0.5) * 1e5; // Assign random x velocity
        electrons[i].vy = ((double)rand() / RAND_MAX - 0.5) * 1e5; // Assign random y velocity
    }
}

// Function to reset and update the grid for spatial partitioning
// Resets the grid.
// Then assigns each electron to a grid cell -> Check neighbouring grid cells not all electrons
void update_grid() {
    // Reset all grid counts to zero
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            grid_count[i][j] = 0;
        }
    }

    // Assign electrons to grid cells based on their positions
    for (int i = 0; i < NUM_ELECTRONS; i++) {
        int gx = (int)(electrons[i].x / CELL_SIZE); // Determine grid cell index in x direction
        int gy = (int)(electrons[i].y / CELL_SIZE); // Determine grid cell index in y direction
        if (gx >= 0 && gx < GRID_SIZE && gy >= 0 && gy < GRID_SIZE) { // Ensure valid grid range
            grid[gx][gy][grid_count[gx][gy]++] = &electrons[i]; // Store electron pointer in the grid cell
        }
    }
}

// Function to compute forces acting on electrons using grid-based partitioning
// Use Coulomb’s law to compute forces
// Avoid self-interaction
// Checks neighbouring cells only
void compute_forces(double ax[], double ay[]) {
    for (int i = 0; i < NUM_ELECTRONS; i++) {
        ax[i] = 0.0; // Initialise acceleration in x direction
        ay[i] = 0.0; // Initialise acceleration in y direction

        int gx = (int)(electrons[i].x / CELL_SIZE); // Get electron's grid cell in x
        int gy = (int)(electrons[i].y / CELL_SIZE); // Get electron's grid cell in y

        // Iterate through the neighboring grid cells to calculate forces
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                int nx = gx + dx; // Neighboring cell in x direction
                int ny = gy + dy; // Neighboring cell in y direction
                if (nx >= 0 && nx < GRID_SIZE && ny >= 0 && ny < GRID_SIZE) { // Ensure valid grid cell
                    for (int j = 0; j < grid_count[nx][ny]; j++) {
                        Electron *neighbor = grid[nx][ny][j]; // Access neighboring electron
                        if (neighbor != &electrons[i]) { // Avoid self-interaction
                            double dx = neighbor->x - electrons[i].x;
                            double dy = neighbor->y - electrons[i].y;
                            double r2 = dx * dx + dy * dy + 1e-18; // Add small value to prevent division by zero
                            double r = sqrt(r2);
                            double force = (E_CHARGE * E_CHARGE) / (4 * M_PI * EPSILON_0 * r2); // Coulomb's law
                            ax[i] += force * dx / (r * E_MASS); // Compute acceleration in x
                            ay[i] += force * dy / (r * E_MASS); // Compute acceleration in y
                        }
                    }
                }
            }
        }
    }
}

// Function to update electron positions and velocities
void update_electrons(double ax[], double ay[]) {
    for (int i = 0; i < NUM_ELECTRONS; i++) {
        electrons[i].vx += ax[i] * DT; // Update velocity in x direction
        electrons[i].vy += ay[i] * DT; // Update velocity in y direction
        electrons[i].x += electrons[i].vx * DT; // Update position in x direction
        electrons[i].y += electrons[i].vy * DT; // Update position in y direction

        // Apply periodic boundary conditions to maintain simulation
        // Electrons wrap around to the start.
        if (electrons[i].x > L) electrons[i].x -= L;
        if (electrons[i].x < 0) electrons[i].x += L;
        if (electrons[i].y > L) electrons[i].y -= L;
        if (electrons[i].y < 0) electrons[i].y += L;
    }
}

// Function to run the full simulation
void run_simulation() {
    FILE *file = fopen("electron_data.csv", "w"); // Open file for output
    fprintf(file, "Time, Electron, X, Y, VX, VY\n"); // Write CSV header

    double ax[NUM_ELECTRONS], ay[NUM_ELECTRONS]; // Arrays for electron accelerations
    for (int t = 0; t < TIME_STEPS; t++) {
        update_grid(); // Update spatial partitioning grid
        compute_forces(ax, ay); // Compute forces using grid partitioning
        update_electrons(ax, ay); // Update electron positions and velocities

        // Log data for each electron in the CSV file
        for (int i = 0; i < NUM_ELECTRONS; i++) {
            fprintf(file, "%d, %d, %.6e, %.6e, %.6e, %.6e\n", t, i, electrons[i].x, electrons[i].y, electrons[i].vx, electrons[i].vy);
        }
    }
    fclose(file); // Close file after writing data
}

// Main function
int main() {
    struct timespec start_time, end_time, time_diff;
    double runtime = 0.0;
    timespec_get(&start_time, TIME_UTC); // Capture start time

    initialize_electrons(); // Initialise electron positions and velocities
    run_simulation(); // Run the simulation

    timespec_get(&end_time, TIME_UTC); // Capture end time

    time_diff.tv_sec = end_time.tv_sec - start_time.tv_sec;
    time_diff.tv_nsec = end_time.tv_nsec - start_time.tv_nsec;
    if (time_diff.tv_nsec < 0) {
        time_diff.tv_sec -= 1;
        time_diff.tv_nsec += 1000000000;
    }
    runtime = time_diff.tv_sec + time_diff.tv_nsec / 1e9;

    printf("Simulation complete. Data saved to electron_data.csv.\n"); // Check if program has ran
    printf("Runtime: %lf seconds\n", runtime);
    return 0;
}
