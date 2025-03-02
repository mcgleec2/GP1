#include <stdio.h>   // Standard input/output functions
#include <stdlib.h>  // Standard library -> memory allocation, random numbers.
#include <math.h>    // Math library -> sqrt
#include <time.h>    // Time library ->seeding random number gen.
#include <mpi.h>     // MPI library

// Define constants
#define NUM_ELECTRONS 10  // Number of electrons in the simulation
#define TIME_STEPS 10     // Number of time steps for the simulation
#define DT 1e-15          // Time step size in seconds
#define L 1e-6            // Length of the simulation area (1 micrometer)
#define GRID_SIZE 10      // Number of divisions per dimension in the spatial partition grid
#define CELL_SIZE (L / GRID_SIZE) // Size of each grid cell
#define E_CHARGE 1.602e-19 // Elementary charge of an electron (Coulombs [C])
#define E_MASS 9.109e-31  // Mass of an electron (kg)
#define EPSILON_0 8.854e-12 // Vacuum permittivity (Farads per meter[F/m])

// Electron structure
// Each electron has an x & y position component -> 2D space.
// Each electron has a vx & vy velocity component -> 2D space.
typedef struct {
    double x, y;   // Position coordinates
    double vx, vy; // Velocity components
} Electron;

// Array to store electrons for the simulation
Electron electrons[NUM_ELECTRONS];

// Function to initialise electrons with random positions and velocities
// Randomly place electrons into the simulation.
// Then assigns random positions and velocities.
void initialize_electrons() {
    srand(time(NULL)); // Seed the random number generator with the current time
    for (int i = 0; i < NUM_ELECTRONS; i++) {
        electrons[i].x = ((double)rand() / RAND_MAX) * L;  // Assign random x position
        electrons[i].y = ((double)rand() / RAND_MAX) * L;  // Assign random y position
        electrons[i].vx = ((double)rand() / RAND_MAX - 0.5) * 1e5; // Random x velocity
        electrons[i].vy = ((double)rand() / RAND_MAX - 0.5) * 1e5; // Random y velocity
    }
}

// Function to compute forces acting on electrons using Coulomb’s law
void compute_forces(Electron *electrons, double ax[], double ay[]) {
    for (int i = 0; i < NUM_ELECTRONS; i++) {
        ax[i] = 0.0; // Initialise acceleration in x direction
        ay[i] = 0.0; // Initialise acceleration in y direction
        for (int j = 0; j < NUM_ELECTRONS; j++) {
            if (i != j) { // Avoid self-interaction
                double dx = electrons[j].x - electrons[i].x;
                double dy = electrons[j].y - electrons[i].y;
                double r2 = dx * dx + dy * dy + 1e-18; // Prevent division by zero
                double r = sqrt(r2);
                double force = (E_CHARGE * E_CHARGE) / (4 * M_PI * EPSILON_0 * r2); // Coulomb's law
                ax[i] += force * dx / (r * E_MASS); // Compute acceleration in x
                ay[i] += force * dy / (r * E_MASS); // Compute acceleration in y
            }
        }
    }
}

// Function to update electron positions and velocities
void update_electrons(Electron *electrons, double ax[], double ay[]) {
    for (int i = 0; i < NUM_ELECTRONS; i++) {
        electrons[i].vx += ax[i] * DT; // Update velocity in x direction
        electrons[i].vy += ay[i] * DT; // Update velocity in y direction
        electrons[i].x += electrons[i].vx * DT; // Update position in x direction
        electrons[i].y += electrons[i].vy * DT; // Update position in y direction

        // Apply periodic boundary conditions (electrons wrap around)
        if (electrons[i].x > L) electrons[i].x -= L;
        if (electrons[i].x < 0) electrons[i].x += L;
        if (electrons[i].y > L) electrons[i].y -= L;
        if (electrons[i].y < 0) electrons[i].y += L;
    }
}

// Convert timespec to seconds as a float
double to_second_float(struct timespec in_time) {
    double out_time = 0.0;
    long int seconds = in_time.tv_sec;
    long int nanoseconds = in_time.tv_nsec;

    out_time = seconds + nanoseconds / 1e9; // Convert to seconds
    return out_time;
}

// Calculate the difference between start and end times
struct timespec calculate_runtime(struct timespec start_time, struct timespec end_time) {
    struct timespec time_diff;
    long int seconds = end_time.tv_sec - start_time.tv_sec;
    long int nanoseconds = end_time.tv_nsec - start_time.tv_nsec;

    if (nanoseconds < 0) {
        seconds -= 1;
        nanoseconds += 1000000000; // Carry the 1
    }

    time_diff.tv_sec = seconds;
    time_diff.tv_nsec = nanoseconds;

    return time_diff;
}

int main(int argc, char **argv) {
    int rank, size; // MPI process rank and total number of processes
    MPI_Init(&argc, &argv); // Initialise MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get current process rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get total number of processes

    double ax[NUM_ELECTRONS], ay[NUM_ELECTRONS]; // Arrays for electron accelerations

    struct timespec start_time, end_time, time_diff; // Time tracking variables
    double runtime = 0.0;

    if (rank == 0) {
        // Start the time tracking for initialization
        timespec_get(&start_time, TIME_UTC); // Capture start time for initialization

        // Root process initialises electrons
        initialize_electrons();

        // End time tracking for initialization
        timespec_get(&end_time, TIME_UTC); // Capture end time for initialization
        time_diff = calculate_runtime(start_time, end_time);
        runtime = to_second_float(time_diff);
        printf("Initialization complete. Time taken: %lf seconds.\n", runtime);
    }

    // Broadcast electrons to all processes
    MPI_Bcast(electrons, NUM_ELECTRONS * sizeof(Electron), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Each process handles a portion of the time steps
    for (int t = rank; t < TIME_STEPS; t += size) {
        // Start time tracking for force computation
        timespec_get(&start_time, TIME_UTC); // Capture start time for force computation

        compute_forces(electrons, ax, ay); // Compute forces on electrons

        // End time tracking for force computation
        timespec_get(&end_time, TIME_UTC); // Capture end time for force computation
        time_diff = calculate_runtime(start_time, end_time);
        runtime = to_second_float(time_diff);
        printf("Time step %d: Force computation time: %lf seconds.\n", t, runtime);

        // Start time tracking for electron update
        timespec_get(&start_time, TIME_UTC); // Capture start time for electron update

        update_electrons(electrons, ax, ay); // Update electron positions and velocities

        // End time tracking for electron update
        timespec_get(&end_time, TIME_UTC); // Capture end time for electron update
        time_diff = calculate_runtime(start_time, end_time);
        runtime = to_second_float(time_diff);
        printf("Time step %d: Electron update time: %lf seconds.\n", t, runtime);
    }

    // Gather updated electron data back to root process
    MPI_Gather(electrons, NUM_ELECTRONS * sizeof(Electron), MPI_BYTE,
               electrons, NUM_ELECTRONS * sizeof(Electron), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (rank == 0) { // Root process writes output to file
        FILE *file = fopen("electron_data.csv", "w"); // Open file
        fprintf(file, "Time, Electron, X, Y, VX, VY\n"); // Write header
        for (int i = 0; i < NUM_ELECTRONS; i++) {
            fprintf(file, "%d, %d, %.6e, %.6e, %.6e, %.6e\n", 0, i, electrons[i].x, electrons[i].y, electrons[i].vx, electrons[i].vy);
        }
        fclose(file); // Close file
        printf("Simulation complete. Data saved to electron_data.csv.\n");
    }

    MPI_Finalize(); // Finalize MPI
    return 0;
}

