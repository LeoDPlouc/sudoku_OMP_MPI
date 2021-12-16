#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#define BLOCK_LOW(id, np, n) ((id) * (n) / (np))
#define BLOCK_HIGH(id, np, n) (BLOCK_LOW((id) + 1, np, n) - 1)
#define BLOCK_SIZE(id, np, n) (BLOCK_LOW((id) + 1, np, n) - BLOCK_LOW(id, np, n))
#define BLOCK_OWNER(i, np, n) ((np * (i + 1) - 1) / n)

#define SUBGRID_N_INDEX(n, subgrid_index, subgrid_size, grid_size) (n % subgrid_size + (n / subgrid_size) * grid_size + (subgrid_index % subgrid_size) * subgrid_size + (subgrid_index / subgrid_size) * (grid_size * subgrid_size))
#define COLUMN_N_INDEX(n, column_index, grid_size) (column_index + (n * grid_size))
#define LINE_N_INDEX(n, line_index, grid_size) (n + (line_index * grid_size))

int parse_grid(int grid_size, int *grid)
{
    for (int i = 0; i < grid_size * grid_size; i++)
    {
        int n;
        scanf("%i", &n);
        grid[i] = n;
    }
}

int print_grid(int grid_size, int *grid)
{
    for (int i = 0; i < grid_size * grid_size; i++)
    {
        printf("%i", grid[i]);
        if ((i + 1) % grid_size == 0)
            printf("\n");
    }
    printf("\n");
}

int is_grid_completed(int grid_size, int *grid)
{
    int count = 0;
    for (int i = 0; i < grid_size * grid_size; i++)
    {
        if (!grid[i])
        {
            count++;
        }
    }
    return count;
}

int compute_only_subgrid(int subgrid_size, int grid_size, int *local_grid)
{
#pragma omp parallel for
    for (int id = 0; id < grid_size; id++)
    {
        int pos = -1;
        for (int n = 0; n < grid_size; n++)
        {
            if (local_grid[SUBGRID_N_INDEX(n, id, subgrid_size, grid_size)])
            {
                if (pos == -1)
                {
                    pos = SUBGRID_N_INDEX(n, id, subgrid_size, grid_size);
                }
                else
                {
                    pos = -1;
                    break;
                }
            }
        }
        if (pos != -1)
        {
            local_grid[pos] = 2;
        }
    }
}
int compute_only_column(int grid_size, int *local_grid)
{
#pragma omp parallel for
    for (int id = 0; id < grid_size; id++)
    {
        int pos = -1;
        for (int n = 0; n < grid_size; n++)
        {
            if (local_grid[COLUMN_N_INDEX(n, id, grid_size)])
            {
                if (pos == -1)
                {
                    pos = COLUMN_N_INDEX(n, id, grid_size);
                }
                else
                {
                    pos = -1;
                    break;
                }
            }
        }
        if (pos != -1)
        {
            local_grid[pos] = 2;
        }
    }
}
int compute_only_line(int grid_size, int *local_grid)
{
#pragma omp parallel for
    for (int id = 0; id < grid_size; id++)
    {
        int pos = -1;
        for (int n = 0; n < grid_size; n++)
        {
            if (local_grid[LINE_N_INDEX(n, id, grid_size)])
            {
                if (pos == -1)
                {
                    pos = LINE_N_INDEX(n, id, grid_size);
                }
                else
                {
                    pos = -1;
                    break;
                }
            }
        }
        if (pos != -1)
        {
            local_grid[pos] = 2;
        }
    }
}

int compute_occupied(int grid_size, int *local_grid, int *grid)
{
#pragma omp parallel for
    for (int i = 0; i < grid_size * grid_size; i++)
    {
        if (grid[i])
        {
            local_grid[i] = 0;
        }
    }
}

int compute_subgrid(int n, int subgrid_size, int grid_size, int *local_grid, int *grid)
{
#pragma omp parallel for
    for (int id = 0; id < grid_size; id++)
    {
        int occupied = 0;
        for (int i = 0; i < grid_size; i++)
        {
            if (grid[SUBGRID_N_INDEX(i, id, subgrid_size, grid_size)] == n + 1)
            {
                occupied = 1;
                break;
            }
        }
        if (occupied)
        {
            for (int i = 0; i < grid_size; i++)
            {
                local_grid[SUBGRID_N_INDEX(i, id, subgrid_size, grid_size)] = 0;
            }
        }
    }
}
int compute_column(int n, int grid_size, int *local_grid, int *grid)
{
#pragma omp parellel for
    for (int id = 0; id < grid_size; id++)
    {
        int occupied = 0;
        for (int i = 0; i < grid_size; i++)
        {
            if (grid[COLUMN_N_INDEX(i, id, grid_size)] == n + 1)
            {
                occupied = 1;
                break;
            }
        }
        if (occupied)
        {
            for (int i = 0; i < grid_size; i++)
            {
                local_grid[COLUMN_N_INDEX(i, id, grid_size)] = 0;
            }
        }
    }
}
int compute_line(int n, int grid_size, int *local_grid, int *grid)
{
#pragma omp parallel for
    for (int id = 0; id < grid_size; id++)
    {
        int occupied = 0;
        for (int i = 0; i < grid_size; i++)
        {
            if (grid[LINE_N_INDEX(i, id, grid_size)] == n + 1)
            {
                occupied = 1;
            }
        }
        if (occupied)
        {
            for (int i = 0; i < grid_size; i++)
            {
                local_grid[LINE_N_INDEX(i, id, grid_size)] = 0;
            }
        }
    }
}

int compute_results(int grid_size, int **local_grids, int *grid)
{
#pragma omp parellel for
    for (int n = 0; n < grid_size; n++)
    {
        for (int i = 0; i < grid_size * grid_size; i++)
        {
            if (local_grids[n][i] == 2)
            {
                grid[i] = n + 1;
            }
        }
    }
}
int compute_only(int grid_size, int **local_grids, int *grid)
{
#pragma omp parellel for
    for (int i = 0; i < grid_size * grid_size; i++)
    {
        int guessed = 0;
        for (int n = 0; n < grid_size; n++)
        {
            if (local_grids[n][i])
            {
                if (guessed)
                {
                    guessed = 0;
                    break;
                }
                else
                {
                    guessed = n + 1;
                }
            }
        }
        if (guessed)
        {
            grid[i] = guessed;
        }
    }
}

int main(int argc, char *argv[])
{
    int subgrid_size, grid_size, id, np;
    int *grid;
    int **local_grids;
    double elapsed_time;

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //parse grid size
    if (!id)
    {
        scanf("%i", &subgrid_size);
        grid_size = subgrid_size * subgrid_size;
    }

    MPI_Bcast(&subgrid_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&grid_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    grid = (int *)malloc(grid_size * grid_size * sizeof(int));
    //parse grid
    if (!id)
    {
        parse_grid(grid_size, grid);
        print_grid(grid_size, grid);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //Inits an array in wich the results of all the blocks will be stored
    local_grids = (int **)(malloc(grid_size * sizeof(int *)));

    //block 0 initialize all the local grids
    if (!id)
    {
#pragma omp parallel for
        for (int n = 0; n < grid_size; n++)
        {
            local_grids[n] = (int *)malloc(grid_size * grid_size * sizeof(int));
        }
    }

    //init local grids for the bloc to work on
#pragma omp parallel for
    for (int n = BLOCK_LOW(id, np, grid_size); n <= BLOCK_HIGH(id, np, grid_size); n++)
    {
        if (id)
        {
            local_grids[n] = (int *)malloc(grid_size * grid_size * sizeof(int));
        }
#pragma omp parallel for
        for (int i = 0; i < grid_size * grid_size; i++)
            local_grids[n][i] = 1;
    }

    while (1)
    //for (int a = 0; a < 1; a++)
    {
        //Block 0 broadcasts the grid
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(grid, grid_size * grid_size, MPI_INT, 0, MPI_COMM_WORLD);

        //for each number this block is responsible of
#pragma omp parallel for
        for (int n = BLOCK_LOW(id, np, grid_size); n <= BLOCK_HIGH(id, np, grid_size); n++)
        {
            //drop positions in the local grid if there's already a number in the global grid
            compute_occupied(grid_size, local_grids[n], grid);
            // if (n == 6)
            //     print_grid(grid_size, local_grids[6]);

            //drop position in the local grid for every subgrid where the number already existes in the global grid
            compute_subgrid(n, subgrid_size, grid_size, local_grids[n], grid);

            //drop position in the local grid for every column where the number already existes in the global grid
            compute_column(n, grid_size, local_grids[n], grid);

            //drop position in the local grid for every line where the number already existes in the global grid
            compute_line(n, grid_size, local_grids[n], grid);

            //mark the position if it's the only one with its number in the subgrid
            compute_only_subgrid(subgrid_size, grid_size, local_grids[n]);

            //mark the position if it's the only one with its number in the column
            compute_only_column(grid_size, local_grids[n]);

            //mark the position if it's the only one with its number in the line
            compute_only_line(grid_size, local_grids[n]);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        //every block will send its local grids to block 0
        if (id)
        {
            for (int n = BLOCK_LOW(id, np, grid_size); n <= BLOCK_HIGH(id, np, grid_size); n++)
            {
                MPI_Send(local_grids[n], grid_size * grid_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }
        if (!id)
        {
            for (int n = BLOCK_LOW(1, np, grid_size); n < grid_size; n++)
            {
                MPI_Status status;
                MPI_Recv(local_grids[n], grid_size * grid_size, MPI_INT, BLOCK_OWNER(n, np, grid_size), MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            }
        }

        if (!id)
        {
            //Block 0 looks through the results of all the blocks and put found numbers in the grid
            compute_results(grid_size, local_grids, grid);
            //print_grid(grid_size, grid);

            //Block 0 looks if only one local grid has a guess for every position in the grid and will consider it a found number
            compute_only(grid_size, local_grids, grid);

            int count = is_grid_completed(grid_size, grid);
            if (!count)
            {
                elapsed_time += MPI_Wtime();
                print_grid(grid_size, grid);
                printf("\nTotal elapsed time: %10.6f\n", elapsed_time);
                MPI_Finalize();
                return 0;
            }
            else
            {
                printf("%i position remaing\n", count);
                //print_grid(grid_size, grid);
            }
        }
    }
}