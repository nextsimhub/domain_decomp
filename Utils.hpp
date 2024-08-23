/*!
 * @file Utils.hpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 23 Aug 2024
 */

#pragma once

#include <mpi.h>

#define CHECK_MPI(func)                                                                            \
    {                                                                                              \
        int mpi_status = (func);                                                                   \
        if (mpi_status != MPI_SUCCESS) {                                                           \
            char mpi_error_string[MPI_MAX_ERROR_STRING];                                           \
            int mpi_error_string_length = 0;                                                       \
            MPI_Error_string(mpi_status, mpi_error_string, &mpi_error_string_length);              \
            if (mpi_error_string != NULL)                                                          \
                fprintf(stderr, "ERROR: MPI call \"%s\" at line %d of file %s with %s (%d)\n",     \
                    #func, __LINE__, __FILE__, mpi_error_string, mpi_status);                      \
            else                                                                                   \
                fprintf(stderr, "ERROR: MPI call \"%s\" at line %d of file %s failed with %d.\n",  \
                    #func, __LINE__, __FILE__, mpi_status);                                        \
            exit(mpi_status);                                                                      \
        }                                                                                          \
    }
