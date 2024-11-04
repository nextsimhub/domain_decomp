/*!
 * @file MainMPI.cpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @author Joe Wallwork <jw2423@cam.ac.uk>
 * @date 04 Nov 2024
 */

#define DOCTEST_CONFIG_IMPLEMENT

#include <doctest/extensions/doctest_mpi.h>

// Global options
int global_argc;
char** global_argv;

int main(int argc, char** argv)
{
    doctest::mpi_init_thread(argc, argv, MPI_THREAD_MULTIPLE);

    // Assign globals, so they can be used within test cases
    global_argc = argc;
    global_argv = argv;

    doctest::Context ctx;
    ctx.setOption("reporters", "MpiConsoleReporter");
    ctx.setOption("reporters", "MpiFileReporter");
    ctx.setOption("force-colors", true);
    ctx.applyCommandLine(argc, argv);

    int test_result = ctx.run();

    doctest::mpi_finalize();

    return test_result;
}
