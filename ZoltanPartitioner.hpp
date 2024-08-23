/*!
 * @file ZoltanPartitioner.hpp
 * @author Athena Elafrou <ae488@cam.ac.uk>
 * @date 25 June 2022
 */

#pragma once

#include <memory>
#include <zoltan_cpp.h>

#include "Grid.hpp"
#include "Partitioner.hpp"

/*!
 * @class ZoltanPartitioner
 * @brief A class that encapsulates a 2D grid partitioner using the Zoltan
 * toolkit.
 *
 * A class that encapsulates a 2D grid partitioner using the Zoltan toolkit. We
 * use the Recursive Coordinate Bisection (RCB) geometric partitioning
 * algorithm.
 */
class ZoltanPartitioner final : public Partitioner {
public:
    // Disallow compiler-generated special functions
    ZoltanPartitioner(const ZoltanPartitioner&) = delete;
    ZoltanPartitioner& operator=(const ZoltanPartitioner&) = delete;

    /*!
     * @brief Destructor.
     */
    ~ZoltanPartitioner() { }

    /*!
     * @brief Create a Zoltan partitioner.
     *
     * @param comm MPI communicator.
     * @param argc The number of arguments.
     * @param argv The argument vector.
     * @return A ZoltanPartitioner object.
     */
    static ZoltanPartitioner* create(MPI_Comm comm, int argc, char** argv);

    /*!
     * @brief Partitions a 2D grid into rectangular boxes, one per process.
     *
     * Partitions a 2D grid into rectangular boxes, one per process, taking into
     * account a land mask, if provided.
     */
    void partition(Grid& grid) override;

protected:
    // Constructor
    ZoltanPartitioner(MPI_Comm comm, int argc, char** argv);

private:
    std::unique_ptr<Zoltan> _zoltan = nullptr; // Zoltan object
};
