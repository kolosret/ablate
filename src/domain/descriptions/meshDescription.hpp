#ifndef ABLATELIBRARY_MESHDESCRIPTION_HPP
#define ABLATELIBRARY_MESHDESCRIPTION_HPP

#include <petsc.h>
#include <domain/region.hpp>
#include <set>

namespace ablate::domain::descriptions {
/**
 * A simple interface that is responsible for determine cell/vertex locations in the mesh
 */
class MeshDescription {
   public:
    virtual ~MeshDescription() = default;
    /**
     * The overall assumed dimension of the mesh
     * @return
     */
    [[nodiscard]] virtual const PetscInt& GetMeshDimension() const = 0;

    /**
     * Total number of cells in the entire mesh
     * @return
     */
    [[nodiscard]] virtual const PetscInt& GetNumberCells() const = 0;

    /**
     * Total number of nodes/vertices in the entire mesh
     * @return
     */
    [[nodiscard]] virtual const PetscInt& GetNumberVertices() const = 0;

    /**
     * Return the cell type for this cell, based off zero offset
     * @return
     */
    [[nodiscard]] virtual DMPolytopeType GetCellType(PetscInt cell) const = 0;

    /**
     * Builds the topology based upon a zero vertex offset.  The cellNodes should be sized to hold cell
     * @return
     */
    virtual void BuildTopology(PetscInt cell, PetscInt* cellNodes) const = 0;

    /**
     * Builds the node coordinate for each vertex
     * @return
     */
    virtual void SetCoordinate(PetscInt node, PetscReal* coordinate) const = 0;

    /**
     * returns the boundary region to label for the entire region
     */
    [[nodiscard]] virtual std::shared_ptr<ablate::domain::Region> GetBoundaryRegion() const = 0;

    /**
     * returns a region if the nodes set (on a face) belongs to a certain region
     */
    [[nodiscard]] virtual std::shared_ptr<ablate::domain::Region> GetRegion(const std::set<PetscInt>& face) const = 0;
};
}  // namespace ablate::domain::descriptions

#endif  // ABLATELIBRARY_MESHDESCRIPTION_HPP
