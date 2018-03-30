// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "OutputMaterial.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <typeinfo> // USES typeid()

// ----------------------------------------------------------------------
const char* pylith::meshio::OutputMaterial::_pyreComponent = "outputmaterial";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputMaterial::OutputMaterial(void)
{ // constructor
    PyreComponent::name(_pyreComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputMaterial::~OutputMaterial(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputMaterial::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    OutputManager::deallocate();

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputMaterial::verifyConfiguration(const pylith::topology::Field& solution,
                                                    const pylith::topology::Field& auxField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::verifyConfiguration(auxField="<<auxField.label()<<")");

    OutputManager::verifyConfiguration(solution, auxField);

    PYLITH_COMPONENT_ERROR("@brad :TODO: Implement verifyConfiguration().");

    PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Write solution at time step.
void
pylith::meshio::OutputMaterial::writeTimeStep(const PylithReal t,
                                              const PylithInt tindex,
                                              const pylith::topology::Field& solution,
                                              const pylith::topology::Field& auxField) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputMaterial::writeTimeStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<", auxField="<<auxField.label()<<")");

    if (!this->shouldWrite(t, tindex)) {
        PYLITH_METHOD_END;
    } // if

    const pylith::string_vector& subfieldNames = (1 == _vertexDataFields.size() && std::string("all") == _vertexDataFields[0]) ? solution.subfieldNames() : _vertexDataFields;

#if 0
    this->openTimeStep(t, solution.mesh());
    const size_t numFields = subfieldNames.size();
    for (size_t iField = 0; iField < numFields; iField++) {
        // :TODO: STUFF GOES HERE

        pylith::topology::Field& fieldBuffer = this->getBuffer(solution, subfieldNames[iField].c_str());
        this->appendVertexField(t, fieldBuffer, fieldBuffer.mesh());
    } // for
    this->closeTimeStep();
#else
    PYLITH_COMPONENT_ERROR("@brad :TODO: Implemenet writeTimeStep().");
#endif

    PYLITH_METHOD_END;
} // writeTimeStep


// End of file
