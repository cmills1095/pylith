// -*- C++ -*-
//
// ---------------------------------------------------------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ---------------------------------------------------------------------------------------------------------------------
//

#include <portinfo>

#include "IntegratorBoundary.hh" // implementation of object methods

#include "pylith/feassemble/DSLabelAccess.hh" // USES DSLabelAccess
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createLowerDimMesh()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/problems/Physics.hh" // USES Physics
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor::optimizeClosure()

#include "spatialdata/spatialdb/GravityField.hh" // HASA GravityField

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "petscds.h" // USES PetscDS

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ---------------------------------------------------------------------------------------------------------------------
// Local "private" functions.
namespace pylith {
    namespace feassemble {
        class _IntegratorBoundary {
public:

            static const char* genericComponent;
        }; // _IntegratorBoundary
        const char* _IntegratorBoundary::genericComponent = "integratorboundary";

    } // feassemble
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorBoundary::IntegratorBoundary(pylith::problems::Physics* const physics) :
    Integrator(physics),
    _boundaryMesh(NULL),
    _boundarySurfaceLabel("") {
    GenericComponent::setName(_IntegratorBoundary::genericComponent);
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::IntegratorBoundary::~IntegratorBoundary(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorBoundary::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    Integrator::deallocate();

    delete _boundaryMesh;_boundaryMesh = NULL;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set label marking boundary associated with boundary condition surface.
void
pylith::feassemble::IntegratorBoundary::setMarkerLabel(const char* value) {
    PYLITH_JOURNAL_DEBUG("setMarkerLabel(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for boundary condition integrator label.");
    } // if

    _boundarySurfaceLabel = value;
} // setMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get label marking boundary associated with boundary condition surface.
const char*
pylith::feassemble::IntegratorBoundary::getMarkerLabel(void) const {
    return _boundarySurfaceLabel.c_str();
} // getMarkerLabel


// ---------------------------------------------------------------------------------------------------------------------
// Get mesh associated with integration domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorBoundary::getPhysicsDomainMesh(void) const {
    assert(_boundaryMesh);
    return *_boundaryMesh;
} // domainMesh


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorBoundary::setKernelsResidual(const std::vector<ResidualKernels>& kernels,
                                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernelsResidual(# kernels="<<kernels.size()<<")");

    PetscErrorCode err;
    DSLabelAccess dsLabel(solution.dmMesh(), _labelName.c_str(), _labelValue);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_field = solution.subfieldInfo(kernels[i].subfield.c_str()).index;
        const PetscInt i_part = kernels[i].part;
        err = PetscWeakFormAddBdResidual(dsLabel.weakForm(), dsLabel.label(), dsLabel.value(), i_field, i_part,
                                         kernels[i].r0, kernels[i].r1);PYLITH_CHECK_ERROR(err);

        switch (kernels[i].part) {
        case RESIDUAL_LHS:
            _hasLHSResidual = true;
            break;
        case RESIDUAL_RHS:
            _hasRHSResidual = true;
            break;
        default:
            PYLITH_JOURNAL_LOGICERROR("Unknown residual part " << kernels[i].part <<".");
        } // switch
    } // for

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        err = PetscDSView(dsLabel.ds(), PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // setKernelsResidual


// ---------------------------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::IntegratorBoundary::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("intialize(solution="<<solution.getLabel()<<")");

    delete _boundaryMesh;
    _boundaryMesh = pylith::topology::MeshOps::createLowerDimMesh(solution.mesh(), _boundarySurfaceLabel.c_str());
    assert(_boundaryMesh);
    pylith::topology::CoordsVisitor::optimizeClosure(_boundaryMesh->dmMesh());

    Integrator::initialize(solution);

    assert(_auxiliaryField);
    PetscErrorCode err;
    PetscDM dmSoln = solution.dmMesh();assert(dmSoln);
    PetscDMLabel dmLabel = NULL;
    err = DMGetLabel(dmSoln, _labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, _auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing auxiliary field.");
        _auxiliaryField->view("Auxiliary field");
    } // if

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Update auxiliary field values to current time.
void
pylith::feassemble::IntegratorBoundary::updateState(const double t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("updateState(t="<<t<<")");

    Integrator::updateState(t);

    assert(_physics);
    _physics->updateAuxiliaryField(_auxiliaryField, t);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        assert(_auxiliaryField);
        PYLITH_JOURNAL_DEBUG("IntegratorInterface component '" << GenericComponent::getName() << "' for '"
                                                               <<_physics->getIdentifier()
                                                               << "': viewing auxiliary field.");
        _auxiliaryField->view("IntegratorInterface auxiliary field", pylith::topology::Field::VIEW_ALL);
    } // if

    PYLITH_METHOD_END;
} // updateState


// ---------------------------------------------------------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorBoundary::computeRHSResidual(pylith::topology::Field* residual,
                                                           const PylithReal t,
                                                           const PylithReal dt,
                                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");
    if (!_hasRHSResidual) { PYLITH_METHOD_END;}
    assert(residual);

    DSLabelAccess dsLabel(solution.dmMesh(), _labelName.c_str(), _labelValue);
    _setKernelConstants(solution, dt);

    PetscFormKey key;
    key.label = dsLabel.label();
    key.value = dsLabel.value();
    key.part = pylith::feassemble::Integrator::RESIDUAL_RHS;

    PetscErrorCode err;
    assert(solution.localVector());
    assert(residual->localVector());
    PetscVec solutionDotVec = NULL;
#if 0
    // I think this is the interface we want.
    err = DMPlexComputeResidual_XXX(dsLabel.dm(), key, solution.localVector(), solutionDotVec, t,
                                    residual->localVector());PYLITH_CHECK_ERROR(err);

    // This doesn't work because we need the key to be able to pass the key as an argument to
    // distinguish between RHS and LHS residuals.
    err = DMPlexComputeBdResidualSingle(dsLabel.dm(), t, dsLabel.weakForm(), dsLabel.label(),
                                        numLabelValues, &labelValue, _fieldIndex, solution.localVector(), solutionDotVec,
                                        residual->localVector());PYLITH_CHECK_ERROR(err);
#endif

    PYLITH_METHOD_END;
} // computeRHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorBoundary::computeLHSResidual(pylith::topology::Field* residual,
                                                           const PylithReal t,
                                                           const PylithReal dt,
                                                           const pylith::topology::Field& solution,
                                                           const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<")");
    if (!_hasLHSResidual) { PYLITH_METHOD_END;}
    assert(residual);

    DSLabelAccess dsLabel(solution.dmMesh(), _labelName.c_str(), _labelValue);
    _setKernelConstants(solution, dt);

    PetscFormKey key;
    key.label = dsLabel.label();
    key.value = dsLabel.value();
    key.part = pylith::feassemble::Integrator::RESIDUAL_LHS;

    PetscErrorCode err;
    assert(solution.localVector());
    assert(residual->localVector());
#if 0
    // I think this is the interface we want.
    err = DMPlexComputeResidual_XXX(dsLabel.dm(), key, solution.localVector(), solutionDotVec, t,
                                    residual->localVector());PYLITH_CHECK_ERROR(err);
#endif

    PYLITH_METHOD_END;
} // computeLHSResidual


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorBoundary::computeLHSJacobian(PetscMat jacobianMat,
                                                           PetscMat precondMat,
                                                           const PylithReal t,
                                                           const PylithReal dt,
                                                           const PylithReal s_tshift,
                                                           const pylith::topology::Field& solution,
                                                           const pylith::topology::Field& solutionDot) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<", solutionDot="<<solutionDot.getLabel()<<") empty method");

    _needNewLHSJacobian = false;
    // No implementation needed for boundary.

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ---------------------------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorBoundary::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                    const PylithReal t,
                                                                    const PylithReal dt,
                                                                    const PylithReal s_tshift,
                                                                    const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", solution="<<solution.getLabel()<<") empty method");

    _needNewLHSJacobianLumped = false;
    // No implementation needed for boundary.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// End of file
