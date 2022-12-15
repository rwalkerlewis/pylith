#if !defined(dynamicspontaneousrupture_hh)
#define dynamicspontaneousrupture_hh

#include <portinfo>

#include "Formulation.hh" // ISA Formulation

class DynamicSpontaneousRupture : public Formulation {
public:

    DynamicSpontaneousRupture(void);
    ~DynamicSpontaneousRupture(void);

protected:

    void _computeLHSResidual(const PetscReal t,
                             const PetscReal dt,
                             const PetscVec solution,
                             const PetscVec solutionDot,
                             PetscVec residual);

    void _computeRHSResidual(const PetscReal t,
                             const PetscReal dt,
                             const PetscVec solution,
                             PetscVec residual);

    void _computeLHSJacobian(const PetscReal t,
                             const PetscReal dt,
                             const PetscVec solution,
                             const PetscVec solutionDot,
                             const PetscReal shift,
                             PetscMat jacobian,
                             PetscMat preconditioner);

private:

    // Not implemented
    DynamicSpontaneousRupture(const DynamicSpontaneousRupture&);
    const DynamicSpontaneousRupture& operator=(const DynamicSpontaneousRupture&);

};

#endif

// End of file
