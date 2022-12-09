from .TestInitialCondition import TestInitialCondition
from .TestInitialConditionDomain import TestInitialConditionDomain
from .TestInitialConditionPatch import TestInitialConditionPatch
from .TestPhysics import TestPhysics
from .TestProblem import TestProblem
from .TestTimeDependent import TestTimeDependent
from .TestProblemDefaults import TestProblemDefaults
from .TestProgressMonitor import TestProgressMonitor
from .TestProgressMonitorTime import TestProgressMonitorTime
from .TestSingleObserver import (TestSinglePhysicsObserver, TestSingleSolnObserver)
from .TestSolution import TestSolution
from .TestSolnDisp import TestSolnDisp
from .TestSolnDispLagrange import (TestSolnDispLagrange, TestSolutionDispLagrange)
from .TestSolnDispPres import (TestSolnDispPres, TestSolutionDispPres)
from .TestSolnDispLagrange import (TestSolnDispLagrange, TestSolutionDispLagrange)
from .TestSolnDispPresTracStrain import (TestSolnDispPresTracStrain, TestSolutionDispPresTracStrain)
from .TestSolnDispPresTracStrainVelPdotTdot import (TestSolnDispPresTracStrainVelPdotTdot, TestSolutionDispPresTracStrainVelPdotTdot)
from .TestSolnDispPresTracStrainLagrange import (TestSolnDispPresTracStrainLagrange, TestSolutionDispPresTracStrainLagrange)
from .TestSolnDispPresTracStrainLagrangeFaultPres import (TestSolnDispPresTracStrainLagrangeFaultPres, TestSolutionDispPresTracStrainLagrangeFaultPres)
from .TestSolnDispVel import (TestSolnDispVel, TestSolutionDispVel)
from .TestSolnDispVelLagrange import (TestSolnDispVelLagrange, TestSolutionDispVelLagrange)
from .TestSolutionSubfield import TestSolutionSubfield
from .TestSubfieldDisplacement import TestSubfieldDisplacement
from .TestSubfieldLagrangeFault import TestSubfieldLagrangeFault
from .TestSubfieldFaultPressure import TestSubfieldFaultPressure
from .TestSubfieldPressure import TestSubfieldPressure
from .TestSubfieldPressureDot import TestSubfieldPressureDot
from .TestSubfieldTemperature import TestSubfieldTemperature
from .TestSubfieldTraceStrain import TestSubfieldTraceStrain
from .TestSubfieldTraceStrainDot import TestSubfieldTraceStrainDot
from .TestSubfieldVelocity import TestSubfieldVelocity


def test_classes():
    classes = [
        TestInitialCondition,
        TestInitialConditionDomain,
        TestInitialConditionPatch,
        TestPhysics,
        TestProblem,
        TestTimeDependent,
        TestProblemDefaults,
        TestProgressMonitorTime,
        TestProgressMonitorTime,
        TestSingleSolnObserver,
        TestSinglePhysicsObserver,
        TestSolution,
        TestSolnDisp,
        TestSolnDispLagrange,
        TestSolutionDispLagrange,
        TestSolnDispPresTracStrain, 
        TestSolutionDispPresTracStrain,
        TestSolnDispPresTracStrainLagrange, 
        TestSolutionDispPresTracStrainLagrange,
        TestSolnDispPresTracStrainLagrangeFaultPres, 
        TestSolutionDispPresTracStrainLagrangeFaultPres,
        TestSolnDispPresTracStrainVelPdotTdot,
        TestSolutionDispPresTracStrainVelPdotTdot,
        TestSolnDispPresTracStrainLagrange, 
        TestSolnDispPresTracStrainLagrange,         
        TestSolnDispVel, 
        TestSolutionDispVel,
        TestSolnDispVelLagrange, 
        TestSolutionDispVelLagrange,
        TestSolutionSubfield,
        TestSubfieldDisplacement,
        TestSubfieldLagrangeFault,
        TestSubfieldFaultPressure,
        TestSubfieldPressure,
        TestSubfieldPressureDot,
        TestSubfieldTemperature,
        TestSubfieldTraceStrain,
        TestSubfieldTraceStrainDot,
        TestSubfieldVelocity,
    ]
    return classes


# End of file
