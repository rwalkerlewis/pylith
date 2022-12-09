from .TestFaultCohesive import TestFaultCohesive
from .TestFaultCohesiveKin import TestFaultCohesiveKin
from .TestFaultCohesiveKinPoro import TestFaultCohesiveKinPoro
from .TestFaultCohesiveImpulses import TestFaultCohesiveImpulses
from .TestKinSrc import (
    TestKinSrc, 
    TestKinSrcConstRate,
    TestKinSrcStep,
    TestKinSrcRamp,
    TestKinSrcBrune,
    TestKinSrcLiuCos,
    TestKinSrcTimeHistory)
from .TestKinSrcPoro import (
    TestKinSrcPoro, 
    TestKinSrcPoroStep)    
from .TestSingleRupture import TestSingleRupture


def test_classes():
    return [
        TestFaultCohesive,
        TestFaultCohesiveKin,
        TestFaultCohesiveKinPoro,
        TestFaultCohesiveImpulses,
        TestKinSrc,
        TestKinSrcConstRate,
        TestKinSrcStep,
        TestKinSrcRamp,
        TestKinSrcBrune,
        TestKinSrcLiuCos,
        TestKinSrcTimeHistory,
        TestKinSrcPoro,
        TestKinSrcPoroStep,
        TestSingleRupture,
    ]


# End of file
