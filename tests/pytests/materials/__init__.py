from . import (
    TestMaterial,
    TestElasticity,
    TestIncompressibleElasticity,
    TestPoroelasticity,
    TestThermoelasticity,
    TestThermoporoelasticity,
    TestRheologies,
    TestAuxiliarySubfields,
    TestDerivedSubfields,
    TestHomogeneous,
)


def test_modules():
    return [
        TestMaterial,
        TestElasticity,
        TestIncompressibleElasticity,
        TestPoroelasticity,
        TestThermoelasticity,
        TestThermoporoelasticity,
        TestRheologies,
        TestAuxiliarySubfields,
        TestDerivedSubfields,
        TestHomogeneous,
    ]


# End of file
