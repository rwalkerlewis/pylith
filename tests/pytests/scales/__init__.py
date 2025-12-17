from . import (
    TestGeneral,
    TestQuasistaticElasticity,
    TestDynamicElasticity,
    TestQuasistaticPoroelasticity,
    TestDynamicPoroelasticity,
    TestQuasistaticThermoelasticity,
    TestQuasistaticThermoporoelasticity,
)

def test_modules():
    return [
        TestGeneral,
        TestQuasistaticElasticity,
        TestDynamicElasticity,
        TestQuasistaticPoroelasticity,
        TestDynamicPoroelasticity,
        TestQuasistaticThermoelasticity,
        TestQuasistaticThermoporoelasticity,
    ]
