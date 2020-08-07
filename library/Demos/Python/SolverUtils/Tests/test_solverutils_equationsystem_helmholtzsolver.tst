<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description> Unit test of the Python interface for the
                  Nektar::SolverUtils::EquationSystem class.
    </description>
    <executable python="true">test_solverutils_equationsystem_helmholtzsolver.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">../../../../../solvers/ADRSolver/Tests/Helmholtz3D_nodal.xml</file>
    </files>
     <metrics>
      <metric type="regex" id="1">
        <regex>^.*testEquationSystemCreateHelmholtzSolver: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
    </metrics>
</test>
