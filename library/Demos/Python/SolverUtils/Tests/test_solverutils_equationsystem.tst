<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description> Unit test of the Python interface for the
                  Nektar::SolverUtils::EquationSystem class.
    </description>
    <executable python="true"> test_solverutils_equationsystem.py </executable>
    <parameters></parameters>
    <files>
        <file description="Session File">../../../../../solvers/ADRSolver/Tests/Helmholtz3D_nodal.xml</file>
        <file description="Session File">../../../../../solvers/ADRSolver/Tests/SteadyAdvDiffReact2D_modal.xml</file>
        <file description="Session File">../../../../../solvers/ADRSolver/Tests/UnsteadyAdvection_FRDG_3DHomo1D_MVM.xml</file>
        <file description="Session File">../../../../../solvers/ADRSolver/Tests/MMFAdvSphere.xml</file>
        <file description="Session File">../../../../../solvers/ADRSolver/Tests/UnsteadyDiffusion_LDG_3DHomo1D_MVM.xml</file>
    </files>
     <metrics>
     <metric type="regex" id="1">
        <regex>^.*testEquationSystemCreateLaplaceSolver: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="2">
        <regex>^.*testEquationSystemCreatePoissonSolver: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="3">
        <regex>^.*testEquationSystemCreateHelmholtzSolver: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="4">
        <regex>^.*testEquationSystemCreateProjectionSolver: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="5">
        <regex>^.*testEquationSystemCreateSteadyAdvectionDiffusionSolver: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="6">
        <regex>^.*testEquationSystemCreateSteadyAdvectionDiffusionReactionSolver: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="7">
        <regex>^.*testEquationSystemCreateUnsteadyAdvectionSolver: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="8">
        <regex>^.*testEquationSystemCreateMMFAdvectionSolver: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PASS</field>
                  </match>
              </matches>
      </metric>
      <metric type="regex" id="9">
        <regex>^.*testEquationSystemCreateUnsteadyDiffusionSolver: (.*)</regex>
              <matches>
                  <match>
                      <field id="0">PAXX</field>
                  </match>
              </matches>
      </metric>

    </metrics>
</test>
