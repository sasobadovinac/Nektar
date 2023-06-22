<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow with Taylor Expansion of P=9,8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_VCSWeakPressure_TaylorHood.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_VCSWeakPressure_TaylorHood.xml</file>
        <file description="Session File">KovaFlow_VCSWeakPressure_TaylorHood.rst</file>
    </files>
    <metrics>
      <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.47795e-09</value>
            <value variable="v" tolerance="1e-12">9.18025e-09</value>
            <value variable="p" tolerance="1e-12">2.74933e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.77977e-09</value>
            <value variable="v" tolerance="1e-12">1.70534e-08</value>
            <value variable="p" tolerance="1e-12">1.41102e-07</value>
        </metric>
    </metrics>
</test>
