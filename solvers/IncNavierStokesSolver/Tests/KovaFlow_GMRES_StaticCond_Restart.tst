<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow with IterativeStaticCond and GMRES solver to test restart procedure</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_GMRES_StaticCond_Restart.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_GMRES_StaticCond_Restart.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.92379e-09</value>
            <value variable="v" tolerance="1e-12">8.34003e-09</value>
            <value variable="p" tolerance="1e-12">1.29442e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.52073e-09</value>
            <value variable="v" tolerance="1e-12">1.63561e-08</value>
            <value variable="p" tolerance="1e-12">5.14353e-08</value>
        </metric>
    </metrics>
</test>
