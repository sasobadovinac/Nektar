<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2.5D Kovasznay flow with non-trivial w-velocity</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P8_8modes_FFTW_VCSWeakPressure.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P8_8modes_FFTW_VCSWeakPressure.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-09">4.5e-10</value>
            <value variable="v" tolerance="1e-09">2.4e-10</value>
            <value variable="w" tolerance="1e-09">2.2e-10</value>
            <value variable="p" tolerance="1e-08">2.2e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">1.3e-09</value>
            <value variable="v" tolerance="1e-08">2.4e-09</value>
            <value variable="w" tolerance="1e-09">7.3e-10</value>
            <value variable="p" tolerance="1e-08">5.8e-09</value>
        </metric>
    </metrics>
</test>
