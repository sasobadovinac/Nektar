<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow 2D with Radiation outflow </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Channel_Flow_3modes_rad.xml</parameters>
    <files>
        <file description="Session File">Channel_Flow_3modes_rad.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">9.62277e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.77556e-16</value>
            <value variable="v" tolerance="1e-12">3.01203e-17</value>
            <value variable="p" tolerance="1e-12">1.77636e-15</value>
        </metric>
    </metrics>
</test>


