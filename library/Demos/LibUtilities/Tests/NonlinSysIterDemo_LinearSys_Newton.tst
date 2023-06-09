<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>NonlinSysIterDemo for a small nonlinear problem of dimension 4</description>
    <executable>NonlinSysIterDemo</executable>
    <parameters>NonlinSysIterDemo_LinearSys_Newton.xml</parameters>
    <files>
        <file description="Session File">NonlinSysIterDemo_LinearSys_Newton.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-06">1.000000e+00</value>
            <value variable="v" tolerance="1e-15">3.52284e-15</value>
        </metric>
    </metrics>
</test>
