<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--Npoints 100 --Ntimesteps 100 --NTimeIntegrationMethod 1</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.224678</value>
        </metric>
    </metrics>
</test>
