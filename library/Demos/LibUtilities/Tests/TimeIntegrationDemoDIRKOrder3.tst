<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 1000 --timesteps 100 --method DIRK --order 3</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">6.11835e-06</value>
        </metric>
    </metrics>
</test>
