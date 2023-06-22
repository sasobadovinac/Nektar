<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 100 --timesteps 100 --method ExplicitSDC --variant Equidistant --order 1 --param 1.0 1</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="5e-03">0.223365</value>
        </metric>
    </metrics>
</test>
