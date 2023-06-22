<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 800 --timesteps 200 --method ExplicitSDC --variant Equidistant --order 8 --param 1.0 8</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-10">6.48390e-05</value>
        </metric>
    </metrics>
</test>
