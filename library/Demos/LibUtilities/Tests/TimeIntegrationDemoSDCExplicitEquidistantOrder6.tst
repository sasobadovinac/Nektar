<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 600 --timesteps 200 --method ExplicitSDC --variant Equidistant --order 6 --param 1.0 6</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-09">0.000115413</value>
        </metric>
    </metrics>
</test>
