<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 400 --timesteps 200 --method ExplicitSDC --variant Equidistant --order 4 --param 1.0 4</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-09">0.000260333</value>
        </metric>
    </metrics>
</test>
