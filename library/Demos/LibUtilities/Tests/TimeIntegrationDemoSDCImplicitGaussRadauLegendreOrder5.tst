<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 10000 --timesteps 10 --method ImplicitSDC --variant GaussRadauLegendre --order 5 --param 1.0 3</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-10">2.26209e-07</value>
        </metric>
    </metrics>
</test>
