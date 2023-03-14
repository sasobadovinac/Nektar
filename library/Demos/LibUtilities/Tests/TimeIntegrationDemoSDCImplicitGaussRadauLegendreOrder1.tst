<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 10000 --timesteps 10 --method ImplicitSDC --variant GaussRadauLegendre --order 1 --param 1.0 1</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-08">0.188226</value>
        </metric>
    </metrics>
</test>
