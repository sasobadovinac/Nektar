<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 10000 --timesteps 10 --method ImplicitSDC --variant GaussGaussLegendre --order 4 --param 1.0 2</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-10">1.83202e-05</value>
        </metric>
    </metrics>
</test>
