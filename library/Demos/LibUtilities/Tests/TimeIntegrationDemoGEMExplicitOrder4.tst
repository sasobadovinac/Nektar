<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 200 --timesteps 100 --method ExtrapolationMethod --variant ExplicitEuler --order 4</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.00104995</value>
        </metric>
    </metrics>
</test>
