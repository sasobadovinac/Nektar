<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 10000 --timesteps 10 --method ExtrapolationMethod --variant ImplicitEuler --order 3</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-08">0.000425624</value>
        </metric>
    </metrics>
</test>
