<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 200 --timesteps 200 --method ExplicitSDC --variant GaussLobattoLegendre --order 4 --param 1.0 3</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-08">0.00104917</value>
        </metric>
    </metrics>
</test>
