<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 100 --timesteps 100 --method ExplicitSDC --variant GaussLobattoLegendre --order 2 --param 1.0 2</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-09">0.000217302</value>
        </metric>
    </metrics>
</test>
