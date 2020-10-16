<test>
    <description>BC test for implicit cfs, fix rho&p on the top</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>BC_EnforceRhoP.xml</parameters>
    <files>
        <file description="Session File">BC_EnforceRhoP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-4">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-3">0</value>
        </metric>
    </metrics>
</test>

