<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Variational optimiser test on all tet cube/sphere</description>
    <executable>NekMesh</executable>
    <parameters>cube-sphere.msh test.xml:xml:test -v -m varopti:hyperelastic:numthreads=2:maxiter=5:nq=4</parameters>
    <files>
        <file description="Input File">cube-sphere.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="0">
            <regex>Worst at end: (-?\d+(?:\.\d*)?(?:[eE][+\-]?\d+)?)</regex>
            <matches>
                <match>
                    <field id="0">6.663593e-01</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
