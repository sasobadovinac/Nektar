<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D internal circle with refinement inside a square</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list 2d_circle_square.mcf 2d_circle_square-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">2d_circle_square.mcf</file>
        <file description="Input File 2">2d_circle_square.stp</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>.*Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
