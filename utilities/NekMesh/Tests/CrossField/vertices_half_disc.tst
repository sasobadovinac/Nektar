<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Export of vertices for the cross field</description>
    <executable>NekMesh</executable>
    <parameters>-v vertices_half_disc.mcf:mcf:vertices=vertices.csv vertices_half_disc.xml</parameters>
    <files>
        <file description="Input File">vertices_half_disc.mcf</file>
        <file description="Input File 2">half_disc.geo</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^\[InputMCF\]\s+Found (\d+) vertices to export for cross field processing:</regex>
            <matches>
                <match>
                    <field id="0">2</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="2">
            <regex>^\[InputMCF\]\s+0\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)</regex>
            <matches>
                <match>
                    <field id="1">-5.84047e-17</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="3">
            <regex>^\[InputMCF\]\s+1\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)</regex>
            <matches>
                <match>
                    <field id="1">3.34524e-17</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
