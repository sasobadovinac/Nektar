<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Export of vertices for the cross field</description>
    <executable>NekMesh</executable>
    <parameters>-v vertices_half_disc.mcf:mcf:vertices=vertices.csv vertices_half_disc.xml</parameters>
    <files>
        <file description="Input File">vertices_half_disc.mcf</file>
        <file description="Input File 2">vertices_half_disc.geo</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Found (\d+) vertices to export for cross field processing:</regex>
            <matches>
                <match>
                    <field id="0">2</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
