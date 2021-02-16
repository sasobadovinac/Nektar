<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Process cross field</description>
    <executable>FieldConvert</executable>
    <parameters>-m crossfield:mergetol=2:vertices=crossfield_half_disc.csv:step=0.05:outcsv=streamlines.csv crossfield_half_disc.xml crossfield_half_disc.fld out.stdout</parameters>
    <files>
        <file description="Input File">crossfield_half_disc.xml</file>
        <file description="Input File 2">crossfield_half_disc.fld</file>
        <file description="Input File 3">crossfield_half_disc.csv</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^([^\s]+) #(\d+) has (\d+) quadrants at \(([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?), ([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\)</regex>
            <matches>
                <match>
                    <field id="0">Vertex</field>
                    <field id="1">1</field>
                    <field id="2">1</field>
                    <field id="3">0</field>
                    <field id="4">-5.84047e-17</field>
                </match>
                <match>
                    <field id="0">Vertex</field>
                    <field id="1">2</field>
                    <field id="2">1</field>
                    <field id="3">1</field>
                    <field id="4">3.34524e-17</field>
                </match>
                <match>
                    <field id="0">Singularity</field>
                    <field id="1">1</field>
                    <field id="2">3</field>
                    <field id="3">0.247445</field>
                    <field id="4">0.252448</field>
                </match>
                <match>
                    <field id="0">Singularity</field>
                    <field id="1">2</field>
                    <field id="2">3</field>
                    <field id="3">0.752555</field>
                    <field id="4">0.252448</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="2">
            <regex>^Final number of streamlines \(after merging\): (\d+)</regex>
            <matches>
                <match>
                    <field id="0">5</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
