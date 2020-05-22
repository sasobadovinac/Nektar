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
        <!-- <metric type="regex" id="1">
            <regex>^Vertex #1 has (\d+) quadrants at</regex>
            <matches>
                <match>
                    <field id="0">1</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="2">
            <regex>^Vertex #2 has (\d+) quadrants at</regex>
            <matches>
                <match>
                    <field id="0">1</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="3">
            <regex>^Singularity #1 has (\d+) quadrants at</regex>
            <matches>
                <match>
                    <field id="0">3</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="4">
            <regex>^Singularity #2 has (\d+) quadrants at</regex>
            <matches>
                <match>
                    <field id="0">3</field>
                </match>
            </matches>
        </metric> -->
        <metric type="regex" id="5">
            <regex>^Final number of streamlines \(after merging\): (\d+)</regex>
            <matches>
                <match>
                    <field id="0">5</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
