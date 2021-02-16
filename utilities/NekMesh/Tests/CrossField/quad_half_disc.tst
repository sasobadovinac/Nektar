<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Split the geometry into quads, export it, mesh it then refine the mesh in all directions</description>
    <executable>NekMesh</executable>
    <parameters>-v -m bl:surf=4,7:surf2=2,10:layers=7:r=1.0:reproject -m bl:surf=4,3,2:surf2=6,11,8:layers=5:r=1.0:reproject -m bl:surf=5,1,9:surf2=7,3,10:layers=3:r=1.0:reproject -m jac:list quad_half_disc.mcf quad_half_disc.xml</parameters>
    <files>
        <file description="Input File">quad_half_disc.mcf</file>
        <file description="Input File 2">half_disc.geo</file>
        <file description="Input File 3">quad_half_disc.csv</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^\[CADSystem\]\s+Found (\d+) wires</regex>
            <matches>
                <match>
                    <field id="0">1</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="2">
            <regex>^\[CADSystem\]\s+There are (\d+) faces after splitting</regex>
            <matches>
                <match>
                    <field id="0">4</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="3">
            <regex>.*Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
