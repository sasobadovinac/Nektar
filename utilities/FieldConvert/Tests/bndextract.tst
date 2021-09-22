<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>extract boundary, 3D channel flow, Hexahedral elements, P=3</description>
    <executable>FieldConvert</executable>
    <parameters>-m extract  Hex_channel_m3.xml Hex_channel_m3_0.chk bnd.fld -e -f</parameters>
    <files>
        <file description="Session File">Hex_channel_m3.xml</file>
        <file description="field File">Hex_channel_m3_0.chk</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-9">0.57735</value>
            <value variable="y" tolerance="1e-9">0.57735</value>
            <value variable="z" tolerance="1e-9">0</value>
            <value variable="u" tolerance="1e-9">0.549417</value>
            <value variable="v" tolerance="1e-9">0.00120599</value>
            <value variable="w" tolerance="1e-9">0.0104854</value>
            <value variable="p" tolerance="1e-9">1.24022e-14</value>
            <value variable="x" tolerance="1e-9">0.57735</value>
            <value variable="y" tolerance="1e-9">0.57735</value>
            <value variable="z" tolerance="1e-9">1</value>
            <value variable="u" tolerance="1e-9">0.296847</value>
            <value variable="v" tolerance="1e-9">0.0837328</value>
            <value variable="w" tolerance="1e-9">0.484771</value>
            <value variable="p" tolerance="1e-9">0.408285</value>
            <value variable="x" tolerance="1e-9">1</value>
            <value variable="y" tolerance="1e-9">0.57735</value>
            <value variable="z" tolerance="1e-9">0.57735</value>
            <value variable="u" tolerance="1e-9">0.00164865</value>
            <value variable="v" tolerance="1e-9">0.00153491</value>
            <value variable="w" tolerance="1e-9">0.531203</value>
            <value variable="p" tolerance="1e-9">0.00082937</value>
            <value variable="x" tolerance="1e-9">0</value>
            <value variable="y" tolerance="1e-9">0.57735</value>
            <value variable="z" tolerance="1e-9">0.57735</value>
            <value variable="u" tolerance="1e-9">0.00164865</value>
            <value variable="v" tolerance="1e-9">0.00153491</value>
            <value variable="w" tolerance="1e-9">0.47056</value>
            <value variable="p" tolerance="1e-9">0.00082937</value>
            <value variable="x" tolerance="1e-9">0.57735</value>
            <value variable="y" tolerance="1e-9">1</value>
            <value variable="z" tolerance="1e-9">0.57735</value>
            <value variable="u" tolerance="1e-9">0.520171</value>
            <value variable="v" tolerance="1e-9">0.276368</value>
            <value variable="w" tolerance="1e-9">0.0723644</value>
            <value variable="p" tolerance="1e-9">0.408285</value>
            <value variable="x" tolerance="1e-9">0.57735</value>
            <value variable="y" tolerance="1e-9">0</value>
            <value variable="z" tolerance="1e-9">0.57735</value>
            <value variable="u" tolerance="1e-9">0.000215851</value>
            <value variable="v" tolerance="1e-9">0.511513</value>
            <value variable="w" tolerance="1e-9">0.00104225</value>
            <value variable="p" tolerance="1e-9">1.32712e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-9">1</value>
            <value variable="y" tolerance="1e-9">1</value>
            <value variable="z" tolerance="1e-9">0</value>
            <value variable="u" tolerance="1e-9">0.971219</value>
            <value variable="v" tolerance="1e-9">4.97516e-06</value>
            <value variable="w" tolerance="1e-9">0.0203256</value>
            <value variable="p" tolerance="1e-9">3.93574e-14</value>
            <value variable="x" tolerance="1e-9">1</value>
            <value variable="y" tolerance="1e-9">1</value>
            <value variable="z" tolerance="1e-9">1</value>
            <value variable="u" tolerance="1e-9">0.524745</value>
            <value variable="v" tolerance="1e-9">0.138487</value>
            <value variable="w" tolerance="1e-9">0.939717</value>
            <value variable="p" tolerance="1e-9">0.997512</value>
            <value variable="x" tolerance="1e-9">1</value>
            <value variable="y" tolerance="1e-9">1</value>
            <value variable="z" tolerance="1e-9">1</value>
            <value variable="u" tolerance="1e-9">8.90468e-07</value>
            <value variable="v" tolerance="1e-9">4.97516e-06</value>
            <value variable="w" tolerance="1e-9">0.939717</value>
            <value variable="p" tolerance="1e-9">2.90323e-14</value>
            <value variable="x" tolerance="1e-9">0</value>
            <value variable="y" tolerance="1e-9">1</value>
            <value variable="z" tolerance="1e-9">1</value>
            <value variable="u" tolerance="1e-9">8.90468e-07</value>
            <value variable="v" tolerance="1e-9">4.97516e-06</value>
            <value variable="w" tolerance="1e-9">0.138374</value>
            <value variable="p" tolerance="1e-9">2.90323e-14</value>
            <value variable="x" tolerance="1e-9">1</value>
            <value variable="y" tolerance="1e-9">1</value>
            <value variable="z" tolerance="1e-9">1</value>
            <value variable="u" tolerance="1e-9">0.860318</value>
            <value variable="v" tolerance="1e-9">0.536351</value>
            <value variable="w" tolerance="1e-9">0.138374</value>
            <value variable="p" tolerance="1e-9">0.997512</value>
            <value variable="x" tolerance="1e-9">1</value>
            <value variable="y" tolerance="1e-9">0</value>
            <value variable="z" tolerance="1e-9">1</value>
            <value variable="u" tolerance="1e-9">8.90468e-07</value>
            <value variable="v" tolerance="1e-9">0.992701</value>
            <value variable="w" tolerance="1e-9">0.00199298</value>
            <value variable="p" tolerance="1e-9">2.57097e-14</value>
        </metric>
    </metrics>
</test>

