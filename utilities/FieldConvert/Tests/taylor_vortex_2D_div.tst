<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D divergence output </description>
    <executable>FieldConvert</executable>
    <parameters> -f -m divergence -e taylor_vortex_2D.xml taylor_vortex_2D.fld taylor_vortex_2D_div.fld</parameters>
    <files>
        <file description="Session File">taylor_vortex_2D.xml</file>
	<file description="Session File">taylor_vortex_2D.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">5.69822</value>
            <value variable="y" tolerance="1e-6">5.69822</value>
            <value variable="u" tolerance="1e-4">1.38207</value>
            <value variable="v" tolerance="1e-4">1.38206</value>
            <value variable="p" tolerance="1e-4">0.608019</value>
            <value variable="divV" tolerance="1e-4">2.55371e-05</value>
        </metric>
    </metrics>
</test>
