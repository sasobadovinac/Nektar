<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process add fields output </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m addfld:fromfld=chanflow_m8_addfld.fld:scale=1 chanflow_m8_addfld.xml chanflow_m8_addfld.fld chanflow_m8_addedfld.fld</parameters>
    <files>
        <file description="Session File">chanflow_m8_addfld.xml</file>
        <file description="Session File">chanflow_m8_addfld.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.57735</value>
            <value variable="y" tolerance="1e-6">0.57735</value>
            <value variable="u" tolerance="1e-6">0.365148</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">2.3094</value>
        </metric>
    </metrics>
</test>

