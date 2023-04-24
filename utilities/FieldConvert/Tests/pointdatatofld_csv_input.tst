<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Project values at quadrature points onto expansion when values are stored in a .csv file </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m pointdatatofld:frompts=flat-plane.csv -m addfld:fromfld=flat-plane.fld:scale=-1 flat-plane.xml flat-plane-new.fld </parameters>
    <files>
        <file description="CSV File">flat-plane.csv</file>
        <file description="Field File">flat-plane.fld</file>
        <file description="Session File">flat-plane.xml</file>
    </files>
    <!--
        This test is designed to test the implementation of the pointdatatofld module

        The pointdatatofld module can accept input in the form of a .pts file or a .csv file,
        this test is designed to test the case when the input is a .csv file.

        The file flat-plane.csv contains the x,y,z coordinates at the Gauss points defined in
        the flat-plane.xml file (exported with the "noequispaced" flag). The .csv file also 
        contains one additional field, defined as u = x + y. I.e., u is a linear function.

        The flat-plane.fld also contains the field u.

        The test is designed to check that the difference between the solution obtained with
        pointdatatofld and the solution stored in flat-plane.fld are equal.
    -->
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-12">0.57735</value>
            <value variable="y" tolerance="1e-12">0.57735</value>
            <value variable="z" tolerance="1e-12">0.</value>
            <value variable="u" tolerance="1e-12">0.</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-12">1.</value>
            <value variable="y" tolerance="1e-12">1.</value>
            <value variable="z" tolerance="1e-12">0.</value>
            <value variable="u" tolerance="1e-12">0.</value>
        </metric>
    </metrics>
</test>

