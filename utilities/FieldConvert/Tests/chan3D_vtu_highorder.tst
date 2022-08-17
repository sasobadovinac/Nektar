<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3D vtu output high order with uncompress option </description>
    <executable>FieldConvert</executable>
    <parameters> -f -n 8 chan3D.xml chan3D.fld chan3D_highorder.vtu:vtu:highorder:uncompress</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3D_highorder.vtu">
                <sha1>7cdde2cd97eec048b6dc50b992ad20150517dfe9</sha1>
             </file>
         </metric>
    </metrics>
</test>
