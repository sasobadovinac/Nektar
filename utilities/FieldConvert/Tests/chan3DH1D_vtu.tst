<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3DH1D vtu output </description>
    <executable>FieldConvert</executable>
    <parameters> -f -n 8 chan3DH1D.xml chan3DH1D.fld chan3DH1D.vtu</parameters>
    <files>
        <file description="Session File">chan3DH1D.xml</file>
	<file description="Session File">chan3DH1D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3DH1D.vtu">
                <sha1>617726b84303d00c41e194c015f206dfe33d1f12</sha1>
             </file>
         </metric>
    </metrics>
</test>
