<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 3DH1D vtu output using legacy writer </description>
    <executable>FieldConvert</executable>
    <parameters> -f -n 8 chan3DH1D.xml chan3DH1D.fld chan3DH1D_legacy.vtu:vtu:legacy</parameters>
    <files>
        <file description="Session File">chan3DH1D.xml</file>
	<file description="Session File">chan3DH1D.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="chan3DH1D_legacy.vtu">
                <sha1>bb645fe6ccdd6b7ef12a4f5aa477ae36d871e988</sha1>
             </file>
         </metric>
    </metrics>
</test>
