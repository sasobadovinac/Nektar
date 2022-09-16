<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D multiblock VTK output </description>
    <executable>FieldConvert</executable>
    <parameters> -f bfs_tg.xml bfs_tg.fld bfs_tg_multiblock.vtu:vtu:multiblock</parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="bfs_tg_multiblock.vtm">
                <sha1>197560b4809c3c63385b2551b94480effdcc1e53</sha1>
             </file>
         </metric>
    </metrics>
</test>
