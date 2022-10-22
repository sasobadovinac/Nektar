<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D vtk output with a range restriction, with hdf5 input files</description>
    <executable>FieldConvert</executable>
    <parameters> -f -r -1,1,-1,1 -e bfs_tg_hdf5.xml bfs_tg_hdf5.fld bfs_tg_hdf5_rng.vtu</parameters>
    <files>
	<file description="Mesh File">bfs_tg_hdf5.xml</file>
	<file description="Mesh File">bfs_tg_hdf5.nekg</file>
	<file description="Session File">bfs_tg_hdf5.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">1.35354</value>
            <value variable="y" tolerance="1e-6">1.09498</value>
            <value variable="u" tolerance="1e-6">1.15931</value>
            <value variable="v" tolerance="1e-8">0.0195245</value>
            <value variable="p" tolerance="1e-7">0.119152</value>
        </metric>
    </metrics>
</test>

