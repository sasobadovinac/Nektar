<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady CG Adv Diff, with QUad & Tri elmts P=5 with Gradient Stabilisation </description>
    <executable>ADRSolver</executable>
    <parameters>AdvDiff_QuadTri_GJP.xml</parameters>
    <files>
        <file description="Session File">AdvDiff_QuadTri_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-7">0.0371413</value>
            <value variable="v" tolerance="1e-7">0.0371413</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">0.0731431</value>
            <value variable="v" tolerance="1e-7">0.0731431</value>
        </metric>
    </metrics>
</test>

