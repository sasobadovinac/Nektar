<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and IP diffusion, Implicit parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_IM.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM.xml</file>
    </files>
    <metrics>
      <metric type="L2" id="1">
        <value variable="rho" tolerance="1e-8">0.00225256</value>
        <value variable="rhou" tolerance="1e-8">0.00217252</value>
        <value variable="rhov" tolerance="1e-8">0.00161124</value>
        <value variable="E" tolerance="1e-8">0.00526196</value>
      </metric>
      <metric type="Linf" id="2">
        <value variable="rho" tolerance="1e-8">0.00249404</value>
        <value variable="rhou" tolerance="1e-8">0.00366045</value>
        <value variable="rhov" tolerance="1e-8">0.00210457</value>
        <value variable="E" tolerance="1e-8">0.00480471</value>
      </metric>
    </metrics>
</test>


