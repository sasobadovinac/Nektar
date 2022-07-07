<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and IP diffusion, Implicit</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</file>
    </files>
    <metrics>
      <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-8">2.23593e-03</value>
            <value variable="rhou" tolerance="1e-8">2.15592e-03</value>
            <value variable="rhov" tolerance="1e-8">1.53450e-03</value>
            <value variable="E" tolerance="1e-8">5.23103e-03</value>
      </metric>
      <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-8">2.39146e-03</value>
            <value variable="rhou" tolerance="1e-8">3.21154e-03</value>
            <value variable="rhov" tolerance="1e-8">2.00662e-03</value>
            <value variable="E" tolerance="1e-8">4.59174e-03</value>
        </metric>
    </metrics>
</test>
