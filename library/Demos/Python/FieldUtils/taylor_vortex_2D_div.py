#-f -m divergence -e taylor_vortex_2D.xml taylor_vortex_2D.fld taylor_vortex_2D_div.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml",  field, "taylor_vortex_2D.xml").Run()
InputModule.Create("fld",  field, "taylor_vortex_2D.fld").Run()
ProcessModule.Create("divergence", field).Run()
OutputModule.Create("dat", field, "taylor_vortex_2D_div.fld").Run()
