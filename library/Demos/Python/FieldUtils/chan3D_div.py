#-f -m divergence -e chan3D.xml chan3D.fld chan3D_div.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, forceoutput=True, error=True)

InputModule.Create("xml",  field, "chan3D.xml").Run()
InputModule.Create("fld",  field, "chan3D.fld").Run()
ProcessModule.Create("divergence", field).Run()
OutputModule.Create("dat", field, "chan3D_div.fld").Run()
