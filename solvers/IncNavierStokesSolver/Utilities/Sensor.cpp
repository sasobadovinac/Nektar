#include <cstdio>
#include <cstdlib>

#include <SolverUtils/Driver.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        fprintf(stderr,"Usage: ./Sensor file.xml \n");
        fprintf(stderr,"\t Method will read intiial conditions section of .xml file for input \n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr session;
    string vDriverModule;
    DriverSharedPtr drv;  
    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);
        
        // Create driver
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        drv = GetDriverFactory().CreateInstance(vDriverModule, session);


        EquationSystemSharedPtr EqSys = drv->GetEqu()[0];
        IncNavierStokesSharedPtr IncNav = EqSys->as<IncNavierStokes>();
        
        IncNav->SetInitialConditions(0.0,false);

        int i,n,nquad,cnt;

        Array<OneD, MultiRegions::ExpListSharedPtr> fields = IncNav->UpdateFields();
        int nfields = fields.num_elements();
        nquad = fields[0]->GetTotPoints();

        Array<OneD, int> vel = IncNav->GetVelocity();
        for(n = 0; n < nfields; ++n)
        {
            if(session->GetVariable(n) == "p")
            {
                break;
            }
        }
        ASSERTL0(n != nfields, "Could not find field named p in m_fields");

        Array<OneD, NekDouble> sensor = fields[n]->UpdatePhys();
        Array<OneD, NekDouble> tmp(nquad);
        Array<OneD, NekDouble> energy(nquad,0.0);

        // Evaluate velocity magnitude
        for(i = 0; i < vel.num_elements(); ++i)
        {
            Vmath::Vmul(nquad,fields[vel[i]]->GetPhys(),1,
                        fields[vel[i]]->GetPhys(),1,tmp,1);
            Vmath::Vadd(nquad,energy,1,tmp,1,energy,1);
        }
        Vmath::Vsqrt(nquad,energy,1,energy,1);
        IncNav->GetSensor(energy,sensor);
        
        fields[n]->FwdTrans_IterPerExp(fields[n]->GetPhys(),fields[n]->UpdateCoeffs());
        
        // Need to reset varibale name for output
        session->SetVariable(n,"Sensor");

        // Reset session name for output file
        std::string outname = IncNav->GetSessionName();
        
        outname += "_Sensor";
        IncNav->ResetSessionName(outname);
        IncNav->Output();

    }
    catch (const std::runtime_error&)
    {
        return 1;
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
    }


    return 0;
}

