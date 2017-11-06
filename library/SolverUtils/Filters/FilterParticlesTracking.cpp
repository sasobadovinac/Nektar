///////////////////////////////////////////////////////////////////////////////
//
// File FilterParticlesTracking.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Particle tracking filter.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <SolverUtils/Filters/FilterParticlesTracking.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <boost/format.hpp>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterParticlesTracking::className =
    GetFilterFactory().RegisterCreatorFunction(
        "ParticlesTracking",FilterParticlesTracking::create);

/**
 *
 */
FilterParticlesTracking::FilterParticlesTracking(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const ParamMap &pParams):
    Filter(pSession)
{
    // Read parameters
        ParamMap::const_iterator it;

    // OutputFile
    it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        m_outputFile = it->second;
    }
    if (!(m_outputFile.length() >= 4
          && m_outputFile.substr(m_outputFile.length() - 4) == ".his"))
    {
        m_outputFile += ".his";
    }

     // Delta T
    it = pParams.find("DeltaT");
    if (it == pParams.end())
    {
        //ASSERTL0(1, "Missing parameter 'DeltaT'.");
        //std::cout << "1" << std::endl;
    }
    else
    {
		//ASSERTL0(it->second.length() > 0, "Missing parameter 'DeltaT'.");
		m_DeltaT = std::stod(it->second);
		//std::cout << "2" << std::endl;
    }
    
    
    // Max Iterations
    it = pParams.find("MaxIter");
    if (it == pParams.end())
    {
        //ASSERTL0(0, "Missing parameter 'MaxIter'.");
        //std::cout << "3" << std::endl;
    }
    else
    {
		//ASSERTL0(it->second.length() > 0, "Missing parameter 'MaxIter'.");
		m_MaxIter = std::stoi(it->second);
		//std::cout << "4" << std::endl;
		
    }
    

    // Points
    it = pParams.find("Points");
    ASSERTL0(it != pParams.end(), "Missing parameter 'Points'.");
    m_historyPointStream.str(it->second);
    m_index = 0;

}


/**
 *
 */
FilterParticlesTracking::~FilterParticlesTracking()
{
}


/**
 *
 */
void FilterParticlesTracking::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
	//std::cout << "Calling FilterParticlesTracking" << std::endl;
	
	// Open output stream
        Array<OneD, NekDouble>  gloCoord(3,0.0);
        int dim = pFields[0]->GetGraph()->GetSpaceDimension();

        bool adaptive;
        m_session->MatchSolverInfo("Driver", "Adaptive",
                                    adaptive, false);
        if (adaptive)
        {
            m_outputStream.open(m_outputFile.c_str(), ofstream::app);
        }
        else
        {
            m_outputStream.open(m_outputFile.c_str());
        }    
        
		m_historyPointStream>> gloCoord[0]
							>> gloCoord[1]
							>> gloCoord[2];


        //m_outputStream << "Positions for particle in X(0) = [";
		//m_outputStream << "  " << gloCoord[0];
		//m_outputStream << ", " << gloCoord[1];
		//m_outputStream << ", " << gloCoord[2];
		//m_outputStream << " ]" ;
		//m_outputStream << endl;
            
        m_outputStream << " X , Y , Z , U , V"<< endl;
        
        m_historyPoints.push_back(MemoryManager<SpatialDomains::PointGeom>
								::AllocateSharedPtr(dim, 0, gloCoord[0],
											gloCoord[1], gloCoord[2]));
								
        
        //}
        
////////Imprimir numeros

		//m_historyPoints[0]->GetCoords(  gloCoord[0],
                                        //gloCoord[1],
                                        //gloCoord[2]);

        //m_outputStream << " "   << boost::format("%25.19e") % gloCoord[0];
        //m_outputStream << " , " << boost::format("%25.19e") % gloCoord[1];
        //m_outputStream << " , " << boost::format("%25.19e") % gloCoord[2];
		//m_outputStream << endl;
        
}


/**
 *
 */
void FilterParticlesTracking::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
//std::cout<<"Now updating FilterParticlesTracking ..."<<time<<std::endl;
	
	//Array<OneD, NekDouble>  locCoord(3,0.0);
	//Array<OneD, NekDouble>  gloCoord(3,0.0);
	//Array<OneD, int>       idList  (1, -1   );
	//Array<OneD, NekDouble> physvals;
	//int dim = pFields[0]->GetGraph()->GetSpaceDimension();
	//NekDouble value;
	//int i, iter;
	//NekDouble dt = 0.001;
	
	//iter=m_historyPoints.size();
	
	////std::cout<<dim<<std::endl;
	
	//m_historyPoints[iter-1]->GetCoords(gloCoord[0],
                                  //gloCoord[1],
                                  //gloCoord[2]);
                                    
	//idList[0] = pFields[0]->GetExpIndex(gloCoord,locCoord,
                                //NekConstants::kNekZeroTol);
	
	//if (idList[0] != -1)
	//{
	
	//m_outputStream << " "   << boost::format("%25.19e") % gloCoord[0];
    //m_outputStream << " , " << boost::format("%25.19e") % gloCoord[1];
    //m_outputStream << " , " << boost::format("%25.19e") % gloCoord[2];
	
	
	
	//for (i = 0; i < dim; ++i)
        //{
			//physvals = pFields[i]->UpdatePhys() + pFields[i]->
									//GetPhys_Offset(idList[0]);
			
			//value = pFields[i]->GetExp(idList[0])->
					//StdPhysEvaluate(locCoord,physvals);
		
			//m_outputStream << " , " << boost::format("%25.19e") % value;
			
			//gloCoord[i]+=value*dt;
		//}
		//m_outputStream << endl;
		
		//m_historyPoints.push_back(MemoryManager<SpatialDomains::PointGeom>
								//::AllocateSharedPtr(dim,0, gloCoord[0],
								//gloCoord[1], gloCoord[2]));

		//////iter=m_historyPoints.size(); std::cout<<iter<<std::endl;
		
		
		////m_historyPoints[0]->GetCoords(gloCoord[0],
									  ////gloCoord[1],
                                      ////gloCoord[2]);
                                  
        //////iter=m_historyPoints.size(); std::cout<<"Last"<<iter<<std::endl;
	//}
	//else
	//{
	//std::cout << "Particle out of domain" << std::endl;
	//m_outputStream << endl;

	//}
	
	
	
}

/**
 *
 */
void FilterParticlesTracking::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
std::cout<<"Now evaluating FilterParticlesTracking ..."<<std::endl;
	
	Array<OneD, NekDouble>  locCoord(3,0.0);
	Array<OneD, NekDouble>  gloCoord(3,0.0);
	Array<OneD, int>        idList  (1, -1 );
	Array<OneD, NekDouble> physvals;
	int dim = pFields[0]->GetGraph()->GetSpaceDimension();
	NekDouble value;
	int i, iter;
	
	//m_MaxIter = 10000;
	//m_DeltaT  = 0.01;
	
	
	
	idList[0] = 0;
	iter = 0;
	
	//std::cout<<"maxiter : "<<m_MaxIter<<std::endl;
	//std::cout<<"dt : "<<m_DeltaT<<std::endl;
	
	while(idList[0] != -1 && iter<m_MaxIter)
	{
		iter=m_historyPoints.size();
		
		std::cout<<"Step  : "<<iter<<std::endl;
		
		m_historyPoints[iter-1]->GetCoords(gloCoord[0],
									       gloCoord[1],
										   gloCoord[2]);
										
		idList[0] = pFields[0]->GetExpIndex(gloCoord,locCoord,
									NekConstants::kNekZeroTol);
		
		if (idList[0] != -1)
		{
		
		m_outputStream << " "   << boost::format("%25.19e") % gloCoord[0];
		m_outputStream << " , " << boost::format("%25.19e") % gloCoord[1];
		m_outputStream << " , " << boost::format("%25.19e") % gloCoord[2];
		
		
		
		for (i = 0; i < dim; ++i)
			{
				physvals = pFields[i]->UpdatePhys() + pFields[i]->
										GetPhys_Offset(idList[0]);
				
				value = pFields[i]->GetExp(idList[0])->
						StdPhysEvaluate(locCoord,physvals);
			
				m_outputStream << " , " << boost::format("%25.19e") % value;
				
				gloCoord[i]+=value*m_DeltaT;
			}
			m_outputStream << endl;
			
			m_historyPoints.push_back(MemoryManager<SpatialDomains::PointGeom>
									::AllocateSharedPtr(dim,0, gloCoord[0],
									gloCoord[1], gloCoord[2]));

			////iter=m_historyPoints.size(); std::cout<<iter<<std::endl;
			
			
			//m_historyPoints[0]->GetCoords(gloCoord[0],
										  //gloCoord[1],
										  //gloCoord[2]);
									  
			////iter=m_historyPoints.size(); std::cout<<"Last"<<iter<<std::endl;
		}
		else
		{
		m_outputStream << endl;
		std::cout << "Particle out of domain" << std::endl;
		}
	}
	std::cout << "Time = " << iter*m_DeltaT<<std::endl;
	m_outputStream.close();
}


/**
 *
 */
bool FilterParticlesTracking::v_IsTimeDependent()
{
    return true;
}
}
}
