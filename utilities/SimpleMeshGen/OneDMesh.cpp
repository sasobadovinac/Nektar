#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

#include<LibUtilities/BasicUtils/ErrorUtil.hpp>

// Use the stl version, primarily for string.
#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <tinyxml.h>


void Header(FILE *, int nel);
void Middle(FILE *);
void End(FILE *);

using namespace std;

void  PrintConditions(void);

int main(int argc, char *argv[])
{
    vector<double> xc; 
    int            nx = 0;
    int            i,j,k;
    
    if(argc != 2)
    {
        fprintf(stderr,"Usage OneDMesh file.xml\n");
        exit(1);
    }
    
    try{
        
        TiXmlDocument doc(argv[argc-1]);
        bool loadOkay = doc.LoadFile();
        
        std::stringstream errstr;
        errstr << "Unable to load file: " << argv[argc-1] << " (";
        errstr << doc.ErrorDesc() << ", line " << doc.ErrorRow()
               << ", column " << doc.ErrorCol() << ")";
        ASSERTL0(loadOkay, errstr.str());
        
        TiXmlHandle docHandle(&doc);
        TiXmlElement* master = NULL;   
        TiXmlElement* block = NULL;
        
        master = doc.FirstChildElement("NEKBLOCK");
        ASSERTL0(master, "Unable to find NEKBLOCK tag in file.");
        
        // Find the Mesh tag and same the dim and space attributes
        block = master->FirstChildElement("XBLOCK");
        
        ASSERTL0(block, "Unable to find XBLOCK tag in file.");
        TiXmlElement *val = block->FirstChildElement("X");    
        while (val)
        {
            TiXmlNode *xval = val->FirstChild();
            
            std::istringstream valDataStrm(xval->ToText()->Value());
            
            try
            {
                while(!valDataStrm.fail())
                {
                    double x_val;
                    valDataStrm >> x_val;
                    
                    if (!valDataStrm.fail())
                    {
                        xc.push_back(x_val);
                    }
                }
            }
            catch(...)
            {
                ASSERTL0(false, "Unable to read Xval data.");
            }
            
            val= val->NextSiblingElement("X");
        }
        
        
        cout << "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" << endl;
        cout << "<NEKTAR>" << endl;
 

        cout << "<EXPANSIONS>" << endl;
        cout << "<E COMPOSITE=\"C[0]\" NUMMODES=\"7\" FIELDS=\"u\" TYPE=\"MODIFIED\" />" <<endl;
        cout << "</EXPANSIONS>\n" << endl;

        PrintConditions();

        //Vertices 
        cout << "<GEOMETRY DIM=\"1\" SPACE=\"1\">" << endl;
        cout << "  <VERTEX>" << endl;
        
        nx = xc.size();
        for(k = 0; k < nx; ++k)
        {
            cout << "    <V ID=\"" << k << "\">\t";
            cout << std::setprecision(8)<< xc[k] << " "   << " 0.0  0.0";
            cout << "  </V>" << endl; 
        }
    
        cout << "  </VERTEX>\n" << endl;
        
        cout << "  <ELEMENT>" << endl;
        int cnt = 0;
        for(i = 0; i < nx-1; ++i)
        {
            cout << "    <S ID=\"" << cnt++ << "\">\t";
            cout << i <<"  " <<  i+1; 
            cout << "  </S>" << endl; 
        }
        cout << "  </ELEMENT>\n" << endl;


        cout << "<COMPOSITE>" << endl;
        cout << "<C ID=\"0\"> S[0-" << nx-2 << "] </C>" << endl;

        cout << "<C ID=\"1\"> V[0] </C>   // left  border" << endl;
        cout << "<C ID=\"2\"> V[" <<nx-1 <<"] </C>   // right border" << endl;

        cout << "</COMPOSITE>\n" << endl;

             
        cout << "<DOMAIN> C[0] </DOMAIN>\n" << endl;
        cout << "</GEOMETRY>\n" << endl;

        cout << "</NEKTAR>" << endl;

    }
    catch(...)
    {
        return 1;
    }

    return 0;

}



void  PrintConditions(void)
{
    cout << "<CONDITIONS>" << endl;
    
    cout << "<SOLVERINFO>" << endl;
    cout << "<I PROPERTY=\"SolverType\"        VALUE=\" \"/>" << endl;
    cout << "</SOLVERINFO>\n" << endl;
            
    cout << "<PARAMETERS>" << endl;
    cout << "<P> TimeStep      = 0.002  </P>" << endl;
    cout << "</PARAMETERS>\n" << endl;
    
    cout << "<VARIABLES>" << endl;
    cout << "  <V ID=\"0\"> u </V>" << endl; 
    cout << "</VARIABLES>\n" << endl;
    
    cout << "<BOUNDARYREGIONS>" << endl;
    cout << "<B ID=\"0\"> C[1] </B>" << endl;
    cout << "<B ID=\"1\"> C[2] </B>" << endl;
     cout << "</BOUNDARYREGIONS>\n" << endl;
    
    cout << "<BOUNDARYCONDITIONS>" << endl;
    cout << "  <REGION REF=\"0\"> // Left border " << endl;
    cout << "     <D VAR=\"u\" VALUE=\"0\" />"  << endl;
    cout << "  </REGION>" << endl;
            
    cout << "  <REGION REF=\"1\"> // Right border " << endl;
    cout << "     <D VAR=\"u\" VALUE=\"0\" />"  << endl;
    cout << "  </REGION>" << endl;
    cout << "</BOUNDARYCONDITIONS>" << endl;
    
    cout << "</CONDITIONS>" << endl;
}

