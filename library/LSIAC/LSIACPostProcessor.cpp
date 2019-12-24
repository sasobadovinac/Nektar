///////////////////////////////////////////////////////////////////////////////
//
// File LSIACPostProcessor.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: LSIACPostProcessor definition
//
///////////////////////////////////////////////////////////////////////////////
#include "LSIACPostProcessor.h"
#include <fstream>
#include <iomanip>

namespace Nektar
{
namespace LSIAC
{

void LSIACPostProcessor::printNekArray(const Array<OneD, NekDouble> &ar) const
{
    cout << "size of array: " << ar.num_elements() << endl;
    for (int i = 0; i < ar.num_elements(); i++)
    {
        cout << ar[i] << "\t";
    }
    cout << endl;
}

void LSIACPostProcessor::printNekArray(const vector<NekDouble> &ar) const
{
    cout << "size of array: " << ar.size() << endl;
    for (int i = 0; i < ar.size(); i++)
    {
        cout << ar[i] << "\t";
    }
    cout << endl;
}

void LSIACPostProcessor::printNekArray(const vector<int> &ar) const
{
    cout << "size of array: " << ar.size() << endl;
    for (int i = 0; i < ar.size(); i++)
    {
        cout << ar[i] << "\t";
    }
    cout << endl;
}

void LSIACPostProcessor::printGraphArray(const Array<OneD, NekDouble> &test,
                                         NekDouble down, NekDouble up,
                                         NekDouble increment) const
{
    int size = test.num_elements();
    for (NekDouble i = up; i >= down; i -= increment)
    {
        cout << fixed << setprecision(3) << i << "\t";
        for (int e = 0; e < size; e++)
        {
            if ((test[e] <= i) && (test[e] > i - increment))
            {
                cout << 'e';
            }
            else
            {
                cout << ' ';
            }
        }
        cout << endl;
    }
}

void LSIACPostProcessor::writeNekArray(vector<int> &ar, string filename) const
{
    ofstream myfile;
    myfile.open(filename.c_str());
    for (int i = 0; i < ar.size(); i++)
    {
        myfile << std::setprecision(19) << ar[i] << "\t";
    }
    myfile << endl;
    myfile.close();
}
void LSIACPostProcessor::writeNekArray(vector<NekDouble> &ar,
                                       string filename) const
{
    ofstream myfile;
    myfile.open(filename.c_str());
    for (int i = 0; i < ar.size(); i++)
    {
        myfile << std::setprecision(19) << ar[i] << "\t";
    }
    myfile << endl;
    myfile.close();
}

void LSIACPostProcessor::writeNekArray(Array<OneD, NekDouble> &ar,
                                       string filename) const
{
    ofstream myfile;
    myfile.open(filename.c_str());
    for (int i = 0; i < ar.num_elements(); i++)
    {
        myfile << std::setprecision(19) << ar[i] << "\t";
    }
    myfile << endl;
    myfile.close();
}

void LSIACPostProcessor::writeNekArrayBin(Array<OneD, NekDouble> &ar,
                                          string filename) const
{
    ofstream myfile(filename.c_str(), ios::binary);
    myfile.write((char *)(&(ar[0])), ar.num_elements() * sizeof(NekDouble));
    myfile.close();
}

void LSIACPostProcessor::readNekArray(Array<OneD, NekDouble> &ar,
                                      string filename) const
{
    ifstream source;
    source.open(filename.c_str(), ios_base::in);

    if (!source)
    {
        NEKERROR(ErrorUtil::efatal, "file does not exist.");
    }

    std::string line;
    std::getline(source, line);
    std::istringstream in(line);

    for (int i = 0; i < ar.num_elements(); i++)
    {
        in >> ar[i];
    }
}

void LSIACPostProcessor::readNekArray(vector<int> &ar, string filename) const
{
    ifstream source;
    source.open(filename.c_str(), ios_base::in);

    if (!source)
    {
        NEKERROR(ErrorUtil::efatal, "file does not exist.");
    }
    std::string line;
    std::getline(source, line);
    std::istringstream inter(line);

    int x;
    while (inter >> x)
    {
        ar.push_back(x);
    }
}

void LSIACPostProcessor::readNekArray(vector<NekDouble> &ar,
                                      string filename) const
{
    ifstream source;
    source.open(filename.c_str(), ios_base::in);

    if (!source)
    {
        NEKERROR(ErrorUtil::efatal, "file does not exist.");
    }
    std::string line;
    while (!source.eof())
    {
        std::getline(source, line);
        std::istringstream inter(line);

        NekDouble x;
        while (inter >> x)
        {
            ar.push_back(x);
        }
    }
}

} // namespace LSIAC
} // namespace Nektar
