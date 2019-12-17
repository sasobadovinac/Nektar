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
        assert(false && "file does not exist.");
    }

    std::string line;
    std::getline(source, line);
    std::istringstream in(line);

    for (int i = 0; i < ar.num_elements(); i++)
    {
        in >> ar[i];
        //		cout << std::setprecision(19)<<ar[i] << endl;
    }
}

void LSIACPostProcessor::readNekArray(vector<int> &ar, string filename) const
{
    ifstream source;
    source.open(filename.c_str(), ios_base::in);

    if (!source)
    {
        assert(false && "file does not exist.");
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
        assert(false && "file does not exist.");
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
