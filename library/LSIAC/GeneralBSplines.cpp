#include "GeneralBSplines.h"

namespace Nektar
{
namespace LSIAC
{

GeneralBSplines::GeneralBSplines(const int order)
{
    this->SetOrder(order);
}

GeneralBSplines::GeneralBSplines(const Array<OneD, NekDouble> &knots,
                                 const int order)
{
    this->SetKnotVector(knots);
    this->SetOrder(order);
}

bool GeneralBSplines::SetKnotVector(const Array<OneD, NekDouble> &knots)
{
    this->m_knotVector = knots;
    return true;
}

bool GeneralBSplines::GetKnotVector(Array<OneD, NekDouble> &knots)
{
    knots = this->m_knotVector;
    return true;
}

int GeneralBSplines::GetOrder() const
{
    return this->m_order;
}

bool GeneralBSplines::SetOrder(const int order)
{
    this->m_order = order;
    return true;
}

bool GeneralBSplines::EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                                       const std::vector<NekDouble> &kvec,
                                       const int j_th,
                                       Array<OneD, NekDouble> &t_values,
                                       const NekDouble shift,
                                       const NekDouble meshScaling) const
{
    return BSplinesBasis(t_pos, kvec, m_order - 1, j_th, t_values, shift,
                         meshScaling);
}

bool GeneralBSplines::EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                                       const int j_th,
                                       Array<OneD, NekDouble> &t_values,
                                       const NekDouble shift,
                                       const NekDouble meshScaling) const
{

    BSplinesBasis(t_pos, m_order - 1, j_th, t_values, shift, meshScaling);
    return true;
}

bool GeneralBSplines::EvaluateBSplines(const Array<OneD, NekDouble> &t_pos,
                                       const Array<OneD, NekDouble> &knots,
                                       int j_th,
                                       Array<OneD, NekDouble> &t_values)
{
    boost::ignore_unused(t_pos, knots, j_th, t_values);
    NEKERROR(ErrorUtil::efatal, "Not yet coded");
    return false;
}

bool GeneralBSplines::BSplinesBasis(const Array<OneD, NekDouble> &t_pos,
                                    const int k, const int j,
                                    Array<OneD, NekDouble> &t_val,
                                    const NekDouble shift,
                                    const NekDouble meshScaling) const
{
    // Note here Order of BSplines are k+1.
    // This is done to follow the paper. Sorry if causes confusion.
    //	cout << "came here" << endl;
    if (0 == k)
    {
        if (j >= 0 && j < m_knotVector.num_elements() - 1)
        {
            for (int i = 0; i < t_pos.num_elements(); i++)
            {
                NekDouble u = (t_pos[i] - shift) / meshScaling;
                if (j == m_knotVector.num_elements() - 1)
                {
                    if ((m_knotVector[j] <= u) && (u <= m_knotVector[j + 1]))
                    {
                        t_val[i] = 1.0;
                    }
                    else
                    {
                        t_val[i] = 0.0;
                    }
                }
                else
                {
                    if ((m_knotVector[j] <= u) && (u < m_knotVector[j + 1]))
                    {
                        t_val[i] = 1.0;
                    }
                    else
                    {
                        t_val[i] = 0.0;
                    }
                }
            }
            /* Working not a bad option
                                            if (j==0)
                                            {
                                                    cout << "came here j=0" <<
               endl;
                                                    //if( (
               (t_pos[i]-shift)/meshScaling>= m_knotVector[j]) & (
               (t_pos[i]-shift)/meshScaling<= m_knotVector[j+1])) if(
               (m_knotVector[j] <= u ) && (u <= m_knotVector[j+1]) )
                                                    {
                                                            t_val[i] = 	1.0;
                                                    }else{
                                                            t_val[i] = 0.0;
                                                    }
                                            }else{
                                            //if( (t_pos[i]/meshScaling-shift >=
               m_knotVector[j]) & (t_pos[i]/meshScaling-shift <
               m_knotVector[j+1]))
                                                    //if( (
               (t_pos[i]-shift)/meshScaling>= m_knotVector[j]) & (
               (t_pos[i]-shift)/meshScaling< m_knotVector[j+1])) cout << "came
               here j!=0" << endl; if(  (m_knotVector[j] < u ) && (u <=
               m_knotVector[j+1]) )
                                                    {
                                                            t_val[i] = 	1.0;
                                                    }else{
                                                            t_val[i] = 0.0;
                                                    }
                                            }
            */
        }
        else
        {
            // Add Assert here. Should not come here.
            printf("Add Asset here");
        }
    }
    else
    {
        NekDouble x_eval, w_jlt, w_jlt1;
        Array<OneD, NekDouble> Bspl_k1_j, Bspl_k1_j1;
        Bspl_k1_j  = Array<OneD, NekDouble>(t_pos.num_elements(), 0.0);
        Bspl_k1_j1 = Array<OneD, NekDouble>(t_pos.num_elements(), 0.0);
        this->BSplinesBasis(t_pos, k - 1, j, Bspl_k1_j, shift, meshScaling);
        this->BSplinesBasis(t_pos, k - 1, j + 1, Bspl_k1_j1, shift,
                            meshScaling);
        for (int i = 0; i < t_pos.num_elements(); i++)
        {
            // x_eval = t_pos[i]/meshScaling-shift;
            x_eval = (t_pos[i] - shift) / meshScaling;
            if (abs(m_knotVector[j + k] - m_knotVector[j]) <= 1e-8)
            {
                w_jlt = 0.0;
            }
            else
            {
                w_jlt = (x_eval - m_knotVector[j]) /
                        (m_knotVector[j + k] - m_knotVector[j]);
            }
            if (abs(m_knotVector[j + k + 1] - m_knotVector[j + 1]) <= 1e-8)
            {
                w_jlt1 = 1.0;
            }
            else
            {
                w_jlt1 = (x_eval - m_knotVector[j + 1]) /
                         (m_knotVector[j + k + 1] - m_knotVector[j + 1]);
            }
            t_val[i] = w_jlt * Bspl_k1_j[i] + (1 - w_jlt1) * Bspl_k1_j1[i];
        }
    }
    return true;
}

bool GeneralBSplines::BSplinesBasis(const Array<OneD, NekDouble> &t_pos,
                                    const vector<NekDouble> &kVec, const int k,
                                    const int j, Array<OneD, NekDouble> &t_val,
                                    const NekDouble shift,
                                    const NekDouble meshScaling) const
{
    // cout << "Praful" << endl;
    // Note here Order of BSplines are k+1.
    // This is done to follow the paper. Sorry if causes confusion.
    if (0 == k)
    {
        for (int i = 0; i < t_pos.num_elements(); i++)
        {
            NekDouble u = (t_pos[i] - shift) / meshScaling;
            if (kVec[j + 1] < kVec.back())
            {
                if ((kVec[j] <= u) && (u < kVec[j + 1]))
                {
                    t_val[i] = 1.0;
                }
                else
                {
                    t_val[i] = 0.0;
                }
            }
            else
            {
                // if( (t_pos[i]/meshScaling-shift >= kVec[j]) &
                // (t_pos[i]/meshScaling-shift < kVec[j+1])) if( (
                // (t_pos[i]-shift)/meshScaling>= kVec[j]) & (
                // (t_pos[i]-shift)/meshScaling< kVec[j+1]))
                if ((kVec[j] <= u && u <= kVec[j + 1]))
                {
                    t_val[i] = 1.0;
                }
                else
                {
                    t_val[i] = 0.0;
                }
            }
            /*
                                            if( ( (t_pos[i]-shift)/meshScaling>=
               kVec[j]) & ( (t_pos[i]-shift)/meshScaling< kVec[j+1]))
                                            //if(  (kVec[j] < u ) && (u <=
               kVec[j+1]) )
                                            {
                                                    t_val[i] = 	1.0;
                                            }else{
                                                    t_val[i] = 0.0;
                                            }
            */
        }
    }
    else
    {
        NekDouble x_eval, w_jlt, w_jlt1;
        Array<OneD, NekDouble> Bspl_k1_j, Bspl_k1_j1;
        Bspl_k1_j  = Array<OneD, NekDouble>(t_pos.num_elements(), 0.0);
        Bspl_k1_j1 = Array<OneD, NekDouble>(t_pos.num_elements(), 0.0);
        this->BSplinesBasis(t_pos, kVec, k - 1, j, Bspl_k1_j, shift,
                            meshScaling);
        this->BSplinesBasis(t_pos, kVec, k - 1, j + 1, Bspl_k1_j1, shift,
                            meshScaling);
        for (int i = 0; i < t_pos.num_elements(); i++)
        {
            // x_eval = t_pos[i]/meshScaling-shift;
            x_eval = (t_pos[i] - shift) / meshScaling;
            if (abs(kVec[j + k] - kVec[j]) <= 1e-8)
            {
                w_jlt = 0.0;
            }
            else
            {
                w_jlt = (x_eval - kVec[j]) / (kVec[j + k] - kVec[j]);
            }
            if (abs(kVec[j + k + 1] - kVec[j + 1]) <= 1e-8)
            {
                w_jlt1 = 1.0;
            }
            else
            {
                w_jlt1 =
                    (x_eval - kVec[j + 1]) / (kVec[j + k + 1] - kVec[j + 1]);
            }
            t_val[i] = w_jlt * Bspl_k1_j[i] + (1 - w_jlt1) * Bspl_k1_j1[i];
        }
    }
    return true;
}

/*

bool GeneralBSplines::BSplinesBasis( const Array<OneD,NekDouble> &t_pos,
Array<OneD,NekDouble> &t_val,const int k,const int j)const {
        // Note here Order of BSplines are k+1.
        // This is done to follow the paper. Sorry if causes confusion.
        if (0==k) {
                if ( j>=0 && m_knotVector.num_elements()-1){
                        for (int i =0; i < t_pos.num_elements(); i++)
                        {
                                if( (t_pos[i] >= m_knotVector[j]) & (t_pos[i] <
m_knotVector[j+1]))
                                {
                                        t_val[i] = 	1.0;
                                }else{
                                        t_val[i] = 0.0;
                                }
                        }
                }else{
                        //Add Assert here. Should not come here.
                        printf("Add Asset here");
                }
        }else{
                NekDouble x_eval, w_jlt,w_jlt1;
                Array<OneD,NekDouble> Bspl_k1_j, Bspl_k1_j1;
                Bspl_k1_j = Array<OneD,NekDouble>(t_pos.num_elements(),0.0);
                Bspl_k1_j1 = Array<OneD,NekDouble>(t_pos.num_elements(),0.0);
                this->BSplinesBasis( t_pos, Bspl_k1_j,  k-1, j);
                this->BSplinesBasis( t_pos, Bspl_k1_j1, k-1, j+1);
                for (int i =0; i < t_pos.num_elements(); i++)
                {
                        x_eval = t_pos[i];
                        if( abs(m_knotVector[j+k] - m_knotVector[j]) <= 1e-8){
                                w_jlt =0.0;
                        }else{
                                w_jlt = ( x_eval - m_knotVector[j]
)/(m_knotVector[j+k]-m_knotVector[j]);
                        }
                        if( abs(m_knotVector[j+k+1] - m_knotVector[j+1]) <=
1e-8){ w_jlt1 =0.0; }else{ w_jlt1 = ( x_eval - m_knotVector[j+1]
)/(m_knotVector[j+k+1]-m_knotVector[j+1]);
                        }
                        t_val[i] = w_jlt*Bspl_k1_j[i] +
(1-w_jlt1)*Bspl_k1_j1[i];
                }
        }
        return true;
}

*/
} // namespace LSIAC
} // namespace Nektar
