#include "CCluster3D.hxx"
#include <vector>

//ClassImp(CCluster3D)

    CCluster3D::CCluster3D() {}
    CCluster3D::~CCluster3D() {}

    void CCluster3D::SetId(int id){fID=id;}

    void CCluster3D::AddConstituent(int hit){
        fConstituents.push_back(hit);
    }

    void CCluster3D::SetStartPoint(int hit){
        fStart = hit;
    }

    void CCluster3D::SetEndPoint(int hit){
        fEnd = hit;
    }

    void CCluster3D::AddBond(int bond){
        fBonds.push_back(bond);
    }

    int CCluster3D::GetId() const{return fID;}

    std::vector<int> CCluster3D::GetConstituents() const{
        return fConstituents;
    }

    int CCluster3D::GetStartPoint() const{
        return fStart;
    }

    int CCluster3D::GetEndPoint() const{
        return fEnd;
    }

    std::vector<int> CCluster3D::GetBonds() const{
        return fBonds;
    }

    bool CCluster3D::operator==(const CCluster3D& rhs) const { return this->GetId() == rhs.GetId();}


#ifdef __CLING__
//#pragma link C++ class CCluster3D+;
//#pragma link C++ class vector<CCluster3D>+;
#endif

