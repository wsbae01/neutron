#include "CBond3D.hxx"
#include <vector>
#include "CHit3D.hxx"

//ClassImp(CBond3D)

    CBond3D::CBond3D() {}
    CBond3D::~CBond3D() {}

    void CBond3D::SetId(int id){fID=id;}

    void CBond3D::AddConstituent(int clusterId){
        fConstituents.push_back(clusterId);
    }

    void CBond3D::SetPoint(int hit){
        fPoint = hit;
    }

    int CBond3D::GetId() const{return fID;}

    std::vector<int> CBond3D::GetConstituents() const{
        return fConstituents;
    }

    int CBond3D::GetPoint() const{
        return fPoint;
    }

    bool CBond3D::operator==(const CBond3D& rhs) const { return this->GetId() == rhs.GetId();}


#ifdef __CLING__
//#pragma link C++ class CBond3D+;
//#pragma link C++ class vector<CBond3D>+;
//#pragma link C++ class CHit3D+;
#endif

