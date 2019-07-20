#include "CHit3D.hxx"
//#include "TVector3.h"
#include <vector>

//ClassImp(CHit3D)


   CHit3D::CHit3D() {};
   CHit3D::~CHit3D() {};

   void CHit3D::SetId(int id){fID=id;}

   void CHit3D::SetTime(double time){
       fTime=time;
    }

   void CHit3D::SetPosition(TVector3 position){
        fPosition.SetX(position.X());
        fPosition.SetY(position.Y());
        fPosition.SetZ(position.Z());
    }

   void CHit3D::SetCharge(double charge){
       fCharge=charge;
    }

    void CHit3D::AddConstituent(int hitsId, int plane){
        switch (plane) {
            case 0:
                fConstituentsYZ.push_back(hitsId);
                break;
            case 1:
                fConstituentsXZ.push_back(hitsId);
                break;
            case 2:
                fConstituentsXY.push_back(hitsId);
                break;
            default:
                std::cout<<"Wrong plane, no constituent added"<<std::endl;
                break;
        }
    }

    void CHit3D::SetFiberCharge(double charge, int plane){
        switch (plane) {
            case 0:
                fChargeYZ=charge;
                break;
            case 1:
                fChargeXZ=charge;
                break;
            case 2:
                fChargeXY=charge;
                break;
            default:
                std::cout<<"Wrong plane, no charge set"<<std::endl;
                break;
        }
    }

    std::vector<int> CHit3D::Get2DConstituents(int plane=0) const{
        if(plane==0){
            if(fConstituentsYZ.size()>0){
                return fConstituentsYZ;
            }else{
                std::cout<<"NO YZ Hits"<<std::endl;
                return {};
            }
        }
        if(plane==1){
            if(fConstituentsXZ.size()>0){
                return fConstituentsXZ;
            }else{
                std::cout<<"NO XZ Hits"<<std::endl;
                return {};
            }
        }
        if(plane==2){
            if(fConstituentsXY.size()>0){
                return fConstituentsXY;
            }else{
                std::cout<<"NO XY Hits"<<std::endl;
                return {};
            }
        }
        return {};
    }

    int CHit3D::GetId() const{return fID;}

    double CHit3D::GetTime() const {return fTime;}

    double CHit3D::GetCharge() const{return fCharge;}

   TVector3 CHit3D::GetPosition() const {return fPosition;}

    double CHit3D::GetFiberCharge(int plane) const{
        if(plane==0)return fChargeYZ;
        if(plane==1)return fChargeXZ;
        if(plane==2)return fChargeXY;
        return -1;
    }


    bool CHit3D::operator==(const CHit3D& rhs) const { return this->GetId() == rhs.GetId();}

    bool CHit3D::operator<(const CHit3D& rhs) const{

        if(this->GetPosition().X()!=rhs.GetPosition().X()){
            return this->GetPosition().X()<rhs.GetPosition().X();
        }else if(this->GetPosition().X()==rhs.GetPosition().X()){
            if(this->GetPosition().Z()!=rhs.GetPosition().Z()){
                return this->GetPosition().Z()<rhs.GetPosition().Z();
            }else if(this->GetPosition().Z()==rhs.GetPosition().Z()){
                return this->GetPosition().Y()<rhs.GetPosition().Y();
            }
        }
        return false;
    }



//#ifdef __CLING__
//#pragma link C++ class CHit3D+;
//#pragma link C++ class vector<CHit3D>+;
//#pragma link C++ class TVector3+;
//#endif

