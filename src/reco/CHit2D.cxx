#include "CHit2D.hxx"
#include <vector>

//ClassImp(CHit2D)

    CHit2D::CHit2D(){};
    CHit2D::~CHit2D(){};

    void CHit2D::SetId(int id){fId=id;}

    void CHit2D::SetConstituents(const std::vector<int>& id){
        for(std::size_t i=0;i<id.size();++i){
            fConstituents.push_back(id[i]);
        }
   }

    void CHit2D::SetRow(double r){fRow=r;}

    void CHit2D::SetColumn(double c){fColumn=c;}

    void CHit2D::SetTime(double time){fTime=time;}

    void CHit2D::SetCharge(double charge){fCharge=charge;}

    void CHit2D::SetPlane(double plane){fPlane=plane;}

    int CHit2D::GetId() const {return fId;}

    std::vector<int> CHit2D::GetConstituents() const {return fConstituents;}

    double CHit2D::GetRow() const {return fRow;}

    double CHit2D::GetColumn() const {return fColumn;}

    double CHit2D::GetTime() const {return fTime;}

    double CHit2D::GetCharge()const {return fCharge;}

    int CHit2D::GetPlane() const {return fPlane;}

    bool CHit2D::operator==(const CHit2D& rhs) const { return this->GetId() == rhs.GetId() && this->GetPlane() == rhs.GetPlane();}


//#ifdef __CLING__
//#pragma link C++ class CHit2D+;
//#pragma link C++ class vector<CHit2D>+;
//#endif
