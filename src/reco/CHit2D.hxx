#ifndef CHit2D_hxx_seen
#define CHit2D_hxx_seen

#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include <iostream>
#include <vector>

using namespace std;
/// The base class for a hit detector element 2D. Works like fiber map.

class CHit2D: public TObject {
public:
    CHit2D() ;
    ~CHit2D();    

    void SetId(int id);
    
    void SetConstituents(const std::vector<int>& id);
    
    void SetRow(double r);
    
    void SetColumn(double c);
    
    void SetTime(double time);
    
    void SetCharge(double charge) ;
    
    void SetPlane(double plane) ;
    
    int GetId() const ;
    
    std::vector<int> GetConstituents() const ;
    
    double GetRow() const ;
    
    double GetColumn() const ;
    
    double GetTime() const ;
    
    double GetCharge()const ;
    
    int GetPlane() const ;
    
    bool operator==(const CHit2D& rhs) const ; 
    
    
   
private:
    
    int fId;
    
    std::vector<int> fConstituents;
    
   
    double fTime;
    
    double fRow;
    
    double fColumn;
    
    double fCharge;
    
    //0=YZ, 1=XZ, 2=XY
    double fPlane;

   //ClassDef(CHit2D,1);
};


#endif
