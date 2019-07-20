#ifndef CHit3D_hxx_seen
#define CHit3D_hxx_seen

#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include "TVector3.h"
#include "CHit2D.hxx"


/// The base class for a hit detector element.

class CHit3D : public TObject {
public:
   CHit3D() ;
   virtual ~CHit3D() ;
    
    void SetId(int id);
    
    void SetTime(double time);
    
    void SetPosition(TVector3 position); 
    
    void SetCharge(double charge); 
    
    void AddConstituent(int hitsId, int plane); 
    
    void SetFiberCharge(double charge, int plane);
    
    std::vector<int> Get2DConstituents(int plane) const;

    int GetId() const;
    
    double GetTime() const;
    
    double GetCharge() const;
    
    TVector3 GetPosition() const;
    
    double GetFiberCharge(int plane) const;
    
    bool operator==(const CHit3D& rhs) const;
    
    bool operator<(const CHit3D& rhs) const; 
   
private:
    
    int fID;
    
    double fTime;
    
    TVector3 fPosition;
    
    double fCharge;
    
    double fChargeXY;
    double fChargeXZ;
    double fChargeYZ;
    
    //0=YZ, 1=XZ, 2=XY
    double fPlane;
    
    std::vector<int> fConstituentsXY;
    std::vector<int> fConstituentsXZ;
    std::vector<int> fConstituentsYZ;
    
//ClassDef(CHit3D,1);

};
#endif
