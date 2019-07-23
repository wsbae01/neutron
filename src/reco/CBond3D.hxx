#ifndef CBond3D_hxx_seen
#define CBond3D_hxx_seen

#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include "TVector3.h"
#include "CHit3D.hxx"

/// The base class for a clustered hits

class CBond3D : public TObject {
public:
    CBond3D();
    ~CBond3D() ;
    
    void SetId(int id);
    
    void AddConstituent(int clusterId);
    
    void SetPoint(int hit);
   
    int GetId() const;
    
    std::vector<int> GetConstituents() const;
    
    int GetPoint() const;
    
    bool operator==(const CBond3D& rhs) const; 

private:
    
    int fID;
    
    std::vector<int> fConstituents;
    
    int fPoint;
  
    //ClassDef(CBond3D,1);

};
#endif
