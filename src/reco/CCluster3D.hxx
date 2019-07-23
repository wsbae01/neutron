#ifndef CCluster3D_hxx_seen
#define CCluster3D_hxx_seen

#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include "TVector3.h"
#include "CHit3D.hxx"
#include "CBond3D.hxx"

/// The base class for a clustered hits

class CCluster3D : public TObject {
public:
    CCluster3D() ;
    ~CCluster3D() ;
    
    void SetId(int id);
    
    void AddConstituent(int hit) ;
    
    void SetStartPoint(int hit);
    
    void SetEndPoint(int hit);
   
    void AddBond(int bond);
   
    int GetId() const;
    
    std::vector<int> GetConstituents() const;
    
    int GetStartPoint() const;
    
    int GetEndPoint() const;
    
    std::vector<int> GetBonds() const;
    
    bool operator==(const CCluster3D& rhs) const ;
    
   
private:
    
    int fID;
    
    std::vector<int> fConstituents;
    
    int fStart;
    
    int fEnd;
    
    std::vector<int> fBonds;

    //ClassDef(CCluster3D,1);

};
#endif
