#ifndef B5DetectorConstruction_h
#define B5DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4TwistedBox.hh"

#include <vector>

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;

namespace twisted
{

  class MagneticField;

  /// Detector construction

  class DetectorConstruction : public G4VUserDetectorConstruction
  {
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    void ConstructMaterials();
    G4VPhysicalVolume* GetWorldVolume() const { return fPhysWorld; }

  private:
    void DefineCommands();

    G4GenericMessenger* fMessenger;

    G4LogicalVolume* fScreen_log = nullptr;
    G4LogicalVolume* fFront_log = nullptr;
    G4LogicalVolume* fBack_log = nullptr;

    std::vector<G4VisAttributes*> fVisAttributes;

    G4VPhysicalVolume* fPhysWorld = nullptr;
    G4Box* fsolidWorld  = nullptr;
    G4TwistedBox* fsolidSegment = nullptr;

    G4double fTotalTwist  = 360.0 * CLHEP::deg;
    G4int fNSegments  = 8;
    G4double fPDx    = 120.0 * CLHEP:: mm;
    G4double fPDy    = 0.4 * CLHEP:: mm;
    G4double fPDz    = 500.0 * CLHEP:: mm;

  };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
