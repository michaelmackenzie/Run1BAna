#include "DetectorConstruction.hh"
#include "ScreenSD.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4TwistedBox.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4Transform3D.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

namespace twisted
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  DetectorConstruction::DetectorConstruction()
  {
    // define commands for this class
    DefineCommands();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  DetectorConstruction::~DetectorConstruction()
  {
    delete fMessenger;

    for (auto visAttributes: fVisAttributes) {
      delete visAttributes;
    }
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4VPhysicalVolume* DetectorConstruction::Construct()
  {
    // Construct materials
    auto nistMan = G4NistManager::Instance();
    ConstructMaterials();
    auto Vacuum = G4Material::GetMaterial("Vacuum");
    //auto TargetMaterial =  G4Material::GetMaterial("Inconel718");
    //auto TargetMaterial = nistMan->FindOrBuildMaterial("G4_C");
    auto TargetMaterial = nistMan->FindOrBuildMaterial("G4_Al");
    //G4cout << "WARNING: using Cu instead of Inconel!!" << G4endl;

    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;

    // geometries --------------------------------------------------------------
    // ****** VOLUMES  ******
    // Define size of Target and volumes around it.
    //
    // - Target (G4Twisted)
    G4double segment_twist = fTotalTwist/fNSegments;
    G4double segment_dx = fPDx;
    G4double segment_dy = fPDy;
    G4double segment_dz = fPDz/fNSegments;

    // - World (Box)
    G4double world_dx = fPDx+1.*cm;
    G4double world_dy = fPDy+1.*cm;
    G4double world_dz = fPDz+5.*cm;

    // ****** SOLID VOLUMES ******
    // - World
    fsolidWorld = new G4Box("World", world_dx/2., world_dy/2., world_dz/2.);
    // - Segment
    fsolidSegment = new G4TwistedBox(
				     "Segment_solid",   // Name
				     segment_twist,             // Twist angle
				     segment_dx / 2.0,         // Half-length along X
				     segment_dy / 2.0,         // Half-length along Y
				     segment_dz / 2.0          // Half-length along Z
				     );
    G4cout << "SEGMENT DZ=" << segment_dz << G4endl;

    // ******* LOGICAL VOLUMES *******
    // - World
    auto world_log =
      new G4LogicalVolume(fsolidWorld,          //its solid
			  Vacuum,           //its material
			  "world_log");            //its name
    // - Segment
    auto segment_log =
      new G4LogicalVolume(fsolidSegment,          //its solid
			  TargetMaterial,           //its material
			  "segment_log");            //its name



    // ******* PHYSICAL VOLUMES *******
    // - World
    //	G4PVPlacement (G4RotationMatrix *pRot, const G4ThreeVector &tlate, G4LogicalVolume *pCurrentLogical, const G4String &pName, G4LogicalVolume *pMotherLogical, G4bool pMany, G4int pCopyNo, G4bool pSurfChk=false)
    fPhysWorld =
      new G4PVPlacement(0,                     //no rotation
			G4ThreeVector(),       //at (0,0,0)
			world_log,             //its logical volume
			"World",               //its name
			0,                     //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking

    // Loop to place the target segments end-to-end with matching rotations

    for (G4int i = 0; i < fNSegments; ++i) {
      // Calculate the Z position for the center of this segment
      // Shifts segments from -totalZ/2 to +totalZ/2
      G4double zPos = -fPDz / 2.0 + (i + 0.5) * segment_dz;

      //Compute the exact progressive alignment angle.
      // Segment 0 must be clocked by half its twist to align its starting face to 0.
      G4double alignmentAngle = (i * segment_twist) + (segment_twist / 2.0);

      //Create the rotation matrix for this specific segment's clocking
      G4RotationMatrix rot;
      rot.rotateZ(alignmentAngle);

      // Combine them into a precise 3D Transformation.
      G4Transform3D transform(rot, G4ThreeVector(0, 0, zPos));

      // Calculate the cumulative rotation for this segment's starting orientation
      G4double currentRotationAngle = i * segment_twist + segment_twist/2.;
      G4RotationMatrix* segRot = new G4RotationMatrix();
      segRot->rotateZ(currentRotationAngle);
      // Dynamic name for tracking each copy
      G4String instanceName = "segment_log_" + std::to_string(i);

      // Place the segment inside your mother volume (worldLogical)
      new G4PVPlacement(transform,                     // rotation
			segment_log,             //its logical volume
			instanceName,               //its name
			world_log,                     //its mother  volume
			false,                 //no boolean operation
			i,                     //copy number
			checkOverlaps);        //overlaps checking
    }
    // visualization attributes ------------------------------------------------

    auto visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    visAttributes->SetVisibility(false);
    world_log->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);

    visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.9));   // LightGray
    segment_log->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);


    // return the world physical volume ----------------------------------------

    return fPhysWorld;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::ConstructSDandField()
  {
    // sensitive detectors -----------------------------------------------------
    auto sdManager = G4SDManager::GetSDMpointer();
    G4String SDname;

  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::ConstructMaterials()
  {
    //    G4double a, z, density,
    G4double density, pressure, temperature;
    // G4String name, symbol;
    density     = universe_mean_density;    //from PhysicalConstants.h
    pressure    = 3.e-18*pascal;
    temperature = 2.73*kelvin;
    new G4Material("Vacuum", 1., 1.01*g/mole, density,
		   kStateGas,temperature,pressure);

    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* Inconel718 = new G4Material( "Inconel718", 8.19*CLHEP::g/CLHEP::cm3, 11);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Al"),0.005);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Cr"),0.190);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Co"),0.010);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Cu"),0.003);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Fe"),0.170);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Mn"),0.003);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Mo"),0.030);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Ni"),0.527);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Nb"),0.052);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Si"),0.003);
    Inconel718->AddMaterial(nistManager->FindOrBuildMaterial("G4_Ti"),0.007);



    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void DetectorConstruction::DefineCommands()
  {
    // Define /inconel/detector command directory using generic messenger class
    fMessenger = new G4GenericMessenger(this,
					"/twisted/detector/",
					"Detector control");

  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
