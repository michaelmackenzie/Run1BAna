void draw_cal_time_peaks() {
    // Create geometry
    TGeoManager *geom = new TGeoManager("cylinders", "Cylinders with Cone");

    TGeoMaterial *matVac = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMedium *medVac = new TGeoMedium("Vacuum", 1, matVac);

    // Create master volume
    TGeoVolume *top = geom->MakeBox("top", medVac, 10, 10, 10);
    geom->SetTopVolume(top);

    // Parameters
    const double dz_target(0.05), offset_trk(1.), dz_trk(2.), dz_calo(0.2), offset_calo(0.2);
    const double z0_trk(dz_target/2. + dz_trk/2. + offset_trk);
    const double z0_calo(z0_trk + dz_trk/2. + dz_calo/2. + offset_calo);
    const double r_target(0.2), r_1_trk(0.35), r_2_trk(0.8);

    // Draw the target
    TGeoTube *target = new TGeoTube("target", 0., r_target, dz_target/2.);
    TGeoVolume *target_vol = new TGeoVolume("target_vol", target, medVac);
    target_vol->SetLineColor(kBlue);
    target_vol->SetTransparency(5);
    top->AddNode(target_vol, 1, new TGeoTranslation(0, 0, 0.));

    // Draw the tracker
    TGeoTube *tracker = new TGeoTube("tracker", r_1_trk, r_2_trk, dz_trk/2.);
    TGeoVolume *trk_vol = new TGeoVolume("tracker_vol", tracker, medVac);
    trk_vol->SetLineColor(kAtlantic);
    trk_vol->SetTransparency(5);
    top->AddNode(trk_vol, 2, new TGeoTranslation(0, 0, z0_trk));

    // Draw the calo
    TGeoTube *calo = new TGeoTube("calo", r_1_trk, r_2_trk, dz_calo/2.);
    TGeoVolume *calo_vol = new TGeoVolume("calo_vol", calo, medVac);
    calo_vol->SetLineColor(kMagenta);
    calo_vol->SetTransparency(5);
    top->AddNode(calo_vol, 3, new TGeoTranslation(0, 0, z0_calo));

    // Create cone from the calo to the target
    TGeoCone *cone = new TGeoCone("cone", z0_calo/2., 0., 0.2, 0., 0.);
    TGeoVolume *cone_vol = new TGeoVolume("conevol", cone, medVac);
    cone_vol->SetLineColor(kRed);
    cone_vol->SetTransparency(0);
    TGeoRotation *rot = new TGeoRotation();
    rot->RotateX(-10.);
    TGeoTranslation *trans = new TGeoTranslation(0., (r_1_trk + r_2_trk)/4., z0_calo/2.);
    TGeoCombiTrans *combi = new TGeoCombiTrans(*trans, *rot);
    top->AddNode(cone_vol, 4, combi);

    //   // Parameters for cylinders
    // Double_t outerR = 2.0;
    // Double_t innerR = 0.8;
    // Double_t thickness = 0.8;
    // Double_t zPos[] = {2, 4, 6, 8};

    // // Create and place four cylinders with holes
    // for(Int_t i = 0; i < 4; i++) {
    //     TGeoTube *outerCyl = new TGeoTube(Form("outer%d", i), 0, outerR, thickness/2);
    //     TGeoTube *hole = new TGeoTube(Form("hole%d", i), 0, innerR, thickness/2);
    //     TGeoCompositeShape *cylinder = new TGeoCompositeShape(
    //         Form("cyl%d", i),
    //         Form("outer%d - hole%d", i, i)
    //     );
    //     TGeoVolume *cylVol = new TGeoVolume(Form("cylvol%d", i), cylinder, medVac);
    //     cylVol->SetLineColor(kBlue + i);

    //     top->AddNode(cylVol, i, new TGeoTranslation(0, 0, zPos[i]));
    // }

    // // Create cone from third cylinder (z=6) to first cylinder (z=2)
    // TGeoCone *cone = new TGeoCone("cone", 2.0, 1.0, 1.5, 1.0, 1.5);
    // TGeoVolume *coneVol = new TGeoVolume("conevol", cone, medVac);
    // coneVol->SetLineColor(kRed);

    // top->AddNode(coneVol, 0, new TGeoTranslation(0, 0, 4));

    geom->CloseGeometry();
    top->Draw("ogl");
}
