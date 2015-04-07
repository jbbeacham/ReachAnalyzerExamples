makeAcceptanceMap(){

  //Double_t delta[11] = {-0.045,-0.04,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.04,0.045};
  Double_t delta[10] = {-0.0425,-0.035,-0.025,-0.015,-0.005,0.005,0.015,0.025,0.035,0.0425};
  Double_t zbeam[5] = {250.0,125.0,0.0,-125.0,-250.0};
  Double_t omega[50] = {2.27,3.26,3.3,3.36,3.41,3.47,3.51,3.4,3.25,2.85,2.9,3.88,3.98,4.0,4.05,4.15,4.17,4.08,3.89,3.38,3.89,4.02,4.15,4.13,4.25,4.32,4.31,4.26,3.99,3.65,3.86,4.14,4.16,4.24,4.18,4.3,4.27,4.33,4.21,3.6,3.64,4.2,4.21,4.24,4.2,4.18,4.19,4.15,4.11,3.35};

  //TGraph2D* omegaSolidAngle = new TGraph2D(50,delta,zbeam,omega);
  TGraph2D* omegaSolidAngle = new TGraph2D(50);

  //for( Int_t jj = 0; jj < 50; ++jj ){
  
  Int_t iOmega = 0;

  for( Int_t mm = 0; mm < 5; ++mm ){
    for( Int_t kk = 0; kk < 10; ++kk ){
      omegaSolidAngle->SetPoint(iOmega,delta[kk],zbeam[mm],omega[iOmega]);
      ++iOmega;
    }
  }


  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  omegaSolidAngle->Draw("lego");
  //std::cout << omegaSolidAngle->Interpolate(0.13,0.03) << std::endl;
  Double_t thetaLo[5] = {4.75,4.5,4.0,3.75,3.5};
  Double_t thetaHi[5] = {7.0,6.75,6.5,6.25,6.25};

  TGraph* thetaHorizLo = new TGraph(5,zbeam,thetaLo);
  TGraph* thetaHorizHi = new TGraph(5,zbeam,thetaHi);

  // Save the graphs

  omegaSolidAngle->SaveAs("Graphs/OmegaSolidAngle.root");

  thetaHorizLo->SaveAs("Graphs/ThetaHorizLo.root");
  thetaHorizHi->SaveAs("Graphs/ThetaHorizHi.root");

  // So, for example, zbeam = 125

  Double_t thisZbeam = 125.0;
  Double_t thisDelta = 0.015;

  Double_t thisOmega = omegaSolidAngle->Interpolate(thisDelta,thisZbeam);

  Double_t xLo = thetaHorizLo->Eval(thisZbeam);
  Double_t xHi = thetaHorizHi->Eval(thisZbeam);

  Double_t vertPlusMinus = thisOmega/(xHi-xLo);

  Double_t yLo = -vertPlusMinus;
  Double_t yHi = vertPlusMinus;

  TBox* accBox = new TBox(xLo,yLo,xHi,yHi);


  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  TH1F* frame = c2->DrawFrame(4.0,-3.0,7.5,3.0);
  accBox->SetFillColor(kRed);
  accBox->Draw();

}
