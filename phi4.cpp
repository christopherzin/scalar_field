//
// File: phi4.cpp
// Chris Zin 
// January 2015
//
// Simulates the equations of motion obtained from 
//  a phi^4 Lagrangian
// Can make a .gif of the motion
// Calculates the field tensor and plots the average
//  value over the grid of some derived quantities
//


#include "TROOT.h" 
#include "TStyle.h" 
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "nr3.h"
#include "zMatOps.h"
#include "zTime.h"
#include "calcTuv.h"
#include "classTensor.h"
#include "prm.h"

// Updates the fields to the next time step
void update(const Mat3DDoub&,const Mat3DDoub&,Mat3DDoub&,
            const Mat3DDoub&,const Mat3DDoub&,Mat3DDoub&,const double[]);
// Makes a .gif if variable makeGIF is set to true
// Delete/Rename previous gif or it will append the new gif to the end            
void MakeGIF(const Mat3DDoub&,TH2D*,TCanvas*);
void MakeGIF4(const Mat3DDoub&,const Mat3DDoub&,const Mat3DDoub&,const Mat3DDoub&,
              TH2D*,TH2D*,TH2D*,TH2D*,TCanvas*);

using namespace std;

int main() 
{
  // This file contains grid output at various times
  // Used to plot in plot_phi4.C
  ofstream outFile, params;
  outFile.open("data/phi4.dat");
  
  // Used to keep track of real time elapsed
  clock_t startTime = clock();
  double secsPassed = 0.;
  
  bool makeGIF = prm::makeGIF,
       makeGIF4 = prm::makeGIF4;

  // Input parameters from prm namespace in file prm.h
  // Note: Any change in the input parameters should be done
  //  in prm.h so it applies to all files  
  const int X_GRID = prm::X_GRID,		
            Y_GRID = prm::Y_GRID,
            ETA_GRID = prm::Y_GRID,
            T_STEPS = prm::T_STEPS;

  double t0 = prm::t0,
         t = prm::t,
         dt = prm::dt,
         dx = prm::dx,
         deta = prm::deta, 
         m = prm::m,
         lambda = prm::lambda;
         
  // Parameter array: p = [dt, dx, deta, m, lambda, tn]
  // Bundles the parameters nicely for passing to functions
  // May be easier/better to use the prm namespace variables
  //  in the functions but I was using this long before 
  //  I thought about namespace constants...
  double p[] = {dt, dx, deta, m, lambda, t};

  // Output parameters to file for reference
  // This file also feeds into plot_phi4.C
  params.open("data/params.dat");
  params << t0 << "\n"
         << dt << "\n"
         << dx << "\n"
         << deta << "\n"
         << m << "\n"
         << lambda << "\n"
         << X_GRID << "\n"
         << Y_GRID << "\n"
         << ETA_GRID << "\n"
         << T_STEPS << "\n";
  params.close();

  // Used to calc avg speed. Array holds:
  // [0][i] = location of wave peak on +x-axis at ith step
  // [1][i] = value of rho at center at ith step
  double wavePeak[2][T_STEPS]; 

  // Grids for different quantities
  Mat3DDoub Rho0(X_GRID, Y_GRID, ETA_GRID),
            Rho1(X_GRID, Y_GRID, ETA_GRID),
            Theta0(X_GRID, Y_GRID, ETA_GRID),	// Values for current time step	
            Theta1(X_GRID, Y_GRID, ETA_GRID),	// Values for next time step
            Phi10(X_GRID, Y_GRID, ETA_GRID),	// Values for time step n
            Phi11(X_GRID, Y_GRID, ETA_GRID),	// Values for time step n+1
            Phi12(X_GRID, Y_GRID, ETA_GRID),	// Values for time step n+2
            Phi20(X_GRID, Y_GRID, ETA_GRID),	// Values for time step n
            Phi21(X_GRID, Y_GRID, ETA_GRID),	// Values for time step n+1
            Phi22(X_GRID, Y_GRID, ETA_GRID);	// Values for time step n+2

  // Declare tensors for Tuv and after boosting          
  Tensor T, B;

  // Set initial conditions                      
//  setInitConst(Theta0,prm::theta0);
//  setInitConst(Theta1,prm::theta1);
  // Set initial Rho^2 and take sqrt
  setInitRhos(Rho0,Rho1,p);
  setInitThetas(Theta0,Theta1,Rho0);
  sqrtMat3D(Rho0);
  sqrtMat3D(Rho1);
  // Set initial values of Phi1 and Phi2 from Rho and Theta
  setInitPhis(Phi10,Phi20,Rho0,Theta0);
  setInitPhis(Phi11,Phi21,Rho1,Theta1);
  

  // Renormalization block
  if (prm::renorm)
  {
    double lambdaRe = lambda;
    double number = 9272.;     //initial TotalT00 value 
    for(int var = 0; var < 10000; var++)
    {
      //cout << "error 1" << endl;
      CalcT00(T.T00,Phi10,Phi11,Phi20,Phi21,p);
      //cout << "error 2" << endl;
      double total = TotalT00(T.T00,p);
      //cout << "error 3" << endl;
      //cout << "loop " << var << endl;
      cout << " ************* " << total << endl;
      if (total/number < 1.1 && total/number > .9)
      {
        cout << "lambdaRe : " << lambdaRe << "\t X_GRID : " << X_GRID << endl;
        return 0;
      }
      else if (total/number < 0.9) 
      {
	cout << "lambdaRe " << lambdaRe << endl;
        lambdaRe = lambdaRe + 1.;
        p[4] = lambdaRe;
        setInitRhos(Rho0,Rho1,p);
        setInitThetas(Theta0,Theta1,Rho0);
        sqrtMat3D(Rho0);
        sqrtMat3D(Rho1);
        setInitPhis(Phi10,Phi20,Rho0,Theta0);
        setInitPhis(Phi11,Phi21,Rho1,Theta1);
      }
      else if (total/number > 0.9) 
      {
	cout << "lambdaRe " << lambdaRe << endl;
        lambdaRe = lambdaRe*0.8;
        p[4] = lambdaRe;
        setInitRhos(Rho0,Rho1,p);
        setInitThetas(Theta0,Theta1,Rho0);
        sqrtMat3D(Rho0);
        sqrtMat3D(Rho1);
        setInitPhis(Phi10,Phi20,Rho0,Theta0);
        setInitPhis(Phi11,Phi21,Rho1,Theta1);
      }
    }
    //cout << "2nd return 0" << endl;
    return 0;
  }  // End Renormalization

  // ROOT stuff for plots
  TFile *TuvFile = new TFile("data/Tuv_plots.root", "RECREATE");
  TCanvas *cLat = new TCanvas("cLat","cLat",600,600);  
  TCanvas *cLat4 = new TCanvas("cLat4","cLat4",800,600);  
  TH2D *hRho = new TH2D("Rho","Rho",X_GRID,0.,dx*X_GRID,Y_GRID,0.,dx*Y_GRID);
  TH2D *hTheta = new TH2D("Theta","Theta",X_GRID,0.,dx*X_GRID,Y_GRID,0.,dx*Y_GRID);
  TH2D *hPhi1 = new TH2D("Phi1","Phi1",X_GRID,0.,dx*X_GRID,Y_GRID,0.,dx*Y_GRID);
  TH2D *hPhi2 = new TH2D("Phi2","Phi2",X_GRID,0.,dx*X_GRID,Y_GRID,0.,dx*Y_GRID);
  gStyle->SetOptStat(0);

  // Declaration of histograms
  double tmax = t0 + dt*(T_STEPS-1);
  TH1D *he = new TH1D("he","he",T_STEPS,t0,tmax);
  TH1D *hP = new TH1D("hP","hP",T_STEPS,t0,tmax);
  TH1D *hPT = new TH1D("hPT","hPT",T_STEPS,t0,tmax);
  TH1D *hPL = new TH1D("hPL","hPL",T_STEPS,t0,tmax);
  TH1D *hPTmPL = new TH1D("hPTmPL","hPTmPL",T_STEPS,t0,tmax);
  TH1D *hEOS = new TH1D("hEOS","hEOS",T_STEPS,t0,tmax);
  TH1D *hxSpeed = new TH1D("hxSpeed","hxSpeed",T_STEPS,t0,tmax);
  TH1D *hySpeed = new TH1D("hySpeed","hySpeed",T_STEPS,t0,tmax);
  TH1D *hRho1D = new TH1D("hRho1D","hRho1D",T_STEPS,t0,tmax);
  TH1D *hTheta1D = new TH1D("hTheta1D","hTheta1D",T_STEPS,t0,tmax);
  TH1D *hBe = new TH1D("hBe","hBe",T_STEPS,t0,tmax);
  TH1D *hBP = new TH1D("hBP","hBP",T_STEPS,t0,tmax);
  TH1D *hBPT = new TH1D("hBPT","hBPT",T_STEPS,t0,tmax);
  TH1D *hBPL = new TH1D("hBPL","hBPL",T_STEPS,t0,tmax);
  TH1D *hBPTmPL = new TH1D("hBPTmPL","hBPTmPL",T_STEPS,t0,tmax);
  TH1D *hBEOS = new TH1D("hBEOS","hBEOS",T_STEPS,t0,tmax);
  TH1D *hBSpeed = new TH1D("hBSpeed","hBSpeed",T_STEPS,t0,tmax);
  TH1D *havgT00 = new TH1D("havgT00","havgT00",T_STEPS,t0,tmax);
  TH1D *havgT11 = new TH1D("havgT11","havgT11",T_STEPS,t0,tmax);
  TH1D *havgT22 = new TH1D("havgT22","havgT22",T_STEPS,t0,tmax);
  TH1D *havgT33 = new TH1D("havgT33","havgT33",T_STEPS,t0,tmax);
  TH1D *havgT01 = new TH1D("havgT01","havgT01",T_STEPS,t0,tmax);
  TH1D *havgT02 = new TH1D("havgT02","havgT02",T_STEPS,t0,tmax);
  TH1D *havgT03 = new TH1D("havgT03","havgT03",T_STEPS,t0,tmax);
  TH1D *havgT12 = new TH1D("havgT12","havgT12",T_STEPS,t0,tmax);
  TH1D *havgT13 = new TH1D("havgT13","havgT13",T_STEPS,t0,tmax);
  TH1D *havgT23 = new TH1D("havgT23","havgT23",T_STEPS,t0,tmax);
  TH1D *hT00p = new TH1D("hT00p","hT00p",T_STEPS,t0,tmax);
  TH1D *hT11p = new TH1D("hT11p","hT11p",T_STEPS,t0,tmax);
  TH1D *hT22p = new TH1D("hT22p","hT22p",T_STEPS,t0,tmax);
  TH1D *hT33p = new TH1D("hT33p","hT33p",T_STEPS,t0,tmax);
  TH1D *hT01p = new TH1D("hT01p","hT01p",T_STEPS,t0,tmax);
  TH1D *hT02p = new TH1D("hT02p","hT02p",T_STEPS,t0,tmax);
  TH1D *hT03p = new TH1D("hT03p","hT03p",T_STEPS,t0,tmax);
  TH1D *hT12p = new TH1D("hT12p","hT12p",T_STEPS,t0,tmax);
  TH1D *hT13p = new TH1D("hT13p","hT13p",T_STEPS,t0,tmax);
  TH1D *hT23p = new TH1D("hT23p","hT23p",T_STEPS,t0,tmax);
  
  // Calculate rho = sqrt(phi1^2 + phi2^2) and Theta = atan(phi2/phi1)
  ModPhi(Rho1,Phi11,Phi21);
  Phase(Theta1,Phi11,Phi21);
// DOING THIS MAKES THE CONDENSATE (looks like it) DEPLETE VERY QUICKLY
  // Rho now becomes Rho^2 to plot densities
  //squareMat3D(Rho0);
  //squareMat3D(Rho1);

  // Print out time step 0
  outFile << 0 << " \n";	
  printMat3D(outFile,Rho1);
/*  if (makeGIF) {MakeGIF(Rho0,hRho,cLat); MakeGIF(Rho1,hRho,cLat);}
  if (makeGIF4) 
  {
    MakeGIF4(Rho0,Theta0,Phi10,Phi20,hRho,hTheta,hPhi1,hPhi2,cLat4);
    MakeGIF4(Rho1,Theta1,Phi11,Phi21,hRho,hTheta,hPhi1,hPhi2,cLat4);
  }
*/  
  // Begin time step loop
  for (int iTime = 1; iTime < T_STEPS+1; iTime++) 
  { 
    // Prints out time
    if (iTime%int(T_STEPS*.1) == 0) 
    {
      secsPassed = (clock() - startTime)/CLOCKS_PER_SEC;
      cout << "Calculating... at time step " << iTime << "\t";
      printTime(secsPassed); 
      cout << "\n";
    }

    // Update values
    update(Phi10,Phi11,Phi12,Phi20,Phi21,Phi22,p);

    ModPhi(Rho0,Phi11,Phi21);
    ModPhi(Rho1,Phi12,Phi22);
    Phase(Theta0,Phi11,Phi21);
    Phase(Theta1,Phi12,Phi22);

    wavePeak[0][iTime] = findMax(Rho1);
    wavePeak[1][iTime] = Rho1[int(X_GRID/2)][int(Y_GRID/2)][int(ETA_GRID/2)];
//    cout << wavePeak[0][iTime] << endl;
//    cout << wavePeak[1][iTime] << endl;

    // Calculates Tuv 
    CalcTuv(T,Phi11,Phi12,Phi21,Phi22,p);

    // Fills histograms
    he->Fill(t,avgMat3D(T.T00));
    hP->Fill(t,avgMat3D(T.T11)+avgMat3D(T.T22)+t*t*avgMat3D(T.T33));
    hPT->Fill(t,avgMat3D(T.T11)+avgMat3D(T.T22));
    hPL->Fill(t,t*t*avgMat3D(T.T33));
    hPTmPL->Fill(t,avgMat3D(T.T11)+avgMat3D(T.T22)-t*t*avgMat3D(T.T33));
    hEOS->Fill(t,avgMat3D(T.T00)-avgMat3D(T.T11)-avgMat3D(T.T22)-t*t*avgMat3D(T.T33));
    havgT00->Fill(t,avgMat3D(T.T00));
    havgT11->Fill(t,avgMat3D(T.T11));
    havgT22->Fill(t,avgMat3D(T.T22));
    havgT33->Fill(t,avgMat3D(T.T33));
    havgT01->Fill(t,avgMat3D(T.T01));
    havgT02->Fill(t,avgMat3D(T.T02));
    havgT03->Fill(t,avgMat3D(T.T03));
    havgT12->Fill(t,avgMat3D(T.T12));
    havgT13->Fill(t,avgMat3D(T.T13));
    havgT23->Fill(t,avgMat3D(T.T23));
    hT00p->Fill(t,T.T00[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    hT11p->Fill(t,T.T11[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    hT22p->Fill(t,T.T22[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    hT33p->Fill(t,T.T33[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    hT01p->Fill(t,T.T01[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    hT02p->Fill(t,T.T02[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    hT03p->Fill(t,T.T03[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    hT12p->Fill(t,T.T12[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    hT13p->Fill(t,T.T13[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    hT23p->Fill(t,T.T23[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);

    double vx = xSpeed(Theta0,Theta1,int(.75*X_GRID),int(.5*Y_GRID),int(.5*ETA_GRID),p);
    hxSpeed->Fill(t,vx);
    double vy = ySpeed(Theta0,Theta1,int(.5*X_GRID),int(.75*Y_GRID),int(.5*ETA_GRID),p);
    hySpeed->Fill(t,vy);
    hRho1D->Fill(t,Rho1[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    hTheta1D->Fill(t,Theta1[int(.75*X_GRID)][int(.5*Y_GRID)][int(.5*ETA_GRID)]);
    
    //cout << avgSpeed(Rho0,Rho1) << endl;

    // Boost Tuv to the rest frame of a fluid element
    BoostTuv(B,T,vx,p);

    hBe->Fill(t,avgMat3D(B.T00));
    hBP->Fill(t,avgMat3D(B.T11)+avgMat3D(B.T22)+t*t*avgMat3D(B.T33));
    hBPT->Fill(t,avgMat3D(B.T11)+avgMat3D(B.T22));
    hBPL->Fill(t,t*t*avgMat3D(B.T33));
    hBPTmPL->Fill(t,avgMat3D(B.T11)+avgMat3D(B.T22)-t*t*avgMat3D(B.T33));
    hBEOS->Fill(t,avgMat3D(B.T00)-avgMat3D(B.T11)-avgMat3D(B.T22)-t*t*avgMat3D(B.T33));

    // Rho now becomes Rho^2 to plot densities
    //squareMat3D(Rho1);

    // Prints out 10 time steps. Note that total output will have time
    //   steps 0, 100, 200, ..., 1000 (for T_STEPS == 1000), but these
    //   are at times 1, 101, 201, ..., 1001 (for t0 == 1 && dt == 1). 
    if ((iTime+2)%int(T_STEPS*.1) == 0) 
    {
      outFile << iTime+2 << " \n";
      printMat3D(outFile,Rho1);
    }
    
    // Make a GIF (200 time steps long, in total)
    if (makeGIF) 
    {
      if (T_STEPS < 200) {MakeGIF(Rho1,hRho,cLat);}
      else if((iTime+2)%int(T_STEPS/200) == 0) {MakeGIF(Rho1,hRho,cLat);}
    }    
    if (makeGIF4) 
    {
      if (T_STEPS < 200) {MakeGIF4(Rho0,Theta0,Phi10,Phi20,hRho,hTheta,hPhi1,hPhi2,cLat4);}
      else if((iTime+2)%int(T_STEPS/200) == 0) 
      {
        MakeGIF4(Rho1,Theta1,Phi12,Phi22,hRho,hTheta,hPhi1,hPhi2,cLat4);
      }
    }    
        
    // Changes the "next" time step into the "current" one
    Assign(Phi10,Phi11);
    Assign(Phi11,Phi12);
    Assign(Phi20,Phi21);
    Assign(Phi21,Phi22);

    // Advance time
    t += dt;
    p[5] = t;
  } 	// End time step loop
  
  avgSpeed(wavePeak,2);

  //===============================================================//
  // ROOT stuff for plots                                          //
  //===============================================================//
  TCanvas *ce = new TCanvas("ce","ce");
  TCanvas *cP = new TCanvas("cP","cP");
  TCanvas *cPT = new TCanvas("cPT","cPT");
  TCanvas *cPL = new TCanvas("cPL","cPL");
  TCanvas *cPTmPL = new TCanvas("cPTmPL","cPTmPL");
  TCanvas *cEOS = new TCanvas("cEOS","cEOS");
  TCanvas *cSpeed = new TCanvas("cSpeed","cSpeed");

  he->SetMarkerStyle(7);
  hP->SetMarkerStyle(7);
  hPT->SetMarkerStyle(7);
  hPL->SetMarkerStyle(7);
  hPTmPL->SetMarkerStyle(7);
  hEOS->SetMarkerStyle(7);
  hxSpeed->SetMarkerStyle(7);
  hySpeed->SetMarkerStyle(8);
  hRho1D->SetMarkerStyle(7);
  hTheta1D->SetMarkerStyle(5);
  hBe->SetMarkerStyle(7);
  hBP->SetMarkerStyle(7);
  hBPT->SetMarkerStyle(7);
  hBPL->SetMarkerStyle(7);
  hBPTmPL->SetMarkerStyle(7);
  hBEOS->SetMarkerStyle(7);
  hBSpeed->SetMarkerStyle(7);
  havgT00->SetMarkerStyle(7);
  havgT11->SetMarkerStyle(7);
  havgT22->SetMarkerStyle(7);
  havgT33->SetMarkerStyle(7);
  havgT01->SetMarkerStyle(7);
  havgT02->SetMarkerStyle(7);
  havgT03->SetMarkerStyle(7);
  havgT12->SetMarkerStyle(7);
  havgT13->SetMarkerStyle(7);
  havgT23->SetMarkerStyle(7);
  hT00p->SetMarkerStyle(7);
  hT11p->SetMarkerStyle(7);
  hT22p->SetMarkerStyle(7);
  hT33p->SetMarkerStyle(7);
  hT01p->SetMarkerStyle(7);
  hT02p->SetMarkerStyle(7);
  hT03p->SetMarkerStyle(7);
  hT12p->SetMarkerStyle(7);
  hT13p->SetMarkerStyle(7);
  hT23p->SetMarkerStyle(7);

  hySpeed->SetMarkerColor(kRed);
  hRho1D->SetMarkerColor(kBlue);
  hTheta1D->SetMarkerColor(kRed);

  hxSpeed->GetYaxis()->SetRangeUser(prm::sMin,prm::sMax);  
  
  ce->cd();
  gPad->SetLogy();  
  he->Draw("p");

  cP->cd();
  gPad->SetLogy();
  hP->Draw("p");

  cPT->cd();
  gPad->SetLogy();
  hPT->Draw("p");

  cPL->cd();
  gPad->SetLogy();
  hPL->Draw("p");
  
  cPTmPL->cd();
  gPad->SetLogy();
  hPTmPL->Draw("p");

  cEOS->cd();
  //gPad->SetLogy();
  hEOS->Draw("p");

  cSpeed->cd();
  hxSpeed->Draw("p");
  cSpeed->Update();
  double rightmax = 1.1*hRho1D->GetMaximum();
  double scale = gPad->GetUymax()/rightmax;
  hRho1D->Scale(scale);
  hRho1D->Draw("pSAME");
//  hTheta1D->Draw("pSAME");
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                            gPad->GetUxmax(),gPad->GetUymax(),
                            -rightmax,rightmax,510,"+L");
  axis->SetLineColor(kBlue);
  axis->SetLabelColor(kBlue);
  axis->SetLabelFont(40);
  axis->Draw();                            
  hySpeed->Draw("pSAME");
  hxSpeed->Draw("pSAME");
  
  TLegend *lSpeed = new TLegend(0.20,0.20,0.50,0.40,"");
  lSpeed->AddEntry(hxSpeed,"v_{x} at (x,0)","p");
  lSpeed->AddEntry(hySpeed,"v_{y} at (0,y)","p");
  lSpeed->AddEntry(hRho1D,"#rho at (x,0)","p");
  lSpeed->Draw();


  // Saves plots to .root file
  he->Write();
  hP->Write();
  hPT->Write();
  hPL->Write();
  hPTmPL->Write();
  hEOS->Write();
  hxSpeed->Write();
  hySpeed->Write();
  hBe->Write();
  hBP->Write();
  hBPT->Write();
  hBPL->Write();
  hBPTmPL->Write();
  hBEOS->Write();
  hBSpeed->Write();
  havgT00->Write();
  havgT11->Write();
  havgT22->Write();
  havgT33->Write();
  havgT01->Write();
  havgT02->Write();
  havgT03->Write();
  havgT12->Write();
  havgT13->Write();
  havgT23->Write();
  hT00p->Write();
  hT11p->Write();
  hT22p->Write();
  hT33p->Write();
  hT01p->Write();
  hT02p->Write();
  hT03p->Write();
  hT12p->Write();
  hT13p->Write();
  hT23p->Write();

  // Makes .png files for plots
  ce->Print("Plots/he.png","png");
  cP->Print("Plots/hP.png","png");
  cPT->Print("Plots/hPT.png","png");
  cPL->Print("Plots/hPL.png","png");
  cPTmPL->Print("Plots/hPTmPL.png","png");
  cEOS->Print("Plots/hEOS.png","png");
  //cSpeed->Print("Plots/hSpeed.png","png"); CAUSING CRASH. TRY COUT RHO NEAR EDGE
  //===============================================================//

  // Close up files
  TuvFile->Close();
  outFile.close();

  return 0;
}


//=========================== Functions ===========================//

void update(const Mat3DDoub &Phi10, const Mat3DDoub &Phi11, Mat3DDoub &Phi12, 
            const Mat3DDoub &Phi20, const Mat3DDoub &Phi21, Mat3DDoub &Phi22, const double p[])
{
  int X_GRID = Phi10.dim1(),
      Y_GRID = Phi10.dim2(),
      ETA_GRID = Phi10.dim3();

  // Parameter array: p = [dt, dx, deta, m, lambda, tn]
  double dt = p[0],
         dx = p[1],
         deta = p[2],
         m = p[3],
         lambda = p[4],
         tn = p[5];

  for (int iETA = 0; iETA < ETA_GRID; iETA++)
  {
    int iETAm1 = (ETA_GRID+iETA-1)%(ETA_GRID),
        iETAp1 = (ETA_GRID+iETA+1)%(ETA_GRID);

    for (int iY = 0; iY < Y_GRID; iY++)
    {
      int iYm1 = (Y_GRID+iY-1)%(Y_GRID),
          iYp1 = (Y_GRID+iY+1)%(Y_GRID);

      for (int iX = 0; iX < X_GRID; iX++)
      { 
        int iXm1 = (X_GRID+iX-1)%(X_GRID),
            iXp1 = (X_GRID+iX+1)%(X_GRID);
            
        Phi12[iX][iY][iETA] = 2*Phi11[iX][iY][iETA] - Phi10[iX][iY][iETA] 
                         - (dt/tn)*(Phi11[iX][iY][iETA] - Phi10[iX][iY][iETA]) 
                         + (dt*dt/(tn*tn*deta*deta))*(Phi11[iX][iY][iETAp1]
                                                     - 2*Phi11[iX][iY][iETA]
                                                     + Phi11[iX][iY][iETAm1])
                         + (dt*dt/(2*dx*dx))*(Phi11[iXp1][iY][iETA] + Phi11[iXm1][iY][iETA]
                                            + Phi11[iX][iYp1][iETA] + Phi11[iX][iYm1][iETA]
                                            + 0.5*(Phi11[iXp1][iYp1][iETA] + Phi11[iXm1][iYp1][iETA]
                                                  + Phi11[iXp1][iYm1][iETA] + Phi11[iXm1][iYm1][iETA])
                                            - 6*Phi11[iX][iY][iETA])
                         - dt*dt*m*m*Phi11[iX][iY][iETA]
                         - 2*dt*dt*lambda*(pow(Phi11[iX][iY][iETA],3)
                                    + Phi11[iX][iY][iETA]*Phi21[iX][iY][iETA]*Phi21[iX][iY][iETA]);

        Phi22[iX][iY][iETA] = 2*Phi21[iX][iY][iETA] - Phi20[iX][iY][iETA] 
                         - (dt/tn)*(Phi21[iX][iY][iETA] - Phi20[iX][iY][iETA]) 
                         + (dt*dt/(tn*tn*deta*deta))*(Phi21[iX][iY][iETAp1]
                                                     - 2*Phi21[iX][iY][iETA]
                                                     + Phi21[iX][iY][iETAm1])
                         + (dt*dt/(2*dx*dx))*(Phi21[iXp1][iY][iETA] + Phi21[iXm1][iY][iETA]
                                            + Phi21[iX][iYp1][iETA] + Phi21[iX][iYm1][iETA]
                                            + 0.5*(Phi21[iXp1][iYp1][iETA] + Phi21[iXm1][iYp1][iETA]
                                                  + Phi21[iXp1][iYm1][iETA] + Phi21[iXm1][iYm1][iETA])
                                            - 6*Phi21[iX][iY][iETA])
                         - dt*dt*m*m*Phi21[iX][iY][iETA]
                         - 2*dt*dt*lambda*(pow(Phi21[iX][iY][iETA],3)
                                    + Phi11[iX][iY][iETA]*Phi11[iX][iY][iETA]*Phi21[iX][iY][iETA]);
      }
    }
  }  
  return;
}

void MakeGIF(const Mat3DDoub &M, TH2D *hLat, TCanvas *cLat)
{
  int X_GRID = M.dim1(),
      Y_GRID = M.dim2(),
      ETA_GRID = M.dim3();
  
  for (int iY = 0; iY < Y_GRID; iY++) 
  {
    for (int iX = 0; iX < X_GRID; iX++) 
    {      
      double temp = M[iX][iY][(ETA_GRID-1)/2];  
      hLat->SetBinContent(iX+1,iY+1,temp);
    }
  }
//  hLat->GetZaxis()->SetRangeUser(prm::zMin,prm::zMax);
  cLat->cd();
//  gPad->SetTheta(60);
  hLat->DrawCopy(prm::plotType);
  cLat->Update();
  cLat->Print("Plots/rho_temp.gif+");

  return;
}

void MakeGIF4(const Mat3DDoub &Rho,const Mat3DDoub &Theta,const Mat3DDoub &Phi1,const Mat3DDoub &Phi2,
              TH2D *hRho,TH2D *hTheta,TH2D *hPhi1,TH2D *hPhi2,TCanvas *cLat4)
{
  int X_GRID = Rho.dim1(),
      Y_GRID = Rho.dim2(),
      ETA_GRID = Rho.dim3();
  
  for (int iY = 0; iY < Y_GRID; iY++) 
  {
    for (int iX = 0; iX < X_GRID; iX++) 
    {      
      double rho = Rho[iX][iY][(ETA_GRID-1)/2];  
      hRho->SetBinContent(iX+1,iY+1,rho);
      double theta = Theta[iX][iY][(ETA_GRID-1)/2];  
      hTheta->SetBinContent(iX+1,iY+1,theta);
      double phi1 = Phi1[iX][iY][(ETA_GRID-1)/2];  
      hPhi1->SetBinContent(iX+1,iY+1,phi1);
      double phi2 = Phi2[iX][iY][(ETA_GRID-1)/2];  
      hPhi2->SetBinContent(iX+1,iY+1,phi2);
    }
  }
  
  hRho->GetZaxis()->SetRangeUser(prm::zMin,prm::zMax);
  hTheta->GetZaxis()->SetRangeUser(-prm::PI,prm::PI);
  hPhi1->GetZaxis()->SetRangeUser(-prm::zMax,prm::zMax);
  hPhi2->GetZaxis()->SetRangeUser(-prm::zMax,prm::zMax);

  cLat4->Clear();
  cLat4->Divide(2,2);
  cLat4->cd(1);
  gPad->SetTheta(60);
  hRho->DrawCopy(prm::plotType);
  cLat4->cd(2);
  gPad->SetTheta(60);
  hTheta->DrawCopy(prm::plotType);
  cLat4->cd(3);
  gPad->SetTheta(60);
  hPhi1->DrawCopy(prm::plotType);
  cLat4->cd(4);
  gPad->SetTheta(60);
  hPhi2->DrawCopy(prm::plotType);

  cLat4->Update();
  cLat4->Print("Plots/phi4_temp.gif+");
}

