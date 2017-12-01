void energy_linearity_zcut_1_mipcut0.5()
{
//=========Macro generated from canvas: c_energy_linearity/c_energy_linearity
//=========  (Thu Aug 10 13:55:10 2017) by ROOT version 6.11/01
   TCanvas *c_energy_linearity = new TCanvas("c_energy_linearity", "c_energy_linearity",757,46,800,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c_energy_linearity->Range(-0.6465822,-182.2785,6.644557,956.962);
   c_energy_linearity->SetFillColor(0);
   c_energy_linearity->SetBorderMode(0);
   c_energy_linearity->SetBorderSize(2);
   c_energy_linearity->SetTickx(1);
   c_energy_linearity->SetTicky(1);
   c_energy_linearity->SetLeftMargin(0.16);
   c_energy_linearity->SetRightMargin(0.05);
   c_energy_linearity->SetTopMargin(0.05);
   c_energy_linearity->SetBottomMargin(0.16);
   c_energy_linearity->SetFrameBorderMode(0);
   c_energy_linearity->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1001[5] = {
   1,
   2,
   3,
   4,
   5.8};
   Double_t Graph0_fy1001[5] = {
   78.05215,
   165.7478,
   242.3306,
   325.5464,
   448.5272};
   Double_t Graph0_fex1001[5] = {
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fey1001[5] = {
   0.478296,
   1.328766,
   2.274384,
   3.519734,
   6.184533};
   TGraphErrors *gre = new TGraphErrors(5,Graph0_fx1001,Graph0_fy1001,Graph0_fex1001,Graph0_fey1001);
   gre->SetName("Graph0");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(4);
   gre->SetMarkerSize(1.2);
   
   TH1F *Graph_Graph1001 = new TH1F("Graph_Graph1001","Graph",100,0.52,6.28);
   Graph_Graph1001->SetMinimum(0);
   Graph_Graph1001->SetMaximum(900);
   Graph_Graph1001->SetDirectory(0);
   Graph_Graph1001->SetStats(0);
   Graph_Graph1001->SetLineWidth(2);
   Graph_Graph1001->SetMarkerStyle(20);
   Graph_Graph1001->SetMarkerSize(1.2);
   Graph_Graph1001->GetXaxis()->SetTitle("E^{beam}/GeV");
   Graph_Graph1001->GetXaxis()->SetLabelFont(42);
   Graph_Graph1001->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1001->GetXaxis()->SetTitleFont(42);
   Graph_Graph1001->GetYaxis()->SetTitle("E^{raw}/MIP");
   Graph_Graph1001->GetYaxis()->SetLabelFont(42);
   Graph_Graph1001->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1001->GetYaxis()->SetTitleFont(42);
   Graph_Graph1001->GetZaxis()->SetLabelFont(42);
   Graph_Graph1001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1001);
   
   gre->Draw("alp");
   
   Double_t Graph1_fx1002[5] = {
   1,
   2,
   4,
   5,
   5.8};
   Double_t Graph1_fy1002[5] = {
   80.73043,
   170.4848,
   334.3313,
   402.386,
   453.8424};
   Double_t Graph1_fex1002[5] = {
   0,
   0,
   0,
   0,
   0};
   Double_t Graph1_fey1002[5] = {
   0.4085916,
   0.8814578,
   2.529535,
   2.043969,
   3.596176};
   gre = new TGraphErrors(5,Graph1_fx1002,Graph1_fy1002,Graph1_fex1002,Graph1_fey1002);
   gre->SetName("Graph1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineStyle(2);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   
   TH1F *Graph_Graph1002 = new TH1F("Graph_Graph1002","Graph",100,0.52,6.28);
   Graph_Graph1002->SetMinimum(42.61016);
   Graph_Graph1002->SetMaximum(495.1503);
   Graph_Graph1002->SetDirectory(0);
   Graph_Graph1002->SetStats(0);
   Graph_Graph1002->SetLineWidth(2);
   Graph_Graph1002->SetMarkerStyle(20);
   Graph_Graph1002->SetMarkerSize(1.2);
   Graph_Graph1002->GetXaxis()->SetLabelFont(42);
   Graph_Graph1002->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1002->GetXaxis()->SetTitleFont(42);
   Graph_Graph1002->GetYaxis()->SetLabelFont(42);
   Graph_Graph1002->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1002->GetYaxis()->SetTitleFont(42);
   Graph_Graph1002->GetZaxis()->SetLabelFont(42);
   Graph_Graph1002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1002);
   
   gre->Draw("lp");
   
   Double_t Graph2_fx1003[5] = {
   1,
   2,
   4,
   5,
   5.8};
   Double_t Graph2_fy1003[5] = {
   95.981,
   205.1596,
   415.7359,
   516.13,
   593.4313};
   Double_t Graph2_fex1003[5] = {
   0,
   0,
   0,
   0,
   0};
   Double_t Graph2_fey1003[5] = {
   0.7923781,
   1.278759,
   3.054484,
   3.382752,
   7.016038};
   gre = new TGraphErrors(5,Graph2_fx1003,Graph2_fy1003,Graph2_fex1003,Graph2_fey1003);
   gre->SetName("Graph2");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineColor(8);
   gre->SetLineWidth(2);
   gre->SetMarkerColor(8);
   gre->SetMarkerStyle(22);
   gre->SetMarkerSize(1.2);
   
   TH1F *Graph_Graph1003 = new TH1F("Graph_Graph1003","Graph",100,0.52,6.28);
   Graph_Graph1003->SetMinimum(44.66275);
   Graph_Graph1003->SetMaximum(650.9732);
   Graph_Graph1003->SetDirectory(0);
   Graph_Graph1003->SetStats(0);
   Graph_Graph1003->SetLineWidth(2);
   Graph_Graph1003->SetMarkerStyle(20);
   Graph_Graph1003->SetMarkerSize(1.2);
   Graph_Graph1003->GetXaxis()->SetLabelFont(42);
   Graph_Graph1003->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1003->GetXaxis()->SetTitleFont(42);
   Graph_Graph1003->GetYaxis()->SetLabelFont(42);
   Graph_Graph1003->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1003->GetYaxis()->SetTitleFont(42);
   Graph_Graph1003->GetZaxis()->SetLabelFont(42);
   Graph_Graph1003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1003);
   
   gre->Draw("lp");
   
   Double_t Graph3_fx1004[5] = {
   1,
   2,
   4,
   5,
   5.8};
   Double_t Graph3_fy1004[5] = {
   112.3267,
   246.7995,
   512.8868,
   641.4205,
   725.0539};
   Double_t Graph3_fex1004[5] = {
   0,
   0,
   0,
   0,
   0};
   Double_t Graph3_fey1004[5] = {
   0.8329203,
   1.800342,
   2.60877,
   4.379014,
   6.916748};
   gre = new TGraphErrors(5,Graph3_fx1004,Graph3_fy1004,Graph3_fex1004,Graph3_fey1004);
   gre->SetName("Graph3");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineColor(4);
   gre->SetLineWidth(2);
   gre->SetMarkerColor(4);
   gre->SetMarkerStyle(23);
   gre->SetMarkerSize(1.2);
   
   TH1F *Graph_Graph1004 = new TH1F("Graph_Graph1004","Graph",100,0.52,6.28);
   Graph_Graph1004->SetMinimum(49.44608);
   Graph_Graph1004->SetMaximum(794.0184);
   Graph_Graph1004->SetDirectory(0);
   Graph_Graph1004->SetStats(0);
   Graph_Graph1004->SetLineWidth(2);
   Graph_Graph1004->SetMarkerStyle(20);
   Graph_Graph1004->SetMarkerSize(1.2);
   Graph_Graph1004->GetXaxis()->SetLabelFont(42);
   Graph_Graph1004->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1004->GetXaxis()->SetTitleFont(42);
   Graph_Graph1004->GetYaxis()->SetLabelFont(42);
   Graph_Graph1004->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1004->GetYaxis()->SetTitleFont(42);
   Graph_Graph1004->GetZaxis()->SetLabelFont(42);
   Graph_Graph1004->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1004);
   
   gre->Draw("lp");
   
   TLegend *leg = new TLegend(0.2,0.7,0.7,0.9,NULL,"brNDC");
   leg->SetTextFont(62);
   leg->SetLineColor(0);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("NULL","SiW-ECAL: reconstructed energy linearity","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph0","Wafer 3, W-configuration 1","lp");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(4);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph1","Wafer 4, W-configuration 1","lp");
   entry->SetLineColor(1);
   entry->SetLineStyle(2);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph2","Wafer 4, W-configuration 2","lp");
   entry->SetLineColor(8);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(8);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph3","Wafer 4, W-configuration 3","lp");
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(62);
   leg->Draw();
   c_energy_linearity->Modified();
   c_energy_linearity->cd();
   c_energy_linearity->SetSelected(c_energy_linearity);
}
