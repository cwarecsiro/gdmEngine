//---------------------------------------------------------------------------
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <FileCtrl.hpp>
#include <vcl.h>
#pragma hdrstop

#include "defines.h"
#include "esri_headers.h"
#include "Splash.h"
#include "CompactorFunctions.h"
#include "BulkCompactor.h"
#include "Transpactor.h"
#include "OneTranspactThread.h"
#pragma package(smart_init)
//---------------------------------------------------------------------------

//   Important: Methods and properties of objects in VCL can only be
//   used in a method called using Synchronize, for example:
//
//      Synchronize(&UpdateCaption);
//
//   where UpdateCaption could look like:
//
//      void __fastcall TOneTranspactThread::UpdateCaption()
//      {
//        Form1->Caption = "Updated in a thread";
//      }
//---------------------------------------------------------------------------

__fastcall TOneTranspactThread::TOneTranspactThread(bool CreateSuspended)
	: TThread(CreateSuspended)
{
}
//---------------------------------------------------------------------------
void __fastcall TOneTranspactThread::Execute()
{
AnsiString s_path,s_prefix,s_mess,s_bpath;
AnsiString s1,s2,s3;
string s_layer_files[50];
int i_layer, n_layers;
TCursor Save_Cursor;
string s_name;
AnsiString s_gdp_name,s_realm_path,s_res_path,s_env_path;
AnsiString s_clim_path,s_subs_path,s_land_path,s_resid_path;
int n_biomes, n_scenarios,n_sgrids,n_bgrids;
int i_biome, i_scenario, i_grid;
int i_stack,n_stacks;
int i_pos;
float f_pos;
int i_length;
int i_calcs;
double d_bare_scale[14];
string s_buffer;
string i_buffer;
ifstream InFileStream;
TTranspactor* Transpactor;

string s_line;   //for csv file read
string s_word;

	//---- Place thread code here ----
i_lock=1; //LOCK
//Synchronize(LockButtons);
	  i_calcs=0;
  //Extract from GUI
  n_biomes=BCMain->BiomeListBox->Items->Count;
  n_scenarios=BCMain->STranListBox->Items->Count;
  s_gdp_name=BCMain->GDMPEdit->Text.c_str();
  s_realm_path= BCMain->FolderEdit->Text;
  s_env_path= BCMain->EnvFolderEdit->Text;
  s1=BCMain->PrefixEdit->Text;
  // Get Realm code
  s_path=s_realm_path;
  i_length=s_path.Length();
  // PARSING LAST TWO LETTERS FOR REALM
  s2= s_path[i_length-2];
  s3= s_path[i_length-1];

  s_prefix=s1+s2+s3;

  n_layers= 18;
  n_stacks=n_biomes*n_scenarios;
  // Fire up a Transpactor (how else would you do this?)
  Transpactor=new TTranspactor(1);
  Transpactor->SetNLayers(18);
  Transpactor->SetNRows(StrToInt(BCMain->NRowsEdit->Text));
  Transpactor->SetNCols(StrToInt(BCMain->NColsEdit->Text));

  //Load in the Realm Level BARE parameters
  s_bpath=s_path+BCMain->BGEdit->Text;
  InFileStream.open(s_bpath.c_str());
  // Check file is present
  if(InFileStream.is_open())
	{
	//Load this here, and we can pass in the bare params as required in the main biome for loop
	//Use the order in the file which numbers 0 to n_biomes, and is not the biome number
	//First the header line
	s_line="";
	getline(InFileStream, s_line);
	//Now the rest
	for(i_biome=0;i_biome<n_biomes;i_biome++)
	  {
	  getline(InFileStream, s_line);
	  stringstream strstr(s_line);
	  string word = "";
	  getline(strstr,word, ',');
	  strstr>>d_bare_scale[i_biome];
	  }  // end for biome
	InFileStream.close();
	}

  //Loop over biomes and scenarios  .. note i_biome IS NOT the biome number..it is the ordered biomes in the realm
  for(i_biome=0;i_biome<n_biomes;i_biome++)
	{
	// Reset flags to zero
	Transpactor->Refresh();
	//Load in the GDM Parameters file
	s_path=s_realm_path+BCMain->BiomeListBox->Items->Strings[i_biome]+"\\"+s_gdp_name;
	Transpactor->LoadGlobalGDMParameters(s_path.c_str());
	// This sets the spline values in Transpactor, but we need to set up the extrapolation
	Transpactor->SetExtrapolation(1); //minimum linear extrapolation
	for(i_layer=1;i_layer<18;i_layer++)
	  {
	  if(!Transpactor->IsZero(i_layer))
		Transpactor->SetUpExtrapolation(i_layer);
	  }

	// And the bare scale
	Transpactor->SetBareScale(d_bare_scale[i_biome]);
	Transpactor->SetBare(0);
	// just set the RESID grids to zero anyway
	Transpactor->SetZero(5);
	Transpactor->SetZero(6);
	Transpactor->SetZero(7);
	// Now we have the correct values for bare and the grids to be transformed
	//Set up all the file paths
	s_land_path=s_env_path+"LAND\\";
	s_subs_path=s_env_path+"SUBS\\";
	s_resid_path=s_env_path+"RESID\\";

	for(i_scenario=0;i_scenario<n_scenarios;i_scenario++)
	  {
	  //Finish file path set up
	  s_clim_path=s_env_path+BCMain->STranListBox->Items->Strings[i_scenario]+"\\";
	  Transpactor->SetUpGlobFilePaths(s_clim_path.c_str(),
									  s_land_path.c_str(),
									  s_subs_path.c_str(),
									  s_resid_path.c_str());
	  i_calcs++;
	  f_pos=100*(float)i_calcs/(float)n_stacks;
	  i_progress2=(int)f_pos;
	  Synchronize(UpdateProgress2);

	  //Auto results name in MCG folder
	  s_res_path=BCMain->OutFolderEdit->Text+s_prefix+"_"+BCMain->STranListBox->Items->Strings[i_scenario]+"_"+BCMain->BiomeListBox->Items->Strings[i_biome]+".mcg";
	   s_caption="Writing to "+s_res_path;
	   Synchronize(UpdateProgressCaption);
	   if(!Transpactor->TransBinaryGridsToFloat(s_res_path.c_str()))
		 {
//		 s_mess= "Faulty inputs: biome "+BCMain->BiomeListBox->Items->Strings[i_biome]+ " scen "+ BCMain->ScenListBox->Items->Strings[i_scenario];
		 MessageDlg(s_mess, mtError, TMsgDlgButtons() << mbOK << mbCancel,0);
		 }

	 } // END FOR SCENARIO
   } // END FOR BIOMES
	i_progress2=100;
	Synchronize(UpdateProgress2);
	s_caption= "Done";
	Synchronize(UpdateProgressCaption);
	      i_lock=0; //UNLOCK
//Synchronize(LockButtons);
}
//---------------------------------------------------------------------------

/**************************************************************************
* GUI STUFF
**************************************************************************/
void __fastcall TOneTranspactThread::UpdateProgress(void)
{
  BCMain->ProgressBar1->Position=i_progress;

}
void __fastcall TOneTranspactThread::UpdateProgress2(void)
{
  BCMain->ProgressBar2->Position=i_progress2;

}
void __fastcall TOneTranspactThread::UpdateProgressCaption(void)
{
  BCMain->StatusBar1->Panels->Items[0]->Text=s_caption;

}

