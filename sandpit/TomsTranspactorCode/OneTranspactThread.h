//---------------------------------------------------------------------------

#ifndef OneTranspactThreadH
#define OneTranspactThreadH
//---------------------------------------------------------------------------
#include <System.Classes.hpp>
//---------------------------------------------------------------------------
class TOneTranspactThread : public TThread
{
private:
protected:
	void __fastcall Execute();
public:
	__fastcall TOneTranspactThread(bool CreateSuspended);
	void __fastcall UpdateProgress(void);
	void __fastcall UpdateProgress2(void);
	void __fastcall UpdateProgressCaption(void) ;
  //	void __fastcall LockButtons(void) ;
 //	void __fastcall UpdateSelected(void);
	int i_scen;
	int i_progress;
	int i_progress2;
	int i_lock;
	AnsiString s_caption;
};
//---------------------------------------------------------------------------
#endif
