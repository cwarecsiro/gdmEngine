//
// clsHashtab.h
//
#ifndef __CLSHASHTAB_H__
#define __CLSHASHTAB_H__

#include <iostream> 
#include <cstdlib> 
#include <cstring> 
#include <iomanip> 
//#define SIZE_KEY       64 
#define SIZE_KEY       128 
//#define SIZE_CODE      64 
#define SIZE_CODE      128 
#define DEFAULT_TABLESIZE    101 

struct NODE 
{ 
   NODE (const char* Key1 = "\0", const char* Code1 = "\0")
   { 
      strcpy(Key, Key1); 
      strcpy(Code, Code1); 
      next = NULL; 
   } 
   char Key[SIZE_KEY]; 
   char Code[SIZE_CODE]; 
   NODE  *next; 
}; 
  


class Hashtable 
{ 
   private: 
      int table_size; 
      NODE** table; 
      int size; 
      long hashString(char* Key); 
      NODE* find(char* Key); 
      NODE* current_entry; 
      int current_index; 
   public: 
      Hashtable(int T = DEFAULT_TABLESIZE);//constructor 
      virtual ~Hashtable();//destructor 
      bool put(NODE *); 
      bool get(NODE *); 
      bool contains(char* Key); 
      bool remove(char* Key); 
      void removeAll(); 
      int getSize(); 
      void initIterator(); 
      bool hasNext(); 
      void getNextKey(char* Key); 
      friend void disp(NODE *); 
}; 



#endif // __CLSHASHTAB_H__