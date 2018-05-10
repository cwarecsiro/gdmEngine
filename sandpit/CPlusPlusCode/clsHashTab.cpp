//
// clsHashTab.cpp
//
#include "stdafx.h"
#include "clsHashTab.h"
#include "Message.h"

  
Hashtable::Hashtable(int T) 
{ 
   size = 0; 
   table_size = T; 
   table = new NODE*[table_size]; 
   for(int i=0; i<table_size; i++) 
   { 
      table[i] = NULL; 
   } 
} 
  
Hashtable::~Hashtable() 
{ 
   removeAll(); 
   delete[] table; 
} 
  
void disp(NODE *N1) 
{ 
	char qq[64];
	sprintf( qq, "Key: <%s>   Code: <%s>", N1->Key, N1->Key);
	Message(qq);
} 
  
bool Hashtable::put(NODE *N) 
{//start put 
   if(find(N->Key) != NULL) 
   { 
      return false;
   } 
   NODE* entry = new NODE(N->Key, N->Code); 
   int bucket = hashString(N->Key); 
   entry->next = table[bucket]; 
   table[bucket] = entry; 
   size++; 
   return true; 
}//end put 
  
  
bool Hashtable::get(NODE* N) 
{//start get 
   NODE* temp = find(N->Key); 
   if(temp == NULL) 
   { 
      N->Code[0] = '\0'; 
      return false; 
   } 
   else
   { 
      strcpy(N->Code, temp->Code); 
      return true; 
   } 
}//end get 
  
bool Hashtable::contains(char* Key) 
{//start contains 
   if(find(Key) == NULL) 
   { 
      return false; 
   } 
   else
   { 
      return true; 
   } 
}//end contains 
  
  
bool Hashtable::remove(char* Key) 
{//start remove 
   int bucket = hashString(Key); 
   NODE* temp = table[bucket]; 
   if(temp == NULL) 
   { 
      return false; 
   } 
   else if(strcmp(Key, temp->Key) == 0) 
   { 
      table[bucket] = temp->next; 
      delete temp; 
      size--; 
      return true; 
   } 
   else
   { 
      NODE* temp_next = temp->next; 
      while(temp_next != NULL) 
      { 
         if(strcmp(Key, temp_next->Key) == 0) 
         { 
            temp->next = temp_next->next; 
            delete temp_next; 
            size--; 
            return true; 
         } 
         temp = temp->next; 
         temp_next = temp_next->next; 
      } 
   } 
   return false; 
}//end remove 
  
  
void Hashtable::removeAll() 
{//start removeAll 
   for(int i=0; i<table_size; i++) 
   { 
      NODE* temp = table[i]; 
      while(temp != NULL) 
      { 
         NODE* next = temp->next; 
         //disp(temp); 
         delete temp; 
         temp = next; 
      } 
   } 
   size = 0; 
}//end removeAll 
  
int Hashtable::getSize() 
{ 
   return size; 
} 
  
NODE* Hashtable::find(char* Key) 
{ //start find 
   int bucket = hashString(Key); 
   NODE* temp = table[bucket]; 
   while(temp != NULL) 
   { 
      if(strcmp(Key, temp->Key) == 0) 
      { 
         return temp; 
      } 
      temp = temp->next; 
   } 
   return NULL; 
}//end find 
  
long Hashtable::hashString(char* Key) 
{//start hashString 
   int n = (int)strlen(Key); 
   long h = 0; 
   for(int i=0; i<n; i++) 
   { 
      //To get almost fair distributions of nodes over the array 
      h = (h << 3) ^ Key[i]; 
   } 
    return abs(h % table_size ); 
}//end hashString 
  
void Hashtable::initIterator() 
{//start initIterator 
   current_entry = NULL; 
   current_index = table_size; 
   for(int i=0; i<table_size; i++) 
   { 
      if(table[i] == NULL) 
      { 
          continue; 
      } 
      else
      { 
         current_entry = table[i]; 
         current_index = i; 
         break; 
      } 
   } 
}//end initIterator 
  
bool Hashtable::hasNext() 
{ 
   if(current_entry == NULL) 
   { 
      return false; 
   } 
   else
   { 
      return true; 
   } 
} 
void Hashtable::getNextKey(char* Key) 
{ 
   if(current_entry == NULL) 
   { 
      Key[0] = '\0'; 
      return; 
   } 
   strcpy(Key, current_entry->Key); 
   if(current_entry->next != NULL) 
   { 
      current_entry = current_entry->next; 
   } 
   else
   { 
     for(int i=current_index+1; i<table_size; i++) 
     { 
        if(table[i] == NULL) 
        { 
           continue; 
        } 
        current_entry = table[i]; 
        current_index = i; 
        return; 
     } 
     current_entry = NULL; 
     current_index = table_size; 
   } 
} 
  




