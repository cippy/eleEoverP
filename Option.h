#ifndef Option_h
#define Option_h

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h                                                                                                                                      
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file                                                                        
#include <string>
#include <vector>
#include <iomanip> //for input/output manipulators  


class Option {

 public:
  Option();
  ~Option();
  Int_t Get_data2016() const { return data2016;}
  Int_t Get_skim1lep1jet80X() const { return skim1lep1jet80X;}
  Int_t Get_fit2sideCB() const { return fit2sideCB;}
  Int_t Get_setScaleOnY() const { return setScaleOnY;}
  Int_t Get_useE() const { return useE;}
  std::string Get_dirName() const { return dirName;}

  void Set_data2016(const Int_t &value) { data2016 = value; }
  void Set_skim1lep1jet80X(const Int_t &value) { skim1lep1jet80X = value; }
  void Set_fit2sideCB(const Int_t &value) { fit2sideCB = value; }
  void Set_setScaleOnY(const Int_t &value) { setScaleOnY = value; }
  void Set_useE(const Int_t &value) { useE = value; }
  void Set_dirName(const std::string &value) { dirName = value; }


 private:
  Int_t data2016;
  Int_t skim1lep1jet80X;
  Int_t fit2sideCB;
  Int_t setScaleOnY;
  Int_t useE;
  std::string dirName;
};

Option::Option() {
  data2016 = 1;
  skim1lep1jet80X = 1;
  fit2sideCB = 0;
  setScaleOnY = 0;
  useE = 0;
  dirName = "default";
}

Option::~Option() {
  std::cout << "Option::~Option() called" << std::endl;
}


#endif
