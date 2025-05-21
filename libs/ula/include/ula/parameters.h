#ifndef WESSEL_ULA_PARAMETERS
#define WESSEL_ULA_PARAMETERS

#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

class parameters {
 public:
  parameters() : key_(), value_() {}
  bool get(string jobfn);
  unsigned int size() const {return key_.size();}
  void clear() {key_.clear();value_.clear();}
  void push_back(string  k, string v);
  double operator[] (string k) const; 
  string key(int i) const;
  string value(int i) const;
  bool defined(string k) const;
  void dump(ostream& ) const;
  void restore(istream& );
 private:
  vector<string> key_;
  vector<string> value_;
};


inline void parameters::push_back(string k, string v) {
  key_.push_back(k);
  value_.push_back(v);
}

inline double parameters::operator[](string k) const {
  for (unsigned int j=0;j<key_.size();++j) 
    if (key_[j]==k) 
      return atof(value_[j].c_str());
  cerr << " Parameter key not defined : " << k << endl;
  exit(1);
  return 0;
}

inline bool parameters::defined(string k) const {
  for (unsigned int j=0;j<key_.size();++j) 
    if (key_[j]==k) 
      return 1;
  return 0;
}

inline string parameters::key(int i) const {
  if (i>=0 && i<(int)key_.size() ) 
    return key_[i];
   cerr << " Parameter key not defined at index: " << i << endl;
  exit(1);
  return "";
}

inline string parameters::value(int i) const {
  if (i>=0 && i<(int)value_.size()) 
    return value_[i];
  cerr << " Parameter value not defined at index: " << i << endl;
  exit(1);
  return "";
}

bool parameters::get(string jobfn) {
  vector<string> key;
  vector<string> value_row;
  vector<vector<string> > values;
  vector<int> times;
  ifstream jobf;
  jobf.open(jobfn.c_str());
  if (!jobf) {
    cerr<<"Jobfile not found: "<<jobfn<<endl;
    exit(1);
  }
  string p;
  jobf >> p;
  while(p!="##") {
    key.push_back(p);
    jobf >> p;
  };
  while (jobf >> p) {
    value_row.push_back(p);
    if (value_row.size()==key.size()) {
      values.push_back(value_row);
      value_row.clear();
      int t;
      jobf >> t;
      times.push_back(t);
    }
  };
  jobf.close();
  ofstream jobfo;
  jobfo.open(jobfn.c_str());
  for (unsigned int i=0;i<key.size();++i)
    jobfo << key[i] <<'\t';
  jobfo <<"##\n";
  bool jobfound=false;
  for (unsigned int i=0;i<values.size();++i) {
    if (times[i]!=0 && !jobfound) {
      jobfound=true;
      this->clear();
      for (unsigned int j=0;j<key.size();j++) 
        this->push_back(key[j],values[i][j]);
      if (times[i]>0) {
        times[i]--;
        this->push_back("separate","0");
      }
      else {
        times[i]++;
        this->push_back("separate","1");
      }
    }
    for (unsigned int j=0;j<key.size();j++)
     jobfo << values[i][j] <<'\t';
    jobfo << times[i] <<'\n';
  };
  jobfo.close();
  return jobfound;
}


void parameters::dump(ostream& dumpfile) const {
  dumpfile<<key_.size()<<" ";
  for (unsigned int j=0;j<key_.size();++j) 
    dumpfile << key_[j] <<" " << value_[j] <<" ";
}


void parameters::restore(istream& dumpfile) {
  this->clear();
  int sz;
  dumpfile >> sz;
  for (;sz;--sz) {
    string hs;
    dumpfile >> hs;
    key_.push_back(hs);
    dumpfile >> hs;
    value_.push_back(hs);
  };
}


#endif
