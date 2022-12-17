#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
typedef double deci;

using namespace std;

deci myRan();

deci myNorm(deci *vec){
  return sqrt(pow(vec[0],2)+pow(vec[1],2)+pow(vec[2],2));
}

deci genWavelen(ifstream* f){
  string fn;
  if((*f).is_open()){
    getline(*f,fn);
    return atof(fn.c_str());
  }
  return 0.0;
}

class record{
  public:
  deci myp[3],myd[2];
  deci myl;
  int plan_index;
  record(deci *pos, deci *dir, deci l, int p_i){
    plan_index=p_i;
    for(int i=0;i<3;i++){
      myp[i]=pos[i];
      if(i==2) continue;
      myd[i]=dir[i];
    }
    myl=l; 
  }
};

class glob_vars{
  public:
  deci r_atm,n_atm,rayleigh_coef,dt,PI,TWO_PI;
  deci r_plan[5];
  bool rayleigh;
  int disap;
  vector<record> records;
  glob_vars():r_plan{3000,3500,4000,5000,5500}{
    //^constructor initialization list
    r_atm=6000;
    n_atm=1;
    rayleigh=true;
    deci s_low=1e-9;
    rayleigh_coef=pow(s_low,4)*1e9;
    dt=1e-9*299792458;
    disap=0;
    PI=acos(-1);
    TWO_PI=2*PI;
    //r_atm[4]={5500,6000,7000,8000,10000};
  }
  void printRecords(ofstream *ofile){
    //*ofile<<"(n, l, pos, dir). Disappeared: "<<disap<<endl;
    cout<<"Records:"<<endl;
    for(int i=0;i<records.size();i++){
      //cout<<"dis: "<<myNorm(records[i].myp)<<", th: "<<atan2(records[i].myp[1], records[i].myp[0])<<", n: "<<records[i].plan_index<<", l: "<<records[i].myl<<endl;
      *ofile<<records[i].plan_index<<"    "<<records[i].myl<<"    "<<atan2(records[i].myp[1], records[i].myp[0])<<"    "<<records[i].myd[0]<<endl;
    }
  }
}gv;

class photon{
  deci myp[3],myd[2];
  deci myl;
  int plan_index;
  public:
  bool exist;
  photon(deci *pos, deci *dir, deci l){
    plan_index=4;
    for(int i=0;i<3;i++){
      myp[i]=pos[i];
      if(i==2) continue;
      myd[i]=dir[i];
    }
    myl=l;
    exist=true;
  }
  void move(){
    if(!exist) return;
    deci pdis=myNorm(myp);
    //cout<<"Mypos: "<<pdis<<endl;

    //Rayleigh scattering
    if(pdis<gv.r_atm && pdis>gv.r_plan[0]){
      //cout<<"Random num:"<<myRan()<<endl;
      deci ran=myRan();
      deci rc=gv.rayleigh_coef*pow(myl,-4);
      if(ran<rc && gv.rayleigh){
        myd[0]=gv.TWO_PI*myRan();
        myd[1]=gv.PI*myRan();
      }
    }

    //Move
    if(pdis>gv.r_plan[0]){
      myp[0]+=gv.dt*sin(myd[0])*cos(myd[1]);
      myp[1]+=gv.dt*sin(myd[0])*sin(myd[1]);
      myp[2]+=gv.dt*cos(myd[0]);
      if(pdis>7000){
        //Be gone
        exist=false;
        gv.disap++;
        //cout<<"Dissappeared!"<<endl;
      }
    }

    //Absorb
    if(pdis<gv.r_plan[plan_index]){
      //Add
      record my_rec(myp,myd,myl,plan_index);
      gv.records.push_back(my_rec);
      plan_index--;
      if(plan_index<0){
        //Be gone
        exist=false;
      }
    }
  }
};

void testing(){
  ifstream fs("rp_wl.txt");
  deci aa=0.0;
  aa+=genWavelen(&fs);
  aa+=genWavelen(&fs);
  cout<<aa<<endl;
  for(int i=0;i<10;i++){
    cout << myRan()<<endl;
  }

  deci z[]={-6200.0,gv.r_atm*(2*myRan()-1),0.0};
  deci zz[]={gv.PI/2,0.0};
  record a_rec0(z,zz,4e-6,3);
  gv.records.push_back(a_rec0);
  record a_rec1(z,zz,5e-6,1);
  gv.records.push_back(a_rec1);

  fs.close();
}

int main(int argc, char** argv) {
  srand (time(NULL));
  char path[150];
  char opath[150];
  if(argc==2){
    strcpy(path,argv[1]);
    strcpy(opath,argv[1]);
    strcat(path,".txt");
    strcat(opath,"_out.txt");
  } else return -1;
  ifstream fs(path);
  ofstream fo(opath);
  ofstream fl("raypla_log.txt",ios::ate | ios::app);//log-file

  if(fs.is_open()){
    cout<<"It opened"<<endl;
    //string fn;
    //getline(fs,fn);
    //cout<<fn<<endl;
    //cout<<(genWavelen(&fs)*3)<<endl;
  }else{
    cout<<"Not opened"<<endl;
    return 1;
  }
  if(fo.is_open()){
    cout<<"Output opened"<<endl;
  }else{
    cout<<"Output not opened"<<endl;
    return 1;
  }

  //cout << "pow haha "<<pow(0.02,-4)<<endl;
  //deci abc[]={1.,1.,1.};
  //cout <<myNorm(abc)<<std::endl;

  //vector<photon> photons;
  //for(int i=0;i<20000;i++){
  int particlecnt = 0;
  int totmovs = 0;
  while(fs.good()){
    deci initpos[]={-6200.0,gv.r_atm*(2*myRan()-1),0.0};
    deci initvel[]={gv.PI/2,0.0};
    photon my_pho(initpos,initvel,genWavelen(&fs));
    //photons.push_back(my_pho);
    while(my_pho.exist){
      my_pho.move();
	  totmovs++;
    }
	particlecnt++;
  }
  //deci z[]={-6200.0,0.0,0.0};
  //for(int j=0;j<128;j++){
    //for(int i=0;i<10000;i++){
  //}

  gv.printRecords(&fo);
  //Disappeared, TotalParticles, TotalMoves
  fl << path<<" disap "<<gv.disap<<" totpa "<<particlecnt<<" totmov "<<totmovs<<endl;

  fo.close();
  fs.close();
  fl.close();
  return 0;
}

deci myRan(){
  deci num;
  num=(deci)rand()/RAND_MAX;
  return num;
}
//,[>+>+<<-]>[>.<-]