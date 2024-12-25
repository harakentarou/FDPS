#include <particle_simulator.hpp>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#define IN_FILE "../mk_particle/dambreak.prof"
#define FLD 0
#define WLL 1
#define GST -1
const double PCL_DST = 0.02;//初期粒子間距離
const double Reff = PCL_DST * 2.1; //影響半径r = 初期粒子間距離の2.1倍
const double Reff2 = Reff*Reff;
const double KNM_VSC = 0.000001;
const int DIM = 3;
double N0, LMD, A1;
int N;
double calc_weight(const double dist){
    double ret = (Reff/dist - 1.0);
    return std::max(ret, 0.0);
}
void Set_Para(void){
	double tn0 = 0.0;
	double tlmd = 0.0;
	for(int ix = -4; ix < 5; ix++){
	for(int iy = -4; iy < 5; iy++){
	for(int iz = -4; iz < 5; iz++){
		double x = PCL_DST * ix;
		double y = PCL_DST * iy;
		double z = PCL_DST * iz;
		double dist2 = x*x + y*y + z*z;
		if(dist2 <= Reff2){
			if(dist2 == 0.0)continue;
			double dist = sqrt(dist2);
			tn0 += calc_weight(dist);
			tlmd += dist2 * calc_weight(dist);
		}
	}}}
	N0 = tn0;
	LMD = tlmd / tn0;
    A1 = (2.0 * DIM * KNM_VSC) / N0 / LMD;
	std::cout << "N0= " << N0 <<std::endl;
	std::cout << "LMD= " << LMD << std::endl;
}
class Dens{
	public:
	PS::F64 dens;
	void clear(){
		dens = 0;
	}
};
class Hydro{
	public:
	PS::F64vec acc;
	PS::F64 dt;
	void clear(){
		acc = 0;
		dt = 0;
	}
};
class Acc{
    public:
    PS::F64vec acc;
    void clear(){
	acc = 0.0;
	}
};

// Full Particle Class
struct FP{
	PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pres;
    PS::F64 pav;//時間平均された粒子の圧力
    PS::F64 dens;
    PS::S32 Typ;
    inline static PS::F64 c = 22.0;//音速
    //inline static PS::F64 nu = 0.000001;//動粘性係数
    inline static PS::F64 grav = -9.8;
    PS::F64 dt = 0.0005;
    PS::S32 id;
    void copyFromForce(const Dens& dens){
		this->dens = dens.dens;
    }
    void copyFromForce(const Hydro& force){
		this->acc = force.acc;
		this->dt = force.dt;
    }
    void copyFromForce(const Acc& acc){
    	this->acc = acc.acc;
    }
    PS::F64vec getPos() const{
		return this->pos;
    }
    void setPos(const PS::F64vec& pos){
		this->pos = pos;
    }
};
struct EP{
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 dt;
    PS::S32 Typ;
    static double cof;
    void copyFromFP(const FP& rp){
 		this->pos = rp.pos;
		this->vel = rp.vel;
		this->dt = rp.dt;
		this->Typ = rp.Typ;
    }
    PS::F64vec getPos() const{
		return this->pos;
    }
    void setPos(const PS::F64vec& pos){
		this->pos = pos;
    }
    PS::F64 getRSearch() const {
    	return Reff;
    }
};
struct CalcAcc{
    void operator() (const EP* const ep_i, const PS::S32 Nip,
		    		 const EP* const ep_j, const PS::S32 Njp,
		    		 Acc* acc){
		std::cout<<"A1= "<<A1<<std::endl;
		for(PS::S32 i = 0 ; i < Nip ; i++){
	    	acc[i].clear();
	    	for(PS::S32 j = 0 ; j < Njp ; ++j){
	    		const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
        		const double dist2 = dr*dr;
        		if(dist2 == 0.0) continue;
				const double dist = sqrt( dist2 );
        		const double w = calc_weight(dist);
				acc[i].acc += (ep_j[j].vel - ep_i[i].vel)*w;
	    	}
        	acc[i].acc = acc[i].acc * A1 + FP::grav;
			if(ep_i[i].Typ == FLD){
        		std::cout << "acc[" << i << "] = " << acc[i].acc << std::endl;
			}
		}
   }
};
void first_UpPtcl(){
	for(int i = 0; i < ;i++){
		if( == FLD){
			FP::vel += 
/*void makeOutputDirectory(char* dir_name){
	struct stat st;
	PS::S32 ret;
	if(PS::Comm::getRank() == 0){
		if(stat(dir_name, &st) != 0){
			ret = mkdir(dir_name, 0777);
		}
		else{
			ret = 0;
		}
	}
	PS::Comm::broadcast(&ret, 1);
	if(ret == 0){
		if(PS::Comm::getRank() == 0)
			fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
	}
	else{
		if(PS::Comm::getRank() == 0)
			fprintf(stderr, "Directory %s is fails to be made.\n", dir_name);
		PS::Abort();
	}
}*/
int main(int argc, char *argv[]){
	PS::Initialize(argc, argv);
	/*char dir_name[1024];
	sprintf(dir_name,"./result");
	makeOutputDirectry(dir_name);*/
	PS::ParticleSystem<FP> emps;
	emps.initialize();
	PS::F64 dt, end_time;
	PS::DomainInfo dinfo;
	PS::TreeForForceShort<Acc, EP, EP>::Scatter acc_tree;
	dinfo.initialize();
	Set_Para();
	FILE*fp;
	fp = fopen(IN_FILE,"r");
	fscanf(fp,"%d\n",&N);
	emps.setNumberOfParticleLocal(N);
	acc_tree.initialize(N);
	for(int i=0; i<N; i++){
        int a[2];
 		double b[8];
        fscanf(fp," %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",&a[0],&a[1],&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7]);
        emps[i].id = a[0];
        emps[i].Typ = a[1];
        emps[i].pos.x = b[0];
        emps[i].pos.y = b[1];
        emps[i].pos.z = b[2];
        emps[i].vel.x = b[3];
        emps[i].vel.y = b[4];
        emps[i].vel.z = b[5];
        emps[i].pres = b[6];
        emps[i].pav = b[7];
	}
	dinfo.decomposeDomainAll(emps);
	emps.exchangeParticle(dinfo);
	acc_tree.calcForceAllAndWriteBack(CalcAcc(), emps, dinfo);
	PS::Finalize();
	return 0;
}
