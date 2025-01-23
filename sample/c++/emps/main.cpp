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
#define SND 22.0
#define DNS_FLD 1000
#define DNS_WLL 1000
#define PCL_DST 0.02 //初期粒子間距離
#define DST_LMT_RAT 0.9 //これ以上の粒子の接近を許さない距離の係数
#define COL_RAT 0.2 //接近した粒子の反発係数
#define MIN (0.0 - PCL_DST * 3)
#define MAX_X (1.0 + PCL_DST * 3)
#define MAX_Y (0.2 + PCL_DST * 3)
#define MAX_Z (0.6 + PCL_DST * 30)
const double Reff = PCL_DST * 2.1; //影響半径r = 初期粒子間距離の2.1倍
const double Reff2 = Reff*Reff;
const double KNM_VSC = 0.000001;
const double rlim = PCL_DST * DST_LMT_RAT; //初期粒子間距離の0.9倍
const double rlim2 = rlim * rlim;
const double COL = 1.0 + COL_RAT; //接近した粒子に対する反発係数+1.0
const int DIM = 3;
double N0, LMD, A1, A2, A3;
int N;//particle_num
double Dns[2],invDns[2];
double calc_weight(const double dist){
    double ret = ((Reff/dist) - 1.0);
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
	Dns[FLD] = DNS_FLD;
	Dns[WLL] = DNS_WLL;
	invDns[FLD] = 1.0 / DNS_FLD;
	invDns[WLL] = 1.0 / DNS_WLL;
    A1 = (2.0 * DIM * KNM_VSC) / N0 / LMD;
    A2 = SND * SND / N0;
    A3 = -DIM / N0;
	//std::cout << "N0= " << N0 <<std::endl;
	//std::cout << "LMD= " << LMD << std::endl;
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
class Pressure{
	public:
	PS::F64 pres;
	PS::F64vec acc;
	PS::F64vec vel;
	void clear(){
		pres = 0;
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
    //inline static PS::F64 c = 22.0;//音速
    //inline static PS::F64 nu = 0.000001;//動粘性係数
    inline static PS::F64 grav = -9.8;
    inline static PS::F64 dt = 0.0005;
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
	void copyFromForce(const Pressure& prs){
    	this->pres = prs.pres;
    	this->acc = prs.acc;
    	this->vel = prs.vel;
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
    PS::F64 pres;
    PS::S32 Typ;
    static double cof;
    void copyFromFP(const FP& rp){
 		this->pos = rp.pos;
		this->vel = rp.vel;
		this->dt = rp.dt;
		this->Typ = rp.Typ;
		this->pres = rp.pres;
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
		for(PS::S32 i = 0; i < Nip; i++){
			acc[i].clear();
			if(ep_i[i].Typ == FLD){
	    	for(PS::S32 j = 0; j < Njp; ++j){
	    		const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
        		const double dist2 = dr*dr;
        		//std::cout << "dist2= " << dist2 << std::endl;
				if(dist2 == 0.0 || dist2 > Reff2)continue;
   	    		if(i != j && ep_j[j].Typ != GST){
					const double dist = sqrt( dist2 );
        			const double w = calc_weight(dist);
					acc[i].acc += (ep_j[j].vel - ep_i[i].vel)*w;
				}
	    	}//forループ抜け
        	acc[i].acc = acc[i].acc * A1; //+ FP::grav;
        	}
		}
   }
};
template<typename Tpi>
void ChkPcl(Tpi & pi){
	if(pi.pos.x > MAX_X || pi.pos.x < MIN || pi.pos.y > MAX_Y || pi.pos.y < MIN || pi.pos.z > MAX_Z || pi.pos.z < MIN){
		pi.Typ = GST;
		pi.pres = 0.0;
		pi.vel = 0.0;
		std::cout << "ChkPcl is success" << std::endl;
	}
}
template<typename Tptcl>
void first_UpPtcl(Tptcl & ptcl){
	for(int i = 0; i < N;i++){
		if( ptcl[i].Typ == FLD){
			ptcl[i].vel += ptcl[i].acc * FP::dt;
			ptcl[i].pos += ptcl[i].vel * FP::dt;
			ptcl[i].acc = 0.0;
			ChkPcl(ptcl[i]);
		}
	}
}
struct Mkprs{
	void operator()(const EP* const ep_i, const PS::S32 Nip,
					const EP* const ep_j, const PS::S32 Njp,
					Pressure* prs){
		for(PS::S32 i = 0; i < Nip; i++){
			prs[i].clear();
			if(ep_i[i].Typ != GST){
				double mi = Dns[ep_i[i].Typ];
				PS::F64vec vec_i = ep_i[i].vel;
				PS::F64 ni = 0.0;
				for(PS::S32 j = 0; j < Njp; ++j){
					const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
					const double dist2 = dr * dr;
					if(dist2 == 0.0 || dist2 > Reff2)continue;
					const double dist = sqrt(dist2);
					const double w = calc_weight(dist);
					ni += w;
					if(ep_i[i].Typ == FLD && dist2 < rlim){
						if(j != i && ep_j[j].Typ != GST){
							double fDT = (ep_i[i].vel - ep_j[j].vel) * dr;
							if(fDT > 0.0){
								const double mj = Dns[ep_j[j].Typ];
								fDT *= COL * mj / (mi + mj) / dist2;
								vec_i -= dr * fDT;
							}
						}
					}
				}
				prs[i].acc = vec_i;
				double pressure = (ni > N0) * (ni - N0) * A2 * mi;
				std::cout << "pressure[" << i << "] = "<< pressure << std::endl;
				prs[i].pres = pressure;
			}
		}
		for(PS::S32 i = 0; i < Nip; i++){
			prs[i].vel = prs[i].acc;
		}
	}
};
struct PrsGrdTrm{
	void operator()(const EP* const ep_i, const PS::S32 Nip,
					const EP* const ep_j, const PS::S32 Njp,
					Acc* acc){
		for(PS::S32 i = 0; i < Nip; i++){
			if(ep_i[i].Typ == FLD){
				acc[i].clear();
				double pres_min = ep_i[i].pres;
				for(PS::S32 j = 0; j < Njp; j++){
					const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
					const double dist2 = dr * dr;
					if(dist2 == 0.0 || dist2 > Reff2)continue;
					if(j != i && ep_j[j].Typ != GST){
						if(pres_min > ep_j[j].pres)pres_min = ep_j[j].pres;
					}
				}
				for(PS::S32 j = 0; j < Njp; j++){
					const PS::F64vec dr = ep_j[j].pos - ep_i[i].pos;
					const double dist2 = dr * dr;
					if(dist2 == 0.0 || dist2 > Reff2)continue;
					if(j != i && ep_j[j].Typ != GST){
						const double dist = sqrt(dist2);
						double w = calc_weight(dist);
						w *= (ep_j[j].pres - pres_min)/dist2;
						acc[i].acc += dr * w;
					}
				}
				acc[i].acc = acc[i].acc * invDns[FLD] * A3;
				//std::cout << "acc[" << i << "]= "<< acc[i].acc << std::endl;
			}
		}
	}
};
template<typename Tptcl>
void second_UpPtcl(Tptcl & ptcl){
	for(int i = 0; i < N; i++){
		if(ptcl[i].Typ == FLD){
			ptcl[i].vel += ptcl[i].acc * FP::dt;
			ptcl[i].pos += ptcl[i].acc * FP::dt * FP::dt;
			ptcl[i].acc = 0.0;
			ChkPcl(ptcl[i]);
		}
	}
}
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
	PS::TreeForForceShort<Pressure, EP, EP>::Scatter pres_tree;
	dinfo.initialize();
	Set_Para();
	FILE*fp;
	fp = fopen(IN_FILE,"r");
	fscanf(fp,"%d\n",&N);
	emps.setNumberOfParticleLocal(N);
	acc_tree.initialize(N);
	pres_tree.initialize(N);
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
	fclose(fp);
	dinfo.decomposeDomainAll(emps);
	emps.exchangeParticle(dinfo);
	fp = fopen("result/result.csv","w");
	acc_tree.calcForceAllAndWriteBack(CalcAcc(), emps, dinfo);
	for(int i = 0; i<N; i++){
		if(emps[i].Typ == FLD)emps[i].acc.z = emps[i].acc.z + FP::grav;
	}
	first_UpPtcl(emps);
	pres_tree.calcForceAllAndWriteBack(Mkprs(), emps, dinfo);
	acc_tree.calcForceAllAndWriteBack(PrsGrdTrm(), emps, dinfo);
	second_UpPtcl(emps);
	pres_tree.calcForceAllAndWriteBack(Mkprs(), emps, dinfo);
	for(int i=0; i<N; i++){
		fprintf(fp, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",emps[i].Typ,emps[i].pos.x,emps[i].pos.y,emps[i].pos.z,emps[i].vel.x,emps[i].vel.y,emps[i].vel.z,emps[i].acc.x,emps[i].acc.y,emps[i].acc.z,emps[i].pres);
	}
	fclose(fp);
	PS::Finalize();
	return 0;
}
