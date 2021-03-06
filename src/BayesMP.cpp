#include "RCbridge.h"

#include <cstring>
#include <fstream>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <random>
#include <algorithm>

using namespace std;

class Para{
	// membership: DP class: 1,2,3,...: positive DE. -1,-2,-3...: negative DE. 0: nonDE.
	// index: gene index.
	// n: total number of genes.
	int membership, n;
	// sumZ: sum of Z statistics
	// postmu: posterior mu
	// postsd: posterior sd
	double sumZ, postmu, postsd, mu0, sigma0, sigma;
	
public:
	std::vector<int> index;
	
	//indexC * index;
		
	Para(double amu0, double asigma0, double asigma, int anum, double aZ, int amembership)
	{
		membership=amembership;
		mu0=amu0;
		sigma0=asigma0;
		sigma=asigma;	    
		index.push_back(anum);
		//index = new indexC(anum);
		sumZ = aZ;
		n = 1;		
		updateTheta();	
	}

	void setMembership(int amember)
	{
		membership = amember;
	}

	int getMembership()
	{
		return(membership);
	}
	
	void addZ(int anum, double aZ)
	{		
		index.push_back(anum);
		sumZ = sumZ + aZ;
		n = n + 1;
		updateTheta();
	}

	int removeZ(int anum, double aZ)
	{		
		// 0: nothing left. 1: successfully. 2; error: no such anum exist.
		index.erase(std::remove(index.begin(), index.end(), anum), index.end());
		
		sumZ = sumZ - aZ;
		n = n - 1;
		if(n==0)
		{					
			return 0;
		}
		updateTheta();
		return 1;		
	}

	void updateTheta(){
		postmu = (sumZ/(sigma*sigma) + mu0/(sigma0*sigma0))/(n/(sigma*sigma) + 1/(sigma0*sigma0));
		postsd = 1/sqrt(n/(sigma*sigma) + 1/(sigma0*sigma0));			
	}
	
	int GetN()
	{
		return(n);
	}

	double GetSumZ()
	{
		return(sumZ);
	}

	double Getpostmu()
	{
		return(postmu);
	}

	double Getpostsd()
	{
		return(postsd);
	}

	double Getmu0()
	{
		return(mu0);
	}

	double Getsigma0()
	{
		return(sigma0);
	}
	
};



class ParaList{
	int length;
	int newMemPlus;
	int newMemMinus;

public:
	std::vector<Para> paraList;	

	ParaList(){
		length = 0;
		newMemPlus = 1;
		newMemMinus = -1;
	}
	
	void addPara(Para apara){
		paraList.push_back(apara);
		length++;
	}
	
	int getLength(){
		return(length);
	}
	
	Para getPara(int l){
		return(paraList[l]);
	}	
	
	void erasePara(int l){
		paraList.erase(paraList.begin() + l);
		length--;
	}	
	
	int getNewMembership(int direction){
		if(direction == 1){
			return newMemPlus;
		} else if(direction == -1){
			return newMemMinus;
		} 
		cout << "error: no such direction!" << endl;
		return 0;
	}

	void updateNewMembership(int direction){
		int lengthAparaList = length;
		int flag = 1;
		
		if(direction == 1){
			int minMem = 1;
			while(1){
				flag = 1;
				for(int l=0;l<lengthAparaList;l++){
					int aMem = paraList[l].getMembership();
					if(aMem == minMem || minMem == newMemPlus){
						minMem++;
						flag = 0;
						break;
					}			
				}
				if(flag){
					newMemPlus = minMem;
					return;
				}			
			}
			
		} else if(direction == -1){
			int minMem = -1;
			while(1){
				flag = 1;
				for(int l=0;l<lengthAparaList;l++){
					int aMem = paraList[l].getMembership();
					if(aMem == minMem || minMem == newMemMinus){
						minMem--;
						flag = 0;
						break;
					}			
				}
				if(flag){
					newMemMinus = minMem;
					return;					
				}			
			}
		} 

	}
	
	int getParaSumNP(){
		int sumN = 0;		
		int lengthAparaList = length;
		for(int l=0;l<lengthAparaList;l++){
			if(paraList[l].getMembership()>0){
				sumN += paraList[l].GetN();
			}			
		}
		return sumN;
	}

	int getParaSumNN(){
		int sumN = 0;		
		int lengthAparaList = length;
		for(int l=0;l<lengthAparaList;l++){
			if(paraList[l].getMembership()<0){
				sumN += paraList[l].GetN();
			}			
		}
		return sumN;
	}
	
};


class bayesMP{
	// dimension variables G: number of genes. S: number of studies.
	int G, S;
	// Z: observed Z statistics
	//double *Z;
	std::vector<double> Z;
	// gamma: prior for pi_g; beta: prior for delta_g; alpha: hyperparameter for DP; mu0: mean para for G0; sigma0: sd para for G0; sigma: sd for emission;
	double gamma, beta, alpha, mu0, sigma0, sigma, trunc;
	int randomGamma;
	std::vector<double> empMu;
	std::vector<double> empSD;
	
	int countAcceptGamma = 0;
	// pi: prob gene g to be DE. delta: prob d DE gene to be positive.
	std::vector<double> pi;
	std::vector<double> delta;
	//double *pi, *delta;
	// Y: latent variable to be infered.
	std::vector<int> Y;
	//int *Y;
	// whether directly output HSind result and HSall result
	int writeY, writePi, writeDelta, writeGamma, writeHSall;
	// number of mcmc iterations.
	int niter, burnin;
	// current iter
	int thisIter;
	char * file_Y;
	char * file_Pi;
	char * file_Delta;
	char * file_Gamma;

	double MHsd = 0.1;
	
	
	std::vector<int> YHSall;
	//int *YHSall;
	ofstream Stream_Y;	
	ofstream Stream_Pi;	
	ofstream Stream_Delta;	
	ofstream Stream_Gamma;	
	
	default_random_engine generator;
	
 	void SetG(int g){
		G = g;
 	}
	
 	void SetS(int s){
		S = s;
 	}

	/*
 	void SetZ(double *aZ){
		int GS = G*S;
		Z = new double [GS];		
		for(int i=0; i<GS; i++){
			Z[i] = *aZ++;
		}
 	}
	*/
 	void SetZ(double *aZ){
		Z = std::vector<double>(aZ, aZ + G*S);		
 	}

 	void SetPi(double *api){
		pi = std::vector<double>(api, api + G);		
 	}

 	void SetDelta(double *adelta){
		delta = std::vector<double>(adelta, adelta + G);		
 	}

 	void SetY(int *aY){
		Y = std::vector<int>(aY, aY + G*S);		
 	}

 	void SetGamma(double agamma){
		gamma = agamma;		
 	}

 	void SetRandomGamma(double arandomGamma){
		randomGamma = arandomGamma;		
 	}

 	void SetEmpMu(double *aempMu){
		empMu = std::vector<double>(aempMu, aempMu + S);
 	}

 	void SetEmpSD(double *aempSD){
		empSD = std::vector<double>(aempSD, aempSD + S);		
 	}

 	void SetBeta(double abeta){
		beta = abeta;		
 	}

 	void SetAlpha(double aalpha){
		alpha = aalpha;
 	}
	
 	void SetMu0(double amu0){
		mu0 = amu0;
 	}

 	void SetSigma0(double asigma0){
		sigma0 = asigma0;
 	}

 	void SetSigma(double asigma){
		sigma = asigma;
 	}
	
 	void SetTrunc(double atrunc){
		trunc = atrunc;
 	}
	
 	void SetnIter(int aniter){
		niter = aniter;
 	}	
	
 	void SetnBurnin(int aburnin){
		burnin = aburnin;
 	}	
		
 	void SetwriteY(int awriteY){
		writeY = awriteY;
 	}

 	void SetwritePi(int awritePi){
		writePi = awritePi;
 	}

 	void SetwriteDelta(int awriteDelta){
		writeDelta = awriteDelta;
 	}

 	void SetwriteGamma(int awriteGamma){
		writeGamma = awriteGamma;
 	}

 	void SetwriteHSall(int awriteHSall){
		writeHSall = awriteHSall;
 	}
	
	void iniYHSall(){
		YHSall = std::vector<int>(G*S, 0);				
	}

	void SetFilename_Y(char *filename){
		file_Y = filename;
	}

	void SetFilename_Pi(char *filename){
		file_Pi = filename;
	}

	void SetFilename_Delta(char *filename){
		file_Delta = filename;
	}

	void SetFilename_Gamma(char *filename){
		file_Gamma = filename;
	}

	void iniBayesMPparaLists(int S){
		for(int s=0;s<S;s++){
			bayesMPparaLists.push_back(ParaList());
		}		
		
	}
	
	
public:
	
	std::vector<ParaList> bayesMPparaLists;
		
	void initialize(int *aG, int *aS, double *aZ, double *agamma, int *randomGamma, double *aempMu, double *aempSD, double *abeta,double *aalpha ,double *amu0, 
			double *asigma0, double *asigma, double *atrunc, double *api, double *adelta, int *aY, int *niter, int *burnin, 
			char *fileName_Y, char *fileName_Pi, char *fileName_Delta, char *fileName_Gamma, 
			int *writeY, int *writePi, int *writeDelta, int *writeGamma, int *writeHSall)
	{
		SetG(*aG);
		SetS(*aS);
		SetZ(aZ);
		//SetZZ(aZ);
		SetGamma(*agamma);
		SetRandomGamma(*randomGamma);
		SetEmpMu(aempMu);
		SetEmpSD(aempSD);
		SetBeta(*abeta);
		SetSigma0(*asigma0);
		SetSigma(*asigma);
		SetAlpha(*aalpha);
		SetTrunc(*atrunc);
		SetMu0(*amu0);
		SetPi(api);
		SetDelta(adelta);
		SetY(aY);
		SetnIter(*niter);
		SetnBurnin(*burnin);
		SetFilename_Y(fileName_Y);
		SetFilename_Pi(fileName_Pi);
		SetFilename_Delta(fileName_Delta);
		SetFilename_Gamma(fileName_Gamma);
		SetwriteY(*writeY);
		SetwritePi(*writePi);
		SetwriteDelta(*writeDelta);
		SetwriteGamma(*writeGamma);
		SetwriteHSall(*writeHSall);
		iniYHSall();
		iniBayesMPparaLists(S);
		thisIter = 0;								
	}

	void iniPara(){
		
		for(int s=0;s<S;s++){
			ParaList aparaList = bayesMPparaLists[s];
			for(int g=0;g<G;g++){
				int findFlag = 0;
				int sGg = s*G+g;
				if(Y[sGg]!=0){
					int lengthAparaList = aparaList.getLength();
					if(lengthAparaList==0){
						aparaList.addPara(Para(mu0,sigma0,sigma,sGg,Z[sGg],Y[sGg]));
					} else {
						for(int l=0;l<lengthAparaList;l++){
							if(Y[sGg]==aparaList.paraList[l].getMembership())
							{
								aparaList.paraList[l].addZ(sGg, Z[sGg]);
								findFlag = 1;
								break;
							}
						}
						if(findFlag==0){
							aparaList.addPara(Para(mu0,sigma0,sigma,sGg,Z[sGg],Y[sGg]));
						}
					} // if else 
				} // end of if Y[sGg]!=0
			} // end of loop for g of G
			/*
			cout << "s = " << s << endl;
			for(int l=0;l<aparaList.getLength();l++){
				cout << "l = " << l << ". n: " << aparaList.paraList[l].GetN() << endl;
			}
			*/
			// cout <<"here 1" << endl;
			aparaList.updateNewMembership(1);
			// cout <<"here 2" << endl;
			aparaList.updateNewMembership(-1);			
			// cout <<"here 3" << endl;
			bayesMPparaLists[s] = aparaList;		
			// cout <<"here 4" << endl;
				
			/*
			
			for(int l=0;l<bayesMPparaLists[s].getLength();l++){
				cout << "l = " << l << ". n: " << bayesMPparaLists[s].paraList[l].GetN()  << endl;
			}
			*/
		} // end of loop for s of S
	}
	
	
	void deletePara1(int g, int s){
		int sGg = s*G+g;
		int aY = Y[sGg];
		if(aY!=0){
			int lengthAparaList = bayesMPparaLists[s].getLength();
			for(int l=0;l<lengthAparaList;l++){
				if(aY==bayesMPparaLists[s].paraList[l].getMembership()){
					int removeZStatus = bayesMPparaLists[s].paraList[l].removeZ(sGg,Z[sGg]);					
					if(removeZStatus==1){
						return;
					} else if(removeZStatus==0){
						bayesMPparaLists[s].erasePara(l);					
						return;
					} else if(removeZStatus==2){
						cout<<"such index doesn't exist, bug 2"<<endl;
						exit(0);
					} else {
						cout<<"unknown error, bug 0!"<<endl;
						exit(0);
					}														
				}
			} // for loop of l for lengthAparaList
		} // if(aY!=0)
	}

	void addPara1(int g, int s){		
		int findFlag = 0;
		int sGg = s*G+g;
		int aY = Y[sGg];

		if(aY==0){
			return;
		}
		
		int lengthAparaList = bayesMPparaLists[s].getLength();
		for(int l=0;l<lengthAparaList;l++){
			if(aY==bayesMPparaLists[s].paraList[l].getMembership()){
				bayesMPparaLists[s].paraList[l].addZ(sGg,Z[sGg]);
				findFlag = 1;
				break;
			}
		} // for loop of l for lengthAparaList
		if(findFlag==0){
			bayesMPparaLists[s].addPara(Para(mu0,sigma0,sigma,sGg,Z[sGg],aY));
		}		
	}	

	void iterateOne() {
		for(int g=0; g<G; g++){
			for(int s=0; s<S; s++){
				updateOne(g,s);
			}
		}		
		
		updatePi();
		if(writeHSall==1 && thisIter>=burnin) updateHSall();	
		if(randomGamma == 1) updateGamma();							
		if(writeY)	 appendFile_Y();
		if(writePi)	 appendFile_Pi();
		if(writeDelta)	 appendFile_Delta();
		if(writeGamma)	 appendFile_Gamma();
		
		thisIter++;
	}		
		

	void updateOne(int g, int s) {
		deletePara1(g, s);
		updateMembership(g ,s);
		addPara1(g, s);
	}
	
	void updateHSall(){
		int gY;
		for(int g=0; g<G; g++){
			gY = -1;
			for(int s=0;s<S; s++){
				if(Y[s*G+g]!=0){
					gY++;
					YHSall[gY*G + g]++;					
				}
			} 
		}
	}

	void outputHSall(char *myfile){
	    ofstream streamHSall;
	    streamHSall.open(myfile);		
		
		for(int g=0;g<G;g++){
			for(int s=0;s<S;s++){
				streamHSall << YHSall[s*G+g] << "\t";
			}
		    streamHSall << endl;						
		}		
		streamHSall.close();
	}
	
	void appendFile_Y(){
		if(thisIter==0){
			Stream_Y.open(file_Y);
		}
		for(int i=0;i<G*S;i++){
			Stream_Y << Y[i] << "\t";
		}
	    Stream_Y << endl;
		if(thisIter==niter-1){Stream_Y.close();}		
	}

	void appendFile_Pi(){
		if(thisIter==0){
			Stream_Pi.open(file_Pi);
		}
		for(int i=0;i<G*S;i++){
			Stream_Pi << pi[i] << "\t";
		}
	    Stream_Pi << endl;
		if(thisIter==niter-1){Stream_Pi.close();}		
	}

	void appendFile_Delta(){
		if(thisIter==0){
			Stream_Delta.open(file_Delta);
		}
		for(int i=0;i<G*S;i++){
			Stream_Delta << delta[i] << "\t";
		}
	    Stream_Delta << endl;
		if(thisIter==niter-1){Stream_Delta.close();}		
	}

	void appendFile_Gamma(){
		if(thisIter==0){
			Stream_Gamma.open(file_Gamma);
		}
		Stream_Gamma << gamma << endl;
		if(thisIter==niter-1){Stream_Gamma.close();}		
	}


	double falp(double x, double mu0, double sigma0, double sigma, double trunc)
	{
		double s0 = sigma*sigma + sigma0*sigma0;
		double s1 = 1/(sigma*sigma) + 1/(sigma0*sigma0);
		double s2 = x/(sigma*sigma) + mu0/(sigma0*sigma0);
		
		double a = dnorm(x, mu0, sqrt(s0), 0);
		double b = 1 - pnorm(trunc, s2/s1, 1/sqrt(s1) , 1, 0);
		double c = 1-pnorm(trunc, mu0, sigma0, 1, 0);
		
		return(a*b/c);
	}

	double faln(double x, double mu0, double sigma0, double sigma, double trunc)
	{
		double s0 = sigma*sigma + sigma0*sigma0;
		double s1 = 1/(sigma*sigma) + 1/(sigma0*sigma0);
		double s2 = x/(sigma*sigma) + mu0/(sigma0*sigma0);
		
		double a = dnorm(x, mu0, sqrt(s0), 0);
		double b = pnorm(-trunc, s2/s1, 1/sqrt(s1) , 1, 0);
		double c = pnorm(-trunc, mu0, sigma0, 1, 0);
		
		return(a*b/c);
	}

	void updateMembership(int g ,int s){	
		
		ParaList aparaList = bayesMPparaLists[s];
		int lengthAparaList = aparaList.getLength();
		int nSumP = aparaList.getParaSumNP();
		int nSumN = aparaList.getParaSumNN();
		int totalLength = lengthAparaList + 1 + 2;

		std::vector<double> poolYPr(totalLength, 0);
		std::vector<int> poolY(totalLength, 0);

		double aZ = Z[s*G + g];		

		for(int l=0;l<lengthAparaList;l++){
			Para apara = aparaList.paraList[l];
			poolY[l] = apara.getMembership();
			
			int n = apara.GetN();
			double postmu = apara.Getpostmu();
			double postsd = apara.Getpostsd();			
			
			if(poolY[l] > 0){				
				poolYPr[l] = falp(aZ, postmu, postsd, sigma, trunc) * n / (nSumP + alpha) * pi[g] * delta[g];												
			} else {
				poolYPr[l] = faln(aZ, postmu, postsd, sigma, trunc) * n / (nSumN + alpha) * pi[g] * (1 - delta[g]);	
			}
			
		} // for loop of l for lengthAparaList
		
		// 0: normal 0,1;
		poolY[lengthAparaList] = 0;
		// here null component is a standard normal distribution.
		poolYPr[lengthAparaList] = dnorm(aZ, empMu[s], empSD[s], 0) * (1 - pi[g]);
		

		int aNewMemPlus = aparaList.getNewMembership(1);
		int aNewMemMinus = aparaList.getNewMembership(-1);
				
		poolY[lengthAparaList + 1] = aNewMemPlus;
		poolYPr[lengthAparaList + 1] = falp(aZ, mu0, sigma0, sigma, trunc) * alpha / (nSumP + alpha) * pi[g] * delta[g];								
		poolY[lengthAparaList + 2] = aNewMemMinus;
		poolYPr[lengthAparaList + 2] = faln(aZ, mu0, sigma0, sigma, trunc) * alpha / (nSumN + alpha) * pi[g] * (1 - delta[g]);	
						
		discrete_distribution<int> distribution(poolYPr.begin(), poolYPr.end());
		int thisInt = distribution(generator);
		Y[s*G + g] = poolY[thisInt];

		if(thisInt == aNewMemPlus){
			bayesMPparaLists[s].updateNewMembership(1);
		} else if(thisInt == aNewMemMinus){
			bayesMPparaLists[s].updateNewMembership(-1);			
		}
		
	}
	
		
	void updatePi(){
		for(int g=0;g<G; g++){
			updateAPI(g);
		}
	}
	
	void updateAPI(int g){
		double Yplus = 0;
		double Yminus = 0;
		for(int s=0; s<S;s++){
			if(Y[s*G+g]>0){
				Yplus++;
			} else if(Y[s*G+g]<0){
				Yminus++;
			}
		}
		
		// for pi parameter;
		double pia = gamma + Yplus + Yminus;
		double pib = 1 - gamma + S - Yplus - Yminus;
		pi[g] = rbeta(pia, pib);	
		delta[g] = rbeta(beta + Yplus, beta + Yminus);
		
		//cout<<"gene "<<g<<". pia:"<<pia<<". pib:"<<pib <<". pig"<< pi[g] <<endl;
		//cout<<"gene "<<g<<". deltaa:"<<deltaa<<". deltab:"<<deltab <<". deltag"<< delta[g] <<endl;				
	}
	
	double binrarySearch(vector<double>& pig)
	{
	  double tol = 1.0/1e8;
	  double gammaLeft = tol;
	  double gammaRight = 1 - tol;
	  double gammaMiddle;
	  double A = 0.0;
  
	  int n = pig.size();
  
	  for(int i=0;i<n;i++){
	    A += (log(pig[i]) - log(1 - pig[i]))/n;
	  }
  
	  int iter = 0;
	  while(iter<=400 && (gammaRight-gammaLeft)>tol){	
		iter++;  
	    gammaMiddle = (gammaLeft + gammaRight)/2;
	    if(A > digamma(gammaMiddle) - digamma(1 - gammaMiddle)){
	      gammaLeft = gammaMiddle;
	    } else {
	      gammaRight = gammaMiddle;
	    }
	  }
	  //cout<<"iter: "<<iter<<endl;
	  //cout<<"gammaRight+gammaLeft)/2: "<<(gammaRight+gammaLeft)/2<<endl;
	  return (gammaRight+gammaLeft)/2;
	}
	
	long double loglikelihood(double gamma,vector<double>& pig)
	{
	  int G = pig.size();
	  long double loglikelihood = 0.0;
	  for(int g=0;g<G;g++){
	    loglikelihood += dbeta(pig[g],gamma,1-gamma,1);
	  }
	  return(loglikelihood);
	}
		
	void updateGamma(){
		double amu = binrarySearch(pi);
		double aprop = rnorm(amu,MHsd);
		//cout<<"loglikelihood(aprop, pi)"<<loglikelihood(aprop, pi)<<". loglikelihood(gamma, pi)"<<loglikelihood(gamma, pi)<<endl;
	    if(runif(0,1) < exp(loglikelihood(aprop, pi) - loglikelihood(gamma, pi))){
		  gamma = aprop;
		  countAcceptGamma++;
		  MHsd *= 1.01;
    	} else {
    	  MHsd /= 1.01;	
    	}	
	}
		
	void printAcceptRate(){
		cout << "mcmc accepted iter: " << countAcceptGamma	 <<endl;		
		cout << "mcmc final gamma: " << gamma	 <<endl;		
	}	
	
	bayesMP(){
		// cout<<"hi, I am constructing a BayesMP obj"<<endl;
	}	
};

void mcmc(int *G, int *S, double *Z, double *gamma, int *randomGamma, double *empMu, double *empSD, double *beta, double *alpha, double *mu0, 
		double *sigma0, double *sigma, double *atrunc, double *pi, double *delta, int *Y, int *niter, int *burnin, 
		int *silence, int *logDotsPerLine, 
		char *fileName_Y, char *fileName_Pi, char *fileName_Delta, char *fileName_Gamma, char *fileName_HSall, 
		int *writeY, int *writePi, int *writeDelta, int *writeGamma, int *writeHSall){
	
	bayesMP  mcmcobj;	
	mcmcobj.initialize(G,S,Z,gamma, randomGamma, empMu, empSD, beta, alpha, mu0, sigma0, sigma, atrunc, pi, delta, Y, niter, burnin, 
			fileName_Y, fileName_Pi, fileName_Delta, fileName_Gamma, 
			writeY, writePi, writeDelta, writeGamma, writeHSall);

	mcmcobj.iniPara();	
	

	
	for(int b=0;b < *niter;b++){
		if(!*silence){
			cout << "." << flush;
			if((b+1) % *logDotsPerLine == 0)
				cout << b+1 <<endl; 			
		}
		mcmcobj.iterateOne();		
		//mcmcobj->paraSPrint();	
		/*
		for(int s=0;s<3;s++){
			cout<<"c print study "<< s << "para: " << mcmcobj.bayesMPparaLists[s].getLength()<<endl;					
			cout<< "nSumP:" << mcmcobj.bayesMPparaLists[s].getParaSumNP() <<endl;
			cout<< "nSumN:" << mcmcobj.bayesMPparaLists[s].getParaSumNN() <<endl;
	    }
		*/
	}

	if(*writeHSall==1){mcmcobj.outputHSall(fileName_HSall);}
	// mcmcobj.printAcceptRate();
	//delete mcmcobj;
}

extern "C" {
	void mcmc_R3(int *G, int *S, double *Z, double *gamma, int *randomGamma, double *empMu, double *empSD, double *beta, double *alpha, 
		double *mu0, double *sigma0, double *sigma, double *atrunc, double *pi, double *delta, int *Y, 
		int *niter, int *burnin, int *silence, int *logDotsPerLine,
		char **fileName_Y, char ** fileName_Pi, char ** fileName_Delta, char ** fileName_Gamma, char ** fileName_HSall,
		int *writeY, int *writePi, int *writeDelta, int *writeGamma, int *writeHSall){
		mcmc(G, S, Z, gamma, randomGamma, empMu, empSD, beta, alpha, mu0, sigma0, sigma, atrunc, pi, delta, Y, niter, burnin, silence, logDotsPerLine, 
		*fileName_Y, *fileName_Pi, *fileName_Delta, *fileName_Gamma, *fileName_HSall, 
		writeY, writePi, writeDelta, writeGamma, writeHSall);
	}
}

