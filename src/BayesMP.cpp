#include "RCbridge.h"

#include <cstring>
#include <fstream>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <random>
#include <ctime>

using namespace std;

class indexC
{
	// a very simple dynamic link
	int anum;
public:
	indexC *left, *right;
	
	indexC()
	{
		left=NULL;
		right=NULL;		
	}
	
	indexC(int bnum)
	{
		anum = bnum;
		left=NULL;
		right=NULL;		
	}
	
	void setNum(int an){
		anum = an;
	}
	
	int getNum(){
		return(anum);
	}
	
};

class para{
	// membership: DP class: 1,2,3,...: positive DE. -1,-2,-3...: negative DE. 0: nonDE.
	// index: gene index.
	// n: total number of genes.
	int membership, n;
	// sumZ: sum of Z statistics
	// postmu: posterior mu
	// postsd: posterior sd
	double sumZ, postmu, postsd, mu0, sigma0, sigma;
	
public:
	// link
	para *left, *right;
	indexC * index;
	
	para()
	{
		left=NULL;
		right=NULL;
	}
	
	para(double amu0, double asigma0, double asigma, int anum, double aZ, int amembership)
	{
		membership=amembership;
		mu0=amu0;
		sigma0=asigma0;
		sigma=asigma;
		index = new indexC(anum);
		sumZ = aZ;
		n = 1;		
		updateTheta();	
		left=NULL;
		right=NULL;
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
		indexC * thisIndex = index;
		while(thisIndex->right != NULL){
			thisIndex = thisIndex->right;
		}
		thisIndex->right = new indexC(anum);
		thisIndex->right->left = thisIndex;
		sumZ = sumZ + aZ;
		n = n + 1;
		updateTheta();
	}

	int removeZ(int anum, double aZ)
	{		
		// 0: nothing left. 1: successfully. 2; error: no such anum exist.
		indexC * thisIndex = index;
		int headPoint = 0;
		
		while(thisIndex != NULL){
			if(thisIndex->getNum()==anum)
			{
				if(thisIndex->left != NULL){
					thisIndex->left->right = thisIndex->right;						
				}
				
				if(thisIndex->right != NULL){
					thisIndex->right->left = thisIndex->left;							
				}
				if(headPoint==0){
					index = thisIndex->right;
				}										
				delete [] thisIndex;
				sumZ = sumZ - aZ;
				n = n - 1;
				if(n==0)
				{					
					return 0;
				}
				updateTheta();
				return 1;
			}
			thisIndex = thisIndex->right;
			headPoint = headPoint + 1;
		}
		return 2;
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

	void freePara()
	{
		indexC * thisIndex;		
		while(index != NULL){
			thisIndex = index->right;
			delete [] index;
			index = thisIndex;
		}		
	}

	
/*	
	~para(){
		cout<<"such index doesn't exist, bug 2"<<endl;		
		indexC * thisIndex;		
		cout<<"such index doesn't exist, bug 3"<<endl;		
		while(index != NULL){
			cout<<"such index doesn't exist, bug 5"<<endl;			
			thisIndex = index->right;
			cout<<"such index doesn't exist, bug 6"<<endl;			
			delete index;
			cout<<"such index doesn't exist, bug 7"<<endl;			
			index = thisIndex;
			cout<<"such index doesn't exist, bug 8"<<endl;			
		}
	}
*/	
};

class bayesMP{
	// dimension variables G: number of genes. S: number of studies.
	int G, S;
	// Z: observed Z statistics
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
	// whether directly output HSind result and HSall result
	int fullRes, HSall;
	// number of mcmc iterations.
	int niter, burnin;
	// current iter
	int thisIter;
	char * fileFullRes;
	char * fileHSall;

	double MHsd = 0.1;
	
	std::vector<int> YHSall;
	ofstream myStream;	
	default_random_engine generator;
	
 	void SetG(int g){
		G = g;
 	}
	
 	void SetS(int s){
		S = s;
 	}

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
	
 	void SetfullRes(int afullRes){
		fullRes = afullRes;
 	}

 	void SetHSall(int aHSall){
		HSall = aHSall;
 	}
		
	void iniYHSall(){
		YHSall = std::vector<int>(G*S, 0);				
	}

	void SetfullFilename(char *filename){
	    fileFullRes = new char[strlen(filename)+8];
	    strcpy(fileFullRes,filename);
	    strcat(fileFullRes,"Full.txt");		
	}

	void SetHSallFileame(char *filename){
	    fileHSall = new char[strlen(filename)+9];
	    strcpy(fileHSall,filename);
	    strcat(fileHSall,"HSall.txt");		
	}
	
public:
	para ** paraObjS;						
	
	void initialize(int *aG, int *aS, double *aZ, double *agamma, int *randomGamma, double *aempMu, double *aempSD, double *abeta,double *aalpha ,double *amu0, double *asigma0, double *asigma, double *atrunc, double *api, double *adelta, int *aY, int *niter, int *burnin, char *filename, int *fullRes, int *aHSall)
	{
		SetG(*aG);
		SetS(*aS);
		SetZ(aZ);
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
		SetfullRes(*fullRes);
		SetHSall(*aHSall);
		SetfullFilename(filename);
		SetHSallFileame(filename);	
		iniYHSall();
		thisIter = 0;								
	}
	
	void updatePara()
	{
		int findFlag;
		paraObjS = new  para * [S];
		for(int s=0;s<S;s++)
		{			
			para * sparaObj = NULL;
			para * sparaPointer;
			para * sparaPointer0;
			for(int g=0;g<G;g++)
			{
				int sGg = s*G+g;
				if(Y[sGg]!=0){
					if(sparaObj==NULL){
						sparaObj = new para(mu0,sigma0,sigma,sGg,Z[sGg],Y[sGg]);	
						paraObjS[s] = sparaObj;
					}
					else
					{
						findFlag = 0;
						sparaPointer = sparaObj;
						while(sparaPointer!=NULL)
						{
							if(Y[sGg]==sparaPointer->getMembership())
							{
								sparaPointer->addZ(sGg, Z[sGg]);
								findFlag = 1;
								break;
							}
							sparaPointer0 = sparaPointer;
							sparaPointer = sparaPointer->right;						
						}
						if(findFlag==0)
						{
							sparaPointer0->right = new para(mu0,sigma0,sigma,sGg,Z[sGg],Y[sGg]);
							sparaPointer0->right->left = sparaPointer0;
						}
					}
				}
			}
		}
	}
	
	void deletePara(int g, int s){
		int sGg = s*G+g;
		int aY = Y[sGg];
		int headPoint = 0;
		if(aY!=0){
			//cout<<"aY"<<aY<<endl;
			
			para * sparaPointer = paraObjS[s];
			while(sparaPointer!=NULL)
			{
				
				if(sparaPointer->getMembership() == aY)
				{
					int removeZStatus = sparaPointer->removeZ(sGg,Z[sGg]);
					if(removeZStatus==1)
					{
						return;
					} else if(removeZStatus==0)
					{
						if(sparaPointer->left != NULL){
							sparaPointer->left->right = sparaPointer->right;						
						}
						
						if(sparaPointer->right != NULL){
							sparaPointer->right->left = sparaPointer->left;							
						}
						if(headPoint==0){
							paraObjS[s] = sparaPointer->right;
						}						
						delete [] sparaPointer;						
						return;
						
					} else if(removeZStatus==2)
					{
						cout<<"such index doesn't exist, bug 2"<<endl;
						exit(0);
					} else {
						cout<<"unknown error, bug 0!"<<endl;
						exit(0);
					}										
				}
				sparaPointer = sparaPointer->right;
				headPoint = headPoint + 1;
			}
			cout<<"no para class match current class label, bug 1!"<<endl;
			exit(0);		
		}
	}

	void addPara(int g, int s){		
		//cout<<"g:"<<g<<". s:"<<s<<endl;
		int findFlag = 0;
		int sGg=s*G+g;
		
		if(Y[sGg]==0){
			return;
		}
		
		para * sparaPointer = paraObjS[s];
		para * sparaPointer0 = sparaPointer;
		
		while(sparaPointer!=NULL)
		{
			if(Y[sGg]==sparaPointer->getMembership())
			{
				sparaPointer->addZ(sGg, Z[sGg]);
				findFlag = 1;
				break;
			}
			sparaPointer0 = sparaPointer;
			sparaPointer = sparaPointer->right;						
		}
		if(findFlag==0)
		{
			sparaPointer0->right = new para(mu0,sigma0,sigma,sGg,Z[sGg],Y[sGg]);
			sparaPointer0->right->left = sparaPointer0;
		}

	}	

	char * GetfullFilename(){
		return(fileFullRes);
	}

	char * GetHSallFileame(){
		return(fileHSall);
	}

	int getNewMembership(int s, int direction){
		para * sparaPointer0 = paraObjS[s];
		para * sparaPointer = sparaPointer0;
		if(direction==1){
			int i = 1;
			while(sparaPointer!=NULL){
				if(sparaPointer->getMembership() == i){
					i++;
					sparaPointer = sparaPointer0;				
				} else {
					sparaPointer = sparaPointer->right;	
				}			
			}
			return i;
		} else if(direction==-1) {
			int i = -1;
			while(sparaPointer!=NULL){
				if(sparaPointer->getMembership() == i){
					i--;
					sparaPointer = sparaPointer0;				
				} else {
					sparaPointer = sparaPointer->right;	
				}			
			}
			return i;
		}
		cout << "error, couldn't create new membership, bug 4"<<endl;
		return 0;
	}

	int getParaLength(int s){
		para * sparaPointer = paraObjS[s];
		int count = 0;
		while(sparaPointer!=NULL)
		{
			count++;
			sparaPointer = sparaPointer->right;
		}
		return count;
	}

	int getParaSumNP(int s){
		para * sparaPointer = paraObjS[s];
		int sumN = 0;
		while(sparaPointer!=NULL)
		{
			if(sparaPointer->getMembership()>0){
				sumN += sparaPointer->GetN();				
			}
			sparaPointer = sparaPointer->right;
		}
		return sumN;
	}

	int getParaSumNN(int s){
		para * sparaPointer = paraObjS[s];
		int sumN = 0;
		while(sparaPointer!=NULL)
		{
			if(sparaPointer->getMembership()<0){
				sumN += sparaPointer->GetN();				
			}
			sparaPointer = sparaPointer->right;
		}
		return sumN;
	}

	void iterateOne() {
		for(int g=0; g<G; g++){
			for(int s=0; s<S; s++){
				updateOne(g,s);
			}
		}
		updatePi();
		updateHSall();	
		if(randomGamma == 1){
			updateGamma();			
		}

		if(fullRes == 1){appendFile(myStream, thisIter);}
		thisIter++;
			
	}		
		
	void updateOne(int g, int s) {
		int start_s;
		int stop_s;
					 
		start_s=clock();
		deletePara(g, s);
		stop_s=clock();
		cout << "time deletePara: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
		start_s=clock();
		updateMembership(g ,s);
		stop_s=clock();
		cout << "time updateMembership: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
		start_s=clock();
		addPara(g, s);
		stop_s=clock();
		cout << "time addPara: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
	}
	
	void updateHSall(){
		if(thisIter<burnin){
			return;
		} 
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

	void appendFile(ofstream& myStream, int thisIter){
		if(thisIter==0){myStream.open(fileFullRes);}
		for(int g=0;g<G;g++){
			myStream << pi[g] << "\t";
		}

		for(int g=0;g<G;g++){
			myStream << delta[g] << "\t";
		}

		for(int i=0;i<G*S;i++){
			myStream << Y[i] << "\t";
		}
				
	    myStream << endl;
		if(thisIter==niter-1){myStream.close();}
		
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
		para * sparaPointer = paraObjS[s];	
		int sparalength = getParaLength(s);	
		int nSumP = getParaSumNP(s);
		int nSumN = getParaSumNN(s);
		int totalLength = sparalength + 1 + 2;

		std::vector<double> poolYPr(totalLength, 0);
		std::vector<int> poolY(totalLength, 0);
		
		double aZ = Z[s*G + g];
		
		// 0: normal 0,1;
		poolY[0] = 0;
		// here null component is a standard normal distribution.
		poolYPr[0] = dnorm(aZ, empMu[s], empSD[s], 0) * (1 - pi[g]);
		
		for(int i=1;i<=sparalength;i++)
		{
			poolY[i] = sparaPointer->getMembership();
			
			int n = sparaPointer->GetN();
			double postmu = sparaPointer->Getpostmu();
			double postsd = sparaPointer->Getpostsd();			
					
			if(poolY[i] > 0){				
				poolYPr[i] = falp(aZ, postmu, postsd, sigma, trunc) * n / (nSumP + alpha) * pi[g] * delta[g];												
			} else {
				poolYPr[i] = faln(aZ, postmu, postsd, sigma, trunc) * n / (nSumN + alpha) * pi[g] * (1 - delta[g]);	
			}
			sparaPointer = sparaPointer->right;
		}
		poolY[sparalength + 1] = getNewMembership(s,1);
		poolYPr[sparalength + 1] = falp(aZ, mu0, sigma0, sigma, trunc) * alpha / (nSumP + alpha) * pi[g] * delta[g];								
		poolY[sparalength + 2] = getNewMembership(s,-1);
		poolYPr[sparalength + 1] = faln(aZ, mu0, sigma0, sigma, trunc) * alpha / (nSumN + alpha) * pi[g] * (1 - delta[g]);	
				
		discrete_distribution<int> distribution(poolYPr.begin(), poolYPr.end());
		int thisInt = distribution(generator);
		Y[s*G + g] = poolY[thisInt];
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
		cout<<"hi, I am constructing a BayesMP obj"<<endl;
	}
	
	void paraSPrint(){
		for(int s=0; s<S; s++){
			para * thisparaobj = paraObjS[s];
			while(thisparaobj!=NULL)
			{
				cout<<"study"<<s<<"GetMember: "<<thisparaobj->getMembership()<<". GetSumZ: "<<thisparaobj->GetSumZ()<<". GetN: "<<thisparaobj->GetN()<<endl;				
				thisparaobj = thisparaobj->right;
			}
		}
		cout<<endl;
	}
	
	void freeBayesMP()
	{
		para * sparaPointer;
		para * sparaPointerD;
		for(int s=0; s<S; s++){
			sparaPointer = paraObjS[s];
			sparaPointerD = sparaPointer;
			while(sparaPointer != NULL){
				sparaPointer = sparaPointerD->right;
				sparaPointerD->freePara();
				delete [] sparaPointerD;
				sparaPointerD = sparaPointer;		
			}					
		}				
	}
	
	~bayesMP(){
		delete [] fileFullRes;
		delete [] fileHSall;
		delete [] paraObjS;		
	}
	
	
};

void mcmc(int *G, int *S, double *Z, double *gamma, int *randomGamma, double *empMu, double *empSD, double *beta, double *alpha, double *mu0, double *sigma0, double *sigma, double *atrunc, double *pi, double *delta, int *Y, int *niter, int *burnin, char *filename , int *fullRes, int *HSall){

		
	bayesMP * mcmcobj = new bayesMP;
	mcmcobj->initialize(G,S,Z,gamma, randomGamma, empMu, empSD, beta, alpha, mu0, sigma0, sigma, atrunc, pi, delta, Y, niter, burnin, filename, fullRes, HSall);
	mcmcobj->updatePara();
	 	
	for(int b=0;b < *niter;b++){
		mcmcobj->iterateOne();		
		//mcmcobj->paraSPrint();	
		cout << "mcmc iter: " << b <<endl;
	}

	if(*HSall==1){mcmcobj->outputHSall(mcmcobj->GetHSallFileame());}		
	mcmcobj->printAcceptRate();
	mcmcobj->freeBayesMP();
	delete mcmcobj;		
}

extern "C" {
	void mcmc_R3(int *G, int *S, double *Z, double *gamma, int *randomGamma, double *empMu, double *empSD, double *beta, double *alpha, double *mu0, double *sigma0, double *sigma, double *atrunc, double *pi, double *delta, int *Y, int *niter, int *burnin, char **filename, int *fullRes,int *HSall){
		mcmc(G, S, Z, gamma, randomGamma, empMu, empSD, beta, alpha, mu0, sigma0, sigma, atrunc, pi, delta, Y, niter, burnin, *filename, fullRes, HSall);
	}	
}



