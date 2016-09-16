#include "RCbridge.h"

#include <cstring>
#include <fstream>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <random>

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
	double *Z;
	// gamma: prior for pi_g; beta: prior for delta_g; alpha: hyperparameter for DP; mu0: mean para for G0; sigma0: sd para for G0; sigma: sd for emission;
	double gamma, beta, alpha, mu0, sigma0, sigma, trunc;
	// pi: prob gene g to be DE. delta: prob d DE gene to be positive.
	double *pi, *delta;
	// record average number of component
	double *aveComponent;	
	// Y: latent variable to be infered.
	int *Y;
	// whether directly output HSind result and HSall result
	int fullRes, HSall, aveCom;
	// number of mcmc iterations.
	int niter, burnin;
	// current iter
	int thisIter;
	char * fileFullRes;
	char * fileHSall;
	char * fileAveCom;
	
	int *YHSall;
	ofstream myStream;	
	default_random_engine generator;
	
 	void SetG(int g){
		G = g;
 	}
	
 	void SetS(int s){
		S = s;
 	}

 	void SetZ(double *aZ){
		int GS = G*S;
		Z = new double [GS];		
		for(int i=0; i<GS; i++){
			Z[i] = *aZ++;
		}
 	}

 	void SetPi(double *api){
		pi = new double [G];		
		for(int g=0; g<G; g++){
			pi[g] = *api++;
		}
 	}

 	void SetDelta(double *adelta){
		delta = new double [G];		
		for(int g=0; g<G; g++){
			delta[g] = *adelta++;
		}
 	}

 	void SetY(int *aY){
		int GS = G*S;		
		Y = new int [GS];		
		for(int i=0; i<GS; i++){
			Y[i] = *aY++;
		}
 	}

 	void SetGamma(double agamma){
		gamma = agamma;		
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
	
 	void SetAveCom(int aaveCom){
		aveCom = aaveCom;
 	}
	
	void iniYHSall(){
		YHSall = new int [G*S];
		for(int i=0;i<G*S;i++){
			YHSall[i] = 0;
		}
		
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
	
	void SetAveComFileame(char *filename){
	    fileAveCom = new char[strlen(filename)+10];
	    strcpy(fileAveCom,filename);
	    strcat(fileAveCom,"AveCom.txt");		
	}
	
public:
	para ** paraObjS;						
	
	void initialize(int *aG, int *aS, double *aZ, double *agamma, double *abeta,double *aalpha ,double *amu0, double *asigma0, double *asigma, double *atrunc, double *api, double *adelta, int *aY, int *niter, int *burnin, char *filename, int *fullRes, int *aHSall, int *aveCom)
	{
		SetG(*aG);
		SetS(*aS);
		SetZ(aZ);
		SetGamma(*agamma);
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
		SetAveCom(*aveCom);		
		SetfullFilename(filename);
		SetHSallFileame(filename);	
		SetAveComFileame(filename);			
		iniYHSall();
		thisIter = 0;								
		aveComponent = new double [*niter];									
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
		para * sparaPointer0;
		
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

	char * GetAveComFileame(){
		return(fileAveCom);
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
		GetAveComponent();		
		if(fullRes == 1){appendFile(myStream, thisIter);}
		thisIter++;
			
	}		
	
	void GetAveComponent(){
		aveComponent[thisIter] = 0;
		for(int s=0;s<S;s++){
			aveComponent[thisIter] += ((double) getParaLength(s))/S;	
			//cout<<"iter: "<< thisIter << ". " << "s: " << s << ". numComponent" <<getParaLength(s) <<endl	;	
		}
	}
	
	void updateOne(int g, int s) {
		deletePara(g, s);
		updateMembership(g ,s);
		addPara(g, s);
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

	void outputAveCom(char *myfile){
	    ofstream streamAveCom;
	    streamAveCom.open(myfile);		
		double count = 0;
		for(int i=0;i<niter;i++){
			if(i<burnin){
				continue;
			}
			count += aveComponent[i];
		}
		count = count / (niter - burnin);			
	    streamAveCom << count;						
	    streamAveCom << endl;						
		streamAveCom.close();		
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
		int tracei;	
		int sparalength = getParaLength(s);	
		int nSumP = getParaSumNP(s);
		int nSumN = getParaSumNN(s);
		int totalLength = sparalength + 1 + 2;
		double * poolYPr = new double [totalLength]();
		int * poolY = new int [totalLength]();
		
		double aZ = Z[s*G + g];
		
		// 0: normal 0,1;
		poolY[0] = 0;
		// here null component is a standard normal distribution.
		poolYPr[0] = dnorm(aZ, 0, 1, 0) * (1 - pi[g]);
		
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
			tracei = i;
			sparaPointer = sparaPointer->right;
		}
		poolY[++tracei] = getNewMembership(s,1);
		poolYPr[tracei] = falp(aZ, mu0, sigma0, sigma, trunc) * alpha / (nSumP + alpha) * pi[g] * delta[g];								
		poolY[++tracei] = getNewMembership(s,-1);
		poolYPr[tracei] = faln(aZ, mu0, sigma0, sigma, trunc) * alpha / (nSumN + alpha) * pi[g] * (1 - delta[g]);	
		
	    vector<double> vectorPr;
		for(int i=0;i<totalLength;i++){
			vectorPr.push_back(poolYPr[i]);
		}
		
		discrete_distribution<int> distribution(vectorPr.begin(), vectorPr.end());
		int thisInt = distribution(generator);
		Y[s*G + g] = poolY[thisInt];

		delete [] poolYPr;
		delete [] poolY;
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
		double pia = gamma/(G - gamma) + Yplus + Yminus;
		//double pia = gamma/G  + Yplus + Yminus;
		double pib = 1 + S - Yplus - Yminus;
	
		gamma_distribution<double> distGammaPia(pia,1);
		gamma_distribution<double> distGammaPib(pib,1);
		
		double piN = distGammaPia(generator);
		double piD = distGammaPib(generator);
		pi[g] = piN/(piN + piD);
		
		// for delta parameter;
		
		double deltaa = beta + Yplus;
		double deltab = beta + Yminus;
	
		gamma_distribution<double> distGammaDeltaa(deltaa,1);
		gamma_distribution<double> distGammaDeltab(deltab,1);
		
		double deltaN = distGammaDeltaa(generator);
		double deltaD = distGammaDeltab(generator);
		delta[g] = deltaN/(deltaN + deltaD);		
		
		//cout<<"gene "<<g<<". pia:"<<pia<<". pib:"<<pib <<". pig"<< pi[g] <<endl;
		//cout<<"gene "<<g<<". deltaa:"<<deltaa<<". deltab:"<<deltab <<". deltag"<< delta[g] <<endl;				
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
		delete [] Z;
		delete [] Y;
		delete [] pi;
		delete [] delta;
		delete [] YHSall;
		delete [] fileFullRes;
		delete [] fileHSall;
		delete [] fileAveCom;
		delete [] aveComponent;		
		delete [] paraObjS;		
	}
	
	
};

void mcmc(int *G, int *S, double *Z, double *gamma, double *beta, double *alpha, double *mu0, double *sigma0, double *sigma, double *atrunc, double *pi, double *delta, int *Y, int *niter, int *burnin, char *filename , int *fullRes, int *HSall, int *aveCom){

		
	bayesMP * mcmcobj = new bayesMP;
	mcmcobj->initialize(G,S,Z,gamma, beta, alpha, mu0, sigma0, sigma, atrunc, pi, delta, Y, niter, burnin, filename, fullRes, HSall, aveCom);
	mcmcobj->updatePara();
	
    
	
	for(int b=0;b < *niter;b++){
		mcmcobj->iterateOne();		
		//mcmcobj->paraSPrint();	
		cout << "mcmc iter: " << b <<endl;
	}

	if(*HSall==1){mcmcobj->outputHSall(mcmcobj->GetHSallFileame());}
	if(*aveCom==1){mcmcobj->outputAveCom(mcmcobj->GetAveComFileame());}
		
	mcmcobj->freeBayesMP();
	delete mcmcobj;
		
}

extern "C" {
	void mcmc_R(int *G, int *S, double *Z, double *gamma, double *beta, double *alpha, double *mu0, double *sigma0, double *sigma, double *atrunc, double *pi, double *delta, int *Y, int *niter, int *burnin, char **filename, int *fullRes,int *HSall, int *aveCom){
		mcmc(G, S, Z, gamma, beta, alpha, mu0, sigma0, sigma, atrunc, pi, delta, Y, niter, burnin, *filename, fullRes, HSall, aveCom);
	}
}



