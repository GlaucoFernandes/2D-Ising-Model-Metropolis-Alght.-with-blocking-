#include <iostream>
#include <math.h>
#include <random>
#include <ostream>
#include <fstream>

const static int NetLenght = 10;
const static int NumLat = pow(NetLenght, 2);
const static int MedPerBlock = 1000, blocks = 10, meds = blocks*MedPerBlock;// No de Medidas por bloco, blocos e medidas 



double J = 1; //Constante de acoplamento dos spins (positiva para ferromagnetismo, negativa para antiferromagnetismo)

int getIndex(int i, int j){
	return j + i*NetLenght;
	}

//////////////// Definições referentes à geradores de números aleatórios
std::mt19937 eng{ 42 };
std::uniform_real_distribution<double> dist{ 0, 1.0 };

double r_rand() {
	return dist(eng);
}
/////////////////////////////////////////


///////////////// Definições relacionadas ao tratamento dos dados

double med(double* dados, int array_lenght){
	double sum = 0;
	
	for(int i = 0; i < array_lenght; i++){
		sum += dados[i];
	}
	
	return sum/array_lenght;
}

double desvpad(double* dados, int array_lenght){
	double sum = 0, media_dados = med(dados, array_lenght);
	
	for(int i = 0; i < array_lenght; i++){
		sum += pow((dados[i] - media_dados), 2);
	}
	
	return sqrt(sum/(array_lenght - 1));
}

class Net
{
private:
	int Rede[NumLat]; //A rede propriamente dita
	int Rede_mp[NetLenght][NetLenght]; //Uma matrix com ponteiros para os números
	
	double RawData[meds][4]; //Matrix para medida da energia, en^2, magnetixação e mag^2
	//double Blocos[blocks][4]; //Matrix com as medidas de En, Mag, Cv, Susc para cada bloco
	double Energy_d[blocks];
	double Magnet_d[blocks];
	double Cv_d[blocks];
	double Susc_d[blocks];
	
	int Viz[NumLat][4]; //Matrix de Vizinhos
	
	double T, Tm;
	int TermPass = 1000; //Passos para a termalização
	
public:

	Net(double Temp, double Tmax)
		:T(Temp), Tm(Tmax)
	{
		this->NetInit();
		this->getViz();
		this->Termalization();
		this->TakeMed();
		this->ProcMed();
	}

	void NetInit(){ //função coloca todos os spins da rede para cima
		for(int i = 0; i < NumLat; i++){
			Rede[i] = 1;
		}
	}
	
	void getViz(){ //Função que preenche a matrix de Vizinhos
	
		for(int i = 0; i < NetLenght; i++){
			for(int j = 0; j < NetLenght; j++){ //Percorre uma matrix imaginária
				int v1 = j + 1;
				int v2 = i - 1;
				int v3 = j - 1;
				int v4 = i + 1;
				
				//condições de contorno periódicas
				
				if(v1 == NetLenght){
					v1 -= NetLenght;
				}
				if(v2 < 0){
					v2 += NetLenght; 
				}
				if(v3 < 0){
					v3 += NetLenght;
				}
				if(v4 == NetLenght){
					v4 -= NetLenght;
				}
				
				Viz[getIndex(i, j)][0] = getIndex(i, v1);
				Viz[getIndex(i, j)][1] = getIndex(v2, j);
				Viz[getIndex(i, j)][2] = getIndex(i, v3);
				Viz[getIndex(i, j)][3] = getIndex(v4, j);
				
			}
			
			
			
		}
	}
	
	double LatEnergy(int i){ //Calcula a energia relacionada às ligações de um sitio
		return -J*Rede[i]*(Rede[Viz[i][0]] + Rede[Viz[i][1]] + Rede[Viz[i][2]] + Rede[Viz[i][3]]);
	}
	
	void MetropolisStep(){ // Define o passo de metropolis
		for(int i = 0; i < NumLat; i++){
			double DE = -2*(this->LatEnergy(i));
			
			if(DE <= 0){
				Rede[i] = -Rede[i];
			}
			
			else{
				double z = exp(-DE/T), r = r_rand();
				
				if(z >= r){
					Rede[i] = -Rede[i];
				}
			}
			
		}
	}
	
	void Termalization(){ // Faz a termalização da rede
		for(int i = 0; i < TermPass; i++){
			this->MetropolisStep();
		}
	}
	
	double NetEnergy(){ //Calcula a energia atual da rede
		double sum = 0;
		
		for(int i = 0; i < NumLat; i++){
			sum += LatEnergy(i);
		}
		
		return sum/2;
		}
	
	double NetMag(){ //calcula a magnetização atual da rede
		double sum = 0;
		
		for(int i = 0; i < NumLat; i++){
			sum += Rede[i];
		}
		
		return abs(sum);
	}
	
	void TakeMed(){ //Toma as medidas
		for(int i = 0; i < meds; i++){
			RawData[i][0] = this->NetEnergy();
			RawData[i][1] = pow((this->NetEnergy()), 2);
			RawData[i][2] = this->NetMag();
			RawData[i][3] = pow((this->NetMag()), 2);
			
			this->MetropolisStep();
		}
	}
	
	void ProcMed(){ //Divide nos blocos
		
		int cont = 0;
		int b_i;
		double MedEn = 0, MedEnSq = 0; 
		double MedMag = 0, MedMagSq = 0;
		
		for(int i = 0; i < meds; i++){
			MedEn += RawData[i][0]/MedPerBlock;
			MedEnSq += RawData[i][1]/MedPerBlock;
			MedMag += RawData[i][2]/MedPerBlock;
			MedMagSq += RawData[i][3]/MedPerBlock;
			
			if(cont == 999){
				cont = -1;
				b_i = floor(i/MedPerBlock);
				Energy_d[b_i] = MedEn;
				Magnet_d[b_i] = MedMag;
				Cv_d[b_i] = (MedEnSq - pow(MedEn, 2))/pow(T, 2);
				Susc_d[b_i] = (MedMagSq - pow(MedMag, 2))/(T);
				
				MedEn = 0;
				MedEnSq = 0;
				MedMag = 0;
				MedMagSq = 0;
				}
			cont += 1;
		}
	}

	double getEnergy(){
		return med(Energy_d, blocks)/NumLat;
	}
	
	double getEnergyError(){
		return desvpad(Energy_d, blocks)/NumLat;
	}
	
	double getMagnet(){
		return med(Magnet_d, blocks)/NumLat;
	}
	
	double getMagnetError(){
		
		return desvpad(Magnet_d, blocks)/NumLat;
	} 

	double getCv(){
		return med(Cv_d, blocks)/NumLat;
	}
	
	double getCvError(){
		return desvpad(Cv_d, blocks)/NumLat;
	} 

	double getSusc(){
		return med(Susc_d, blocks)/NumLat;
	}
	
	double getSuscError(){
		return desvpad(Susc_d, blocks)/NumLat;
	} 
};

void ExportData(){
	
	std::ofstream myfile;
	myfile.open("Results/Data/Data.csv");
	
	for(double i = 0.1; i < 2; i+= 0.3){
		Net Ising(i, 5);
		
		myfile << i << "," << Ising.getEnergy() << "," << Ising.getEnergyError() << ",";
		myfile << Ising.getMagnet() << "," << Ising.getMagnetError() << ",";
		myfile << Ising.getCv() << "," << Ising.getCvError() << ",";
		myfile << Ising.getSusc() << "," << Ising.getSuscError() << std::endl;
		
		
	}
	
		for(double i = 2; i < 3; i+= 0.05){
		Net Ising(i, 5);
		
		myfile << i << "," << Ising.getEnergy() << "," << Ising.getEnergyError() << ",";
		myfile << Ising.getMagnet() << "," << Ising.getMagnetError() << ",";
		myfile << Ising.getCv() << "," << Ising.getCvError() << ",";
		myfile << Ising.getSusc() << "," << Ising.getSuscError() << std::endl;
		
		
	}
	
		for(double i = 3; i < 5; i+= 0.3){
		Net Ising(i, 5);
		
		myfile << i << "," << Ising.getEnergy() << "," << Ising.getEnergyError() << ",";
		myfile << Ising.getMagnet() << "," << Ising.getMagnetError() << ",";
		myfile << Ising.getCv() << "," << Ising.getCvError() << ",";
		myfile << Ising.getSusc() << "," << Ising.getSuscError() << std::endl;
	}
	
	myfile.close();
}

int main() 
{
	ExportData();
	
	return 0;
}