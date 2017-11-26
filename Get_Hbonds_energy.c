#include <stdio.h>
#include <stdlib.h>
#include <string.h>


double getEnergy(char *,char *);

int main(int argc, char *argv[]) {


char filename[70];
char filename2[70];
char hbondFile[]="hbonds.out";

sprintf(filename,"%s",argv[1]);
sprintf(filename2,"%s",argv[1]);
int steps;
int startconf;
int totalconf;
char sys_command[170];
char froda_files[170];


int dot;
int x,y;
int filename_length;

double Lowest_Energy=100;
double hbond_Energy;

///////aromatic energy varaiables///
FILE *arom_energy;
double arom_energy_temp;
double Lowest_arom_energy=100;
char Arom_Energy_cmd[170];

///////PD_DOCK energy varaiables///
FILE *PD_DOCK_energy;
double PD_DOCK_energy_temp;
double Lowest_PD_DOCK_energy=100;
char PD_DOCK_Energy_cmd[170];

////////////////////////////////////

steps=atoi(argv[2]);
totalconf=atoi(argv[3]);

filename_length=strlen(filename);

puts(filename);

for(x=filename_length-1;x>=0;x--)
	if (filename[x]=='.')
		{dot=x;x=-1;}

filename[dot]='\0';
puts(filename);


system("rm hbonds.out");
for(x=steps;x<=totalconf;x+=steps)
{
	if(x<10)sprintf(froda_files,"%s_froda_0000000%d.pdb",filename,x);
	if((x<100)&&(x>=10))sprintf(froda_files,"%s_froda_000000%d.pdb",filename,x);
	if((x<1000)&&(x>=100))sprintf(froda_files,"%s_froda_00000%d.pdb",filename,x);
	if((x<10000)&&(x>=1000))sprintf(froda_files,"%s_froda_0000%d.pdb",filename,x);

	sprintf(sys_command,"FIRST -L $FIRST_LIB -E -0.1 -non -hbout %s",froda_files);

printf("froda files:%s\n",froda_files);
///////////////////////////////////////////////////////
//////////////// - Aromatic Energy - /////////////////
/////////////////////////////////////////////////////

Arom_Energy_cmd[0]='\0';
sprintf(Arom_Energy_cmd,"./NEW_Aromatic_Energy12.o %s > arom_energy.txt",froda_files);
system(Arom_Energy_cmd);


arom_energy = fopen("arom_energy.txt","r");
fscanf(arom_energy,"%lf",&arom_energy_temp);

if (arom_energy_temp < Lowest_arom_energy)
	Lowest_arom_energy = arom_energy_temp;

fclose(arom_energy);

/////////////////////////////////////////////////////////////
///////////////////// - PD_DOCK Energy - ///////////////////
///////////////////////////////////////////////////////////

PD_DOCK_Energy_cmd[0]='\0';



sprintf(PD_DOCK_Energy_cmd,"reduce.3.14.080821.linuxi386 -trim %s > noH%s",froda_files,froda_files); 
printf("%s",PD_DOCK_Energy_cmd);
system(PD_DOCK_Energy_cmd);



PD_DOCK_Energy_cmd[0]='\0';

sprintf(PD_DOCK_Energy_cmd,"./CalculateBindingEnergies.o noH%s",froda_files);
system(PD_DOCK_Energy_cmd);
printf("%s",PD_DOCK_Energy_cmd);


PD_DOCK_energy = fopen("PDDEnergy.txt","r");

fscanf(PD_DOCK_energy,"%lf",&PD_DOCK_energy_temp);

if (PD_DOCK_energy_temp < Lowest_PD_DOCK_energy)
        Lowest_PD_DOCK_energy = PD_DOCK_energy_temp;

fclose(PD_DOCK_energy);


//////////////////////////////////////////////////////////
/////////////////// - Hbonds - //////////////////////////
////////////////////////////////////////////////////////

system(sys_command);
sys_command[0]='\0';


sprintf(sys_command,"cp %s TempPdbfile.pdb",froda_files);
system(sys_command);
sys_command[0]='\0';


hbond_Energy=getEnergy(froda_files,"hbonds.out");
printf("\n%lf\n",hbond_Energy);



sprintf(sys_command,"rm %s_*",froda_files);
system(sys_command);
sys_command[0]='\0';

if(hbond_Energy<Lowest_Energy)
	Lowest_Energy=hbond_Energy;


}


FILE *LowestENERGY;
LowestENERGY = fopen("LowEnergy.dat","w");


char sequence[10];
fprintf(LowestENERGY,"%lf\t%lf\t%lf\n",Lowest_Energy,Lowest_arom_energy,Lowest_PD_DOCK_energy);
fclose(LowestENERGY);


return EXIT_SUCCESS;
}


double getEnergy(char *pdbFile,char *HBondFile)    {


	FILE *PDB;
	FILE *HBond_out;

	char Line[190];
	char atom[6];
	float aNUM;
	char aNAME[4];
	char rNAME[4];
	float rNUM;
	char chain;

	int DNAMax=0,DNAMin=10000000, ProtMax=0,ProtMin=10000000;

	PDB =fopen(pdbFile,"r");
	HBond_out=fopen(HBondFile,"r");


///////	//new stuff//
	int Prot;	
	int flexRes;
	int NRes;
	int CARes;
	int CRes;
	double flex;
	double flex1;
	double flex2;
	int atm1,atm2;
	FILE *BondFlexFile;
	char bond_flex[150];
	char bondflex[150];

	sprintf(bond_flex,"%s",pdbFile);
	int x;
	for(x=strlen(bond_flex);x>0;x--)
	{
		if(bond_flex[x]=='.')
		{
		bond_flex[x]='\0';
		x=-1;
		}
	}
	sprintf(bondflex,"%s_bond.txt",bond_flex);
	
///////


	while (!feof(PDB))
	{

	fgets(Line,150,PDB);

	sscanf(Line,"%s %f %s %s",atom, &aNUM,aNAME,rNAME);

	if ((atom[0]=='A')&&(rNAME[0]=='D')&&(aNUM<DNAMin))DNAMin=aNUM;
	if ((atom[0]=='A')&&(rNAME[0]=='D')&&(aNUM>DNAMax))DNAMax=aNUM;
	if ((atom[0]=='A')&&(rNAME[0]!='D')&&(aNUM<ProtMin))ProtMin=aNUM;
	if ((atom[0]=='A')&&(rNAME[0]!='D')&&(aNUM>ProtMax))ProtMax=aNUM;

	}
	fclose(PDB);
	

	printf("DNA:%d-%d; Prot:%d-%d",DNAMin,DNAMax,ProtMin,ProtMax);


	int Hydrogen=0,Acceptor=0;
	double Total_HBEnergy=0;
	double HBondEnergy=0;


	while (!feof(HBond_out))
	{


		fgets(Line,150,HBond_out);
		sscanf(Line," %d %d %lf",&Hydrogen,&Acceptor,&HBondEnergy);


		if(((Hydrogen>=DNAMin)&&(Hydrogen<=DNAMax)&&(Acceptor>=ProtMin)&&(Acceptor<=ProtMax))
				||
				((Acceptor>=DNAMin)&&(Acceptor<=DNAMax)&&(Hydrogen>=ProtMin)&&(Hydrogen<=ProtMax)))
		{

		/////////New stuff//////
			
			
			if((Acceptor>=ProtMin)&&(Acceptor<=ProtMax))Prot=Acceptor;
			if((Hydrogen>=ProtMin)&&(Hydrogen<=ProtMax))Prot=Hydrogen;
			
			aNUM=0;
			PDB =fopen(pdbFile,"r");
			while(Prot!=aNUM)
			{
				fgets(Line,150,PDB);
        			sscanf(Line,"%s %f %s %s %c %f",atom, &aNUM,aNAME,rNAME,&chain,&rNUM);
			
				//if(Prot==aNUM)printf("yay %f",rNUM);
			}
			fclose(PDB);

			flexRes=rNUM;

			PDB =fopen(pdbFile,"r");
			while((!feof(PDB))&&(rNUM<=flexRes))
			{
				fgets(Line,150,PDB);
                                sscanf(Line,"%s %f %s %s %c %f",atom, &aNUM,aNAME,rNAME,&chain,&rNUM);

				if((rNUM==flexRes)&&(strcmp(aNAME,"N")==0)) NRes=aNUM;
                                if((rNUM==flexRes)&&(strcmp(aNAME,"CA")==0)) CARes=aNUM;
                                if((rNUM==flexRes)&&(strcmp(aNAME,"C")==0)) CRes=aNUM;

			}
			fclose(PDB);		
			//printf("\n%d %d %d %d %d\n%lf\n",Prot,NRes,CARes,CRes,HBondEnergy,flexRes);
			
			BondFlexFile=fopen(bondflex,"r");
	
			flex1=1;
			flex2=1;
			for(x=0;x<4;x++) fgets(Line,150,BondFlexFile);
			while(!feof(BondFlexFile))
			{
				fgets(Line,150,BondFlexFile);
                               	sscanf(Line,"%d %d %lf",&atm1,&atm2,&flex);

				if((atm1==NRes)&&(atm2==CARes)) flex1=flex;
                                if((atm1==CARes)&&(atm2==CRes)) flex2=flex;
			
			}
			fclose(BondFlexFile);
		/////////

		
			//if(HBondEnergy<-6.2)
			//if ((flex1<0.2)&&(flex2<0.2))
			Total_HBEnergy+=HBondEnergy;
		}

	}

	//printf("\n%lf",Total_HBEnergy);

	//fclose(PDB);
	
	if (Total_HBEnergy>0)Total_HBEnergy=0;

	fclose(HBond_out);
	return Total_HBEnergy;
}

