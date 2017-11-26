/*
 ============================================================================
 Name        : CalculateBindingEnergies.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
		puts("!!!Hello World!!!"); /* prints !!!Hello World!!! */

		
		char PdbFile[70];
		char line[100];


	FILE *EnergyList;

	sprintf(PdbFile,"%s",argv[1]);



	char command[100];
	sprintf(command,"pd_dock_setup %s_out -f %s",PdbFile,PdbFile);

	system("rm structure.list");
	system(command);


	FILE *Input_Pdb;
	FILE *Input_Temp;
	char filename[60];

	sprintf(filename,"input.%s_out",PdbFile);

	Input_Pdb=fopen(filename,"r");
	Input_Temp=fopen("input.Temp","w");

	while(!feof(Input_Pdb))
		{
		fgets(line,100,Input_Pdb);

		if ((line[0]=='D')&&(line[1]=='O')&&(line[2]=='_')&&(line[3]=='C'))
			sprintf(line,"DO_COMPUTE_ENERGY\t1\n");

		if ((line[0]=='D')&&(line[1]=='O')&&(line[2]=='_')&&(line[3]=='D'))
			sprintf(line,"DO_DOCKING\t\t0\n");

		fprintf(Input_Temp,"%s",line);
		}
	fclose(Input_Pdb);
	fclose(Input_Temp);

	system("pd_dock input.Temp > TempEnergy.txt");

	FILE *TempEnergy;

	char word1[30],word2[30],word3[30];
	float energy=0.0;


	TempEnergy=fopen("TempEnergy.txt","r");
	EnergyList=fopen("PDDEnergy.txt","w");


	fscanf(TempEnergy,"%s %s %s %f\n%s %f \t %s %f\n %s %f %s %f\n%s %f \t %s %f\n%s %s %f",word1,word2,word3,&energy,word3,&energy,word3,&energy,word3,&energy,word3,&energy,word3,&energy,word3,&energy,word1,word3,&energy);
	
fprintf(EnergyList,"%f",energy);



	fclose(TempEnergy);
	fclose(EnergyList);




	return EXIT_SUCCESS;
}
