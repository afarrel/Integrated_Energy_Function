/*
Name: generate_DNA_struct.c
Author: Alvin Farrel
Description: Pipeline to prepare PDB file for energy calculations using MB, HB, and pi E-functions
Usage: ./generate_DNA_struct.o Your_File.pdb
Dependencies: Reduce, FIRST (FRODA), 3DNA, Get_Hbonds_energy.o, PDDOCK, NEW_Aromatic_Energy12.o
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char *argv[]) {

FILE *sequences;    //Sequence List File
FILE *complements;  //COmplement Sequence List FIle
FILE *Energy_List;  //Energy Results File
FILE *Low_Energy;   //Lowest Energy File

char PDBFILE[25];   //PDB Input File

sprintf(PDBFILE,"%s",argv[1]);


sequences = fopen("GeneratedSequences.txt","r");
complements = fopen("ComplementSequences.txt","r");

char EnergyListFile[70];
sprintf(EnergyListFile,"Energy_List_%s.txt",argv[1]);
Energy_List = fopen(EnergyListFile,"w");
fprintf(Energy_List,"Motif\tHBond_Energy\tAromatic\tPD_DOCK\n");  //Header For Result File
fclose(Energy_List);

system("bash exports.sh");

/////////////////////////////////////////////////////////
///////////////////////Reducd-Trim//////////////////////
///////////////////////////////////////////////////////
char ReduceCommand1[70];
sprintf(ReduceCommand1,"reduce.3.14.080821.linuxi386 -trim %s > 01_%s ",PDBFILE,PDBFILE);
puts(ReduceCommand1);
system(ReduceCommand1);
///////////////////////////////////////////////////
//////////////////////////////////////////////////



while(!feof(sequences))
{

char DNAseq[12];
char DNAseqComp[12];

fgets(DNAseq,12,sequences);
sscanf(DNAseq,"%s\n",DNAseq);   //Saves Sequences tp DNASeq
fgets(DNAseqComp,12,complements);
sscanf(DNAseqComp,"%s\n",DNAseqComp); //Saves complement dequences to DNAseqComp

int chainA_5;
int chainA_3;
int chainB_5;
int chainB_3;

char chainA;
char chainB;

chainA = 'C';
chainB = 'D';
chainA_5 = 1;
chainA_3 = 8;
chainB_5 = 12;
chainB_3 = 19;

int Sequence_Length=strlen(DNAseq); //the files have '\n' which is counted as a character
int x;
char sys_command[250];



///////////////////////////////////////////////////////
//////////////////////// - 3DNA - /////////////////////
///////////////////////////////////////////////////////
// Mutate the DNA to the DNASeq/DNAseqComp permuation//

sprintf(sys_command,"mutate_bases '");
for(x=0;x<Sequence_Length;x++)
{
	sprintf(sys_command,"%sc=%c s=%d m=D%c;c=%c s=%d m=D%c;",sys_command,chainA,chainA_5+x,DNAseq[x],chainB,chainB_3-x,DNAseqComp[x]);
}
sprintf(sys_command,"%s' 01_%s %s_%s.pdb",sys_command,PDBFILE,PDBFILE,DNAseq);

puts(sys_command);
system(sys_command);

/////////////////////////////////////////////////////////
//////////////////// - Reduce - ////////////////////////
///////////////////////////////////////////////////////
////////// Add Hydrogens and flip residues ///////////
char ReduceCommand2[100];
sprintf(ReduceCommand2,"reduce.3.14.080821.linuxi386 -flip -build %s_%s.pdb > %s_%s-reduced.pdb",PDBFILE,DNAseq,PDBFILE,DNAseq);
puts(ReduceCommand2);
system(ReduceCommand2);




//////////////////////////////////////////////////////
//////////////////////// - 3DNA - ////////////////////
//////////////////////////////////////////////////////
////////////// Repeated to Renumber Atoms ////////////

sprintf(sys_command,"mutate_bases '");
for(x=0;x<Sequence_Length;x++)
{
        sprintf(sys_command,"%sc=%c s=%d m=D%c;c=%c s=%d m=D%c;",sys_command,chainA,chainA_5+x,DNAseq[x],chainB,chainB_3-x,DNAseqComp[x]);
}
sprintf(sys_command,"%s' %s_%s-reduced.pdb  %s_%s.pdb",sys_command,PDBFILE,DNAseq,PDBFILE,DNAseq);


puts(sys_command);
system(sys_command);


/////////////////////////////////////////////////////////////////
/////////////////////// - FIRST/FRODA - ////////////////////////
///////////////////////////////////////////////////////////////
/////////////////// Determine Hbond Data//////////////////////
char FIRST_conformations[200];
sprintf(FIRST_conformations,"FIRST -L $FIRST_LIB -non -FRODA -maxfitc 500 -totconf 1 -freq 1 -dtol 0.3 %s_%s.pdb",PDBFILE,DNAseq);

puts(FIRST_conformations);
system(FIRST_conformations);
printf("%s:%s\n", DNAseq,DNAseqComp);




//////////////////////////////////////////////////////////////////////////
/////////////////////// - Get Energy - //////////////////////////////////
////////////////////////////////////////////////////////////////////////
/// Runs Get_Hbonds_energy.o -> Gets PDDOCk, Arom, and HBOND Energy /// 

char file_checks[60];
sprintf(file_checks,"%s_%s_froda_00000001.pdb",PDBFILE,DNAseq);

if(access(file_checks,F_OK)==0)
{
char Get_Energy[140];
sprintf(Get_Energy,"./Get_Hbonds_energy.o %s_%s.pdb 1 1",PDBFILE,DNAseq);

puts(Get_Energy);
system(Get_Energy);



char EnergyLine[60];
Low_Energy = fopen("LowEnergy.dat","r");
fgets(EnergyLine,50,Low_Energy);
fclose(Low_Energy);

//Add Energies to Results file//
Energy_List = fopen(EnergyListFile,"a");
fprintf(Energy_List,"%s\t%s",DNAseq,EnergyLine);
fclose(Energy_List);
}
else
{
Energy_List = fopen(EnergyListFile,"a");
fprintf(Energy_List,"%s\t0.000000\t0.000000\t0.000000\n",DNAseq);
fclose(Energy_List);
}


system("rm *.pdb_*"); //remove files produced by programs
}

fclose(sequences);
fclose(complements);

return EXIT_SUCCESS;
}
