/*
 ============================================================================
 Name        : Aromatic_Energy12.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glob.h>

double aromatic_energy(char pdb[]);

int main(int argc, char *argv[])
{
    double test;
    char pdb[60];
    
    sprintf(pdb,"%s",argv[1]);
    
    test=aromatic_energy(pdb);
    
    printf("%lf",test);
    
    
    
    
    return EXIT_SUCCESS;
    
}


double aromatic_energy(char pdb[]) {

	FILE *pdbfile;

	char Line[150];


	int Bases=150;
	char atom[6];
	int ANum;
	char AName[5];
	char RName[4];
	char Chain;
	int res_num;
	double xcor;
	double ycor;
	double zcor;
	double C_Array[Bases][40]; //25 and 50
	double A_Array[Bases][40];
	double Tyr_Array[Bases][40];
	char Base_order[Bases];
	int Base_num[Bases];
	int max_base_num=0;
	int base_count=0;
	char chain_check='X';

	int C_count=0;
	int T_count=0;
	int A_count=0;
	int other_atoms;


	//Nearest Neighbor Parameters for  - SantaLucia et al 2006
	double CA = -1.38;
	double AC = -1.43;
	double CT = -1.16;
	double TC = -1.46;
	double CG = -2.09;
	double GC = -2.28;
	double CC = -1.77;
	double AA = -1.02;
	double TT = -1.02;
	double AT = -0.73;
	double TA = -0.60;
	double TG = -1.38;
	double GT = -1.43;
	double AG = -1.16;
	double GA = -1.46;
	double GG = -1.77;




	pdbfile=fopen(pdb,"r");


	while(!feof(pdbfile))
	{
		fgets(Line,150,pdbfile);
		if(Line[0]=='A')
		{
			sscanf(Line,"%s %d %s %s %c %d %lf %lf %lf",atom,&ANum,AName,RName,&Chain,&res_num,&xcor,&ycor,&zcor);

			if(RName[0]=='D')
			if(Chain!=chain_check)
			  {
				max_base_num=0;
				Base_order[base_count]='X';
				Base_num[base_count]=0;
				//for()
				base_count++;
				chain_check=Chain;
			  }

			if(RName[0]=='D')
				if(res_num>max_base_num)
			{
				Base_order[base_count]=RName[1];
				Base_num[base_count]=res_num;
				max_base_num=res_num;
				base_count++;
				other_atoms=25;

			}

			if ((strcmp(RName,"DC")==0)||(strcmp(RName,"C")==0)||(strcmp(RName,"DT")==0)||(strcmp(RName,"T")==0))
			if(RName[0]=='D')
			{
				if(strcmp(AName,"N1")==0)
				{
					C_Array[C_count][0]=xcor;
					C_Array[C_count][1]=ycor;
					C_Array[C_count][2]=zcor;

				}
				else
				if(strcmp(AName,"C2")==0)
				{
					C_Array[C_count][3]=xcor;
					C_Array[C_count][4]=ycor;
					C_Array[C_count][5]=zcor;
				}
				else
				if(strcmp(AName,"N3")==0)
				{
					C_Array[C_count][6]=xcor;
					C_Array[C_count][7]=ycor;
					C_Array[C_count][8]=zcor;
				}
				else
				if(strcmp(AName,"C4")==0)
				{
					C_Array[C_count][9]=xcor;
					C_Array[C_count][10]=ycor;
					C_Array[C_count][11]=zcor;
				}
				else
				if(strcmp(AName,"C5")==0)
				{
					C_Array[C_count][12]=xcor;
					C_Array[C_count][13]=ycor;
					C_Array[C_count][14]=zcor;

					if ((strcmp(RName,"DC")==0)||(strcmp(RName,"C")==0))
					C_Array[C_count][27]=1;
				if ((strcmp(RName,"DT")==0)||(strcmp(RName,"T")==0))
					C_Array[C_count][27]=2;
				C_Array[C_count][28]=base_count-1;   //base_count-1 because after base_count is assigned its incremented.

				}
				else
				if(strcmp(AName,"C6")==0)
				{
					C_Array[C_count][15]=xcor;
					C_Array[C_count][16]=ycor;
					C_Array[C_count][17]=zcor;
					if ((strcmp(RName,"DC")==0)||(strcmp(RName,"C")==0)||(strcmp(RName,"DT")==0)||(strcmp(RName,"T")==0))C_count++;

				}
				else
				if(strcmp(AName,"N4")==0)
				{
					C_Array[C_count][18]=xcor;
					C_Array[C_count][19]=ycor;
					C_Array[C_count][20]=zcor;
				}
				else
				if(strcmp(AName,"O4")==0)
				{
					C_Array[C_count][21]=xcor;
					C_Array[C_count][22]=ycor;
					C_Array[C_count][23]=zcor;
				}

			}

			if ((strcmp(RName,"DA")==0)||(strcmp(RName,"A")==0)||(strcmp(RName,"DG")==0)||(strcmp(RName,"G")==0))
						//if(RName[0]=='D')
						{
							if(strcmp(AName,"C4")==0)
							{
								A_Array[A_count][0]=xcor;
								A_Array[A_count][1]=ycor;
								A_Array[A_count][2]=zcor;
								A_count++;
							}
							else
							if(strcmp(AName,"C5")==0)
							{
								A_Array[A_count][3]=xcor;
								A_Array[A_count][4]=ycor;
								A_Array[A_count][5]=zcor;
							}
							else
							if(strcmp(AName,"N7")==0)
							{
								A_Array[A_count][6]=xcor;
								A_Array[A_count][7]=ycor;
								A_Array[A_count][8]=zcor;
							}
							else
							if(strcmp(AName,"C8")==0)
							{
								A_Array[A_count][9]=xcor;
								A_Array[A_count][10]=ycor;
								A_Array[A_count][11]=zcor;

							}
							else
							if(strcmp(AName,"N9")==0)
							{
								A_Array[A_count][12]=xcor;
								A_Array[A_count][13]=ycor;
								A_Array[A_count][14]=zcor;
							}
							else
							if(strcmp(AName,"N1")==0)
							{
								A_Array[A_count][15]=xcor;
								A_Array[A_count][16]=ycor;
								A_Array[A_count][17]=zcor;

								if ((strcmp(RName,"DA")==0)||(strcmp(RName,"A")==0))
									A_Array[A_count][27]=1;

								if ((strcmp(RName,"DG")==0)||(strcmp(RName,"G")==0))
									A_Array[A_count][27]=2;

								A_Array[A_count][28]=base_count-1;

							}
							else
							if(strcmp(AName,"C2")==0)
							{
								A_Array[A_count][18]=xcor;
								A_Array[A_count][19]=ycor;
								A_Array[A_count][20]=zcor;
							}

							else
							if(strcmp(AName,"N6")==0)
							{
								A_Array[A_count][21]=xcor;
								A_Array[A_count][22]=ycor;
								A_Array[A_count][23]=zcor;
							}
							else
							if(strcmp(AName,"O6")==0)
							{
								A_Array[A_count][24]=xcor;
								A_Array[A_count][25]=ycor;
								A_Array[A_count][26]=zcor;

							}



						}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Aromatic Residues~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~Tyrosine~~~~~*/
			if (strcmp(RName,"TYR")==0)
			{
				if(strcmp(AName,"CG")==0)
				{
					Tyr_Array[T_count][0]=xcor;
					Tyr_Array[T_count][1]=ycor;
					Tyr_Array[T_count][2]=zcor;
					Tyr_Array[T_count][22]=1;
				}
				if(strcmp(AName,"CD1")==0)
				{
					Tyr_Array[T_count][3]=xcor;
					Tyr_Array[T_count][4]=ycor;
					Tyr_Array[T_count][5]=zcor;
				}
				if(strcmp(AName,"CE1")==0)
				{
					Tyr_Array[T_count][6]=xcor;
					Tyr_Array[T_count][7]=ycor;
					Tyr_Array[T_count][8]=zcor;
				}
				if(strcmp(AName,"CZ")==0)
				{
					Tyr_Array[T_count][9]=xcor;
					Tyr_Array[T_count][10]=ycor;
					Tyr_Array[T_count][11]=zcor;
				}
				if(strcmp(AName,"CE2")==0)
				{
					Tyr_Array[T_count][12]=xcor;
					Tyr_Array[T_count][13]=ycor;
					Tyr_Array[T_count][14]=zcor;
				}
				if(strcmp(AName,"CD2")==0)
				{
					Tyr_Array[T_count][15]=xcor;
					Tyr_Array[T_count][16]=ycor;
					Tyr_Array[T_count][17]=zcor;
				}
				if(strcmp(AName,"OH")==0)
				{
					Tyr_Array[T_count][18]=xcor;
					Tyr_Array[T_count][19]=ycor;
					Tyr_Array[T_count][20]=zcor;
					Tyr_Array[T_count][21]=res_num;
					T_count++;
				}

			}

			/*~~~~Phenylalanine~~~~~*/
			if (strcmp(RName,"PHE")==0)
			{
				if(strcmp(AName,"CG")==0)
				{
					Tyr_Array[T_count][0]=xcor;
					Tyr_Array[T_count][1]=ycor;
					Tyr_Array[T_count][2]=zcor;
					Tyr_Array[T_count][22]=2;
				}
				if(strcmp(AName,"CD1")==0)
				{
					Tyr_Array[T_count][3]=xcor;
					Tyr_Array[T_count][4]=ycor;
					Tyr_Array[T_count][5]=zcor;
				}
				if(strcmp(AName,"CE1")==0)
				{
					Tyr_Array[T_count][6]=xcor;
					Tyr_Array[T_count][7]=ycor;
					Tyr_Array[T_count][8]=zcor;
				}
				if(strcmp(AName,"CZ")==0)
				{
					Tyr_Array[T_count][9]=xcor;
					Tyr_Array[T_count][10]=ycor;
					Tyr_Array[T_count][11]=zcor;
					T_count++;
				}
				if(strcmp(AName,"CE2")==0)
				{
					Tyr_Array[T_count][12]=xcor;
					Tyr_Array[T_count][13]=ycor;
					Tyr_Array[T_count][14]=zcor;
				}
				if(strcmp(AName,"CD2")==0)
				{
					Tyr_Array[T_count][15]=xcor;
					Tyr_Array[T_count][16]=ycor;
					Tyr_Array[T_count][17]=zcor;
				}


			}

			/*~~~~Tryptophan~~~~~*/
			if (strcmp(RName,"TRP")==0)
			{
				if(strcmp(AName,"CG")==0)
				{
					Tyr_Array[T_count][0]=xcor;
					Tyr_Array[T_count][1]=ycor;
					Tyr_Array[T_count][2]=zcor;
					Tyr_Array[T_count][22]=3;
				}
				if(strcmp(AName,"NE1")==0)
				{
					Tyr_Array[T_count][3]=xcor;
					Tyr_Array[T_count][4]=ycor;
					Tyr_Array[T_count][5]=zcor;
				}
				if(strcmp(AName,"CD2")==0)
				{
					Tyr_Array[T_count][6]=xcor;
					Tyr_Array[T_count][7]=ycor;
					Tyr_Array[T_count][8]=zcor;
				}
				if(strcmp(AName,"CE2")==0)
				{
					Tyr_Array[T_count][9]=xcor;
					Tyr_Array[T_count][10]=ycor;
					Tyr_Array[T_count][11]=zcor;
				}
				if(strcmp(AName,"CZ3")==0)
				{
					Tyr_Array[T_count][12]=xcor;
					Tyr_Array[T_count][13]=ycor;
					Tyr_Array[T_count][14]=zcor;
				}
				if(strcmp(AName,"CH2")==0)
				{
					Tyr_Array[T_count][15]=xcor;
					Tyr_Array[T_count][16]=ycor;
					Tyr_Array[T_count][17]=zcor;
					T_count++;
				}
				if(strcmp(AName,"CD1")==0)
				{
					Tyr_Array[T_count][18]=xcor;
					Tyr_Array[T_count][19]=ycor;
					Tyr_Array[T_count][20]=zcor;
					Tyr_Array[T_count][21]=res_num;
				}


			}
			/*~~~~Histidine~~~~~*/
			if (strcmp(RName,"HIS")==0)
			{
				if(strcmp(AName,"CG")==0)
				{
					Tyr_Array[T_count][0]=xcor;
					Tyr_Array[T_count][1]=ycor;
					Tyr_Array[T_count][2]=zcor;
					Tyr_Array[T_count][22]=4;
				}
				if(strcmp(AName,"ND1")==0)
				{
					Tyr_Array[T_count][3]=xcor;
					Tyr_Array[T_count][4]=ycor;
					Tyr_Array[T_count][5]=zcor;
				}
				if(strcmp(AName,"CE1")==0)
				{
					Tyr_Array[T_count][6]=xcor;
					Tyr_Array[T_count][7]=ycor;
					Tyr_Array[T_count][8]=zcor;
				}
				if(strcmp(AName,"NE2")==0)
				{
					Tyr_Array[T_count][9]=xcor;
					Tyr_Array[T_count][10]=ycor;
					Tyr_Array[T_count][11]=zcor;
					T_count++;
				}
				if(strcmp(AName,"CD2")==0)
				{
					Tyr_Array[T_count][12]=xcor;
					Tyr_Array[T_count][13]=ycor;
					Tyr_Array[T_count][14]=zcor;
				}

				 if(strcmp(AName,"CG")==0)
				{
					Tyr_Array[T_count][15]=xcor;
					Tyr_Array[T_count][16]=ycor;
					Tyr_Array[T_count][17]=zcor;

				}



			}

/*~~~~~~~~~~~~~~~~~~~~~~~END AROMATIC RESIDUES~~~~~~~~~~~~~~~~~~~~~~*/


		}

	}
	Base_order[base_count]='X';
	fclose(pdbfile);






int x,y,z,w;
double D[7], TempD[7],closest;

double Total_Energy=0;
int num;



double XMax,XMin;
double YMax,YMin;
double ZMax,ZMin;

double TXMax,TXMin;
double TYMax,TYMin;
double TZMax,TZMin;



int Pass=1;
int Fail=0;
int a;


double T_central[5][3];

double T_centroid[3];
int l,m;



	for(x=0;x<T_count;x++)
	{
		for(y=0;y<C_count;y++)

		{

			num=C_Array[y][28];
		if((Base_order[num]!='X')&&(Base_order[num+1]!='X')&&(Base_order[num-1]!='X'))
		{
			for(z=0;z<6;z++)
			{
				closest=1000;

				for(w=0;w<6;w++)
				{
					TempD[w]=(Tyr_Array[x][z*3+0]-C_Array[y][w*3+0])*(Tyr_Array[x][z*3+0]-C_Array[y][w*3+0]) + (Tyr_Array[x][z*3+1]-C_Array[y][w*3+1])*(Tyr_Array[x][z*3+1]-C_Array[y][w*3+1]) + (Tyr_Array[x][z*3+2]-C_Array[y][w*3+2])*(Tyr_Array[x][z*3+2]-C_Array[y][w*3+2]);
					TempD[w]=sqrt(TempD[w]);

					if (TempD[w]<closest)
					{
							closest=TempD[w];
					}

				}

				D[z]=closest;
			}

	///////////////
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			/*~~~~~~ Calculate Centroid of Aromatic Residue ~~~~~~*/
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


if((Tyr_Array[x][22]==1)||(Tyr_Array[x][22]==2))
{
			for(l=0;l<3;l++)
				for(m=0;m<3;m++)
				{
					T_central[l][m]=Tyr_Array[x][l*3+m]+(Tyr_Array[x][(3+l)*3+m]-Tyr_Array[x][l*3+m])/2;
				}

			for(l=0;l<3;l++)
			{
				T_centroid[l]=(T_central[0][l]+T_central[1][l]+T_central[2][l])/3;
			}

}//if((Tyr_Array[T_count][22]==1)||(Tyr_Array[T_count][22]==2))

if(Tyr_Array[x][22]==3)
{

				for(m=0;m<3;m++)
				{
					T_central[0][m]=Tyr_Array[y][4*3+m]+(Tyr_Array[y][5*3+m]-Tyr_Array[y][4*3+m])/2;
				}
				for(m=0;m<3;m++)


			for(l=0;l<3;l++)
			{
				T_centroid[l]=Tyr_Array[y][6*3+l]+(T_central[0][l]-Tyr_Array[y][6*3+l])/2;
			}

}//if(Tyr_Array[T_count][22]==3)


if(Tyr_Array[x][22]==4)
{

	for(l=0;l<3;l++)
		for(m=0;m<3;m++)
		{
		T_central[l][m]=Tyr_Array[x][l*3+m] +
				((Tyr_Array[x][(3+l)*3+m]+(Tyr_Array[x][(3+l)*3+m]-Tyr_Array[x][(2+l)*3+m])/2) -
						Tyr_Array[x][l*3+m])/2;
		}
	for(l=3;l<5;l++)
			for(m=0;m<3;m++)
			{
			T_central[l][m]=Tyr_Array[x][l*3+m] +
					((Tyr_Array[x][(-3+l)*3+m]+(Tyr_Array[x][(-3+l)*3+m]-Tyr_Array[x][(-2+l)*3+m])/2) -
							Tyr_Array[x][l*3+m])/2;
			}

	for(l=0;l<3;l++)
	{
		T_centroid[l]=(T_central[0][l]+T_central[1][l]+T_central[2][l]+T_central[3][l]+T_central[4][l])/5;
	}


}//if(Tyr_Array[T_count][22]==4)

////////////
////////////


double distance,closestAtom[3]={1000,1000,1000};
	int closestC[3]={100,100,100};

	for(m=0;m<=24;m+=3)
		{
			distance=(T_centroid[0]-C_Array[y][m])*(T_centroid[0]-C_Array[y][m])+
					(T_centroid[1]-C_Array[y][m+1])*(T_centroid[1]-C_Array[y][m+1])+
					(T_centroid[2]-C_Array[y][m+2])*(T_centroid[2]-C_Array[y][m+2]);
			distance=sqrt(distance);

			if(distance<closestAtom[0])
			{
				closestAtom[2]=closestAtom[1];
				closestAtom[1]=closestAtom[0];
				closestAtom[0]=distance;
				closestC[2]=closestC[1];
				closestC[1]=closestC[0];
				closestC[0]=m;
			}

			if((distance>closestAtom[0])&&(distance<closestAtom[1]))
			{
				closestAtom[2]=closestAtom[1];
				closestAtom[1]=distance;
				closestC[2]=closestC[1];

				closestC[1]=m;
			}
			if((distance>closestAtom[1])&&(distance<closestAtom[2]))
			{
				closestAtom[2]=distance;
				closestC[2]=m;
			}

		}




for(a=0;a<=15;a+=3)
{
	if(a==0)
	{
		XMax=Tyr_Array[x][a];
		XMin=Tyr_Array[x][a];
		YMax=Tyr_Array[x][a+1];
		YMin=Tyr_Array[x][a+1];
		ZMax=Tyr_Array[x][a+2];
		ZMin=Tyr_Array[x][a+2];
	}

	{
		if(Tyr_Array[x][a]>XMax)XMax=Tyr_Array[x][a];
		if(Tyr_Array[x][a]<XMin)XMin=Tyr_Array[x][a];
		if(Tyr_Array[x][a+1]>YMax)YMax=Tyr_Array[x][a+1];
		if(Tyr_Array[x][a+1]<YMin)YMin=Tyr_Array[x][a+1];
		if(Tyr_Array[x][a+2]>ZMax)ZMax=Tyr_Array[x][a+2];
		if(Tyr_Array[x][a+2]<ZMin)ZMin=Tyr_Array[x][a+2];

		if(C_Array[y][a]>XMax)XMax=C_Array[y][a];
		if(C_Array[y][a]<XMin)XMin=C_Array[y][a];
		if(C_Array[y][a+1]>YMax)YMax=C_Array[y][a+1];
		if(C_Array[y][a+1]<YMin)YMin=C_Array[y][a+1];
		if(C_Array[y][a+2]>ZMax)ZMax=C_Array[y][a+2];
		if(C_Array[y][a+2]<ZMin)ZMin=C_Array[y][a+2];
	}
	if(a==0)
	{


		TXMax=C_Array[y][closestC[0]];
		TXMin=C_Array[y][closestC[0]];
		TYMax=C_Array[y][closestC[0]+1];
		TYMin=C_Array[y][closestC[0]+1];
		TZMax=C_Array[y][closestC[0]+2];
		TZMin=C_Array[y][closestC[0]+2];


		for(l=1;l<3;l++)
		{
			if(C_Array[y][closestC[l]]>TXMax)TXMax=C_Array[y][closestC[l]];
			if(C_Array[y][closestC[l]]<TXMin)TXMin=C_Array[y][closestC[l]];
			if(C_Array[y][closestC[l]]>TYMax)TYMax=C_Array[y][closestC[l]+1];
			if(C_Array[y][closestC[l]]<TYMin)TYMin=C_Array[y][closestC[l]+1];
			if(C_Array[y][closestC[l]]>TZMax)TZMax=C_Array[y][closestC[l]+2];
			if(C_Array[y][closestC[l]]<TZMin)TZMin=C_Array[y][closestC[l]+2];
		}
	}

	{
		if(Tyr_Array[x][a]>TXMax)TXMax=Tyr_Array[x][a];
		if(Tyr_Array[x][a]<TXMin)TXMin=Tyr_Array[x][a];
		if(Tyr_Array[x][a+1]>TYMax)TYMax=Tyr_Array[x][a+1];
		if(Tyr_Array[x][a+1]<TYMin)TYMin=Tyr_Array[x][a+1];
		if(Tyr_Array[x][a+2]>TZMax)TZMax=Tyr_Array[x][a+2];
		if(Tyr_Array[x][a+2]<TZMin)TZMin=Tyr_Array[x][a+2];
	}
}






			Pass=1;
			Fail=0;

			pdbfile=fopen(pdb,"r");
			while((!feof(pdbfile))&&(Pass==1))
			{
				fgets(Line,150,pdbfile);
				if(Line[0]=='A')
				{
					sscanf(Line,"%s %d %s %s %c %d %lf %lf %lf",atom,&ANum,AName,RName,&Chain,&res_num,&xcor,&ycor,&zcor);

					if(AName[0]!='H')
					{

						if((TXMax>xcor)&&(xcor>TXMin)&&(TYMax>ycor)&&(ycor>TYMin)&&(TZMax>zcor)&&(zcor>TZMin))
						{
							Pass=0;
							for(a=0;a<=15;a+=3)
							{
								if((Tyr_Array[x][a]==xcor)&&(Tyr_Array[x][a+1]==ycor)&&(Tyr_Array[x][a+2]==zcor))
								{
									Pass=1;
									Fail=0;
								}
								if((C_Array[y][a]==xcor)&&(C_Array[y][a+1]==ycor)&&(C_Array[y][a+2]==zcor))
								{
									Pass=1;
									Fail=0;
								}
							}
						}
					}

				}
			}
			fclose(pdbfile);


			double avogadros=6.02214129*pow(10,23);
			double kconstant=8.9875517873681764*pow(10,9);  //Coloumbs Constant
			double EEnergy;
			double dielectric=4; //dielectric constant in protein environment
			double electron = 1.602176565*pow(10,-19);
			double Dist_Check;
			double coulomb_meter=3.336*pow(10,-30);
			double tyr_charge;

			if(Tyr_Array[x][22]==1) tyr_charge = -0.372 * 26.9/27.1;
			if(Tyr_Array[x][22]==2) tyr_charge = -0.372 * 27.1/27.1;
			if(Tyr_Array[x][22]==3) tyr_charge = -0.372 * 32.6/27.1;
			if(Tyr_Array[x][22]==4) tyr_charge = -0.372 * 21.0/27.1;


			if(C_Array[y][27]==1)
			{

				double N1,C2,O2,N3,C4,N4,C5,C6;
				double C_ENeg;
				double N4_dist;
				double C5_dist;
				double C6_dist;
				double N4_charge;
				double C5_charge;
				double C6_charge;
				double C,N,O;


				O=3.44;
				N=3.04;
				C=2.55;

				N1=5.697;
				C6=6.157;
				C5=4.904;
				C4=6.738;
				N4=3.867;
				N3=5.939;
				C2=7.002;
				O2=5.322;

				N4_charge = 0.157 + 0.157;
				C5_charge = 0.066;
				C6_charge = 0.085;
	



				N4_dist=sqrt((C_Array[y][18]-T_centroid[0])*(C_Array[y][18]-T_centroid[0])+
						(C_Array[y][19]-T_centroid[1])*(C_Array[y][19]-T_centroid[1])+
						(C_Array[y][20]-T_centroid[2])*(C_Array[y][20]-T_centroid[2]));

				C5_dist=sqrt((C_Array[y][12]-T_centroid[0])*(C_Array[y][12]-T_centroid[0])+
						(C_Array[y][13]-T_centroid[1])*(C_Array[y][13]-T_centroid[1])+
						(C_Array[y][14]-T_centroid[2])*(C_Array[y][14]-T_centroid[2]));

				C6_dist=sqrt((C_Array[y][15]-T_centroid[0])*(C_Array[y][15]-T_centroid[0])+
						(C_Array[y][16]-T_centroid[1])*(C_Array[y][16]-T_centroid[1])+
						(C_Array[y][17]-T_centroid[2])*(C_Array[y][17]-T_centroid[2]));


				//Potential energy = coloumbs constant x avogrado's number x charge1 x charge2/dielectric constant x distance =J/mol

				EEnergy=0;
				EEnergy+= kconstant*avogadros*(N4_charge*electron)*(tyr_charge*electron)/(dielectric*(N4_dist*pow(10,-10)));
				EEnergy+= kconstant*avogadros*(C5_charge*electron)*(tyr_charge*electron)/(dielectric*(C5_dist*pow(10,-10)));
				EEnergy+= kconstant*avogadros*(C6_charge*electron)*(tyr_charge*electron)/(dielectric*(C6_dist*pow(10,-10)));
				//used be -0.088

				EEnergy = EEnergy/1000;  // (J/mol to kJ/mol)
				EEnergy= EEnergy/4.18400; // kJ/mol to kcal/mol)


				Dist_Check=N4_dist;
				if(Dist_Check>C5_dist) Dist_Check=C5_dist;
				if(Dist_Check>C6_dist) Dist_Check=C6_dist;


			}//if(C_Array[C_count][27]==1)

			if(C_Array[y][27]==2)
			{

				double N1,C2,O2,N3,C4,O4,C5,C6;
				double C_ENeg;
				double O4_dist;
				double C5_dist;
				double C6_dist;
				double O4_charge;
				double C5_charge;
				double C6_charge;
				double C,N,O;

				O=3.44;
				N=3.04;
				C=2.55;

				N1=6.156;
				C2=7.506;
				O2=6.116;
				N3=5.839;
				C4=6.895;
				O4=5.143;
				C5=5.866;
				C6=6.383;



				O4_charge= -0.478;
				C6_charge= 0.096;
	


				O4_dist=sqrt((C_Array[y][18]-T_centroid[0])*(C_Array[y][18]-T_centroid[0])+
						(C_Array[y][19]-T_centroid[1])*(C_Array[y][19]-T_centroid[1])+
						(C_Array[y][20]-T_centroid[2])*(C_Array[y][20]-T_centroid[2]));

				C5_dist=sqrt((C_Array[y][12]-T_centroid[0])*(C_Array[y][12]-T_centroid[0])+
						(C_Array[y][13]-T_centroid[1])*(C_Array[y][13]-T_centroid[1])+
						(C_Array[y][14]-T_centroid[2])*(C_Array[y][14]-T_centroid[2]));

				C6_dist=sqrt((C_Array[y][15]-T_centroid[0])*(C_Array[y][15]-T_centroid[0])+
						(C_Array[y][16]-T_centroid[1])*(C_Array[y][16]-T_centroid[1])+
						(C_Array[y][17]-T_centroid[2])*(C_Array[y][17]-T_centroid[2]));


				EEnergy=0;
				EEnergy+= kconstant*avogadros*(O4_charge*electron)*(tyr_charge*electron)/(dielectric*(O4_dist*pow(10,-10)));
				EEnergy+= kconstant*avogadros*(C6_charge*electron)*(tyr_charge*electron)/(dielectric*(C6_dist*pow(10,-10)));

				EEnergy = EEnergy/1000;  // (J/mol to kJ/mol)
				EEnergy= EEnergy/4.18400; // kJ/mol to kcal/mol)

				Dist_Check=O4_dist;
				if(Dist_Check>C6_dist) Dist_Check=C6_dist;


			}//if(C_Array[C_count][27]==2)



if (Pass==1)
	if(Dist_Check<4.5)//8.5
			{
double NNParameter1, NNParameter2;
	if(Base_order[num]=='C')
	{
	if(Base_order[num-1]=='A')NNParameter1=AC;
	if(Base_order[num-1]=='T')NNParameter1=TC;
	if(Base_order[num-1]=='C')NNParameter1=CC;
	if(Base_order[num-1]=='G')NNParameter1=GC;
	if(Base_order[num-1]=='X')NNParameter1=-5;
	if(Base_order[num+1]=='A')NNParameter2=CA;
	if(Base_order[num+1]=='T')NNParameter2=CT;
	if(Base_order[num+1]=='C')NNParameter2=CC;
	if(Base_order[num+1]=='G')NNParameter2=CG;
	if(Base_order[num-1]=='X')NNParameter2=-5;
	}
	if(Base_order[num]=='A')
	{
			if(Base_order[num-1]=='A')NNParameter1=AA;
			if(Base_order[num-1]=='T')NNParameter1=TA;
			if(Base_order[num-1]=='C')NNParameter1=CA;
			if(Base_order[num-1]=='G')NNParameter1=GA;
			if(Base_order[num-1]=='X')NNParameter1=-5;
			if(Base_order[num+1]=='A')NNParameter2=AA;
			if(Base_order[num+1]=='T')NNParameter2=AT;
			if(Base_order[num+1]=='C')NNParameter2=AC;
			if(Base_order[num+1]=='G')NNParameter2=AG;
			if(Base_order[num-1]=='X')NNParameter2=-5;
	}
	if(Base_order[num]=='T')
	{
			if(Base_order[num-1]=='A')NNParameter1=AT;
			if(Base_order[num-1]=='T')NNParameter1=TT;
			if(Base_order[num-1]=='C')NNParameter1=CT;
			if(Base_order[num-1]=='G')NNParameter1=GT;
			if(Base_order[num-1]=='X')NNParameter1=-5;
			if(Base_order[num+1]=='A')NNParameter2=TA;
			if(Base_order[num+1]=='T')NNParameter2=TT;
			if(Base_order[num+1]=='C')NNParameter2=TC;
			if(Base_order[num+1]=='G')NNParameter2=TG;
			if(Base_order[num-1]=='X')NNParameter2=-5;
	}
	if(Base_order[num]=='G')
	{
			if(Base_order[num-1]=='A')NNParameter1=AG;
			if(Base_order[num-1]=='T')NNParameter1=TG;
			if(Base_order[num-1]=='C')NNParameter1=CG;
			if(Base_order[num-1]=='G')NNParameter1=GG;
			if(Base_order[num-1]=='X')NNParameter1=-5;
			if(Base_order[num+1]=='A')NNParameter2=GA;
			if(Base_order[num+1]=='T')NNParameter2=GT;
			if(Base_order[num+1]=='C')NNParameter2=GC;
			if(Base_order[num+1]=='G')NNParameter2=GG;
			if(Base_order[num-1]=='X')NNParameter2=-5;
	}





EEnergy = EEnergy * ( fabs(EEnergy) / ( fabs(EEnergy) + fabs(NNParameter1)+fabs(NNParameter2) ) );



		Total_Energy+=EEnergy;



			}//if Pass=1

		}
		}
	}

///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////
	for(x=0;x<T_count;x++)
			for(y=0;y<A_count;y++)
			{num=A_Array[y][28];
			if((Base_order[num]!='X')&&(Base_order[num+1]!='X')&&(Base_order[num-1]!='X'))
			{

				for(z=0;z<6;z++)
				{
					closest=1000;

					for(w=0;w<9;w++)
					{
						TempD[w]=(Tyr_Array[x][z*3+0]-A_Array[y][w*3+0])*(Tyr_Array[x][z*3+0]-A_Array[y][w*3+0]) + (Tyr_Array[x][z*3+1]-A_Array[y][w*3+1])*(Tyr_Array[x][z*3+1]-A_Array[y][w*3+1]) + (Tyr_Array[x][z*3+2]-A_Array[y][w*3+2])*(Tyr_Array[x][z*3+2]-A_Array[y][w*3+2]);
						TempD[w]=sqrt(TempD[w]);

						if (TempD[w]<closest)
						{
								closest=TempD[w];
						}

					}


					D[z]=closest;
				}

		///////////////
				///////////////
						/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
						/*~~~~~~ Calculate Centroid of Aromatic Residue ~~~~~~*/
						/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


			if((Tyr_Array[x][22]==1)||(Tyr_Array[x][22]==2))
			{
						for(l=0;l<3;l++)
							for(m=0;m<3;m++)
							{
								T_central[l][m]=Tyr_Array[x][l*3+m]+(Tyr_Array[x][(3+l)*3+m]-Tyr_Array[x][l*3+m])/2;
							}

						for(l=0;l<3;l++)
						{
							T_centroid[l]=(T_central[0][l]+T_central[1][l]+T_central[2][l])/3;
						}

			}//if((Tyr_Array[T_count][22]==1)||(Tyr_Array[T_count][22]==2))

			if(Tyr_Array[x][22]==3)
			{

							for(m=0;m<3;m++)
							{
								T_central[0][m]=Tyr_Array[y][4*3+m]+(Tyr_Array[y][5*3+m]-Tyr_Array[y][4*3+m])/2;
							}
							for(m=0;m<3;m++)


						for(l=0;l<3;l++)
						{
							T_centroid[l]=Tyr_Array[y][6*3+l]+(T_central[0][l]-Tyr_Array[y][6*3+l])/2;
						}

			}//if(Tyr_Array[T_count][22]==3)


			if(Tyr_Array[x][22]==4)
			{


				for(l=0;l<3;l++)
					for(m=0;m<3;m++)
					{
					T_central[l][m]=Tyr_Array[x][l*3+m] +
							((Tyr_Array[x][(3+l)*3+m]+(Tyr_Array[x][(3+l)*3+m]-Tyr_Array[x][(2+l)*3+m])/2) -
									Tyr_Array[x][l*3+m])/2;
					}
				for(l=3;l<5;l++)
						for(m=0;m<3;m++)
						{
						T_central[l][m]=Tyr_Array[x][l*3+m] +
								((Tyr_Array[x][(-3+l)*3+m]+(Tyr_Array[x][(-3+l)*3+m]-Tyr_Array[x][(-2+l)*3+m])/2) -
										Tyr_Array[x][l*3+m])/2;
						}

				for(l=0;l<3;l++)
				{
					T_centroid[l]=(T_central[0][l]+T_central[1][l]+T_central[2][l]+T_central[3][l]+T_central[4][l])/5;
				}


			}//if(Tyr_Array[T_count][22]==4)

			////////////
			////////////


				             for(a=0;a<=24;a+=3)
							{
								if(a==0)
								{
									XMax=A_Array[y][a];
									XMin=A_Array[y][a];
									YMax=A_Array[y][a+1];
									YMin=A_Array[y][a+1];
									ZMax=A_Array[y][a+2];
									ZMin=A_Array[y][a+2];
								}
								else
								{
									if(A_Array[y][a]>XMax)XMax=A_Array[y][a];
									if(A_Array[y][a]<XMin)XMin=A_Array[y][a];
									if(A_Array[y][a+1]>YMax)YMax=A_Array[y][a+1];
									if(A_Array[y][a+1]<YMin)YMin=A_Array[y][a+1];
									if(A_Array[y][a+2]>ZMax)ZMax=A_Array[y][a+2];
									if(A_Array[y][a+2]<ZMin)ZMin=A_Array[y][a+2];
								}
							}





								double distance,closestAtom[3]={1000,1000,1000};
									int closestC[3]={100,100,100};

									for(m=0;m<=24;m+=3)
										{
											distance=(T_centroid[0]-A_Array[y][m])*(T_centroid[0]-A_Array[y][m])+
													(T_centroid[1]-A_Array[y][m+1])*(T_centroid[1]-A_Array[y][m+1])+
													(T_centroid[2]-A_Array[y][m+2])*(T_centroid[2]-A_Array[y][m+2]);
											distance=sqrt(distance);

											if(distance<closestAtom[0])
											{
												closestAtom[2]=closestAtom[1];
												closestAtom[1]=closestAtom[0];
												closestAtom[0]=distance;
												closestC[2]=closestC[1];
												closestC[1]=closestC[0];
												closestC[0]=m;
											}

											if((distance>closestAtom[0])&&(distance<closestAtom[1]))
											{
												closestAtom[2]=closestAtom[1];
												closestAtom[1]=distance;
												closestC[2]=closestC[1];

												closestC[1]=m;
											}
											if((distance>closestAtom[1])&&(distance<closestAtom[2]))
											{
												closestAtom[2]=distance;
												closestC[2]=m;
											}

										}



				for(a=0;a<=15;a+=3)
				{
					if(a==0)
					{
						XMax=Tyr_Array[x][a];
						XMin=Tyr_Array[x][a];
						YMax=Tyr_Array[x][a+1];
						YMin=Tyr_Array[x][a+1];
						ZMax=Tyr_Array[x][a+2];
						ZMin=Tyr_Array[x][a+2];
					}

					{
						if(Tyr_Array[x][a]>XMax)XMax=Tyr_Array[x][a];
						if(Tyr_Array[x][a]<XMin)XMin=Tyr_Array[x][a];
						if(Tyr_Array[x][a+1]>YMax)YMax=Tyr_Array[x][a+1];
						if(Tyr_Array[x][a+1]<YMin)YMin=Tyr_Array[x][a+1];
						if(Tyr_Array[x][a+2]>ZMax)ZMax=Tyr_Array[x][a+2];
						if(Tyr_Array[x][a+2]<ZMin)ZMin=Tyr_Array[x][a+2];

						if(a<15)
						{
						if(A_Array[y][a]>XMax)XMax=A_Array[y][a];
						if(A_Array[y][a]<XMin)XMin=A_Array[y][a];
						if(A_Array[y][a+1]>YMax)YMax=A_Array[y][a+1];
						if(A_Array[y][a+1]<YMin)YMin=A_Array[y][a+1];
						if(A_Array[y][a+2]>ZMax)ZMax=A_Array[y][a+2];
						if(A_Array[y][a+2]<ZMin)ZMin=A_Array[y][a+2];
						}
					}
					if(a==0)
					{


						TXMax=A_Array[y][closestC[0]];
						TXMin=A_Array[y][closestC[0]];
						TYMax=A_Array[y][closestC[0]+1];
						TYMin=A_Array[y][closestC[0]+1];
						TZMax=A_Array[y][closestC[0]+2];
						TZMin=A_Array[y][closestC[0]+2];


						for(l=1;l<3;l++)
						{
							if(A_Array[y][closestC[l]]>TXMax)TXMax=A_Array[y][closestC[l]];
							if(A_Array[y][closestC[l]]<TXMin)TXMin=A_Array[y][closestC[l]];
							if(A_Array[y][closestC[l]]>TYMax)TYMax=A_Array[y][closestC[l]+1];
							if(A_Array[y][closestC[l]]<TYMin)TYMin=A_Array[y][closestC[l]+1];
							if(A_Array[y][closestC[l]]>TZMax)TZMax=A_Array[y][closestC[l]+2];
							if(A_Array[y][closestC[l]]<TZMin)TZMin=A_Array[y][closestC[l]+2];
						}
					}

					{
						if(Tyr_Array[x][a]>TXMax)TXMax=Tyr_Array[x][a];
						if(Tyr_Array[x][a]<TXMin)TXMin=Tyr_Array[x][a];
						if(Tyr_Array[x][a+1]>TYMax)TYMax=Tyr_Array[x][a+1];
						if(Tyr_Array[x][a+1]<TYMin)TYMin=Tyr_Array[x][a+1];
						if(Tyr_Array[x][a+2]>TZMax)TZMax=Tyr_Array[x][a+2];
						if(Tyr_Array[x][a+2]<TZMin)TZMin=Tyr_Array[x][a+2];
					}
				}





				Pass=1;
				Fail=0;

				pdbfile=fopen(pdb,"r");
				while((!feof(pdbfile))&&(Pass==1))
				{
					fgets(Line,150,pdbfile);
					if(Line[0]=='A')
					{
						sscanf(Line,"%s %d %s %s %c %d %lf %lf %lf",atom,&ANum,AName,RName,&Chain,&res_num,&xcor,&ycor,&zcor);

						if(AName[0]!='H')
						{

							if((TXMax>xcor)&&(xcor>TXMin)&&(TYMax>ycor)&&(ycor>TYMin)&&(TZMax>zcor)&&(zcor>TZMin))
							{
								Pass=0;
								for(a=0;a<=15;a+=3)
								{
									if((Tyr_Array[x][a]==xcor)&&(Tyr_Array[x][a+1]==ycor)&&(Tyr_Array[x][a+2]==zcor))
									{
										Pass=1;
										Fail=0;
									}
									if((A_Array[y][a]==xcor)&&(A_Array[y][a+1]==ycor)&&(A_Array[y][a+2]==zcor))
									{
										Pass=1;
										Fail=0;
									}
								}
							}
						}

					}
				}
				fclose(pdbfile);




				double avogadros=6.02214129*pow(10,23);
				double kconstant=8.9875517873681764*pow(10,9);  //Coloumbs Constant
				double EEnergy;
				double dielectric=4; //dielectric constant in protein environment
				double electron = 1.602176565*pow(10,-19);
				double Dist_Check;
				double coulomb_meter=3.336*pow(10,-30);
				double tyr_charge;

				if(Tyr_Array[x][22]==1) tyr_charge = -0.372 * 26.9/27.1;
				if(Tyr_Array[x][22]==2) tyr_charge = -0.372 * 27.1/27.1;
				if(Tyr_Array[x][22]==3) tyr_charge = -0.372 * 32.6/27.1;
				if(Tyr_Array[x][22]==4) tyr_charge = -0.372 * 21.0/27.1;


				if(A_Array[y][27]==1)
				{

					double N1,C2,N3,C4,C5,C6,N6,N7,C8,N9;
					double C_ENeg;
					double N6_dist;
					double N7_dist;
					double C8_dist;
					double C5_dist;
					double N6_charge;
					double N7_charge;
					double C8_charge;
					double C5_charge; //addition
					double C,N,O;

					O=3.44;
					N=3.04;
					C=2.55;

					N1=5.760;
					C2=5.527;
					N3=5.776;
					C4=5.959;
					C5=5.469;
					C6=6.018;
					N6=3.688;
					N7=6.025;
					C8=5.921;
					N9=3.967;


					N6_charge = 0.157 + 0.157;
					C5_charge = -0.015;
					N7_charge = -0.21;
					C8_charge = 0.115;




					N6_dist=sqrt((A_Array[y][21]-T_centroid[0])*(A_Array[y][21]-T_centroid[0])+
							(A_Array[y][22]-T_centroid[1])*(A_Array[y][22]-T_centroid[1])+
							(A_Array[y][23]-T_centroid[2])*(A_Array[y][23]-T_centroid[2]));

					N7_dist=sqrt((A_Array[y][6]-T_centroid[0])*(A_Array[y][6]-T_centroid[0])+
							(A_Array[y][7]-T_centroid[1])*(A_Array[y][7]-T_centroid[1])+
							(A_Array[y][8]-T_centroid[2])*(A_Array[y][8]-T_centroid[2]));

					C8_dist=sqrt((A_Array[y][9]-T_centroid[0])*(A_Array[y][9]-T_centroid[0])+
							(A_Array[y][10]-T_centroid[1])*(A_Array[y][10]-T_centroid[1])+
							(A_Array[y][11]-T_centroid[2])*(A_Array[y][11]-T_centroid[2]));

					C5_dist=sqrt((A_Array[y][3]-T_centroid[0])*(A_Array[y][3]-T_centroid[0])+
							(A_Array[y][4]-T_centroid[1])*(A_Array[y][4]-T_centroid[1])+
							(A_Array[y][5]-T_centroid[2])*(A_Array[y][5]-T_centroid[2]));


					EEnergy=0;
					EEnergy+= kconstant*avogadros*(N6_charge*electron)*(tyr_charge*electron)/(dielectric*(N6_dist*pow(10,-10)));
					EEnergy+= kconstant*avogadros*(N7_charge*electron)*(tyr_charge*electron)/(dielectric*(N7_dist*pow(10,-10)));
					EEnergy+= kconstant*avogadros*(C8_charge*electron)*(tyr_charge*electron)/(dielectric*(C8_dist*pow(10,-10)));
					EEnergy+= kconstant*avogadros*(C5_charge*electron)*(tyr_charge*electron)/(dielectric*(C5_dist*pow(10,-10)));


					EEnergy = EEnergy/1000;  // (J/mol to kJ/mol)
					EEnergy= EEnergy/4.18400; // kJ/mol to kcal/mol)

					Dist_Check=N6_dist;
					if(Dist_Check>N7_dist) Dist_Check=N7_dist;
					if(Dist_Check>C8_dist) Dist_Check=C8_dist;
					if(Dist_Check>C5_dist) Dist_Check=C5_dist;

				}//if(A_Array[C_count][27]==1)

				if(A_Array[y][27]==2)
				{

					double N1,C2,N2,N3,C4,C5,C6,O6,N7,C8,N9;
					double C_ENeg;
					double C5_dist;
					double O6_dist;
					double N7_dist;
					double C8_dist;
					double C5_charge;
					double O6_charge;
					double N7_charge;
					double C8_charge;
					double C,N,O;

					O=3.44;
					N=3.04;
					C=2.55;

					N1=5.467;
					C2=6.757;
					N2=3.928;
					N3=6.037;
					C4=6.321;
					C5=5.667;
					C6=6.844;
					O6=5.472;
					N7=5.984;
					C8=5.775;
					N9=3.817;

					C5_charge = 0.007;
					O6_charge = -0.441;
					N7_charge = -0.215;
					C8_charge = 0.107;
	



					O6_dist=sqrt((A_Array[y][24]-T_centroid[0])*(A_Array[y][24]-T_centroid[0])+
							(A_Array[y][25]-T_centroid[1])*(A_Array[y][25]-T_centroid[1])+
							(A_Array[y][26]-T_centroid[2])*(A_Array[y][26]-T_centroid[2]));

					N7_dist=sqrt((A_Array[y][6]-T_centroid[0])*(A_Array[y][6]-T_centroid[0])+
							(A_Array[y][7]-T_centroid[1])*(A_Array[y][7]-T_centroid[1])+
							(A_Array[y][8]-T_centroid[2])*(A_Array[y][8]-T_centroid[2]));

					C8_dist=sqrt((A_Array[y][9]-T_centroid[0])*(A_Array[y][9]-T_centroid[0])+
							(A_Array[y][10]-T_centroid[1])*(A_Array[y][10]-T_centroid[1])+
							(A_Array[y][11]-T_centroid[2])*(A_Array[y][11]-T_centroid[2]));


                                        C5_dist=sqrt((A_Array[y][3]-T_centroid[0])*(A_Array[y][3]-T_centroid[0])+
                                                        (A_Array[y][4]-T_centroid[1])*(A_Array[y][4]-T_centroid[1])+
                                                        (A_Array[y][5]-T_centroid[2])*(A_Array[y][5]-T_centroid[2]));


					EEnergy=0;
					EEnergy+= kconstant*avogadros*(C5_charge*electron)*(tyr_charge*electron)/(dielectric*(C5_dist*pow(10,-10)));
					EEnergy+= kconstant*avogadros*(O6_charge*electron)*(tyr_charge*electron)/(dielectric*(O6_dist*pow(10,-10)));
					EEnergy+= kconstant*avogadros*(N7_charge*electron)*(tyr_charge*electron)/(dielectric*(N7_dist*pow(10,-10)));
					EEnergy+= kconstant*avogadros*(C8_charge*electron)*(tyr_charge*electron)/(dielectric*(C8_dist*pow(10,-10)));
					

					EEnergy = EEnergy/1000;  // (J/mol to kJ/mol)
					EEnergy= EEnergy/4.18400; // kJ/mol to kcal/mol)


					Dist_Check=O6_dist;
					if(Dist_Check>N7_dist) Dist_Check=N7_dist;
					if(Dist_Check>C8_dist) Dist_Check=C8_dist;
					if(Dist_Check>C5_dist) Dist_Check=C5_dist;
				}//if(C_Array[A_count][27]==2)



	if (Pass==1)
		if(Dist_Check<4.5)
				{

	double NNParameter1, NNParameter2;

	if(Base_order[num]=='C')
	{
	if(Base_order[num-1]=='A')NNParameter1=AC;
	if(Base_order[num-1]=='T')NNParameter1=TC;
	if(Base_order[num-1]=='C')NNParameter1=CC;
	if(Base_order[num-1]=='G')NNParameter1=GC;
	if(Base_order[num-1]=='X')NNParameter1=-5;
	if(Base_order[num+1]=='A')NNParameter2=CA;
	if(Base_order[num+1]=='T')NNParameter2=CT;
	if(Base_order[num+1]=='C')NNParameter2=CC;
	if(Base_order[num+1]=='G')NNParameter2=CG;
	if(Base_order[num-1]=='X')NNParameter2=-5;
	}
	if(Base_order[num]=='A')
	{
			if(Base_order[num-1]=='A')NNParameter1=AA;
			if(Base_order[num-1]=='T')NNParameter1=TA;
			if(Base_order[num-1]=='C')NNParameter1=CA;
			if(Base_order[num-1]=='G')NNParameter1=GA;
			if(Base_order[num-1]=='X')NNParameter1=-5;
			if(Base_order[num+1]=='A')NNParameter2=AA;
			if(Base_order[num+1]=='T')NNParameter2=AT;
			if(Base_order[num+1]=='C')NNParameter2=AC;
			if(Base_order[num+1]=='G')NNParameter2=AG;
			if(Base_order[num-1]=='X')NNParameter2=-5;
	}
	if(Base_order[num]=='T')
	{
			if(Base_order[num-1]=='A')NNParameter1=AT;
			if(Base_order[num-1]=='T')NNParameter1=TT;
			if(Base_order[num-1]=='C')NNParameter1=CT;
			if(Base_order[num-1]=='G')NNParameter1=GT;
			if(Base_order[num-1]=='X')NNParameter1=-5;
			if(Base_order[num+1]=='A')NNParameter2=TA;
			if(Base_order[num+1]=='T')NNParameter2=TT;
			if(Base_order[num+1]=='C')NNParameter2=TC;
			if(Base_order[num+1]=='G')NNParameter2=TG;
			if(Base_order[num-1]=='X')NNParameter2=-5;
	}
	if(Base_order[num]=='G')
	{
			if(Base_order[num-1]=='A')NNParameter1=AG;
			if(Base_order[num-1]=='T')NNParameter1=TG;
			if(Base_order[num-1]=='C')NNParameter1=CG;
			if(Base_order[num-1]=='G')NNParameter1=GG;
			if(Base_order[num-1]=='X')NNParameter1=-5;
			if(Base_order[num+1]=='A')NNParameter2=GA;
			if(Base_order[num+1]=='T')NNParameter2=GT;
			if(Base_order[num+1]=='C')NNParameter2=GC;
			if(Base_order[num+1]=='G')NNParameter2=GG;
			if(Base_order[num-1]=='X')NNParameter2=-5;
	}






EEnergy = EEnergy * ( fabs(EEnergy) / ( fabs(EEnergy) + fabs(NNParameter1)+fabs(NNParameter2) ) );


		Total_Energy+=EEnergy;

				}//if Pass=1


		}
		}



return Total_Energy;
}


///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////




