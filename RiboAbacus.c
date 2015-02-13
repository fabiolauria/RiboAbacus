#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>



// user definition part
#define PROTEIN_INPUT   1	//assumes only values 0 and 1: 0 if the input file does not contain data related to protein amounts, 1 otherwise
#define RAMP_HYOUTHESIS	1	//assumes only values 0 and 1: 0 if the ramp hypothesis is not considered, 1 otherwise
#define RAMP_LENGTH     50	//specifies the ramp length (measured in codons)
#define SLOWDOWN_RATE   70	//specifies the ribosome slowdown rate on the ramp (it is a percentage: for example, 70 means that the speed of the ribosomes is reduced by 70%)
#define WIDTH_BIN       150	//specifies the width of the bin of the mRNA lengths distribution





//function to count the number of lines of a file
int nline(FILE *fin)
	{
		char ch;
		int nl=0;
		while ((ch=fgetc(fin))!=EOF)
			{
 				if (ch=='\n')
 					{
 			       			nl++;
 			       		}
 		 	}
		return nl;
	}




//function to find the maximum of 6 values
float maximum(float a, float b, float c, float d, float e, float f)
	{
		if(a>=b && a>=c && a>=d && a>=e && a>=f) return a;
		if(b>=a && b>=c && b>=d && b>=e && b>=f) return b;
		if(c>=a && c>=b && c>=d && c>=e && c>=f) return c;
		if(d>=a && d>=b && d>=c && d>=e && d>=f) return d;
		if(e>=a && e>=b && e>=c && e>=d && e>=f) return e;
		if(f>=a && f>=b && f>=c && f>=d && f>=e) return f;
	}




//functions to solve the system: 
float S1(float kS, float k1)
	{
		return kS/k1;			
	}

float S2(float kS, float k2, float k9, float ctilde, float dvk)
	{
		return kS*(ctilde/(k2*k9*dvk));			
	}

float S3(float kS, float chat, float k9)
	{
		return kS*(chat/k9);			
	}

float S4(float kS, float k4)
	{
		return kS/k4;			
	}

float S5(float kS, float k5)
	{
		return kS/k5;			
	}

float S6(float kS, float k6)
	{
		return kS/k6;			
	}

float S7(float kS, float k7)
	{
		return kS/k7;			
	}

float S8(float kS, float U, float k8, float k9, float G)
	{
		return (kS/(k8*G))*((k8/(k9*U))+1);
	}

float S9(float kS, float U, float k9)
	{
		return kS/(k9*U);
	}

float S21(float kS, float k2, float k2r, float k9, float chat, float dv1)
	{
		return (kS+((k2r*chat*kS)/k9))/(k2*dv1);
	}

float sumS3S7(float kS, float k4, float k5, float k6, float k7, float k9, float chat)
	{
		return S3(kS, chat, k9)+S4(kS, k4)+S5(kS, k5)+S6(kS, k6)+S7(kS, k7);
	}

float sum_rib_downstream(float sum_rib_onecodon[], int n, int L)
	{
		int i;
		float grsum=0;
		for(i=n+1; i<=n+L; i++)
			{
				grsum=grsum+sum_rib_onecodon[i];
			}
		return grsum;
	}

float U(float sum_rib_onecodon[], float Mr, int n, int L)
	{	
		return (Mr-sum_rib_downstream(sum_rib_onecodon, n, L)-sum_rib_onecodon[n+L])/(Mr-sum_rib_downstream(sum_rib_onecodon, n, L));
	}

float sum_rib_one_codon(float kS, float U, float k1, float k2, float k4, float k5, float k6, float k7, float k8, float k9, float G, float ctilde, float chat, float dvk)
	{
		return S1(kS, k1)+S2(kS, k2, k9, ctilde, dvk)+sumS3S7(kS, k4, k5, k6, k7, k9, chat)+S8(kS, U, k8, k9, G)+S9(kS, U, k9);
	}




//function to compute the number of ribosomes on the transcript 
float sum_rib_transcript(float sum_rib_onecodon[], float length)
	{
		int i;
		float totsum=0;
		for(i=0; i<length; i++)
			{
				totsum=totsum+sum_rib_onecodon[i];
			}
		return totsum;
	}




//function to compute the number of ribosomes on the ramp
float sum_rib_ramp(float sum_rib_onecodon[], int len_trip, float ramp)
	{
		int i;
		float sum_rib_ramp=0;
		for(i=0; i<ramp; i++)
			{
				sum_rib_ramp=sum_rib_ramp+sum_rib_onecodon[i];
			}
		return sum_rib_ramp;
	}




//principal function: to obtain the number of ribosomes per trancript and other information
void output(FILE *fin, FILE *fout, int nline, float stax, int ramp, int ramp_presence, int protein)
{
int i, k, j, z, wl, read;
int n_tr=nline, l_max=100000, l_max_tr=l_max/3, len, len_trip;
int b, taxmin=1, taxmax=8192, min_interval=8, m[min_interval+1], B[min_interval+1], tmn, tmx, tmd, mtmd, h;
char divide_v[l_max_tr][4], tr[l_max], info[30];
float division_values[l_max_tr], sum_rib_onecodon[l_max_tr], Uv[l_max_tr], chat, ctilde, control, k1=20, k2=100, k2r=79, k3=207, k3r=3.45, k4=100, k5=638, k6=15, k7=20, k8=150, k8r=140, k9=250, Mr, G=30, L=10, cu_vec[64]={17.6, 15.2, 12.2, 10.6, 20.3, 17.7, 15.3, 12.6, 7.7, 12.2, 1.0, 1.6, 12.9, 4.4, 0.8, 13.2, 13.2, 17.5, 10.9, 4.5, 19.6, 19.8, 15.1, 10.4, 7.2, 16.9, 12.3, 6.2, 39.6, 6.9, 34.2, 11.4, 16.0, 13.1, 17.0, 12.1, 20.8, 18.9, 19.1, 19.5, 7.5, 15.1, 24.4, 12.2, 22.0, 6.1, 31.9, 12.0, 11.0, 18.4, 21.8, 10.8, 14.5, 27.7, 25.1, 22.2, 7.1, 15.8, 29.0, 16.5, 28.1, 7.4, 39.6, 16.5}, cai_vec[l_max_tr], CAIr, CAIt ,n_rib=0, GCr, GCt, n_rib_ramp, sdr=((100-stax)/100);

chat=((k3r+k4)*k9)/(k3*k4);
ctilde=((k2r/k3)+((k2r*k3r)/(k3*k4))+1)*k9;

if(protein==0)
	{
		read=5;
	}
else
	{
		read=6;
	}
for(i=0; i<(n_tr/2); i++)
{	

	//scanning of the input file to acquire the trnascripts copy number Mr and the mRNA sequences tr[] 
	for(z=0; z<read; z++)
		{
		if(z%read!=1 && z%read!=read-1){fscanf(fin, "%s", &info); fprintf(fout, "%s\t", info);}
		else
		    {
			if(z%read==1){fscanf(fin, "%f", &Mr); fprintf(fout, "%f\t", Mr); Mr=Mr;}
			if(z%read==read-1)
			{
			for(k=0; k<l_max; k++)
				{
					tr[k]=' ';
				}
			fscanf(fin, "%s", &tr);
			len=strlen(tr);
			
			
			
			
			//detection of problematic transcripts (unavailable, wrong or too long sequences, unusual start codons)
			wl=0;
			for(k=0; k<len; k++)
				{
					if(tr[k]!='A' && tr[k]!='T' && tr[k]!='C' && tr[k]!='G')
						{
							wl=wl+1;
						}
				}
			if(len%3!=0 || tr[0]!='A' || tr[1]!='T' || tr[2]!='G' || wl!=0 || len>30000)
				{
					fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "na", "na", "na", "na", "na", "na", "na", "na", "na", "na");
				}
			
			
			
			
			else
				{
					
					//initialization of CAI (codon adaptation index) values and percentage of GC content (of the ramp (r) and the entire transcript (t))
					CAIr=0;
					CAIt=0;
					GCr=0;
					GCt=0;
					
					
					
					
					//division of the transcript in triplets and association of the correspondent codon usage bias value and of a number suitable for the computation of the CAI
					for(k=0; k<len; k=k+3)
						{
							divide_v[k/3][0]=tr[k];
							divide_v[k/3][1]=tr[k+1];
							divide_v[k/3][2]=tr[k+2];
						}
					len_trip=len/3;
					for(k=0; k<len_trip; k++)
						{
							if(divide_v[k][0]=='G')
								{	
									if(divide_v[k][1]=='G')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[63]; cai_vec[k]=(division_values[k]/maximum(cu_vec[63],cu_vec[59],cu_vec[55],cu_vec[51],0,0));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[59]; cai_vec[k]=(division_values[k]/maximum(cu_vec[63],cu_vec[59],cu_vec[55],cu_vec[51],0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[55]; cai_vec[k]=(division_values[k]/maximum(cu_vec[63],cu_vec[59],cu_vec[55],cu_vec[51],0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[51]; cai_vec[k]=(division_values[k]/maximum(cu_vec[63],cu_vec[59],cu_vec[55],cu_vec[51],0,0));}
										}else{
									if(divide_v[k][1]=='A')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[62]; cai_vec[k]=(division_values[k]/maximum(cu_vec[62],cu_vec[58],0,0,0,0));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[58]; cai_vec[k]=(division_values[k]/maximum(cu_vec[62],cu_vec[58],0,0,0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[54]; cai_vec[k]=(division_values[k]/maximum(cu_vec[54],cu_vec[50],0,0,0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[50]; cai_vec[k]=(division_values[k]/maximum(cu_vec[54],cu_vec[50],0,0,0,0));}
										}else{
									if(divide_v[k][1]=='C')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[61]; cai_vec[k]=(division_values[k]/maximum(cu_vec[61],cu_vec[57],cu_vec[53],cu_vec[49],0,0));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[57]; cai_vec[k]=(division_values[k]/maximum(cu_vec[61],cu_vec[57],cu_vec[53],cu_vec[49],0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[53]; cai_vec[k]=(division_values[k]/maximum(cu_vec[61],cu_vec[57],cu_vec[53],cu_vec[49],0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[49]; cai_vec[k]=(division_values[k]/maximum(cu_vec[61],cu_vec[57],cu_vec[53],cu_vec[49],0,0));}
										}else{
									if(divide_v[k][1]=='T')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[60]; cai_vec[k]=(division_values[k]/maximum(cu_vec[60],cu_vec[56],cu_vec[52],cu_vec[48],0,0));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[56]; cai_vec[k]=(division_values[k]/maximum(cu_vec[60],cu_vec[56],cu_vec[52],cu_vec[48],0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[52]; cai_vec[k]=(division_values[k]/maximum(cu_vec[60],cu_vec[56],cu_vec[52],cu_vec[48],0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[48]; cai_vec[k]=(division_values[k]/maximum(cu_vec[60],cu_vec[56],cu_vec[52],cu_vec[48],0,0));}
										}
								}}}}else{
							if(divide_v[k][0]=='A')
								{	
									if(divide_v[k][1]=='G')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[47]; cai_vec[k]=(division_values[k]/maximum(cu_vec[47],cu_vec[43],cu_vec[31],cu_vec[27],cu_vec[23],cu_vec[19]));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[43]; cai_vec[k]=(division_values[k]/maximum(cu_vec[47],cu_vec[43],cu_vec[31],cu_vec[27],cu_vec[23],cu_vec[19]));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[39]; cai_vec[k]=(division_values[k]/maximum(cu_vec[39],cu_vec[35],cu_vec[13],cu_vec[9],cu_vec[5],cu_vec[1]));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[35]; cai_vec[k]=(division_values[k]/maximum(cu_vec[39],cu_vec[35],cu_vec[13],cu_vec[9],cu_vec[5],cu_vec[1]));}
										}else{
									if(divide_v[k][1]=='A')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[46]; cai_vec[k]=(division_values[k]/maximum(cu_vec[46],cu_vec[42],0,0,0,0));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[42]; cai_vec[k]=(division_values[k]/maximum(cu_vec[46],cu_vec[42],0,0,0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[38]; cai_vec[k]=(division_values[k]/maximum(cu_vec[38],cu_vec[34],0,0,0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[34]; cai_vec[k]=(division_values[k]/maximum(cu_vec[38],cu_vec[34],0,0,0,0));}
										}else{
									if(divide_v[k][1]=='C')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[45]; cai_vec[k]=(division_values[k]/maximum(cu_vec[45],cu_vec[41],cu_vec[37],cu_vec[33],0,0));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[41]; cai_vec[k]=(division_values[k]/maximum(cu_vec[45],cu_vec[41],cu_vec[37],cu_vec[33],0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[37]; cai_vec[k]=(division_values[k]/maximum(cu_vec[45],cu_vec[41],cu_vec[37],cu_vec[33],0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[33]; cai_vec[k]=(division_values[k]/maximum(cu_vec[45],cu_vec[41],cu_vec[37],cu_vec[33],0,0));}
										}else{
									if(divide_v[k][1]=='T')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[44]; cai_vec[k]=1;}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[40]; cai_vec[k]=(division_values[k]/maximum(cu_vec[40],cu_vec[36],cu_vec[32],0,0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[36]; cai_vec[k]=(division_values[k]/maximum(cu_vec[40],cu_vec[36],cu_vec[32],0,0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[32]; cai_vec[k]=(division_values[k]/maximum(cu_vec[40],cu_vec[36],cu_vec[32],0,0,0));}
										}
								}}}}else{
							if(divide_v[k][0]=='C')
								{	
									if(divide_v[k][1]=='G')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[31]; cai_vec[k]=(division_values[k]/maximum(cu_vec[47],cu_vec[43],cu_vec[31],cu_vec[27],cu_vec[23],cu_vec[19]));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[27]; cai_vec[k]=(division_values[k]/maximum(cu_vec[47],cu_vec[43],cu_vec[31],cu_vec[27],cu_vec[23],cu_vec[19]));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[23]; cai_vec[k]=(division_values[k]/maximum(cu_vec[47],cu_vec[43],cu_vec[31],cu_vec[27],cu_vec[23],cu_vec[19]));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[19]; cai_vec[k]=(division_values[k]/maximum(cu_vec[47],cu_vec[43],cu_vec[31],cu_vec[27],cu_vec[23],cu_vec[19]));}
										}else{
									if(divide_v[k][1]=='A')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[30]; cai_vec[k]=(division_values[k]/maximum(cu_vec[30],cu_vec[26],0,0,0,0));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[26]; cai_vec[k]=(division_values[k]/maximum(cu_vec[30],cu_vec[26],0,0,0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[22]; cai_vec[k]=(division_values[k]/maximum(cu_vec[22],cu_vec[18],0,0,0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[18]; cai_vec[k]=(division_values[k]/maximum(cu_vec[22],cu_vec[18],0,0,0,0));}
										}else{
									if(divide_v[k][1]=='C')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[29]; cai_vec[k]=(division_values[k]/maximum(cu_vec[29],cu_vec[25],cu_vec[21],cu_vec[17],0,0));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[25]; cai_vec[k]=(division_values[k]/maximum(cu_vec[29],cu_vec[25],cu_vec[21],cu_vec[17],0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[21]; cai_vec[k]=(division_values[k]/maximum(cu_vec[29],cu_vec[25],cu_vec[21],cu_vec[17],0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[17]; cai_vec[k]=(division_values[k]/maximum(cu_vec[29],cu_vec[25],cu_vec[21],cu_vec[17],0,0));}
										}else{
									if(divide_v[k][1]=='T')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[28]; cai_vec[k]=(division_values[k]/maximum(cu_vec[28],cu_vec[24],cu_vec[20],cu_vec[16],cu_vec[12],cu_vec[8]));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[24]; cai_vec[k]=(division_values[k]/maximum(cu_vec[28],cu_vec[24],cu_vec[20],cu_vec[16],cu_vec[12],cu_vec[8]));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[20]; cai_vec[k]=(division_values[k]/maximum(cu_vec[28],cu_vec[24],cu_vec[20],cu_vec[16],cu_vec[12],cu_vec[8]));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[16]; cai_vec[k]=(division_values[k]/maximum(cu_vec[28],cu_vec[24],cu_vec[20],cu_vec[16],cu_vec[12],cu_vec[8]));}
										}
								}}}}else{
							if(divide_v[k][0]=='T')
								{	
									if(divide_v[k][1]=='G')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[15]; cai_vec[k]=1;}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[11]; cai_vec[k]=(division_values[k]/maximum(cu_vec[11],cu_vec[14],cu_vec[10],0,0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[7]; cai_vec[k]=(division_values[k]/maximum(cu_vec[7],cu_vec[3],0,0,0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[3]; cai_vec[k]=(division_values[k]/maximum(cu_vec[7],cu_vec[3],0,0,0,0));}
										}else{
									if(divide_v[k][1]=='A')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[14]; cai_vec[k]=(division_values[k]/maximum(cu_vec[11],cu_vec[14],cu_vec[10],0,0,0));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[10]; cai_vec[k]=(division_values[k]/maximum(cu_vec[11],cu_vec[14],cu_vec[10],0,0,0));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[6]; cai_vec[k]=(division_values[k]/maximum(cu_vec[6],cu_vec[2],0,0,0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[2]; cai_vec[k]=(division_values[k]/maximum(cu_vec[6],cu_vec[2],0,0,0,0));}
										}else{
									if(divide_v[k][1]=='C')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[13]; cai_vec[k]=(division_values[k]/maximum(cu_vec[39],cu_vec[35],cu_vec[13],cu_vec[9],cu_vec[5],cu_vec[1]));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[9]; cai_vec[k]=(division_values[k]/maximum(cu_vec[39],cu_vec[35],cu_vec[13],cu_vec[9],cu_vec[5],cu_vec[1]));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[5]; cai_vec[k]=(division_values[k]/maximum(cu_vec[39],cu_vec[35],cu_vec[13],cu_vec[9],cu_vec[5],cu_vec[1]));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[1]; cai_vec[k]=(division_values[k]/maximum(cu_vec[39],cu_vec[35],cu_vec[13],cu_vec[9],cu_vec[5],cu_vec[1]));}
										}else{
									if(divide_v[k][1]=='T')
										{
											if(divide_v[k][2]=='G'){division_values[k]=cu_vec[12]; cai_vec[k]=(division_values[k]/maximum(cu_vec[28],cu_vec[24],cu_vec[20],cu_vec[16],cu_vec[12],cu_vec[8]));}else
											if(divide_v[k][2]=='A'){division_values[k]=cu_vec[8]; cai_vec[k]=(division_values[k]/maximum(cu_vec[28],cu_vec[24],cu_vec[20],cu_vec[16],cu_vec[12],cu_vec[8]));}else
											if(divide_v[k][2]=='C'){division_values[k]=cu_vec[4]; cai_vec[k]=(division_values[k]/maximum(cu_vec[4],cu_vec[0],0,0,0,0));}else
											if(divide_v[k][2]=='T'){division_values[k]=cu_vec[0]; cai_vec[k]=(division_values[k]/maximum(cu_vec[4],cu_vec[0],0,0,0,0));}
										}
								}}}}}}}
								
								
								division_values[k]=division_values[k]/1000;
									
						}
						
							
						
						//computation of the transcript CAI 
						for(k=0; k<len_trip; k++)
							{
								CAIt=CAIt+log(cai_vec[k]);
							}
						CAIt=CAIt/len_trip;
						
						
						
						
						//computation of the ramp CAI 
						if(ramp_presence==1)
							{
								for(k=0; k<ramp; k++)
									{
										CAIr=CAIr+log(cai_vec[k]);
									}
								CAIr=CAIr/ramp;
							}
											
						
						
						
						//computation of the transcript percentage of GC content
						for(k=0; k<len; k++)
							{
								if(tr[k]=='C' || tr[k]=='G')
									{
										GCt=GCt+1;
									}
							}
						GCt=(GCt*100)/len;
						
						
						

						//computation of the ramp percentage of GC content 
						if(ramp_presence==1)
							{
								for(k=0; k<(ramp*3); k++)
									{
										if(tr[k]=='C' || tr[k]=='G')
											{
												GCr=GCr+1;
											}
									}
								GCr=(GCr*100)/(ramp*3);
							}								
							
							
							
							
					//initialization of the number of ribosomes at each codon and definition of the range of values for the ribosome exit fluxes used in the following procedure				
					for(k=0; k<len_trip+L-1; k++)
						{
							sum_rib_onecodon[k]=0;
						}
					tmx=taxmax;
					tmn=taxmin-1;
					tmd=(tmx-tmn)/2;
								
								
							
								
					//computation of both the number of ribosomes and the translocation probability for each codon for a set of exit flux values (the one in the middle of the range). At each step the exit flux range is halved (unitll it reaches a minimum number of elements) depending on the value that satisfied the imposed conditions. The computations are differentiated according to the presence of the ramp: if the ramp hypothesis is taken into account some kinetic constants are multiplied by the fixed slowdown rate (sdr) for the first (chosen ramp length) position in order to decrease the speed of all the reactions.		
					while(tmx-tmn>min_interval)
						{
							mtmd=0;	
							
							
							if(ramp_presence==0)				
							{							
								for(j=len_trip-2; j>0; j--)
									{	
										Uv[j]=U(sum_rib_onecodon, Mr, j, L);
										sum_rib_onecodon[j]=sum_rib_one_codon(tmd, Uv[j], k1, k2, k4, k5, k6, k7, k8, k9, G, ctilde, chat, division_values[j]);
									}
								Uv[0]=U(sum_rib_onecodon, Mr, 0, L);
								sum_rib_onecodon[0]=S21(tmd, k2, k2r, k9, chat, division_values[0])+sumS3S7(tmd, k4, k5, k6, k7, k9, chat)+S8(tmd, Uv[0], k8, k9, G)+S9(tmd, Uv[0], k9);
							}
							
							
							else
							{
								for(j=len_trip-2; j>ramp-1; j--)
									{	
										Uv[j]=U(sum_rib_onecodon, Mr, j, L);
										sum_rib_onecodon[j]=sum_rib_one_codon(tmd, Uv[j], k1, k2, k4, k5, k6, k7, k8, k9, G, ctilde, chat, division_values[j]);
									}
								for(j=ramp-1; j>0; j--)
									{	
										Uv[j]=U(sum_rib_onecodon, Mr, j, L);
										sum_rib_onecodon[j]=sum_rib_one_codon(tmd, Uv[j], k1*sdr, k2, k4*sdr, k5*sdr, k6*sdr, k7*sdr, k8*sdr, k9*sdr, G, ctilde, chat, division_values[j]);
									}
								Uv[0]=U(sum_rib_onecodon, Mr, 0, L);
								sum_rib_onecodon[0]=S21(tmd, k2*sdr, k2r, k9, chat, division_values[0])+sumS3S7(tmd, k4*sdr, k5*sdr, k6*sdr, k7*sdr, k9*sdr, chat)+S8(tmd, Uv[0], k8*sdr, k9, G)+S9(tmd, Uv[0], k9*sdr);	
							}
										
																				
							control=(Mr*len_trip)/L;
							if(sum_rib_transcript(sum_rib_onecodon, len_trip-1)<=control)
								{
									for(h=0; h<=len_trip-2; h++)
										{
											if(Uv[h]>=0 && Uv[h]<=1 && sum_rib_onecodon[h]>=0 && sum_rib_onecodon[h]<=Mr)
												{
													mtmd=mtmd+1;
												}
											else
												{
													mtmd=mtmd+0;		
												}
										}
								}
							else
								{
									mtmd=mtmd+0;
								}
							if(mtmd==len_trip-1)
								{
									tmn=tmd;
									tmx=tmx;
									tmd=tmn+(tmx-tmn)/2;
									
								}
							else
								{
									tmn=tmn;
									tmx=tmd;
									tmd=(tmx-tmn)/2;
								}
						}
						
						
						
						
					//computation of both the number of ribosomes and the translocation probability for each codon for each value of the left range of exit fluxes. The computations are differentiated according to the presence of the ramp: if the ramp hypothesis is taken into account some kinetic constants are multiplied by the fixed slowdown rate (sdr) for the first (chosen ramp length) position in order to decrease the speed of all the reactions.	
					for(k=0; k<min_interval+1; k++)
						{
							m[k]=0;
						}
					for(k=0; k<min_interval+1; k++)
						{
							B[k]=0;
						}
					b=0;
					if(tmn==0)
						{
							tmn=tmn+1;
							tmx=tmx+1;
						}							
					for(k=tmn; k<=tmx; k++)
						{
							
							if(ramp_presence==0)				
							{							
								for(j=len_trip-2; j>0; j--)
									{	
										Uv[j]=U(sum_rib_onecodon, Mr, j, L);
										sum_rib_onecodon[j]=sum_rib_one_codon(k, Uv[j], k1, k2, k4, k5, k6, k7, k8, k9, G, ctilde, chat, division_values[j]);
									}
								Uv[0]=U(sum_rib_onecodon, Mr, 0, L);
								sum_rib_onecodon[0]=S21(k, k2, k2r, k9, chat, division_values[0])+sumS3S7(k, k4, k5, k6, k7, k9, chat)+S8(k, Uv[0], k8, k9, G)+S9(k, Uv[0], k9);
							}
							
							else
							{
								for(j=len_trip-2; j>ramp-1; j--)
									{	
										Uv[j]=U(sum_rib_onecodon, Mr, j, L);
										sum_rib_onecodon[j]=sum_rib_one_codon(k, Uv[j], k1, k2, k4, k5, k6, k7, k8, k9, G, ctilde, chat, division_values[j]);
									}
								for(j=ramp-1; j>0; j--)
									{	
										Uv[j]=U(sum_rib_onecodon, Mr, j, L);
										sum_rib_onecodon[j]=sum_rib_one_codon(k, Uv[j], k1*sdr, k2, k4*sdr, k5*sdr, k6*sdr, k7*sdr, k8*sdr, k9*sdr, G, ctilde, chat, division_values[j]);
									}
								Uv[0]=U(sum_rib_onecodon, Mr, 0, L);
								sum_rib_onecodon[0]=S21(k, k2*sdr, k2r, k9, chat, division_values[0])+sumS3S7(k, k4*sdr, k5*sdr, k6*sdr, k7*sdr, k9*sdr, chat)+S8(k, Uv[0], k8*sdr, k9, G)+S9(k, Uv[0], k9*sdr);
							}

													
							control=(Mr*len_trip)/L;
							if(sum_rib_transcript(sum_rib_onecodon, len_trip-1)<=control)
								{
									for(h=0; h<=len_trip-2; h++)
										{
											if(Uv[h]>=0 && Uv[h]<=1 && sum_rib_onecodon[h]>=0 && sum_rib_onecodon[h]<=Mr)
												{
													m[b]=m[b]+1;
												}
											else
												{
													m[b]=m[b]+0;		
												}
										}
								}
							else
								{
									m[b]=m[b]+0;
								}
							if(m[b]==len_trip-1)
								{
									B[b]=k;
								}
							else
								{
									B[b]=0;
								}
							b=b+1;
						}
						
						
						
						
					//selection of the maximum ribosome exit flux that satisfies all the conditions and computation of both the number of ribosomes and the translocation probability for each codon using that value. The computations are differentiated according to the presence of the ramp: if the ramp hypothesis is taken into account some kinetic constants are multiplied by the fixed slowdown rate (sdr) for the first (chosen ramp length) position in order to decrease the speed of all the reactions.		
					h=0;
					for(k=0; k<=min_interval; k++)
						{
							if(B[k]==0)
								{
									if(k==0)
										{
											h=1;
										}
									else
										{
											h=B[k-1];
											

											if(ramp_presence==0)				
											{											
												for(j=len_trip-2; j>0; j--)
													{	
														Uv[j]=U(sum_rib_onecodon, Mr, j, L);
														sum_rib_onecodon[j]=sum_rib_one_codon(h, Uv[j], k1, k2, k4, k5, k6, k7, k8, k9, G, ctilde, chat, division_values[j]);
													}
												Uv[0]=U(sum_rib_onecodon, Mr, 0, L);
												sum_rib_onecodon[0]=S21(h, k2, k2r, k9, chat, division_values[0])+sumS3S7(k, k4, k5, k6, k7, k9, chat)+S8(k, Uv[0], k8, k9, G)+S9(k, Uv[0], k9);
											}
											else
											{
												for(j=len_trip-2; j>ramp-1; j--)
													{	
														Uv[j]=U(sum_rib_onecodon, Mr, j, L);
														sum_rib_onecodon[j]=sum_rib_one_codon(h, Uv[j], k1, k2, k4, k5, k6, k7, k8, k9, G, ctilde, chat, division_values[j]);
													}
												for(j=ramp-1; j>0; j--)
													{	
														Uv[j]=U(sum_rib_onecodon, Mr, j, L);
														sum_rib_onecodon[j]=sum_rib_one_codon(h, Uv[j], k1*sdr, k2, k4*sdr, k5*sdr, k6*sdr, k7*sdr, k8*sdr, k9*sdr, G, ctilde, chat, division_values[j]);
													}
												Uv[0]=U(sum_rib_onecodon, Mr, 0, L);
												sum_rib_onecodon[0]=S21(h, k2*sdr, k2r, k9, chat, division_values[0])+sumS3S7(k, k4*sdr, k5*sdr, k6*sdr, k7*sdr, k9*sdr, chat)+S8(k, Uv[0], k8*sdr, k9, G)+S9(k, Uv[0], k9*sdr);
											}
									
									
											n_rib=(sum_rib_transcript(sum_rib_onecodon, len_trip-1))/Mr;
											n_rib_ramp=(sum_rib_ramp(sum_rib_onecodon, len_trip, ramp))/Mr;
											k=min_interval+1;
										}
								}
						}
						
						
						
						
					//writing on the output file
					if(h>1)
						{
							fprintf(fout, "%i\t%f\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%i\n", (int)len, n_rib, (int)n_rib_ramp, GCt, GCr, exp(CAIt), exp(CAIr), (n_rib*3000)/(len), (((n_rib*3000)/(len))/100)*Mr, (int)h);
						}
					else
						{
							fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "na", "na", "na", "na", "na", "na", "na", "na", "na", "na");
						}
		 		}
			}
		    }
		}		
}
}




//function to eliminate results associated with unavailable or wrong sequences
void correctline (FILE *fin, FILE *fout, int nline, int protein)
{
	int i, j, n_tr=nline, ninfo;
	
	if(protein==0) {ninfo=14;}
	else {ninfo=15;}
	
	char info[30], info_v[ninfo][30];
	
	

	if(ninfo==14)
	{
		fprintf(fout, "geneID\ttranscript_level\tproteinID\ttranscriptID\tlength\tnumber_rib\tnumber_rib_ramp\tGC%_tot\tGC%_ramp\tCAI_tot\tCAI_ramp\tcoverage\tTE\texit_flux\n");
	}
	else
	{
		fprintf(fout, "geneID\ttranscript_level\tprotein_amount\tproteinID\ttranscriptID\tlength\tnumber_rib\tnumber_rib_ramp\tGC%_tot\tGC%_ramp\tCAI_tot\tCAI_ramp\tcoverage\tTE\texit_flux\n");
	}
	
	
	for(i=0; i<n_tr; i++)
	{
		for(j=0; j<ninfo; j++)
		{
			fscanf(fin, "%s", &info);
			strcpy(info_v[j], info);
		}
		if(strcmp(info_v[ninfo-1], "na")!=0)
		{
			for(j=0; j<ninfo; j++)
			{
				if(j==ninfo-9)
				{
					fprintf(fout, "%i\t",(int) (atof(info_v[j])+0.5));
				}
				else
				{
					fprintf(fout, "%s\t", info_v[j]);
				}
			}
			fprintf(fout, "\n");
		}
	}	
}




//function to generate the data for number of ribosome per transcript distributions
void graph(FILE *fin, FILE *out, int nline, int protein)
{
	int i, j, ntrans=nline, ncount=51, ninfo;
	char si[20];
	float intens, intens_v[ntrans], nrib, nrib_v[ntrans], count[ncount], freq, tot=0;
	
	if(protein==0) {ninfo=14;}
	else {ninfo=15;}
	
	if(fin==NULL)
	{
		printf("\nError\n\n");
	}
	else
	{	
		for(i=0; i<ninfo; i++)
		{
			fscanf(fin, "%s", &si);
		}
		for(i=1; i<ntrans; i++)		
		{
			for(j=0; j<ninfo; j++)
				{
					if(j!=1 && j!=ninfo-9)
					{
						fscanf(fin, "%s", &si);
					}
					else
					{
						if(j==1){fscanf(fin, "%f", &intens);}
						if(j==ninfo-9){fscanf(fin, "%f", &nrib);}
					}
				}
				
			intens_v[i-1]=intens/100;
			nrib_v[i-1]=nrib;
		}
		
		for(i=0; i<ncount; i++)
			{
				count[i]=0;
			}
		fprintf(out, "rib/tr\tevents\tfreq_int\n");		
		for(i=0; i<ncount; i++)
			{
				for(j=0; j<ntrans; j++)
					{
						if(nrib_v[j]==i)
							{
								count[i]=count[i]+intens_v[j];
							}
					}
				tot=tot+count[i];
			}
		freq=0;
		for(i=0; i<ncount; i++)
			{	
				fprintf(out, "%i\t%i\t%f\n", i, (int)count[i], count[i]/tot);
				freq=freq+(count[i]/tot);
			}
	}
}




//function to generate the data for length distribution with specific width of the bins
void graph_length(FILE *fin, FILE *out, int nline, int protein, int fixbin)
{
	printf("%i\n", nline);
	int i, j, ntrans=nline, bin=fixbin, ncount=(15000/bin)+1, ninfo;
	char si[20];
	float intens, intens_v[ntrans], freq, tot=0, len, len_v[ntrans], count[ncount];
	
	if(protein==0) {ninfo=14;}
	else {ninfo=15;}
	
	if(fin==NULL)
	{
		printf("\nError\n\n");
	}
	else
	{	
		for(i=0; i<ninfo; i++)
		{
			fscanf(fin, "%s", &si);
		}
		for(i=1; i<ntrans; i++)		
		{
			for(j=0; j<ninfo; j++)
				{
					if(j!=1 && j!=ninfo-10)
					{
						fscanf(fin, "%s", &si);
					}
					else
					{
						if(j==1){fscanf(fin, "%f", &intens);}
						if(j==ninfo-10){fscanf(fin, "%f", &len);}
					}
				}
				
			intens_v[i-1]=intens/100;
			len_v[i-1]=len;
		}
		
		for(i=0; i<ncount; i++)
			{
				count[i]=0;
			}
			
		fprintf(out, "len/tr\tevents\tfreq_int\n");		
		for(i=0; i<ncount; i++)
			{	
				for(j=0; j<ntrans; j++)
					{
						if(len_v[j]<-50+(i+1)*bin && len_v[j]>=-50+i*bin)
							{
								count[i]=count[i]+intens_v[j];
							}
					}
				tot=tot+count[i];
			}
		freq=0;
		for(i=0; i<ncount; i++)
			{	
				fprintf(out, "%i\t%i\t%f\n", i*bin, (int)count[i], count[i]/tot);
				freq=freq+(count[i]/tot);
			}
	}
}







int main()
{

int n;


FILE* fin;
FILE* out;

fin = fopen ("input.txt", "r");



if(fin==NULL)
{
printf("\nError\n\n");
}
else
{



//to count the number of lines of the input file (i.e. twice the number of transcripts under study)
n=nline(fin);
fclose(fin);



//main function: to obtain the number of ribosome per transcript and other information (also unavailable and wrong sequences are listed) with or without the ramp hypothesis
fin = fopen ("input.txt", "r");
out = fopen ("temporary.txt", "w");
output(fin, out, n, SLOWDOWN_RATE, RAMP_LENGTH, RAMP_HYOUTHESIS, PROTEIN_INPUT);
fclose(out);
fclose(fin);



//to eliminate results associated with unavailable or wrong sequences
fin = fopen ("temporary.txt", "r"); 
n=nline(fin);
fclose(fin);

fin = fopen ("temporary.txt", "r");
out = fopen ("results.txt", "w");
correctline(fin, out, n, PROTEIN_INPUT);
fclose(out);
fclose(fin);



//to generate the data for number of ribosome per transcript distributions
fin = fopen ("results.txt", "r"); 
n=nline(fin);
fclose(fin);

fin = fopen ("results.txt", "r"); 
out = fopen ("ribosome_per_transcript_distribution.txt", "w");
graph(fin, out, n, PROTEIN_INPUT);
fclose(out);
fclose(fin);



//to generate the data for length distribution with specific width of the bins
fin = fopen ("results.txt", "r"); 
n=nline(fin);
fclose(fin);

fin = fopen ("results.txt", "r"); 
out = fopen ("length_distribution.txt", "w");
graph_length(fin, out, n, PROTEIN_INPUT, WIDTH_BIN);
fclose(out);
fclose(fin);

}


return 0;
}


































