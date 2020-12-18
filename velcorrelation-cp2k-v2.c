//BY FILIPE MATUSALEM, DEC 2020     filipematus@gmail.com 
//Program to compute normalized velocity correlation function from CP2K xyz velocity file
//Compilation: g++ -o velcorrelation-cp2k.x velcorrelation-cp2k-v2.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>




int main(int argc, char *argv[])
{
FILE *velocity,*output;
float a,b;
int i,j,k,l,m,natoms,nspecies,nsteps,ntype[10],M;
char str1[150],ch,species[10][10],lixo[150];


if( argc < 2 ){
printf("\n\n");
printf("EEEEEEE   RRRRRRRR   RRRRRRRR   OOOOOOO   RRRRRRRR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EEEE      RRRRRRRR   RRRRRRRR   OO   OO   RRRRRRRR\n");
printf("EE        RRRR       RRRR       OO   OO   RRRR\n");
printf("EE        RR  RR     RR  RR     OO   OO   RR  RR\n");
printf("EEEEEEE   RR    RR   RR    RR   OOOOOOO   RR    RR\n\n");

printf("Enter the name of the cp2k velocity xyz file \n\n");

 exit(0);}

strcpy(str1,argv[1]);


velocity = fopen(str1,"r"); /* Arquivo ASCII, para leitura */
if(!velocity)
{
printf( "Error opening argument 1 file\n");
printf("Enter the name of the cp2k velocity xyz file \n");

exit(0);
}

output = fopen("velcorrelation.dat","w"); /* Arquivo ASCII, para escrita */
if(!output)
{
printf( "Error creating velcorrelation.dat file\n");
exit(0);
}

fscanf(velocity,"%d",&natoms);
printf("No. atoms = %d\n",natoms);

do
fscanf(velocity,"%s",str1);                                      /*posiciona o  após a palavra xxxxx*/
while(strcmp(str1,"time")!=0);

do
ch = getc(velocity);              /*chega ao fim da linha*/
while(ch!='\n');

fscanf(velocity,"%s",species[0]);

do
ch = getc(velocity);              /*chega ao fim da linha*/
while(ch!='\n');

j=k=1;
for(i=0;i<natoms-1;i++){
fscanf(velocity,"%s",species[j]);
k++;
if(strcmp(species[j-1],species[j])!=0){ntype[j-1]=k-1;k=1;j++;}

do
ch = getc(velocity);              /*chega ao fim da linha*/
while(ch!='\n');
}
ntype[j-1]=k;
nspecies=j;

printf("No. Species %d \n\n",nspecies);
printf("Specie  Number\n");
for(i=0;i<nspecies;i++){
printf("   %s      %d \n",species[i],ntype[i]);
}

nsteps=1;
while (fscanf(velocity,"%s",str1) != EOF){            /*conta steps*/
if(strcmp(str1,"time")==0)nsteps++;                      
}

printf("\nNo. steps = %d\n",nsteps);

rewind(velocity);


printf("\n---------------The maximum time lag (in units of MD steps) can be entered as a second argument (Default = No. steps/2)-----------------------\n");

if( argc > 2 ){
M=atoi(argv[2]);
if(M>nsteps/2){
printf("\n\n");
printf("EEEEEEE   RRRRRRRR   RRRRRRRR   OOOOOOO   RRRRRRRR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EEEE      RRRRRRRR   RRRRRRRR   OO   OO   RRRRRRRR\n");
printf("EE        RRRR       RRRR       OO   OO   RRRR\n");
printf("EE        RR  RR     RR  RR     OO   OO   RR  RR\n");
printf("EEEEEEE   RR    RR   RR    RR   OOOOOOO   RR    RR\n\n");

printf("Time lag should be <= No. steps/2 \n\n");

 exit(0);}
}
else {M=nsteps/2;}

printf("\nMaximum time lag is = %d MD steps\n\n",M);


float* velx = new float[natoms*nsteps];
float* vely = new float[natoms*nsteps];
float* velz = new float[natoms*nsteps];

for(j=0;j<nsteps;j++){
do
fscanf(velocity,"%s",str1);                                      /*posiciona o  após a palavra "xxx"*/
while(strcmp(str1,"time")!=0);

do
ch = getc(velocity);              /*chega ao fim da linha*/
while(ch!='\n');

for(i=0;i<natoms;i++){

fscanf(velocity,"%s",lixo);
fscanf(velocity,"%f",&velx[i+j*natoms]); 
fscanf(velocity,"%f",&vely[i+j*natoms]); 
fscanf(velocity,"%f",&velz[i+j*natoms]); 
}

}
fclose(velocity);

float vel0[natoms][3],vel1[natoms][3],vaf[natoms],vaf_per_type[nspecies],norm[nspecies];

//header of output
fprintf(output,"#Step   ");for(i=0;i<nspecies;i++)fprintf(output,"  %s       ",species[i]);fprintf(output,"  No. steps %d   max time lag %d   ",nsteps,M);fprintf(output,"\n");


printf("working (several minutes are required depending on the number of steps)\n\n");
printf("%%\n");
printf("0 ");
b=0;

//loop on vcf(j)
for(j=0;j<M;j++){

//to print the evolution on screen
a=j*100/M;
if(a!=b) fprintf(stderr,"%.0f ",a);   //fprintf(stderr is used  to print imediately on the screen
b=a;


for(i=0;i<nspecies;i++)vaf_per_type[i]=0;
for(i=0;i<natoms;i++)vaf[i]=0;

//loop \sum_{l=0}^{n-1-M}
for(l=0;l<=nsteps-M-1;l++){

//debug
//printf("%d%d ",l,l+j);

//compute vaf
for(i=0;i<natoms;i++){
vaf[i]=vaf[i]+velx[i+(l)*natoms]*velx[i+(l+j)*natoms] + vely[i+(l)*natoms]*vely[i+(l+j)*natoms] + velz[i+(l)*natoms]*velz[i+(l+j)*natoms];
//debug
//printf("%f\n",vaf[i]);
}
}

//debug
//printf("\n");

fprintf(output,"%d ",j);
//average for each specie
i=0;
for(k=0;k<nspecies;k++){for(m=0;m<ntype[k];m++)vaf_per_type[k]=vaf_per_type[k]+vaf[m+i];
i=i+m;

//average by natoms
a=vaf_per_type[k]/ntype[k];

//average by time lag
a=a/(nsteps-M);

if(j==0)norm[k]=a;

fprintf(output,"%e  ",a/norm[k]);
}
fprintf(output,"\n");
fflush(output);        //fflush is used to printo to the file imediate
}

printf("\n\nDONE!! VAF written to velcorrelation.dat file!!\n\n");


fclose(output);
}
