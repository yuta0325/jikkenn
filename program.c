#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 
#define MAX_SEQ_NUM 30 
#define MAX_GENE_NUM 8
#define BACK_A 7519429.0
#define BACK_C 4637676.0
#define BACK_G 4637676.0
#define BACK_T 7519429.0
#define ATGC 4
#define large 10000
int B=0,i=0,j=0;
int hindo[ATGC][15];
float kakuritu[ATGC][15];
float ozzu[ATGC][15];
int property(int w);
int wariai(int x);
int ozz(int y,int u);

char g_motif[MAX_SEQ_NUM][BUFSIZE]; 

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; 

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("motif_region_file open error.\n");
    exit(1); 
  }

  while(fscanf(fp, "%s", buffer) != EOF){ 
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; 
    }
    strcpy(g_motif[seq_num],buffer); 
    seq_num++;
  }
  return seq_num;
}

int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1);
  }

  while(fscanf(fp, "%s", buffer) != EOF){
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0';
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); 
    }else{
      strcpy(g_pro[gene_num].seq,buffer);
      gene_num++;
    }    
  }
  return gene_num;
}

int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); 
  int gene_num = read_promoter(argv[2]);  

  while(g_motif[0][i]!='\0')
  {
    B++;
    i++;
  }
  property(B);
  wariai(B);
  ozz(B,gene_num);

  char random[large];
  int num_random[large];
  char randomhairetu[10][large];
  int number_random[large-500];
  for(i=0;i<large;i++)
  {
    num_random[i]=rand()%10000;
    if(num_random[i]<3093)
    {
      random[i]='A';
    }
    else if(num_random[i]<5000)
    {
      random[i]='C';
    }
    else if(num_random[i]<6907)
    {
      random[i]='G';
    }
    else if(num_random[i]<10000)
    {
      random[i]='T';
    }
  }
  int z=0;
  for(i=0;i<10;i++)
  {
     number_random[i]=rand()%9500;
    for(j=0;j<large;j++)
    {
      z=number_random[i]+j;
      randomhairetu[i][j]=random[z];
    }
  }
}

int property(int w)
{
  for(i=0;i<ATGC;i++)
  {
    for(j=0;j<w;j++)
    {
      hindo[i][j]=0;
    }
  }
  
  for(i=0;i<15;i++)
  {
    for(j=0;j<w;j++)
    {
      if(g_motif[i][j]=='A'){
        hindo[0][j]++;
      }
      else if(g_motif[i][j]=='C')
      {
        hindo[1][j]++;
      }
      else if(g_motif[i][j]=='G')
      {
        hindo[2][j]++;
      }
      else if(g_motif[i][j]=='T')
      {
        hindo[3][j]++;
      }
    }
  }

  for(i=0;i<ATGC;i++)
  {
    for(j=0;j<w;j++)
    {
      printf("%2d,",hindo[i][j]);
      hindo[i][j]+=1;
    }
    printf("\n");
  }
  return 0;
}

int wariai(int x)
{
  int sum=hindo[0][0]+hindo[1][0]+hindo[2][0]+hindo[3][0];
  for(i=0;i<4;i++)
  {
    for(j=0;j<x;j++)
    {
      kakuritu[i][j]=(float)hindo[i][j]/sum;
    }
  }
  for(i=0;i<4;i++)
  {
    for(j=0;j<x;j++)
    {
      printf("%2f,",kakuritu[i][j]);
    }
    printf("\n");
  }
  return 0;
}

int ozz(int y,int u)
{
  float back_sum=BACK_A+BACK_C+BACK_G+BACK_T;
  float qA= BACK_A/back_sum;
  float qC= BACK_C/back_sum;
  float qG= BACK_G/back_sum;
  float qT= BACK_T/back_sum;
  float q[ATGC]={qA,qC,qG,qT};
  printf("%f, %f",qA,qC);
  for(i=0;i<ATGC;i++)
  {
    for(j=0;j<y;j++)
    {
      ozzu[i][j]=log(kakuritu[i][j]/q[i]);
    }
  }
  for(i=0;i<ATGC;i++)
  {
    for(j=0;j<y;j++)
    {
      printf("%10f,",ozzu[i][j]);
    }
    printf("\n");
  }

  int genelen=0;
  int l=0;
  while(g_pro[0].seq[l]!='\0')
  {
    genelen++;
    l++;
  }
  for(int d=0;d<u;d++)
  {
    float t=0.0,score=0.0;
    int m=0,h=0;
    while(m<genelen-y+1)
    {
      score=0.0;
      for(i=0;i<y;i++)
      {
        if(g_pro[d].seq[m+i]=='A')
        {
          score=score+ozzu[0][i];
        }
        else if(g_pro[d].seq[m+i]=='C')
        {
          score=score+ozzu[1][i];
        }
        else if(g_pro[d].seq[m+i]=='G')
        {
          score=score+ozzu[2][i];
        }
        else if(g_pro[d].seq[m+i]=='T')
        {
          score=score+ozzu[3][i];
        }
      }
      if(t<score)
      {
        t=score;
        h=m;
      }
      m++;
    }
    printf("%s\n",g_pro[d].name);
    printf("score: %f\n",t);
    printf("gene : ");
    for(j=h;j<h+y;j++)
    {
      printf("%c",g_pro[d].seq[j]);
    }
    printf("\n\n");
  }
}