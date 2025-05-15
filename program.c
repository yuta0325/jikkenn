#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 
#define MAX_SEQ_NUM 30 
#define MAX_GENE_NUM 8
#define back_A 7519429.0
#define back_C 4637676.0
#define back_G 4637676.0
#define back_T 7519429.0


char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("motif_region_file open error.\n");
    exit(1); //ファイルが開けなかった場合プログラムを終了
  }

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
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
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む

  printf("motif region:\n");
  for(int i = 0; i < seq_num; i++){
    printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
  }
  printf("\n");

  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  
  printf("promoter_sequence:\n");
  for(int i = 0; i < gene_num; i++){
    printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
    printf("%s\n", g_pro[i].seq);
  }
  int B=0;
  int i=0;
  int j=0;
  while(g_motif[0][i]!='\0')
  {
    B++;
    i++;
  }
  int hindo[4][B];
  for(i=0;i<4;i++)
  {
    for(j=0;j<B;j++)
    {
      hindo[i][j]=0;
    }
  }
  
  for(i=0;i<15;i++)
  {
    for(j=0;j<B;j++)
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

  for(i=0;i<4;i++)
  {
    for(j=0;j<B;j++)
    {
      printf("%2d,",hindo[i][j]);
      hindo[i][j]+=1;
    }
    printf("\n");
  }

  float kakuritu[4][B];
  int sum=hindo[0][0]+hindo[1][0]+hindo[2][0]+hindo[3][0];
  for(i=0;i<4;i++)
  {
    for(j=0;j<B;j++)
    {
      kakuritu[i][j]=(float)hindo[i][j]/sum;
    }
  }
  for(i=0;i<4;i++)
  {
    for(j=0;j<B;j++)
    {
      printf("%2f,",kakuritu[i][j]);
    }
    printf("\n");
  }

  float back_sum=back_A+back_C+back_G+back_T;
  float qA= back_A/back_sum;
  float qC= back_C/back_sum;
  float qG= back_G/back_sum;
  float qT= back_T/back_sum;
  float q[4]={qA,qC,qG,qT};
  
  float ozzu[4][B];
  for(i=0;i<4;i++)
  {
    for(j=0;j<B;j++)
    {
      ozzu[i][j]=log(kakuritu[i][j]/q[i]);
    }
  }
  for(i=0;i<4;i++)
  {
    for(j=0;j<B;j++)
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
  int m=0;
  float score=0;
  float t=0;
  int h;
  for(int d=0;d<8;d++)
  {
    while(m<genelen-B+1)
    {
      for(i=0;i<B;i++)
      {
        if(g_pro[m+i].seq=='A')
        {
          score+=ozzu[0][i];
        }
        else if(g_pro[m+i].seq=='C')
        {
          score+=ozzu[1][i];
        }
        else if(g_pro[m+i].seq=='G')
        {
          score+=ozzu[2][i];
        }
        else if(g_pro[m+i].seq=='T')
        {
          score+=ozzu[3][i];
        }
      }
      if(t<score)
      {
        t=score;
        h=m;
      }
      m++;
    }
    printf("name=%s\n",g_pro[d].name);
    printf("score=%f\n",score);
    for(j=h;j<B;j++)
    {
      printf("%c",g_pro[j].seq);
    }
  }
  return 0;
}