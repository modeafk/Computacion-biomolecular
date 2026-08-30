#include <iostream>
#include <string>
#include <vector>
#include <algorithm>


using namespace std;


void algor_wush(const string& seq1, const string& seq2, int igual=1, int diferente=-1, int gap=-2){

int n=seq1.length();
int m=seq2.length();

vector<vector<int>> score(n+1,vector<int>(m+1,0));

//para los bordes
for(int i=0; i <=n; i++){

    score[i][0]=i*gap;
}
for(int j=0; j<=m;j++){

    score[0][j]=j*gap;
}
//llenar la matriz

for(int i=1; i<=n;i++){
    for(int j=1; j<=m;j++){

     int diagonal1=diferente;
     if(seq1[i-1]==seq2[j-1]){
        diagonal1=igual;

     }
     int diagonal=score[i-1][j-1]+diagonal1;

     int arriba=score[i-1][j]+gap;
     int izquierda =score[i][j-1]+gap;
     
     score[i][j]=max({diagonal,arriba,izquierda});
    }
}

string primero="";
string segundo="";
int i=n,j=m;

while(i>0||j>0){

    int score_tmp;
    if (i > 0 && j > 0 && seq1[i - 1] == seq2[j - 1]) {
    score_tmp = igual;
    } else {
    score_tmp = diferente;
    }
    if (i > 0 && j > 0 && score[i][j] == score[i - 1][j - 1] + score_tmp) {
            primero += seq1[i - 1];
            segundo += seq2[j - 1];
            --i; --j;
        } else if (i > 0 && score[i][j] == score[i - 1][j] + gap) {
            primero += seq1[i - 1];
            segundo += '-';
            --i;
        } else if( j>0) {
            primero += '-';
            segundo += seq2[j - 1];
            --j;
        }
}

reverse(primero.begin(),primero.end());
reverse(segundo.begin(),segundo.end());
cout<<"Total: "<<score[n][m]<<endl;
cout<<"Primera cadena :" <<primero<<endl;
cout<<"Segunda cadena: "<<segundo<<endl;


}


int main(){

    string seq1="attaaaggtt";
    string seq2="tataccttcc";
    algor_wush(seq1,seq2);
    return 0;

}
