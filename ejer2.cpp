#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cctype>

using namespace std;

//preprocesamiento para eliminar los numeros que contaban el txt
string limpiar_secuencia(const string& cruda) {
    string limpia = "";
    for (char c : cruda) {
        if (isalpha(c)) {
            limpia += toupper(c);
        }
    }
    return limpia;
}

//backtracking para las soluciones de las cadenas pequeñas
void encontrar_todos_alineamientos(const vector<vector<int>>& score, const string& seq1, const string& seq2, 
                                   int i, int j, string align1, string align2, int igual, int diferente, int gap) {
    if (i == 0 && j == 0) {
        reverse(align1.begin(), align1.end());
        reverse(align2.begin(), align2.end());
        cout << "Opcion de alineamiento:\n" << align1 << "\n" << align2 << "\n\n";
        return;
    }

    int score_tmp = (i > 0 && j > 0 && seq1[i - 1] == seq2[j - 1]) ? igual : diferente;

    //diagonal
    if (i > 0 && j > 0 && score[i][j] == score[i - 1][j - 1] + score_tmp) {
        encontrar_todos_alineamientos(score, seq1, seq2, i - 1, j - 1, align1 + seq1[i - 1], align2 + seq2[j - 1], igual, diferente, gap);
    }
    //arriba
    if (i > 0 && score[i][j] == score[i - 1][j] + gap) {
        encontrar_todos_alineamientos(score, seq1, seq2, i - 1, j, align1 + seq1[i - 1], align2 + "-", igual, diferente, gap);
    }
    //izquierda
    if (j > 0 && score[i][j] == score[i][j - 1] + gap) {
        encontrar_todos_alineamientos(score, seq1, seq2, i, j - 1, align1 + "-", align2 + seq2[j - 1], igual, diferente, gap);
    }
}

void algor_wush_avanzado(const string& seq1, const string& seq2, bool mostrar_todas, const string& nombre_dotplot, int igual=1, int diferente=-1, int gap=-2) {
    int n = seq1.length();
    int m = seq2.length();
    
    vector<vector<int>> score(n + 1, vector<int>(m + 1, 0));

    for (int i = 0; i <= n; i++) score[i][0] = i * gap;
    for (int j = 0; j <= m; j++) score[0][j] = j * gap;

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int diagonal1 = (seq1[i - 1] == seq2[j - 1]) ? igual : diferente;
            int diagonal = score[i - 1][j - 1] + diagonal1;
            int arriba = score[i - 1][j] + gap;
            int izquierda = score[i][j - 1] + gap;
            score[i][j] = max({diagonal, arriba, izquierda});
        }
    }

    cout << "Score Total: " << score[n][m] << endl;

    if (mostrar_todas) {
        encontrar_todos_alineamientos(score, seq1, seq2, n, m, "", "", igual, diferente, gap);
    } else {
        string primero = "", segundo = "";
        int i = n, j = m;

        while (i > 0 || j > 0) {
            int score_tmp = (i > 0 && j > 0 && seq1[i - 1] == seq2[j - 1]) ? igual : diferente;
            
            //evaluamos la diagonal primero
            if (i > 0 && j > 0 && score[i][j] == score[i - 1][j - 1] + score_tmp) {
                primero += seq1[i - 1];
                segundo += seq2[j - 1];
                --i; --j;
            } else if (i > 0 && score[i][j] == score[i - 1][j] + gap) {
                primero += seq1[i - 1];
                segundo += '-';
                --i;
            } else {
                primero += '-';
                segundo += seq2[j - 1];
                --j;
            }
        }
        reverse(primero.begin(), primero.end());
        reverse(segundo.begin(), segundo.end());
        cout << "Alineamiento optimo \n";
    }

    // Exportar Dot Plot a archivo de texto
    if (!nombre_dotplot.empty()) {
        ofstream out(nombre_dotplot);
        if (out.is_open()) {
            out << "  " << seq2 << "\n";
            for (int i = 0; i < n; i++) {
                out << seq1[i] << " ";
                for (int j = 0; j < m; j++) {
                    if (seq1[i] == seq2[j]) out << "*";
                    else out << " ";
                }
                out << "\n";
            }
            out.close();
            cout << "grafico guarado en : " << nombre_dotplot << "\n";
        }
    }
}

int main() {
    cout << "pueba con AAAC vs AGC ---\n";
    algor_wush_avanzado("AAAC", "AGC", true, "");

    ifstream archivo("Sequencias.txt");
    string linea, bac_raw, sars_raw, inf_raw;
    int estado = 0;

    if (archivo.is_open()) {
        while (getline(archivo, linea)) {
            if (linea.find("Bacteria") != string::npos) estado = 1;
            else if (linea.find("Sars-Cov") != string::npos) estado = 2;
            else if (linea.find("Influenza") != string::npos) estado = 3;
            else {
                if (estado == 1) bac_raw += linea;
                else if (estado == 2) sars_raw += linea;
                else if (estado == 3) inf_raw += linea;
            }
        }
        archivo.close();
    } else {
        cout << "\nError: No se encontro el archivo Sequencias.txt\n";
        return 1;
    }

    string bacteria = limpiar_secuencia(bac_raw);
    string sars = limpiar_secuencia(sars_raw);
    string influenza = limpiar_secuencia(inf_raw);

    cout << "\n BACTERIA vs SARS-COV-2 \n";
    algor_wush_avanzado(bacteria, sars, false, "Bac_Sars.txt");

    cout << "\nBACTERIA vs INFLUENZA \n";
    algor_wush_avanzado(bacteria, influenza, false, "Bac_Inf.txt");

    cout << "\nSARS-COV-2 vs INFLUENZA\n";
    algor_wush_avanzado(sars, influenza, false, "Sars_Inf.txt");

    return 0;
}