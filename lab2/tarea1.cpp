#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// Parametros de puntuacion
const int MATCH_SCORE = 1;
const int MISMATCH_SCORE = -1;
const int GAP_PENALTY = -1;

struct AlignmentPair {
    int score;
    string seq1;
    string seq2;
};

// Algoritmo Needleman-Wunsch (Alineamiento Pareado Global)
AlignmentPair needlemanWunsch(const string& s1, const string& s2) {
    int n = s1.length();
    int m = s2.length();
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));

    for (int i = 0; i <= n; ++i) dp[i][0] = i * GAP_PENALTY;
    for (int j = 0; j <= m; ++j) dp[0][j] = j * GAP_PENALTY;

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int match = dp[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_SCORE);
            int del = dp[i - 1][j] + GAP_PENALTY;
            int ins = dp[i][j - 1] + GAP_PENALTY;
            dp[i][j] = max({match, del, ins});
        }
    }

    string a1 = "", a2 = "";
    int i = n, j = m;
    while (i > 0 || j > 0) {
        int current = dp[i][j];
        if (i > 0 && j > 0 && current == dp[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? MATCH_SCORE : MISMATCH_SCORE)) {
            a1 += s1[i - 1]; a2 += s2[j - 1]; i--; j--;
        } else if (i > 0 && current == dp[i - 1][j] + GAP_PENALTY) {
            a1 += s1[i - 1]; a2 += '-'; i--;
        } else {
            a1 += '-'; a2 += s2[j - 1]; j--;
        }
    }
    reverse(a1.begin(), a1.end());
    reverse(a2.begin(), a2.end());

    return {dp[n][m], a1, a2};
}

// Calculo del Score SP (Sum-of-Pairs)
int calculateSPScore(const vector<string>& msa) {
    int totalScore = 0;
    int k = msa.size();
    if (k == 0) return 0;
    int len = msa[0].length();

    for (int i = 0; i < k; ++i) {
        for (int j = i + 1; j < k; ++j) {
            for (int p = 0; p < len; ++p) {
                char c1 = msa[i][p];
                char c2 = msa[j][p];
                if (c1 == '-' && c2 == '-') continue;
                if (c1 == '-' || c2 == '-') totalScore += GAP_PENALTY;
                else if (c1 == c2) totalScore += MATCH_SCORE;
                else totalScore += MISMATCH_SCORE;
            }
        }
    }
    return totalScore;
}

// Algoritmo MSA Estrella con ensamblado mediante Gaps Maximos
void runStarMSA(const vector<string>& sequences, const string& label) {
    int k = sequences.size();
    vector<vector<int>> scoreMatrix(k, vector<int>(k, 0));
    vector<int> totalScores(k, 0);

    // 1. Matriz de alineamiento pareado
    for (int i = 0; i < k; ++i) {
        for (int j = i + 1; j < k; ++j) {
            AlignmentPair res = needlemanWunsch(sequences[i], sequences[j]);
            scoreMatrix[i][j] = res.score;
            scoreMatrix[j][i] = res.score;
        }
    }

    // 2. Suma de puntuaciones para hallar la secuencia centro
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            if (i != j) totalScores[i] += scoreMatrix[i][j];
        }
    }

    int centerIdx = 0;
    int maxScore = totalScores[0];
    for (int i = 1; i < k; ++i) {
        if (totalScores[i] > maxScore) {
            maxScore = totalScores[i];
            centerIdx = i;
        }
    }

    cout << "=== " << label << " ===" << endl;
    cout << "Secuencia Centro elegida: S" << (centerIdx + 1) << " (Score acumulado: " << maxScore << ")" << endl;

    string Sc = sequences[centerIdx];
    int L = Sc.length();

    // Estructuras de descomposición por posición
    vector<vector<string>> ins(k, vector<string>(L + 1, ""));
    vector<vector<char>> ch(k, vector<char>(L, '-'));

    // Configurar secuencia centro
    for (int p = 0; p < L; ++p) {
        ch[centerIdx][p] = Sc[p];
    }

    // Alinear cada secuencia contra el centro
    for (int i = 0; i < k; ++i) {
        if (i == centerIdx) continue;

        AlignmentPair res = needlemanWunsch(Sc, sequences[i]);
        int p = 0;
        for (size_t pos = 0; pos < res.seq1.length(); ++pos) {
            if (res.seq1[pos] == '-') {
                ins[i][p] += res.seq2[pos];
            } else {
                ch[i][p] = res.seq2[pos];
                p++;
            }
        }
    }

    // Calcular gaps maximos necesarios entre caracteres del centro
    vector<int> maxGaps(L + 1, 0);
    for (int p = 0; p <= L; ++p) {
        for (int i = 0; i < k; ++i) {
            maxGaps[p] = max(maxGaps[p], (int)ins[i][p].length());
        }
    }

    // Reconstruccion final alineada con longitud uniforme
    vector<string> alignedSeqs(k, "");
    for (int p = 0; p <= L; ++p) {
        int g = maxGaps[p];
        for (int i = 0; i < k; ++i) {
            string currentIns = ins[i][p];
            currentIns.append(g - currentIns.length(), '-');
            alignedSeqs[i] += currentIns;
            if (p < L) {
                alignedSeqs[i] += ch[i][p];
            }
        }
    }

    cout << "Alineamiento Múltiple Resultante:" << endl;
    for (int i = 0; i < k; ++i) {
        cout << "S" << (i + 1) << ": " << alignedSeqs[i] << endl;
    }

    int spScore = calculateSPScore(alignedSeqs);
    cout << "Score Final (Sum-of-Pairs): " << spScore << "\n\n";
}

int main() {
    // 1. Ejemplo de clase
    vector<string> classExample = {
        "ATTGCCATT",
        "ATGGCCATT",
        "ATCCAATTTT",
        "ATCTTCTT",
        "ACTGACC"
    };
    runStarMSA(classExample, "Ejemplo de Clase (S1 - S5)");

    // 2. Muestras BRCA1 Forward (F)
    vector<string> brca1_F = {
        "TGCCGGCAGGGATGTGCTTG",
        "GTTTAGGTTTTTGCTTATGCAGCATCCA",
        "GGAAAAGCACAGAACTGGCCAACA",
        "GCCAGTTGGTTGATTTCCACCTCCA",
        "ACCCCCGACATGCAGAAGCTG",
        "TGACGTGTCTGCTCCACTTCCA"
    };
    runStarMSA(brca1_F, "BRCA1 - Cadenas Forward (F)");

    // 3. Muestras BRCA1 Reverse (R)
    vector<string> brca1_R = {
        "TGCTTGCAGTTTGCTTTCACTGATGGA",
        "TCAGGTACCCTGACCTTCTCTGAAC",
        "GTGGGTTGTAAAGGTCCCAAATGGT",
        "TGCCTTGGGTCCCTCTGACTGG",
        "GTGGTGCATTGATGGAAGGAAGCA",
        "AGTGAGAGGAGCTCCCAGGGC"
    };
    runStarMSA(brca1_R, "BRCA1 - Cadenas Reverse (R)");

    return 0;
}