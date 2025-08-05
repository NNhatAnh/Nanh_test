#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iomanip>
using namespace std;
using Matrix = vector<vector<int>>;
const int M = 512;
const int K = 1024;
const int N = 2048;
const int MAX_ITER = 50;
const int NUM_FRAMES = 100;
const int theta[8] = {3, 0, 1, 2, 2, 3, 0, 1};
const int phi[8][4] = {
    {16, 0, 0, 0},
    {103, 53, 8, 35},
    {105, 74, 119, 97},
    {0, 45, 89, 112},
    {50, 47, 31, 64},
    {29, 0, 122, 93},
    {115, 59, 1, 99},
    {30, 102, 69, 94}
};

default_random_engine rng(random_device{}());
uniform_real_distribution<double> dist(0.0, 1.0);

// === Matrix Utilities ===

Matrix identity_matrix(int M) {
    cout<<"M: "<<M<<endl;
    Matrix I(M, vector<int>(M, 0));
    for (int i = 0; i < M; ++i) I[i][i] = 1;
    return I;
}

Matrix zero_matrix(int M) {
    return vector<vector<int>>(M, vector<int>(M, 0));
}

Matrix xor_matrix(const Matrix& A, const Matrix& B) {
    Matrix C = A;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            C[i][j] ^= B[i][j];
    return C;
}

Matrix generate_permutation_matrix(int M, int k) {
    Matrix P(M, vector<int>(M, 0));
    for (int i = 0; i < M; ++i) {
        int floor_4i_M = (4 * i) / M;
        int t = ((theta[k] + floor_4i_M) % 4) * (M / 4);
        int f = (phi[k][floor_4i_M] + i) % (M / 4);
        int j = t + f;
        P[i][j] = 1;
    }
    return P;
}

Matrix mat_mult(const Matrix& A, const Matrix& B) {
    int n = A.size(), m = B[0].size(), p = B.size();
    Matrix C(n, vector<int>(m, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            for (int k = 0; k < p; ++k)
                C[i][j] ^= A[i][k] & B[k][j];
    return C;
}

// Sai ở đây
Matrix mat_inv(Matrix A) {
    int n = A.size();
    Matrix I = identity_matrix(n);
    for (int i = 0; i < n; ++i) {
        if (A[i][i] == 0) {
            bool found = false;
            for (int j = i + 1; j < n; ++j) {
                // cout << "A["<<j<<"]["<<i<<"] = "<<A[j][i]<<endl;
                if (A[j][i]) {
                    swap(A[i], A[j]);
                    swap(I[i], I[j]);
                    found = true;
                    break;
                }
            }
            if (!found) {
                cout<< "Matrix is not invertible";
                exit(1);
            } else {
                cout<< "Matrix invertible success";
            }
        } else {
            cout<<"A[i][i] invalid";
            exit(1);
        }
        for (int j = 0; j < n; ++j) {
            if (i != j && A[j][i]) {
                for (int k = 0; k < n; ++k) {
                    A[j][k] ^= A[i][k];
                    I[j][k] ^= I[i][k];
                }
            }
        }
    }
    return I;
}

// === Build H and G ===

void build_H_irig106(Matrix& H) {
    cout<<"Matrix H" <<endl;
    H.assign(2 * M, vector<int>(4 * M, 0));
    auto Pi = vector<Matrix>(8);
    for (int k = 0; k < 8; ++k)
        Pi[k] = generate_permutation_matrix(M, k);
    auto I = identity_matrix(M);
    auto Z = zero_matrix(M);
        
    vector<Matrix> row0 = {Z, Z, I, Z, xor_matrix(I, Pi[0]), Z};
    vector<Matrix> row1 = {I, I, Z, I, xor_matrix(xor_matrix(Pi[1], Pi[2]), Pi[3]), Z};
    // vector<Matrix> row2 = {I, xor_matrix(Pi[4], Pi[5]), Z, xor_matrix(Pi[6], Pi[7]), I, Z};
    
    for (int r = 0; r < 2; ++r) {
        const auto& row = (r == 0) ? row0 : row1;
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < M; ++k)
                    H[r * M + i][j * M + k] = row[j][i][k];;
                    // cout << "row["<<j<<"]["<<i<<"]["<<"k] = "<<row[j][i][k]<<endl;
    }

    cout<<H.size()<<"x"<<H[0].size()<<endl;
}

void build_G_from_H(const Matrix& H, Matrix& G) {
    cout<<"Matrix G" <<endl;
    Matrix P(2 * M, vector<int>(2 * M)), Q(2 * M, vector<int>(2 * M));

    for (int i = 0; i < 2 * M; ++i)
        for (int j = 0; j < 2 * M; ++j) {
            Q[i][j] = H[i][j];
            P[i][j] = H[i][j + 2 * M];
        }

    Matrix P_inv = mat_inv(P);
    Matrix W_T = mat_mult(P_inv, Q);
    Matrix W(2 * M, vector<int>(2 * M));
    for (int i = 0; i < 2 * M; ++i)
        for (int j = 0; j < 2 * M; ++j)
            W[j][i] = W_T[i][j];
    G.assign(K, vector<int>(N, 0));
    for (int i = 0; i < K; ++i) {
        G[i][i] = 1;
        for (int j = 0; j < N - K; ++j)
            G[i][j + K] = W[i][j];
    }
}

// === Channel & Decoder ===

vector<int> random_bits(int len) {
    vector<int> bits(len);
    for (int i = 0; i < len; ++i)
        bits[i] = dist(rng) > 0.5 ? 1 : 0;
    return bits;
}

vector<int> encode(const vector<int>& u, const Matrix& G) {
    vector<int> x(N, 0);
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < K; ++i)
            x[j] ^= u[i] & G[i][j];
    return x;
}

vector<double> add_awgn_noise(const vector<int>& x, double snr_db) {
    double snr_linear = pow(10.0, snr_db / 10.0);
    double sigma = sqrt(1.0 / (2 * snr_linear));
    normal_distribution<double> noise(0.0, sigma);
    vector<double> y(N);
    for (int i = 0; i < N; ++i)
        y[i] = (1 - 2 * x[i]) + noise(rng);
    return y;
}

vector<double> calc_llr(const vector<double>& y, double snr_db) {
    double snr_linear = pow(10.0, snr_db / 10.0);
    double sigma2 = 1.0 / (2 * snr_linear);
    vector<double> llr(N);
    for (int i = 0; i < N; ++i)
        llr[i] = 2 * y[i] / sigma2;
    return llr;
}

vector<int> sum_product_decode(const vector<double>& llr, const Matrix& H, int& iterations) {
    int m = H.size(), n = H[0].size();
    vector<vector<double>> v_to_c(m, vector<double>(n, 0));
    vector<vector<double>> c_to_v(m, vector<double>(n, 0));
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; ++i)
            if (H[i][j]) v_to_c[i][j] = llr[j];

    for (iterations = 0; iterations < MAX_ITER; ++iterations) {
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                if (H[i][j]) {
                    double prod = 1.0;
                    for (int k = 0; k < n; ++k)
                        if (k != j && H[i][k])
                            prod *= tanh(v_to_c[i][k] / 2.0);
                    const double eps = 1e-12;
prod = max(min(prod, 1.0 - eps), -1.0 + eps);
c_to_v[i][j] = 2 * atanh(prod);
                }

        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i)
                if (H[i][j]) {
                    double sum = llr[j];
                    for (int k = 0; k < m; ++k)
                        if (k != i && H[k][j])
                            sum += c_to_v[k][j];
                    v_to_c[i][j] = sum;
                }

        vector<int> hard(n);
        for (int j = 0; j < n; ++j) {
            double total = llr[j];
            for (int i = 0; i < m; ++i)
                if (H[i][j]) total += c_to_v[i][j];
            hard[j] = total < 0;
        }

        bool pass = true;
        for (int i = 0; i < m; ++i) {
            int sum = 0;
            for (int j = 0; j < n; ++j)
                sum ^= H[i][j] & hard[j];
            if (sum) { pass = false; break; }
        }
        if (pass) return hard;
    }

    vector<int> hard(n);
    for (int j = 0; j < n; ++j) {
        double total = llr[j];
        for (int i = 0; i < m; ++i)
            if (H[i][j]) total += c_to_v[i][j];
        hard[j] = total < 0;
    }
    return hard;
}

// === MAIN ===

int main() {
    Matrix H, G;
    build_H_irig106(H);
    build_G_from_H(H, G);

    cout << "\n### Performance\n";
    cout << "SNR(dB) | Success | Rate(%) | Avg_Iter | BER\n";
    cout << "-----------------------------------------------\n";

    vector<double> snr_list = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};
    const int MAX_ERRORS = 10;

    for (double snr : snr_list) {
        int success = 0, total_iter = 0, bit_errors = 0;
        vector<int> error_stats(MAX_ERRORS + 1, 0);  // Đếm số frame với 0,1,...,10 lỗi

        for (int f = 0; f < NUM_FRAMES; ++f) {
            auto u = random_bits(K);
            auto x = encode(u, G);
            auto y = add_awgn_noise(x, snr);
            auto llr = calc_llr(y, snr);
            int iter = 0;
            auto decoded = sum_product_decode(llr, H, iter);
            total_iter += iter;

            vector<int> u_hat(decoded.begin(), decoded.begin() + K);

            int frame_errors = 0;
            for (int i = 0; i < K; ++i) {
                if (u[i] != u_hat[i]) {
                    frame_errors++;
                    bit_errors++;
                }
            }

            if (frame_errors == 0)
                success++;

            if (frame_errors <= MAX_ERRORS)
                error_stats[frame_errors]++;
        }

        double rate = 100.0 * success / NUM_FRAMES;
        double avg_iter = (double)total_iter / NUM_FRAMES;
        double ber = (double)bit_errors / (NUM_FRAMES * K);

        printf("  %4.1f   |   %4d   |  %6.1f |   %.1f   | %.5f\n",
               snr, success, rate, avg_iter, ber);

        // In bảng Error Correction
        cout << "\n### Error Correction (SNR = " << snr << " dB)\n";
        cout << "Errors | Success | Rate(%)\n";
        cout << "----------------------------\n";
        for (int i = 0; i <= MAX_ERRORS; ++i) {
            double err_rate = 100.0 * error_stats[i] / NUM_FRAMES;
            printf("  %2d    |  %4d    |  %6.1f\n", i, error_stats[i], err_rate);
        }

        cout << endl;
    }

    return 0;
}
