#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
using namespace std;

using Matrix = vector<vector<int>>;
using Vector = vector<int>;

// XOR hai vector nhị phân
Vector xor_vector(const Vector &a, const Vector &b)
{
    Vector res(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        res[i] = a[i] ^ b[i];
    return res;
}

// Nhân vector với ma trận nhị phân modulo 2
Vector vector_times_matrix(const Vector &vec, const Matrix &mat)
{
    Vector res(mat[0].size(), 0);
    for (size_t j = 0; j < mat[0].size(); ++j)
        for (size_t i = 0; i < vec.size(); ++i)
            res[j] ^= (vec[i] & mat[i][j]);
    return res;
}

// Tạo ma trận sinh G từ H theo dạng hệ thống H = [P | I]
Matrix build_generator_matrix(const Matrix &H)
{
    size_t r = H.size();    // số hàng
    size_t n = H[0].size(); // số cột
    size_t k = n - r;

    Matrix P(r, Vector(k));
    for (size_t i = 0; i < r; ++i)
        for (size_t j = 0; j < k; ++j)
            P[i][j] = H[i][j];

    Matrix G(k, Vector(n));
    for (size_t i = 0; i < k; ++i)
    {
        G[i][i] = 1; // Ma trận đơn vị
        for (size_t j = 0; j < r; ++j)
            G[i][k + j] = P[j][i]; // P^T
    }

    return G;
}

// Sinh vector ngẫu nhiên nhị phân
Vector generate_random_vector(size_t length)
{
    Vector v(length);
    for (size_t i = 0; i < length; ++i)
        v[i] = rand() % 2;
    return v;
}

// Syndrome H * c^T
Vector compute_syndrome(const Matrix &H, const Vector &c)
{
    Vector syn(H.size());
    for (size_t i = 0; i < H.size(); ++i)
        for (size_t j = 0; j < H[0].size(); ++j)
            syn[i] ^= H[i][j] & c[j];
    return syn;
}

// Bit-Flipping Decoder đơn giản
Vector bit_flipping(const Matrix &H, Vector c_received, int max_iter = 10)
{
    size_t n = c_received.size();
    size_t m = H.size();

    for (int iter = 0; iter < max_iter; ++iter)
    {
        Vector syndrome = compute_syndrome(H, c_received);
        if (all_of(syndrome.begin(), syndrome.end(), [](int v)
                   { return v == 0; }))
            break;

        // Đếm số lần bit vi phạm
        vector<int> error_count(n, 0);
        for (size_t i = 0; i < m; ++i)
        {
            if (syndrome[i])
            {
                for (size_t j = 0; j < n; ++j)
                    if (H[i][j])
                        error_count[j]++;
            }
        }

        // Lật những bit sai nhiều
        for (size_t i = 0; i < n; ++i)
            if (error_count[i] > 2)
                c_received[i] ^= 1;
    }

    return c_received;
}

void print_vector(const Vector &v, const string &label)
{
    cout << label << ": ";
    for (int bit : v)
        cout << bit;
    cout << endl;
}

int main()
{
    srand(time(nullptr));

    // Ma trận H (4x7) – ví dụ toy
    Matrix H = {
        {1, 1, 0, 1, 0, 0, 0},
        {0, 1, 1, 0, 1, 0, 0},
        {1, 0, 1, 0, 0, 1, 0},
        {1, 1, 1, 0, 0, 0, 1}};

    // Sinh ma trận G
    Matrix G = build_generator_matrix(H);
    size_t K = G.size();
    size_t N = G[0].size();

    // Tạo message ngẫu nhiên
    Vector u = generate_random_vector(K);
    Vector c = vector_times_matrix(u, G);

    // Giả lập nhiễu
    Vector noisy = c;
    noisy[1] ^= 1; // flip 1 bit

    // Giải mã
    Vector decoded = bit_flipping(H, noisy);

    // In kết quả
    print_vector(u, "Original data");
    print_vector(c, "Encoded codeword");
    print_vector(noisy, "Noisy received");
    print_vector(decoded, "Decoded codeword");

    Vector decoded_info(decoded.begin(), decoded.begin() + K);
    print_vector(decoded_info, "Decoded info bits");

    return 0;
}
