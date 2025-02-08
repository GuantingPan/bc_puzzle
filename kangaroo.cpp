#include <gmp.h>
#include <gmpxx.h>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <set>

using namespace std;

typedef array<mpz_class, 2> Point;

const mpz_class modulo("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
const mpz_class Gx("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
const mpz_class Gy("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16);
const Point PG = {Gx, Gy};
const Point Z = {0, 0};

auto starttime = chrono::high_resolution_clock::now();

Point add(const Point& P, const Point& Q, const mpz_class& p = modulo) {
    if (P == Z) return Q;
    if (Q == Z) return P;
    const mpz_class& Px = P[0];
    const mpz_class& Py = P[1];
    const mpz_class& Qx = Q[0];
    const mpz_class& Qy = Q[1];
    if (Px == Qx && (Py != Qy || Py == 0)) return Z;
    mpz_class m, num, denom, inv;
    if (Px == Qx) {
        num = 3 * Px * Px;
        denom = 2 * Py;
    } else {
        num = Qy - Py;
        denom = Qx - Px;
    }
    mpz_invert(inv.get_mpz_t(), denom.get_mpz_t(), p.get_mpz_t());
    m = (num * inv) % p;
    mpz_class x = (m * m - Px - Qx) % p;
    if (x < 0) x += p;
    mpz_class y = (m * (Px - x) - Py) % p;
    if (y < 0) y += p;
    return {x, y};
}

Point mul(const mpz_class& k, const Point& P = PG, const mpz_class& p = modulo) {
    Point R = Z;
    Point current = P;
    mpz_class k_copy = k;
    while (k_copy > 0) {
        if (k_copy % 2 == 1) {
            R = add(R, current, p);
        }
        current = add(current, current, p);
        k_copy >>= 1;
    }
    return R;
}

mpz_class X2Y(const mpz_class& X, int y_parity, const mpz_class& p = modulo) {
    mpz_class X_cubed = (X * X * X) % p;
    mpz_class tmp = (X_cubed + mpz_class(7)) % p;
    mpz_class Y;
    mpz_class exp = (p + mpz_class(1)) / mpz_class(4);
    mpz_powm(Y.get_mpz_t(), tmp.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    if ((Y % 2) != y_parity) {
        Y = p - Y;
    }
    return Y;
}

bool comparator(const Point& P, const mpz_class& Pindex, const mpz_class& DP_rarity,
                std::vector<Point>& T, std::vector<mpz_class>& t, const std::vector<Point>& W,
                const std::vector<mpz_class>& w) {
    if (P[0] % DP_rarity == 0) {
        T.push_back(P);
        t.push_back(Pindex);
        std::set<mpz_class> T_set;
        for (const auto& tp : T) T_set.insert(tp[0]);
        for (const auto& match : W) {
            auto it = find(T.begin(), T.end(), match);
            if (it != T.end()) {
                int index = distance(T.begin(), it);
                mpz_class tP = t[index];
                mpz_class wP = w[distance(W.begin(), find(W.begin(), W.end(), match))];
                mpz_class dec = abs(tP - wP);
                auto end = chrono::system_clock::now();
                time_t end_time = chrono::system_clock::to_time_t(end);
                cout << "\n\033[01;33m[+]\033[32m PUZZLE SOLVED: \033[32m" << ctime(&end_time)
                     << "\r";
                cout << "\r\033[01;33m[+]\033[32m Private key (dec): \033[32m" << dec << "\033[0m"
                     << endl;
                ofstream file("KEYFOUNDKEYFOUND.txt", ios::app);
                file << "\n" << string(140, '-') << endl;
                file << "SOLVED " << ctime(&end_time);
                file << "Private Key (decimal): " << dec << endl;
                file << "Private Key (hex): " << dec.get_str(16) << endl;
                file << string(140, '-') << endl;
                file.close();
                return true;
            }
        }
    }
    return false;
}

vector<mpz_class> generate_powers_of_two(int hop_modulo) {
    vector<mpz_class> powers(hop_modulo);
    for (int pw = 0; pw < hop_modulo; ++pw) {
        powers[pw] = mpz_class(1) << pw;
    }
    return powers;
}

string search(const vector<Point>& P, const Point& W0, const mpz_class& DP_rarity, int Nw, int Nt, int hop_modulo,
              const mpz_class& upper_range_limit, const mpz_class& lower_range_limit, const vector<mpz_class>& powers_of_two) {
    vector<Point> T(Nt, Z), W(Nw, Z);
    vector<mpz_class> t(Nt), w(Nw);

    gmp_randclass rand(gmp_randinit_default);

    for (int k = 0; k < Nt; ++k) {
        t[k] = lower_range_limit + rand.get_z_range(upper_range_limit - lower_range_limit);
        T[k] = mul(t[k]);
    }

    for (int k = 0; k < Nw; ++k) {
        w[k] = rand.get_z_range(upper_range_limit - lower_range_limit);
        W[k] = add(W0, mul(w[k]));
    }

    long long Hops = 0, Hops_old = 0;
    auto t0 = chrono::high_resolution_clock::now();
    vector<mpz_class> memo(hop_modulo);

    for (int pw = 0; pw < hop_modulo; ++pw) {
        memo[pw] = powers_of_two[pw];
    }

    bool solved = false;
    while (!solved) {
        for (int k = 0; k < (Nt + Nw); ++k) {
            ++Hops;
            if (k < Nt) {
                mpz_class pw = T[k][0] % hop_modulo;
                solved = comparator(T[k], t[k], DP_rarity, T, t, W, w);
                if (solved) break;
                t[k] += memo[pw.get_ui()];
                T[k] = add(P[pw.get_ui()], T[k], modulo);
            } else {
                int n = k - Nw;
                mpz_class pw = W[n][0] % hop_modulo;
                solved = comparator(W[n], w[n], DP_rarity, W, w, T, t);
                if (solved) break;
                w[n] += memo[pw.get_ui()];
                W[n] = add(P[pw.get_ui()], W[n], modulo);
            }
        }

        auto t1 = chrono::high_resolution_clock::now();
        double elapsed_seconds = chrono::duration_cast<chrono::duration<double>>(t1 - t0).count();
        if (elapsed_seconds > 1.0) {
            double hops_per_second = static_cast<double>(Hops - Hops_old) / elapsed_seconds;
            auto elapsed_time = chrono::duration_cast<chrono::seconds>(t1 - starttime);
            int hours = chrono::duration_cast<chrono::hours>(elapsed_time).count();
            int minutes = chrono::duration_cast<chrono::minutes>(elapsed_time % chrono::hours(1)).count();
            int seconds = (elapsed_time % chrono::minutes(1)).count();
            stringstream elapsed_time_str;
            elapsed_time_str << setfill('0') << setw(2) << hours << ":"
                             << setfill('0') << setw(2) << minutes << ":"
                             << setfill('0') << setw(2) << seconds;
            double p_2 = log2(Hops);
            cout << "\r[+] [Hops: 2^" << fixed << setprecision(2) << p_2 << " <-> " << fixed << setprecision(0) 
                 << hops_per_second << " h/s] ["
                 << elapsed_time_str.str() << "]" << flush;
            t0 = t1;
            Hops_old = Hops;
        }
    }

    cout << "\r[+] Hops: " << Hops << endl;
    auto end = chrono::high_resolution_clock::now();
    double elapsed_seconds = chrono::duration_cast<chrono::duration<double>>(end - starttime).count();
    return "\r[+] Solution time: " + to_string(elapsed_seconds) + " sec";
}

int main()
{
    int puzzle = 64;
    string compressed_public_key =
    "03100611c54dfef604163b8358f7b7fac13ce478e02cb224ae16d45526b25d9d4d";
    int kangaroo_power = 6;
    mpz_class lower_range_limit = mpz_class(1) << (puzzle - 1);
    mpz_class upper_range_limit = (mpz_class(1) << puzzle) - 1;

    mpz_class DP_rarity = mpz_class(1) << ((puzzle - 2 * kangaroo_power) / 2 - 2);
    int hop_modulo = ((puzzle - 1) / 2) + kangaroo_power;

    int Nt = 1 << kangaroo_power;
    int Nw = 1 << kangaroo_power;

    vector<mpz_class> powers_of_two = generate_powers_of_two(hop_modulo);

    mpz_class X, Y;
    if (compressed_public_key.length() == 66) {
        X = mpz_class(compressed_public_key.substr(2), 16);
        Y = X2Y(X, stoi(compressed_public_key.substr(0, 2)) - 2);
    } else {
        cout << "[error] pubkey len(66/130) invalid!" << endl;
        return 1;
    }

    Point W0 = {X, Y};
    auto starttime = chrono::high_resolution_clock::now();
    time_t currentTime = time(nullptr);
    cout << "\r\033[01;33m[+]\033[32m KANGAROO: \033[01;33m" << ctime(&currentTime) << "\033[0m" << "\r";
    cout << "[+] [Puzzle]: " << puzzle << endl;
    cout << "[+] [Lower range limit]: " << lower_range_limit << endl;
    cout << "[+] [Upper range limit]: " << upper_range_limit << endl;
    cout << "[+] [EC Point Coordinate X]: " << X << endl;
    cout << "[+] [EC Point Coordinate Y]: " << Y << endl;
    mpz_class expected_hops = 2.2 * sqrt(mpz_class(1) << (puzzle - 1));
    double log_expected_hops = log2(expected_hops.get_d());
    cout << "[+] [Expected Hops: 2^" << fixed << setprecision(2) 
            << log_expected_hops << " (" << expected_hops << ")]" << endl;

    vector<Point> P = {PG};
    P.reserve(256);
    for (int k = 0; k < 255; ++k) {
        P.push_back(add(P[k], P[k]));
    }

    unsigned long seed = static_cast<unsigned long>(time(nullptr));
    gmp_randclass rand(gmp_randinit_default);
    rand.seed(seed);

    search(P, W0, DP_rarity, Nw, Nt, hop_modulo, upper_range_limit, lower_range_limit, powers_of_two);

    cout << "\r[+] Average time to solve: " << chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - starttime).count() << " sec" << endl;
    return 0;
}

