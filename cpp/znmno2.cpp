// Zn/MnO2 cell model: Zn anode | probe region | separator | porous MnO2 cathode.
//
// A self-contained C++17 program: edit the input file, run, and it writes the output table. It reads
// the same namelist input as the Fortran program and the Python package, and contains its own copy of
// the BAND block-tridiagonal solver (adapted from bandsolver v0.1.2,
// https://github.com/jcbernard87/bandsolver, BSD-3-Clause), so it needs no library.
//
//   build:  c++ -std=c++17 -O2 -ffp-contract=off znmno2.cpp -o znmno2_cpp
//   usage:  znmno2_cpp [input.nml]          (default input file: znmno2.nml)
//
// &run mode = 'corrected' (default), 'faithful_charge' or 'faithful_phcell'; data_dir = the folder of
// the spline tables. Equations: docs/model.md. Parameters: docs/parameters.md.
//
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ============================== BAND solver ==============================
// Solves A[j] dc[j-1] + B[j] dc[j] + D[j] dc[j+1] = G[j], j = 0..nj-1 (row-major n x n blocks).
namespace band {

enum class Pivot { partial, legacy };

// Solve Bm * S = R in place (Bm n x n, R n x m, row-major). Returns false if singular.
bool solve_partial(int n, int m, double* Bm, double* R) {
    const double eps = std::numeric_limits<double>::epsilon();
    double bscale = 0;
    for (int i = 0; i < n * n; ++i) bscale = std::max(bscale, std::abs(Bm[i]));
    if (bscale == 0) return false;
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i)
            if (std::abs(Bm[i * n + k]) > std::abs(Bm[p * n + k])) p = i;
        if (std::abs(Bm[p * n + k]) <= n * eps * bscale) return false;
        if (p != k) {
            std::swap_ranges(Bm + p * n, Bm + p * n + n, Bm + k * n);
            std::swap_ranges(R + p * m, R + p * m + m, R + k * m);
        }
        double f = 1.0 / Bm[k * n + k];
        for (int c = k; c < n; ++c) Bm[k * n + c] *= f;
        for (int c = 0; c < m; ++c) R[k * m + c] *= f;
        for (int i = 0; i < n; ++i) {
            if (i == k) continue;
            f = Bm[i * n + k];
            if (f == 0) continue;
            for (int c = k; c < n; ++c) Bm[i * n + c] -= f * Bm[k * n + c];
            for (int c = 0; c < m; ++c) R[i * m + c] -= f * R[k * m + c];
        }
    }
    return true;
}

// Pivot heuristic of the archival MATINV (Newman, Appendix C), with the operation order kept.
// exact: only an exactly zero pivot counts as singular (as the archival code), else |pivot| <= n eps max|B|.
bool solve_legacy(int n, int m, double* Bm, double* R, bool exact) {
    const double eps = std::numeric_limits<double>::epsilon();
    const double bmax0 = static_cast<double>(1.1f);
    double bscale = 0;
    for (int i = 0; i < n * n; ++i) bscale = std::max(bscale, std::abs(Bm[i]));
    const double tol = exact ? 0.0 : n * eps * bscale;
    std::vector<char> used(n, 0);
    int irow = 0, jcol = 0, jc = 0;
    for (int nn = 0; nn < n; ++nn) {
        double bmax = bmax0;
        bool found = false;
        for (int i = 0; i < n; ++i) {
            if (used[i]) continue;
            double bnext = 0, btry = 0;
            for (int j = 0; j < n; ++j) {
                if (used[j]) continue;
                const double a = std::abs(Bm[i * n + j]);
                if (a <= bnext) continue;
                bnext = a;
                if (bnext <= btry) continue;
                bnext = btry;
                btry = a;
                jc = j;
                found = true;
            }
            if (bnext >= bmax * btry) continue;
            bmax = bnext / btry;
            irow = i;
            jcol = jc;
        }
        if (!found) return false;
        used[jcol] = 1;
        if (jcol != irow) {
            std::swap_ranges(Bm + irow * n, Bm + irow * n + n, Bm + jcol * n);
            std::swap_ranges(R + irow * m, R + irow * m + m, R + jcol * m);
        }
        if (std::abs(Bm[jcol * n + jcol]) <= tol) return false;
        double f = 1.0 / Bm[jcol * n + jcol];
        for (int j = 0; j < n; ++j) Bm[jcol * n + j] *= f;
        for (int k = 0; k < m; ++k) R[jcol * m + k] *= f;
        for (int i = 0; i < n; ++i) {
            if (i == jcol) continue;
            f = Bm[i * n + jcol];
            for (int j = 0; j < n; ++j) Bm[i * n + j] -= f * Bm[jcol * n + j];
            for (int k = 0; k < m; ++k) R[i * m + k] -= f * R[jcol * m + k];
        }
    }
    return true;
}

// Block elimination and back substitution. Returns false on a singular or non-finite system.
bool solve(int n, int nj, const std::vector<double>& A, const std::vector<double>& B,
           const std::vector<double>& D, const std::vector<double>& G, std::vector<double>& dc, Pivot pivot,
           bool exact = false) {
    const std::size_t nn = static_cast<std::size_t>(n) * n;
    dc.assign(static_cast<std::size_t>(n) * nj, 0.0);
    for (const auto* v : {&A, &B, &D, &G})
        for (double x : *v)
            if (!std::isfinite(x)) return false;
    const int np1 = n + 1;
    std::vector<double> E(static_cast<std::size_t>(nj) * n * np1), Am(nn), Bm(nn), R(static_cast<std::size_t>(n) * np1);
    auto Ej = [&](int j) { return E.data() + static_cast<std::size_t>(j) * n * np1; };
    auto blk = [&](const std::vector<double>& M, int j) { return M.data() + static_cast<std::size_t>(j) * nn; };
    auto block_solve = [&](int m) {
        return pivot == Pivot::legacy ? solve_legacy(n, m, Bm.data(), R.data(), exact)
                                      : solve_partial(n, m, Bm.data(), R.data());
    };
    std::copy(blk(B, 0), blk(B, 0) + nn, Bm.begin());
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) R[i * np1 + k] = blk(D, 0)[i * n + k];
        R[i * np1 + n] = G[i];
    }
    if (!block_solve(np1)) return false;
    for (int k = 0; k < n; ++k) {
        Ej(0)[k * np1 + n] = R[k * np1 + n];
        for (int l = 0; l < n; ++l) Ej(0)[k * np1 + l] = -R[k * np1 + l];
    }
    for (int j = 1; j < nj; ++j) {
        std::copy(blk(A, j), blk(A, j) + nn, Am.begin());
        std::copy(blk(B, j), blk(B, j) + nn, Bm.begin());
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k) R[i * np1 + k] = blk(D, j)[i * n + k];
        const double* E1 = Ej(j - 1);
        for (int i = 0; i < n; ++i) {
            R[i * np1 + n] = -G[static_cast<std::size_t>(j) * n + i];
            for (int l = 0; l < n; ++l) {
                R[i * np1 + n] += Am[i * n + l] * E1[l * np1 + n];
                for (int k = 0; k < n; ++k) Bm[i * n + k] += Am[i * n + l] * E1[l * np1 + k];
            }
        }
        if (!block_solve(np1)) return false;
        for (int k = 0; k < n; ++k)
            for (int q = 0; q < np1; ++q) Ej(j)[k * np1 + q] = -R[k * np1 + q];
    }
    for (int k = 0; k < n; ++k) dc[static_cast<std::size_t>(nj - 1) * n + k] = Ej(nj - 1)[k * np1 + n];
    for (int j = nj - 2; j >= 0; --j) {
        const double* Ec = Ej(j);
        double* x = dc.data() + static_cast<std::size_t>(j) * n;
        const double* xn = x + n;
        for (int k = 0; k < n; ++k) {
            x[k] = Ec[k * np1 + n];
            for (int l = 0; l < n; ++l) x[k] += Ec[k * np1 + l] * xn[l];
        }
    }
    return true;
}

}  // namespace band

// ============================== input ==============================
std::string lower(std::string s) {
    for (auto& ch : s) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return s;
}

// Minimal Fortran-namelist reader: &group name = value, ... / ; '!' comments; d exponents.
// Returns {group: {name: raw value}} with lower-case names; quotes are removed from strings.
using Namelist = std::map<std::string, std::map<std::string, std::string>>;

Namelist read_namelist(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("input file not found: " + path);
    std::string text, line;
    while (std::getline(in, line)) {
        std::string out;
        char quote = 0;
        for (char ch : line) {
            if (quote) { if (ch == quote) quote = 0; }
            else if (ch == '\'' || ch == '"') quote = ch;
            else if (ch == '!') break;
            out += ch;
        }
        text += out + "\n";
    }
    Namelist nl;
    const std::regex group(R"(&(\w+)((?:'[^']*'|"[^"]*"|[^/'"])*)/)");
    const std::regex entry(R"((\w+)\s*=\s*('[^']*'|"[^"]*"|[^,\s/]+))");
    for (std::sregex_iterator g(text.begin(), text.end(), group), end; g != end; ++g) {
        const std::string gname = lower((*g)[1]);
        const std::string body = (*g)[2];
        if (nl.count(gname)) throw std::runtime_error("namelist group &" + gname + " appears twice");
        auto& entries = nl[gname];
        for (std::sregex_iterator e(body.begin(), body.end(), entry); e != end; ++e) {
            std::string val = (*e)[2];
            if (!val.empty() && (val.front() == '\'' || val.front() == '"')) val = val.substr(1, val.size() - 2);
            entries[lower((*e)[1])] = val;
        }
    }
    return nl;
}

double to_real(std::string v) {
    std::replace(v.begin(), v.end(), 'd', 'e');
    std::replace(v.begin(), v.end(), 'D', 'e');
    return std::stod(v);
}

bool to_bool(const std::string& v) {
    const std::string s = lower(v);
    return s == ".true." || s == "t" || s == ".t." || s == "true";
}

// ============================== gfortran list-directed output ==============================
// REAL(8): G25.17E3 (F editing for 0.1 <= |x| < 1e17, else 1P E25.17E3); INTEGER: I11; text as is;
// each item after one separator blank (docs/validation.md; python/znmno2_model/fortran_io.py).
std::string rjust(const std::string& s, std::size_t w) { return s.size() >= w ? s : std::string(w - s.size(), ' ') + s; }

std::string f_real(double x) {
    const int W = 25, Dd = 17;
    if (std::isnan(x)) return rjust("NaN", W);
    if (std::isinf(x)) return rjust(x > 0 ? "Infinity" : "-Infinity", W);
    if (x == 0.0) {
        std::string body = std::string(std::signbit(x) ? "-" : "") + "0." + std::string(Dd - 1, '0');
        return rjust(body, W - 5) + "     ";
    }
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.16e", std::abs(x));
    std::string s(buf);
    const auto epos = s.find('e');
    const int e = std::stoi(s.substr(epos + 1));
    std::string digits = s.substr(0, 1) + s.substr(2, epos - 2);
    const std::string sign = x < 0 ? "-" : "";
    const int k = e + 1;
    if (k >= 0 && k <= Dd) {
        std::string body = k == 0 ? "0." + digits : digits.substr(0, k) + "." + digits.substr(k);
        return rjust(sign + body, W - 5) + "     ";
    }
    char ex[16];
    std::snprintf(ex, sizeof ex, "E%c%03d", e < 0 ? '-' : '+', std::abs(e));
    return rjust(sign + digits.substr(0, 1) + "." + digits.substr(1) + ex, W);
}

struct Rec {
    std::string s;
    Rec& r(double x) { s += " " + f_real(x); return *this; }
    Rec& i(int v) { char b[32]; std::snprintf(b, sizeof b, "%11d", v); s += " " + std::string(b); return *this; }
    Rec& t(const std::string& v) { s += " " + v; return *this; }
};

// ============================== spline tables ==============================
std::vector<float> read_floats(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("data file not found: " + path);
    std::vector<float> v;
    std::string line;
    while (std::getline(in, line)) {
        auto p = line.find_first_not_of(" \t");
        if (p == std::string::npos || line[p] == '#') continue;
        std::istringstream ss(line);
        std::string w;
        while (ss >> w) v.push_back(std::strtof(w.c_str(), nullptr));
    }
    return v;
}

struct Tables {
    std::vector<float> zn_pts, mn_pts, coef, ocp_t, ocp_c;      // coef: 149 x 149 row-major; ocp_c: 50 x 4
    void load(const std::string& dir) {
        auto ph = read_floats(dir + "/ph_spline_phreeqc.txt");
        const int nz = static_cast<int>(ph[0]), nm = static_cast<int>(ph[1]), cz = static_cast<int>(ph[2]),
                  cm = static_cast<int>(ph[3]);
        zn_pts.assign(ph.begin() + 4, ph.begin() + 4 + nz);
        mn_pts.assign(ph.begin() + 4 + nz, ph.begin() + 4 + nz + nm);
        coef.assign(ph.begin() + 4 + nz + nm, ph.begin() + 4 + nz + nm + cz * cm);
        auto oc = read_floats(dir + "/r3_ocp_spline.txt");
        const int n = static_cast<int>(oc[0]);
        ocp_t.assign(oc.begin() + 1, oc.begin() + 1 + n);
        ocp_c.assign(oc.begin() + 1 + n, oc.begin() + 1 + n + (n - 1) * 4);
    }
};

// pow() exactly as the C library computes it. With a known base or exponent, compilers may rewrite
// pow(x, 0.5) as sqrt(x) or pow(10, y) as exp10(y) (clang does, without fast-math), which can change
// the last bit; gfortran calls pow, so both arguments are hidden from the optimizer.
double lpow(double x, double y) {
    volatile double xv = x, yv = y;
    return std::pow(static_cast<double>(xv), static_cast<double>(yv));
}

// ============================== faithful ports ==============================
// The charge line (ZnMn02_v3) and the pH-cell line (ZnMn02_v2.1_GITT_no_probe), operation by operation in
// the original's order and precision; mirrors fortran/faithful.f90. float = the original's REAL(4).
struct Faithful {
    bool ph = false;
    int N = 5, NJ = 122, S = 51;
    // template values (single-precision literals)
    float rxnk_2 = -8.5f, rxnk_3 = -7.5f, frac_zmcx = 0.4f, frac_zmcmax = 0.03f, xmax_t = 0.022f,
          applied_current = 0.000121f, porosity = 0.8f, volfrac_mno2 = 0.07f, stated_mass_loading = 0.00121f,
          zhs_ksp = 30.0f, rxnk_5 = -9.0f, fraction_kmno2 = 0.7f;
    Tables tab;
    double Rigc, Temp, Fconst, diff_Zn, diff_Mn, diff_SO4, diff_H, pH_init;
    double molar_mass_KMn8O16, density_KMn8O16, molar_mass_ZHS, density_ZHS, density_ZMC;
    double molar_mass_MnO2, density_MnO2, xmax_c, V_at_Zmin = 0, V_at_Zmax = 0, BL_thickness, sigma_sep;
    double Rxn2_K = 0, Rxn3_K = 0, Rxn5_K = 0, K_sp, Zmin, Zmax, ratio_initial, cbulk_Zn, cbulk_Mn, cbulk_H = 0;
    double sigma, xmax, Area_CS, len_sep, eps_sep, eps, volfrac_inert, AM_Grams, AM_Grams_sim;
    double KMn8O16_init, ZHS_init, ZMCx_init, ZMC_max_init, MnO2_init, Phi_1_init, molar_mass_ZMC_max;
    double resevoir_scaling, current_target;
    std::array<double, 7> z_ion{}, c_initial{}, diff_ion{};       // 1-based
    static constexpr int RAMP_ITERS = 1000;
    static constexpr double RAMP_DELT = 1.0e-9, DELT_NOMINAL = 1.0, WRITE_DENSITY = 10.0;
    // state, 1-based node index j (vectors of size NJ+1), unknown k (1..N)
    std::vector<double> c, delC, xx, delx, diff_term, mig_term;
    std::vector<double> KMn, ZHS, ZMCx, ZMCm, MnO2, por, tort, a_K, a_ZHS, a_ZMCx, a_ZMCm, a_MnO2, znm, mnm, ratio, MW;
    std::vector<double> pH_phreeqc, pH_ZHS, pH_standard;
    std::vector<double> Ab, Bb, Db, Gb;
    double time = 0, delT = DELT_NOMINAL, current = 0, c_density = 0, c_specific = 0, mAhg = 0, anode_pot = 0,
           ramp_initial = 0;
    int ramp_count = 0, last_write_time = 0;
    bool ramp_on = true, last_nan = false;
    char state = 'D';
    std::vector<std::string> rows;

    double& C(int k, int j) { return c[static_cast<std::size_t>(j) * 7 + k]; }
    double& DC(int k, int j) { return delC[static_cast<std::size_t>(j) * 7 + k]; }
    double& DT(int k, int j) { return diff_term[static_cast<std::size_t>(j) * 7 + k]; }
    double& MT(int k, int j) { return mig_term[static_cast<std::size_t>(j) * 7 + k]; }

    static float pow10_sp(float x) {
        if (x == std::trunc(x)) return static_cast<float>(lpow(10.0, static_cast<int>(std::lround(x))));
        return static_cast<float>(lpow(10.0, static_cast<double>(x)));
    }

    void constants() {
        Rigc = 8.314f; Temp = 298; Fconst = 96485;
        diff_Zn = 7.15e-6; diff_Mn = 6.88e-6; diff_SO4 = 1.07e-5; diff_H = 9.0e-5;
        pH_init = 5.5;
        molar_mass_KMn8O16 = 734.59f; density_KMn8O16 = 5.03f;
        molar_mass_ZHS = 549.819f; density_ZHS = 2.67f; density_ZMC = 5.0;
        molar_mass_MnO2 = 86.9368f; density_MnO2 = 5.03f;
        xmax_c = 200.0e-5;
        BL_thickness = 0.5 * 1.0e-4;
        sigma_sep = 1.0e-20;
        Area_CS = 0.178134094f;
        len_sep = 600.0 * 1.0e-4;
        eps_sep = 0.9f;
        z_ion.fill(0); c_initial.fill(0); diff_ion.fill(0);
        if (!ph) {
            N = 5; NJ = 122; S = 51;
            V_at_Zmin = 1.75; V_at_Zmax = 1.45f;
            Rxn2_K = pow10_sp(rxnk_2); Rxn3_K = pow10_sp(rxnk_3);
            K_sp = (zhs_ksp == std::trunc(zhs_ksp)) ? static_cast<float>(lpow(10.0, static_cast<int>(std::lround(zhs_ksp))))
                                                    : pow10_sp(zhs_ksp);
            Zmin = 0.2f; Zmax = 0.5;
            ratio_initial = Zmin * static_cast<double>(1.001f);
            cbulk_Zn = 0.002f; cbulk_Mn = 0.00005f;
            sigma = 0.1f;
            xmax = xmax_t;
            current_target = applied_current;
            eps = porosity;
            double vZMCx = static_cast<double>(volfrac_mno2 * frac_zmcx);
            double vZMCm = static_cast<double>(volfrac_mno2 * frac_zmcmax);
            double tmp = 0.0001f;
            volfrac_inert = 1.0 - eps - tmp - static_cast<double>(0.0001f) - vZMCx - vZMCm;
            AM_Grams = stated_mass_loading;
            AM_Grams_sim = vZMCx * Area_CS * xmax * density_ZMC;
            KMn8O16_init = tmp * density_KMn8O16;
            ZHS_init = static_cast<double>(0.0001f) * density_ZHS;
            ZMCx_init = vZMCx * density_ZMC;
            ZMC_max_init = vZMCm * density_ZMC;
            MnO2_init = 0;
            Phi_1_init = 1.7f;
            molar_mass_ZMC_max = static_cast<double>(65.38f) * Zmax + static_cast<double>(54.93f) + static_cast<double>(15.999f * 2.0f);
            z_ion[3] = 2.0; z_ion[4] = 0.0; z_ion[5] = -2.0;      // z_ion(4) = z_Mn0, zero (M-1)
            c_initial[3] = cbulk_Zn; c_initial[4] = cbulk_Mn;
            c_initial[5] = -1 * ((z_ion[3] * c_initial[3]) + (z_ion[4] * c_initial[4])) / z_ion[5];
            diff_ion[3] = diff_Zn; diff_ion[4] = diff_Mn; diff_ion[5] = diff_SO4;
            resevoir_scaling = (static_cast<double>(0.0001f) + (len_sep * Area_CS * eps_sep) + (xmax * Area_CS * eps)) /
                               (len_sep * Area_CS * eps_sep);
            len_sep = len_sep * resevoir_scaling;
        } else {
            N = 6; NJ = 134; S = 63;
            Rxn5_K = pow10_sp(rxnk_5);
            K_sp = 7.0e-26;
            Zmin = 0.35f; Zmax = 0.65f;
            ratio_initial = Zmax * static_cast<double>(0.9999f);
            cbulk_Zn = 0.002f; cbulk_Mn = 0.00005f; cbulk_H = 0.000543f;
            if ((-1.0 * std::log10(1000 * cbulk_H)) <= pH_init) pH_init = -1.0 * std::log10(1000 * cbulk_H);
            sigma = 0.01f;
            xmax = 0.0218f;
            current_target = 0.107f / 1000.0f;
            eps = 0.815f;
            double vMnO2 = static_cast<double>(0.0616f * fraction_kmno2);
            double tmp = 0.00001f;
            volfrac_inert = 1.0 - eps - static_cast<double>(0.00001f) - static_cast<double>(0.000001f) -
                            tmp * static_cast<double>(0.01f) - tmp - vMnO2;
            AM_Grams = 0.00108f;
            AM_Grams_sim = tmp * static_cast<double>(0.01f) * Area_CS * xmax * density_ZMC;
            KMn8O16_init = static_cast<double>(0.00001f) * density_KMn8O16;
            MnO2_init = vMnO2 * density_MnO2;
            ZHS_init = static_cast<double>(0.000001f) * density_ZHS;
            ZMCx_init = tmp * static_cast<double>(0.01f) * density_ZMC;
            ZMC_max_init = tmp * density_ZMC;
            Phi_1_init = 1.1f;
            molar_mass_ZMC_max = static_cast<double>(65.38f) * Zmax + static_cast<double>(54.93f) + static_cast<double>(15.999f * 2.0f);
            z_ion[3] = 2.0; z_ion[4] = 2.0; z_ion[5] = -2.0; z_ion[6] = 1.0;
            c_initial[3] = cbulk_Zn; c_initial[4] = cbulk_Mn; c_initial[6] = cbulk_H;
            c_initial[5] = -1 * ((z_ion[3] * c_initial[3]) + (z_ion[4] * c_initial[4]) + (z_ion[6] * c_initial[6])) / z_ion[5];
            diff_ion[3] = diff_Zn; diff_ion[4] = diff_Mn; diff_ion[5] = diff_SO4; diff_ion[6] = diff_H;
            resevoir_scaling = 1.0;
        }
    }

    void transport_terms(int j) {
        for (int ic = 3; ic <= N; ++ic) {
            if (j < S && !ph) {
                DT(ic, j) = (resevoir_scaling * resevoir_scaling) * por[j] * diff_ion[ic] / tort[j];
                MT(ic, j) = (resevoir_scaling * resevoir_scaling) * por[j] * z_ion[ic] * diff_ion[ic] * Fconst / (Rigc * Temp * tort[j]);
            } else {
                DT(ic, j) = por[j] * diff_ion[ic] / tort[j];
                MT(ic, j) = por[j] * z_ion[ic] * diff_ion[ic] * Fconst / (Rigc * Temp * tort[j]);
            }
        }
    }

    void initial_condition() {
        const std::size_t n1 = NJ + 1;
        c.assign(n1 * 7, 0); delC.assign(n1 * 7, 0); xx.assign(n1, 0); delx.assign(n1, 0);
        diff_term.assign(n1 * 7, 0); mig_term.assign(n1 * 7, 0);
        for (auto* v : {&KMn, &ZHS, &ZMCx, &ZMCm, &MnO2, &por, &tort, &a_K, &a_ZHS, &a_ZMCx, &a_ZMCm, &a_MnO2, &znm, &mnm,
                        &ratio, &MW})
            v->assign(n1, 0.0);
        pH_phreeqc.assign(n1, pH_init); pH_ZHS.assign(n1, pH_init); pH_standard.assign(n1, pH_init);
        const double h_sep = len_sep / static_cast<double>(static_cast<float>(S - 2));
        const double h_cath = xmax / static_cast<double>(static_cast<float>(NJ - S - 1));
        for (int j = 1; j <= NJ; ++j) {
            if (j == 1) xx[j] = 0.0;
            else if (j < S) xx[j] = h_sep * static_cast<float>(j - 1) - h_sep / 2.0;
            else if (j == S) xx[j] = len_sep;
            else if (j == NJ) xx[j] = xmax + len_sep;
            else xx[j] = len_sep + h_cath * static_cast<float>(j - S) - h_cath / 2.0;
        }
        for (int j = 2; j <= NJ - 1; ++j) {
            if (j < S) delx[j] = h_sep;
            else if (j > S) delx[j] = h_cath;
        }
        delx[1] = 0; delx[S] = 0; delx[NJ] = 0;
        for (int j = 1; j <= NJ; ++j) {
            C(1, j) = Phi_1_init;
            C(2, j) = 0;
            for (int ic = 3; ic <= N; ++ic) C(ic, j) = c_initial[ic];
            if (j < S) {
                por[j] = eps_sep;
                tort[j] = 2 * lpow(por[j], -0.5);
            } else {
                KMn[j] = KMn8O16_init; MnO2[j] = MnO2_init; ZHS[j] = ZHS_init; ZMCx[j] = ZMCx_init; ZMCm[j] = ZMC_max_init;
                por[j] = eps;
                tort[j] = 2 * lpow(por[j], -0.5);
                a_K[j] = 3.0 * KMn[j] / (density_KMn8O16 * (xmax_c));
                if (ph) a_MnO2[j] = 3.0 * MnO2[j] / (density_MnO2 * (xmax_c));
                a_ZHS[j] = 3.0 * ZHS[j] / (density_ZHS * (xmax_c));
                a_ZMCx[j] = 3.0 * ZMCx[j] / (density_ZMC * (xmax_c));
                a_ZMCm[j] = 3.0 * ZMCm[j] / (density_ZMC * (xmax_c));
                ratio[j] = ratio_initial;
                MW[j] = static_cast<double>(65.38f) * ratio[j] + static_cast<double>(54.938f) + static_cast<double>(15.999f * 2.0f);
                znm[j] = ZMCx[j] * ratio[j] / MW[j];
                mnm[j] = ZMCx[j] / MW[j];
            }
            transport_terms(j);
        }
        state = current_target >= 0.0 ? 'D' : 'C';
    }

    // ---------------------------------------------------------------- pH
    float basis(int i, int k, float x, const std::vector<float>& t) const {   // 1-based knot index i
        if (i + k > static_cast<int>(t.size())) {
            if (std::isnan(x)) return std::numeric_limits<float>::quiet_NaN();
            throw std::runtime_error("pH spline knot read out of bounds: not reproducible");
        }
        if (k == 1) return (t[i - 1] <= x && x < t[i]) ? 1.0f : 0.0f;
        float a = (t[i + k - 2] != t[i - 1]) ? (x - t[i - 1]) / (t[i + k - 2] - t[i - 1]) : 0.0f;
        float b = (t[i + k - 1] != t[i]) ? (t[i + k - 1] - x) / (t[i + k - 1] - t[i]) : 0.0f;
        return a * basis(i, k - 1, x, t) + b * basis(i + 1, k - 1, x, t);
    }

    static int find_interval(float x, const std::vector<float>& knots) {
        int low = 1, high = static_cast<int>(knots.size());
        while (low < high) {
            int mid = (low + high) / 2;
            if (x < knots[mid - 1]) high = mid; else low = mid + 1;
        }
        return std::max(std::min(low - 1, static_cast<int>(knots.size()) - 2), 0);
    }

    float manual_spline(float conc1, float conc2) const {
        const float zn = conc2, mn = conc1;
        const float lg_zn = std::log10(zn), lg_mn = std::log10(mn);
        const int kx = 3, ky = 3;
        const int ix = find_interval(lg_zn, tab.zn_pts), iy = find_interval(lg_mn, tab.mn_pts);
        float bx[4], by[4];
        for (int i = 0; i <= kx; ++i)
            bx[i] = ((ix - kx + i) >= 1 && (ix - kx + i) <= 153 - kx) ? basis(ix - kx + i, kx + 1, lg_zn, tab.zn_pts) : 0.0f;
        for (int i = 0; i <= ky; ++i)
            by[i] = ((iy - ky + i) >= 1 && (iy - ky + i) <= 153 - ky) ? basis(iy - ky + i, ky + 1, lg_mn, tab.mn_pts) : 0.0f;
        float ph_ = 0.0f;
        for (int i = 1; i <= 4; ++i)
            for (int j = 1; j <= 4; ++j) {
                const int ci = ix - kx + i - 1, cj = iy - ky + j - 1;
                if (ci >= 1 && ci <= 149 && cj >= 1 && cj <= 149)
                    ph_ = ph_ + tab.coef[static_cast<std::size_t>(ci - 1) * 149 + (cj - 1)] * bx[i - 1] * by[j - 1];
            }
        return ph_;
    }

    float eval_pH(float zn, float mn, int j) const {
        const float zn_m = zn * 1000.0f, mn_m = mn * 1000.0f;
        const float c_so4 = zn_m + mn_m;
        const float z2 = zn_m * zn_m;
        const float ph_ksp = static_cast<float>(-std::log10(lpow(static_cast<double>(z2 * z2 * c_so4) / K_sp,
                                                                     static_cast<double>(1.0f / 6.0f))));
        float r = manual_spline(zn_m, mn_m);
        if (j > S && ph_ksp < r) r = ph_ksp;
        return r;
    }

    double zhs_pH(double c_zn, double c_so4) const {
        const double zn = c_zn * 1000.0, so4 = c_so4 * 1000.0;
        const double z2 = zn * zn;
        const double c_h = lpow((z2 * z2) * so4 / K_sp, ph ? (-1.0 / 6.0) : (1.0 / 6.0));
        return -std::log10(c_h);
    }

    double zn_anode_pot(double c_zn) const {
        const double zn_ref = 0.001f, n_e = 2.0;
        return static_cast<double>(-0.762f) + (Rigc * Temp / (n_e * Fconst)) * std::log(c_zn / zn_ref);
    }

    double r3_ocp_spline(double theta) const {
        int ind = 1;
        while (ind <= 51) {
            if (theta <= tab.ocp_t[ind - 1]) break;
            ++ind;
        }
        ind = ind - 1;
        if (std::isnan(theta)) return theta;
        if (ind < 1 || ind > 50) throw std::runtime_error("R3 OCP spline read out of bounds (M-17): not reproducible");
        const double p1 = tab.ocp_c[(ind - 1) * 4], p2 = tab.ocp_c[(ind - 1) * 4 + 1], p3 = tab.ocp_c[(ind - 1) * 4 + 2],
                     p4 = tab.ocp_c[(ind - 1) * 4 + 3];
        const double d = theta - tab.ocp_t[ind - 1];
        double v = p1 * (d * d * d) + p2 * (d * d) + p3 * d + p4;
        return ((V_at_Zmin - V_at_Zmax) * v) + V_at_Zmax;
    }

    // ---------------------------------------------------------------- reactions
    void reaction_2(double p1, double p2, double c3, double c4, double c5, int j, double* out) const {
        const double c03_ref = 0.002f, c04_ref = 0.0001f, c05_ref = 0.0021f, lim = 0.1f;
        std::fill(out, out + 14, 0.0);
        const double xreact = Zmax;
        const double pH_func = zhs_pH(c3, c5);
        const double alpha_a = 0.5, alpha_c = 1.0 - alpha_a, rk = Rxn2_K;
        const double n_e = 2.0 * (1.0 - xreact);
        const double OCP = static_cast<double>(1.78f + 0.76f) + ((Rigc * Temp / (2 * Fconst)) * (1.0 * (-1) * std::log(c3 / c03_ref))) +
                           ((Rigc * Temp / (2 * Fconst)) * (1.0 * (-2.0) * std::log(c4 / c04_ref))) -
                           static_cast<double>(0.0592f * (8.0f / 2.0f)) * pH_func;
        const double exi = Fconst * rk * lpow(c3 / c03_ref, 1.0 * (static_cast<double>(8.0f / 3) - xreact) * alpha_a / n_e) *
                           lpow(c4 / c04_ref, 1.0 * (-1.0) * alpha_c / n_e) *
                           lpow(c5 / c05_ref, 1.0 * static_cast<double>(2.0f / 3.0f) * alpha_a / n_e);
        const double eta = p1 - p2 - OCP;
        const double area = (state == 'C') ? a_ZHS[j] : a_ZMCm[j];
        double irxn = area * exi * (std::exp(alpha_a * Fconst * eta / (Rigc * Temp)) - std::exp(-alpha_c * Fconst * eta / (Rigc * Temp)));
        const double M = molar_mass_ZMC_max, dt = delT;
        double dZn = 0, dMn = 0, dH = 0, dZMC = 0, dSO4 = 0, dZHS = 0;
        auto rates = [&](bool with_h) {
            dZn = (static_cast<double>(8.0f / 3.0f) - xreact) * irxn / (n_e * Fconst);
            dMn = -irxn / (n_e * Fconst);
            if (with_h) dH = 4 * irxn / (n_e * Fconst);
            dZMC = irxn / (n_e * Fconst);
            dSO4 = static_cast<double>(2.0f / 3.0f) * irxn / (n_e * Fconst);
            dZHS = static_cast<double>(-2.0f / 3.0f) * irxn / (n_e * Fconst);
        };
        rates(true);
        if ((-1.0 * dMn / area) > (c4 * diff_Mn / BL_thickness)) { irxn = ((c4 * diff_Mn / BL_thickness) / ((-1.0 * dMn) / area)) * irxn; rates(true); }
        if ((-1.0 * dZn / area) > (c3 * diff_Zn / BL_thickness)) { irxn = ((c3 * diff_Zn / BL_thickness) / ((-1.0 * dZn) / area)) * irxn; rates(true); }
        if ((-1.0 * dSO4 / area) > (c5 * diff_SO4 / BL_thickness)) { irxn = ((c5 * diff_SO4 / BL_thickness) / ((-1.0 * dSO4) / area)) * irxn; rates(true); }
        if ((-1.0 * (dZMC * M * dt) / ZMCm[j]) >= lim) { irxn = irxn * std::abs(ZMCm[j] / (dZMC * M * dt)) * lim; rates(false); }
        if ((-1.0 * (dZHS * M * dt) / ZHS[j]) >= lim) { irxn = irxn * std::abs(ZHS[j] / (dZHS * M * dt)) * lim; rates(false); }
        if ((ZMCm[j] + (dZMC * M * dt)) <= 0.0) { irxn = irxn * std::abs(ZMCm[j] / (dZMC * M * dt)) * static_cast<double>(0.1f); rates(true); }
        if ((ZHS[j] + (dZHS * M * dt)) <= 0.0) { irxn = irxn * std::abs(ZHS[j] / (dZHS * M * dt)) * static_cast<double>(0.1f); rates(true); }
        out[1] = OCP; out[2] = eta; out[3] = exi; out[4] = irxn;
        out[5] = dZn; out[6] = dMn; out[7] = dSO4; out[8] = dH; out[10] = dZHS; out[11] = dZMC;
    }

    void reaction_3(double p1, double p2, double c3, double c4, double c5, int j, double* out) const {
        (void)c4;
        const double c03_ref = 0.002f;
        std::fill(out, out + 14, 0.0);
        (void)zhs_pH(c3, c5);                     // computed and unused in the original
        const double alpha_a = 0.5, alpha_c = 1.0 - alpha_a, rk = Rxn3_K, n_e = 2.0;
        const double theta = (ratio[j] - Zmin) / (Zmax - Zmin);
        float Vint;                               // implicitly REAL(4) in the original
        if (theta <= static_cast<double>(-0.1f)) Vint = 1.8f;
        else Vint = static_cast<float>(r3_ocp_spline(theta) - static_cast<double>(0.762f));
        const double nernst = ((Rigc * Temp / (n_e * Fconst)) * (1.0 * std::log(c3 / c03_ref)));
        const double OCP = nernst + Vint - anode_pot;
        double exi;
        if (ratio[j] < Zmax) exi = Fconst * rk * lpow(c3 / c03_ref, 0.5 * alpha_c) * lpow(ratio[j], alpha_a) * lpow(Zmax - ratio[j], alpha_c);
        else exi = 1.0e-15;
        const double eta = p1 - p2 - OCP;
        const double area = a_ZMCx[j];
        double irxn = area * exi * (std::exp(alpha_a * Fconst * eta / (Rigc * Temp)) - std::exp(-alpha_c * Fconst * eta / (Rigc * Temp)));
        double dZn = irxn / (n_e * Fconst);
        if ((-1.0 * dZn / area) > (c3 * diff_Zn / BL_thickness)) {
            irxn = ((c5 * diff_Zn / BL_thickness) / std::abs(dZn)) * irxn;       // c5: the original's typo (M-6)
            dZn = irxn / (n_e * Fconst);
        }
        if ((znm[j] + (-1.0 * dZn * delT)) <= 0.0) {
            irxn = irxn * std::abs(znm[j] / (-1.0 * dZn * delT)) * static_cast<double>(0.1f);
            dZn = irxn / (n_e * Fconst);
        }
        out[1] = OCP; out[2] = eta; out[3] = exi; out[4] = irxn; out[5] = dZn; out[12] = -1.0 * dZn;
    }

    void reaction_5(double p1, double p2, double c3, double c4, double c5, double c6, int j, double* out) const {
        const double c04_ref = 0.0001, lim = 0.1;
        std::fill(out, out + 14, 0.0);
        const double pH_std = -std::log10(1000.0 * c6);
        const double pH_zhs_r = zhs_pH(c3, c5);
        const double alpha_a = 0.5, alpha_c = 1.0 - alpha_a, Uref = 1.2225, Uanode = anode_pot, rk = Rxn5_K, n_e = 2.0;
        const double area = a_MnO2[j];
        double OCP, exi, eta, irxn, dMn = 0, dH = 0, dZn = 0, dSO4 = 0, dZHS = 0, dMnO2 = 0;
        auto ocp_irxn = [&](double pH_func) {
            OCP = (Uref - Uanode) + (Rigc * Temp / (n_e * Fconst)) * (-std::log(c4 / c04_ref)) - 0.0592 * (4.0 / n_e) * pH_func;
            exi = Fconst * rk * lpow(c4 / c04_ref, -alpha_c * 1.0 / n_e);
            eta = p1 - p2 - OCP;
            irxn = area * exi * (std::exp(alpha_a * Fconst * eta / (Rigc * Temp)) - std::exp(-alpha_c * Fconst * eta / (Rigc * Temp)));
        };
        if (pH_std < pH_zhs_r) {
            ocp_irxn(pH_std);
            dMn = -irxn / (n_e * Fconst); dH = 4.0 * irxn / (n_e * Fconst); dMnO2 = irxn / (n_e * Fconst);
            dZHS = 0.0; dZn = 0.0; dSO4 = 0.0;
            if ((-dMn / area) > (c4 * diff_Mn / BL_thickness)) {
                irxn = ((c4 * diff_Mn / BL_thickness) / ((-dMn) / area)) * irxn;
                dMn = -irxn / (n_e * Fconst); dH = 4.0 * irxn / (n_e * Fconst); dMnO2 = irxn / (n_e * Fconst);
            }
            if ((-dMnO2 * molar_mass_MnO2 * delT) / MnO2[j] >= lim) {
                irxn = irxn * std::abs(MnO2[j] / (dMnO2 * molar_mass_MnO2 * delT)) * lim;
                dMn = -irxn / (n_e * Fconst); dH = 4.0 * irxn / (n_e * Fconst); dMnO2 = irxn / (n_e * Fconst);
            }
            const double H_post = (c6 * por[j] * delx[j] * Area_CS) + (dH * delT);
            const double ZHS_H = lpow(10.0, -pH_zhs_r) / 1000.0;
            if ((H_post < ZHS_H) && (static_cast<double>(1.1f) * pH_std >= pH_zhs_r)) {
                const double H_needed = ZHS_H - H_post;
                const double H_rate = H_needed / delT;
                const double r_ZHS = H_rate / 6.0;
                dH = dH + H_rate;
                dZn = dZn - 4.0 * r_ZHS;
                dSO4 = dSO4 - 1.0 * r_ZHS;
                dZHS = dZHS + r_ZHS;
            }
        } else {
            ocp_irxn(pH_zhs_r);
            auto rates5 = [&]() {
                dMn = -irxn / (n_e * Fconst); dMnO2 = irxn / (n_e * Fconst); dH = 0.0;
                dZHS = -(4.0 * irxn / (n_e * Fconst)) / 6.0;
                dZn = (8.0 / 3.0) * irxn / (n_e * Fconst);
                dSO4 = (2.0 / 3.0) * irxn / (n_e * Fconst);
            };
            rates5();
            if ((-dZn / area) > (c3 * diff_Zn / BL_thickness)) { irxn = ((c3 * diff_Zn / BL_thickness) / ((-dZn) / area)) * irxn; rates5(); }
            if ((-dSO4 / area) > (c5 * diff_SO4 / BL_thickness)) { irxn = ((c5 * diff_SO4 / BL_thickness) / ((-dSO4) / area)) * irxn; rates5(); }
            if ((-dMnO2 * molar_mass_MnO2 * delT) / MnO2[j] >= lim) { irxn = irxn * std::abs(MnO2[j] / (dMnO2 * molar_mass_MnO2 * delT)) * lim; rates5(); }
            if ((-dZHS * molar_mass_MnO2 * delT) / ZHS[j] >= lim) { irxn = irxn * std::abs(ZHS[j] / (dZHS * molar_mass_MnO2 * delT)) * lim; rates5(); }
        }
        out[1] = OCP; out[2] = eta; out[3] = exi; out[4] = irxn;
        out[5] = dZn; out[6] = dMn; out[7] = dSO4; out[8] = dH; out[10] = dZHS; out[13] = dMnO2;
    }

    void react_tot(const double* v, int j, double* tot) const {   // v 1-based, tot[1..13]
        std::fill(tot, tot + 14, 0.0);
        double r2[14], r3[14], r5[14];
        if (!ph) {
            reaction_2(v[1], v[2], v[3], v[4], v[5], j, r2);
            reaction_3(v[1], v[2], v[3], v[4], v[5], j, r3);
            for (int n = 4; n <= 12; ++n) {
                double acc = 0.0;
                acc = 0.0 + acc;
                acc = r2[n] + acc;
                acc = r3[n] + acc;
                tot[n] = acc;
            }
        } else {
            reaction_5(v[1], v[2], v[3], v[4], v[5], v[6], j, r5);
            for (int n = 4; n <= 13; ++n) {
                double acc = 0.0;
                for (int k = 0; k < 4; ++k) acc = 0.0 + acc;
                acc = r5[n] + acc;
                for (int k = 0; k < 4; ++k) acc = 0.0 + acc;
                tot[n] = acc;
            }
        }
    }

    // d[m][q]: derivative of output q (1 = irxn, 2.. = species) with respect to unknown m (1-based)
    void drx_dc(const double* v, int j, double d[7][7]) const {
        double st[7], a[7], b[7], t1[14], t2[14];
        const int nout = N - 1;
        for (int m = 1; m <= N; ++m) st[m] = v[m] * static_cast<double>(0.001f);
        if (!ph) st[2] = v[1] * static_cast<double>(0.001f);       // the phi2 step uses phi1 (M-7)
        else for (int m = 1; m <= N; ++m) if (st[m] == 0.0) st[m] = 1.0e-6;
        for (int m = 1; m <= N; ++m) {
            for (int q = 1; q <= N; ++q) a[q] = v[q];
            a[m] = v[m] + st[m];
            if (m >= 3 && v[m] <= st[m]) {
                react_tot(a, j, t1);
                react_tot(v, j, t2);
                for (int q = 1; q <= nout; ++q) d[m][q] = (t1[q + 3] - t2[q + 3]) / (st[m]);
            } else {
                for (int q = 1; q <= N; ++q) b[q] = v[q];
                b[m] = v[m] - st[m];
                react_tot(a, j, t1);
                react_tot(b, j, t2);
                for (int q = 1; q <= nout; ++q) d[m][q] = (t1[q + 3] - t2[q + 3]) / (2.0 * st[m]);
            }
        }
    }

    double en_row(const double* v) const {
        if (ph) return -z_ion[3] * v[3] - z_ion[4] * v[4] - z_ion[5] * v[5] - z_ion[6] * v[6];
        return -z_ion[3] * v[3] - z_ion[4] * v[4] - z_ion[5] * v[5];
    }

    // ---------------------------------------------------------------- fillmat + ABDGXY
    void assemble() {
        const std::size_t nn = static_cast<std::size_t>(N) * N;
        Ab.assign(nn * NJ, 0); Bb.assign(nn * NJ, 0); Db.assign(nn * NJ, 0); Gb.assign(static_cast<std::size_t>(N) * NJ, 0);
        double alphaE = 0, alphaW = 0, betaE = 0, betaW = 0;
        last_nan = false;
        for (int j = 1; j <= NJ; ++j) {
            double dE[7][7] = {}, dW[7][7] = {}, fE[7][7] = {}, fW[7][7] = {}, rj[7][7] = {}, smG[7] = {};
            double cE[7] = {}, cW[7] = {}, dcdxE[7] = {}, dcdxW[7] = {}, v[7] = {};
            for (int k = 1; k <= N; ++k) v[k] = C(k, j);
            auto store = [&](bool hasA, bool hasD) {
                for (int i = 1; i <= N; ++i) {
                    for (int k = 1; k <= N; ++k) {
                        const std::size_t ix = static_cast<std::size_t>(j - 1) * nn + (i - 1) * N + (k - 1);
                        if (hasA) Ab[ix] = (1.0 - alphaW) * fW[i][k] - betaW * dW[i][k];
                        if (!hasA) Bb[ix] = rj[i][k] - (1.0 - alphaE) * fE[i][k] + betaE * dE[i][k];
                        else if (!hasD) Bb[ix] = rj[i][k] + betaW * dW[i][k] + alphaW * fW[i][k];
                        else Bb[ix] = rj[i][k] + betaW * dW[i][k] + alphaW * fW[i][k] - (1.0 - alphaE) * fE[i][k] + betaE * dE[i][k];
                        if (hasD) Db[ix] = -alphaE * fE[i][k] - betaE * dE[i][k];
                    }
                    Gb[static_cast<std::size_t>(j - 1) * N + (i - 1)] = smG[i];
                }
            };
            if (j == 1) {
                alphaE = delx[j] / (delx[j + 1] + delx[j]);
                betaE = 2.0 / (delx[j] + delx[j + 1]);
                for (int ic = 1; ic <= N; ++ic) {
                    cE[ic] = alphaE * C(ic, j + 1) + (1.0 - alphaE) * C(ic, j);
                    dcdxE[ic] = betaE * (C(ic, j + 1) - C(ic, j));
                }
                dE[1][1] = -1.0;
                smG[1] = -(dE[1][1] * dcdxE[1]);
                smG[2] = en_row(v);
                for (int ic = 3; ic <= N; ++ic) rj[2][ic] = z_ion[ic];
                dE[3][3] = -1.0 * DT(3, j);
                dE[3][2] = -1.0 * MT(3, j) * cE[3];
                fE[3][3] = -1.0 * MT(3, j) * dcdxE[2];
                smG[3] = -1.0 * c_density / (z_ion[3] * Fconst) + (dE[3][3] * dcdxE[3] + fE[3][3] * cE[3]);
                for (int ic = 4; ic <= N - 1; ++ic) {
                    dE[ic][ic] = -1.0 * DT(ic, j);
                    dE[ic][2] = -1.0 * MT(ic, j) * cE[ic];
                    fE[ic][ic] = -1.0 * MT(ic, j) * dcdxE[2];
                    smG[ic] = 0.0 + (dE[ic][ic] * dcdxE[ic] + fE[ic][ic] * cE[ic]);
                }
                smG[N] = 0.0 - C(2, j);
                rj[N][2] = 1.0;
                store(false, true);
                continue;
            }
            alphaW = delx[j - 1] / (delx[j - 1] + delx[j]);
            betaW = 2.0 / (delx[j - 1] + delx[j]);
            if (j < NJ) {
                alphaE = delx[j] / (delx[j + 1] + delx[j]);
                betaE = 2.0 / (delx[j] + delx[j + 1]);
            }
            for (int ic = 1; ic <= N; ++ic) {
                cW[ic] = alphaW * C(ic, j) + (1.0 - alphaW) * C(ic, j - 1);
                dcdxW[ic] = betaW * (C(ic, j) - C(ic, j - 1));
                if (j < NJ) {
                    cE[ic] = alphaE * C(ic, j + 1) + (1.0 - alphaE) * C(ic, j);
                    dcdxE[ic] = betaE * (C(ic, j + 1) - C(ic, j));
                }
            }
            if (j == S) {
                dE[1][1] = -1.0;
                smG[1] = -(dE[1][1] * dcdxE[1]);
                smG[2] = en_row(v);
                for (int ic = 3; ic <= N; ++ic) rj[2][ic] = z_ion[ic];
                for (int ic = 3; ic <= N; ++ic) {
                    dW[ic][ic] = -1 * DT(ic, j - 1);
                    dE[ic][ic] = -1 * DT(ic, j + 1);
                    fW[ic][ic] = -1.0 * MT(ic, j - 1) * dcdxW[2];
                    fE[ic][ic] = -1 * MT(ic, j + 1) * dcdxE[2];
                    dW[ic][2] = -1 * MT(ic, j - 1) * cW[ic];
                    dE[ic][2] = -1 * MT(ic, j + 1) * cE[ic];
                    smG[ic] = 0.0 - (fW[ic][ic] * cW[ic] + dW[ic][ic] * dcdxW[ic]) + (fE[ic][ic] * cE[ic] + dE[ic][ic] * dcdxE[ic]);
                }
            } else if (j == NJ) {
                dW[1][1] = -(1.0 - por[j]) * sigma;
                smG[1] = (c_density) - dW[1][1] * dcdxW[1];
                smG[2] = en_row(v);
                for (int ic = 3; ic <= N; ++ic) rj[2][ic] = z_ion[ic];
                for (int ic = 3; ic <= N; ++ic) {
                    dW[ic][ic] = -1.0 * DT(ic, j);
                    dW[ic][2] = -1.0 * MT(ic, j) * cW[ic];
                    fW[ic][ic] = -1.0 * MT(ic, j) * dcdxW[2];
                    smG[ic] = 0.0 - (dW[ic][ic] * dcdxW[ic] + fW[ic][ic] * cW[ic]);
                }
                store(true, false);
                for (int i = 1; i <= N; ++i) {
                    if (std::isnan(smG[i])) last_nan = true;
                    for (int k = 1; k <= N; ++k)
                        if (std::isnan(dE[i][k]) || std::isnan(dW[i][k]) || std::isnan(fE[i][k]) || std::isnan(fW[i][k]) ||
                            std::isnan(rj[i][k]))
                            last_nan = true;
                }
                continue;
            } else if (j < S) {
                const double p = por[j];
                dE[1][1] = -(1.0 - p) * sigma_sep;
                dW[1][1] = -(1.0 - p) * sigma_sep;
                smG[1] = 0.0 - (fW[1][1] * cW[1] + dW[1][1] * dcdxW[1]) + (fE[1][1] * cE[1] + dE[1][1] * dcdxE[1]);
                smG[2] = en_row(v);
                for (int ic = 3; ic <= N; ++ic) rj[2][ic] = z_ion[ic];
                for (int ic = 3; ic <= N; ++ic) {
                    dW[ic][ic] = -1 * DT(ic, j);
                    dE[ic][ic] = -1 * DT(ic, j);
                    fW[ic][ic] = -1 * MT(ic, j) * dcdxW[2];
                    fE[ic][ic] = -1 * MT(ic, j) * dcdxE[2];
                    dW[ic][2] = -1 * MT(ic, j) * cW[ic];
                    dE[ic][2] = -1 * MT(ic, j) * cE[ic];
                    smG[ic] = -(fW[ic][ic] * cW[ic] + dW[ic][ic] * dcdxW[ic]) + (fE[ic][ic] * cE[ic] + dE[ic][ic] * dcdxE[ic]);
                }
                for (int iq = 3; iq <= N; ++iq)
                    for (int iv = 1; iv <= N; ++iv) rj[iq][iv] = (iq == iv) ? -(p / delT) * delx[j] : 0.0;
            } else {
                double rxn[14], drx[7][7] = {};
                react_tot(v, j, rxn);
                drx_dc(v, j, drx);
                double pW, pE, mW[7], mE[7], dfW[7], dfE[7];
                if (ph) {
                    pW = alphaW * por[j] + (1.0 - alphaW) * por[j - 1];
                    pE = alphaE * por[j + 1] + (1.0 - alphaE) * por[j];
                    for (int k = 1; k <= N; ++k) {
                        mW[k] = alphaW * MT(k, j) + (1.0 - alphaW) * MT(k, j - 1);
                        mE[k] = alphaE * MT(k, j + 1) + (1.0 - alphaE) * MT(k, j);
                        dfW[k] = alphaW * DT(k, j) + (1.0 - alphaW) * DT(k, j - 1);
                        dfE[k] = alphaE * DT(k, j + 1) + (1.0 - alphaE) * DT(k, j);
                    }
                } else {
                    pW = por[j]; pE = por[j];
                    for (int k = 1; k <= N; ++k) { mW[k] = MT(k, j); mE[k] = MT(k, j); dfW[k] = DT(k, j); dfE[k] = DT(k, j); }
                }
                dE[1][1] = -(1.0 - pE) * sigma;
                dW[1][1] = -(1.0 - pW) * sigma;
                for (int ic = 1; ic <= N; ++ic) rj[1][ic] = -drx[ic][1] * delx[j];
                smG[1] = rxn[4] * delx[j] - (fW[1][1] * cW[1] + dW[1][1] * dcdxW[1]) + (fE[1][1] * cE[1] + dE[1][1] * dcdxE[1]);
                smG[2] = en_row(v);
                for (int ic = 3; ic <= N; ++ic) rj[2][ic] = z_ion[ic];
                for (int ic = 3; ic <= N; ++ic) {
                    dW[ic][ic] = -1 * dfW[ic];
                    dE[ic][ic] = -1 * dfE[ic];
                    fW[ic][ic] = -1 * mW[ic] * dcdxW[2];
                    fE[ic][ic] = -1 * mE[ic] * dcdxE[2];
                    dW[ic][2] = -1 * mW[ic] * cW[ic];
                    dE[ic][2] = -1 * mE[ic] * cE[ic];
                    smG[ic] = -(rxn[ic + 2]) * delx[j] - (fW[ic][ic] * cW[ic] + dW[ic][ic] * dcdxW[ic]) +
                              (fE[ic][ic] * cE[ic] + dE[ic][ic] * dcdxE[ic]);
                }
                const double p = por[j];
                for (int iq = 3; iq <= N; ++iq)
                    for (int iv = 1; iv <= N; ++iv)
                        rj[iq][iv] = (iq == iv) ? drx[iv][iq - 1] * delx[j] - (p / delT) * delx[j] : drx[iv][iq - 1] * delx[j];
            }
            store(true, true);
        }
    }

    void solve() {
        std::vector<double> dc;
        const bool ok = band::solve(N, NJ, Ab, Bb, Db, Gb, dc, band::Pivot::legacy, true);
        for (int j = 1; j <= NJ; ++j)
            for (int k = 1; k <= N; ++k)
                DC(k, j) = ok ? dc[static_cast<std::size_t>(j - 1) * N + (k - 1)] : std::numeric_limits<double>::quiet_NaN();
    }

    // ---------------------------------------------------------------- after the solve
    void update_band_variables() {
        for (int j = 1; j <= NJ; ++j) {
            for (int k = 1; k <= N; ++k) C(k, j) = C(k, j) + DC(k, j);
            if (ph)
                for (int k = 3; k <= N; ++k) C(k, j) = std::max(C(k, j), 1.0e-20);
        }
    }

    void copy_boundary(int to, int from) {
        for (auto* v : {&KMn, &ZMCx, &ZMCm, &ZHS, &MnO2, &znm, &mnm, &ratio, &MW, &por, &tort, &a_K, &a_ZMCx, &a_ZMCm, &a_ZHS,
                        &a_MnO2})
            (*v)[to] = (*v)[from];
    }

    bool update_other_variables() {             // returns false where the original stops
        const double cutoff_theta = 0.985f;
        for (int j = S + 1; j <= NJ - 1; ++j) {
            double v[7] = {}, r[14];
            for (int q = 1; q <= N; ++q) v[q] = C(q, j) - (DC(q, j) / 2);
            react_tot(v, j, r);
            KMn[j] = KMn[j] + (molar_mass_KMn8O16 * r[9] * delT);
            ZHS[j] = ZHS[j] + (molar_mass_ZHS * r[10] * delT);
            ZMCm[j] = ZMCm[j] + (molar_mass_ZMC_max * r[11] * delT);
            if (ph) MnO2[j] = MnO2[j] + (molar_mass_MnO2 * r[13] * delT);
            znm[j] = znm[j] + (r[12] * delT);
            ratio[j] = znm[j] / mnm[j];
            if (ratio[j] >= (cutoff_theta * Zmax)) {
                const double zn_t = znm[j], mn_t = mnm[j];
                const double mn1 = (zn_t - Zmax * mn_t) / (cutoff_theta * Zmax - Zmax);
                const double mn2 = mn_t - mn1;
                const double zn1 = cutoff_theta * Zmax * mn1;
                znm[j] = zn1; mnm[j] = mn1;
                ratio[j] = znm[j] / mnm[j];
                ZMCm[j] = ZMCm[j] + (mn2 * molar_mass_ZMC_max);
            }
            MW[j] = static_cast<double>(65.38f) * ratio[j] + static_cast<double>(54.93f) + static_cast<double>(15.999f * 2.0f);
            ZMCx[j] = MW[j] * mnm[j];
            if (KMn[j] < 0.0) return false;
            if (ZHS[j] < 0.0) ZHS[j] = 1.0e-50;
            if (ZMCx[j] < 0.0) ZMCx[j] = 1.0e-50;
            if (ZMCm[j] < 0.0) ZMCm[j] = 1.0e-50;
            if (ph && MnO2[j] < 0.0) MnO2[j] = 1.0e-50;
            if (znm[j] < 0.0) znm[j] = 1.7e-50;
            if (mnm[j] < 0.0) mnm[j] = 1.0e-50;
            a_K[j] = 3.0 * KMn[j] / (density_KMn8O16 * (xmax_c));
            a_ZMCx[j] = 3.0 * ZMCx[j] / (density_ZMC * (xmax_c));
            a_ZMCm[j] = 3.0 * ZMCm[j] / (density_ZMC * (xmax_c));
            a_ZHS[j] = 3.0 * ZHS[j] / (density_ZHS * (xmax_c));
            if (ph) a_MnO2[j] = 3.0 * MnO2[j] / (density_MnO2 * (xmax_c));
            const double por_past = por[j];
            if (ph) {
                por[j] = 1.0 - volfrac_inert - (KMn[j] / density_KMn8O16) - (ZMCx[j] / density_ZMC) - (ZHS[j] / density_ZHS) -
                         (ZMCm[j] / density_ZMC) - (MnO2[j] / density_MnO2);
                tort[j] = 2.0 * lpow(por[j], -0.5);
                transport_terms(j);
            }
            for (int q = 3; q <= N; ++q) C(q, j) = C(q, j) * por_past / por[j];
            pH_phreeqc[j] = eval_pH(static_cast<float>(C(3, j)), static_cast<float>(C(4, j)), j);
            pH_ZHS[j] = zhs_pH(C(3, j), C(5, j));
            if (ph) pH_standard[j] = static_cast<float>(-1.0f * std::log10(static_cast<float>(1000 * C(6, j))));
            else pH_standard[j] = static_cast<float>(-1.0f * std::log10(static_cast<float>(1000 * C(1, j + 1))));   // M-14
        }
        copy_boundary(S, S + 1);
        pH_phreeqc[NJ] = pH_phreeqc[S + 1]; pH_ZHS[NJ] = pH_ZHS[S + 1]; pH_standard[NJ] = pH_standard[S + 1];   // M-16
        copy_boundary(NJ, NJ - 1);
        pH_phreeqc[NJ] = pH_phreeqc[NJ - 1]; pH_ZHS[NJ] = pH_ZHS[NJ - 1]; pH_standard[NJ] = pH_standard[NJ - 1];
        anode_pot = zn_anode_pot(C(3, 1));
        return true;
    }

    // ---------------------------------------------------------------- output
    std::string header() const {
        if (ph)
            return " time Voltage Current mAh/g State  C_Zn(avg) C_Mn(avg) C_SO4(avg) C_H(avg)  pH_phreeq(avg) "
                   "pH_ZHS(avg) pH_standard(avg)  KMn8O16(avg) ZMC(avg) ZMC_max(avg) ZHS(avg) MnO2(avg) "
                   "intercalation_Zn_to_Mn_ratio(avg)  R1_OCP(avg) R2_OCP(avg) R3_OCP(avg) R4_OCP(avg) R5_OCP(avg)  "
                   "R1_Eta(avg) R2_Eta(avg) R3_Eta(avg) R4_Eta(avg) R5_Eta(avg)  R1_Exi(avg) R2_Exi(avg) R3_Exi(avg) "
                   "R4_Exi(avg) R5_Exi(avg)  R1_I(avg) R2_I(avg) R3_I(Avg) R4_I(avg) R5_I(avg)  Area_KMn8O16 Area_ZMC "
                   "Area_ZMC_max Area_ZHS Area_MnO2  mAh/g_sim dom_react";
        return " time Voltage Current mAh/g State C_Zn(avg) C_Mn(avg) C_SO4(avg) pH_phreeq(avg) pH_ZHS(avg) "
               "pH_standard(avg)  KMn8O16(avg) ZMC(avg) ZMC_max(avg) ZHS(avg) intercalation_Zn_to_Mn_ratio(avg)  "
               "R1_OCP(avg) R2_OCP(avg) R3_OCP(avg)  R1_Eta(avg) R2_Eta(avg) R3_Eta(avg)  R1_Exi(avg) R2_Exi(avg) "
               "R3_Exi(avg)  R1_I(avg) R2_I(avg) R3_I(Avg)  Area_KMn8O16 Area_ZMC Area_ZMC_max Area_ZHS  "
               "mAh/g_sim dom_react";
    }

    void write_row() {
        const int nn = NJ - S - 1;
        std::vector<std::array<double, 5>> r2a(NJ + 1), r3a(NJ + 1);
        for (int j = S + 1; j <= NJ - 1; ++j) {
            double out[14];
            if (!ph) {
                reaction_2(C(1, j), C(2, j), C(3, j), C(4, j), C(5, j), j, out);
                for (int q = 1; q <= 4; ++q) r2a[j][q] = out[q];
                reaction_3(C(1, j), C(2, j), C(3, j), C(4, j), C(5, j), j, out);
                for (int q = 1; q <= 4; ++q) r3a[j][q] = out[q];
            } else {
                reaction_5(C(1, j), C(2, j), C(3, j), C(4, j), C(5, j), C(6, j), j, out);
                for (int q = 1; q <= 4; ++q) r2a[j][q] = out[q];
            }
        }
        auto mean_div = [&](auto get) { double acc = 0; for (int j = S + 1; j <= NJ - 1; ++j) acc = acc + get(j) / nn; return acc; };
        auto sum_div = [&](auto get) { double acc = 0; for (int j = S + 1; j <= NJ - 1; ++j) acc = acc + get(j); return acc / nn; };
        double zr = sum_div([](int) { return 0.0; }), r2[5], r3[5];
        for (int q = 1; q <= 4; ++q) {
            r2[q] = sum_div([&](int j) { return r2a[j][q]; });
            r3[q] = sum_div([&](int j) { return r3a[j][q]; });
        }
        Rec rec;
        rec.r(time).r(C(1, NJ)).r(current).r(mAhg).t(std::string(1, state));
        if (!ph) {
            const double temp[4] = {std::abs(zr), std::abs(r2[4]), std::abs(r3[4]), std::abs(zr)};
            int dom = 1;
            for (int q = 1; q < 4; ++q) if (temp[q] > temp[dom - 1]) dom = q + 1;
            rec.r(mean_div([&](int j) { return C(3, j); })).r(mean_div([&](int j) { return C(4, j); })).r(mean_div([&](int j) { return C(5, j); }));
            rec.r(mean_div([&](int j) { return pH_phreeqc[j]; })).r(mean_div([&](int j) { return pH_ZHS[j]; })).r(mean_div([&](int j) { return pH_standard[j]; }));
            rec.r(mean_div([&](int j) { return KMn[j]; })).r(mean_div([&](int j) { return ZMCx[j]; })).r(mean_div([&](int j) { return ZMCm[j]; }))
               .r(mean_div([&](int j) { return ZHS[j]; })).r(mean_div([&](int j) { return ratio[j]; }));
            for (int q = 1; q <= 4; ++q) rec.r(zr).r(r2[q]).r(r3[q]);
            rec.r(sum_div([&](int j) { return a_K[j]; })).r(sum_div([&](int j) { return a_ZMCx[j]; }))
               .r(sum_div([&](int j) { return a_ZMCm[j]; })).r(sum_div([&](int j) { return a_ZHS[j]; }));
            rec.r(mAhg * AM_Grams / AM_Grams_sim).i(dom);
        } else {
            for (int k = 3; k <= 6; ++k) rec.r(mean_div([&](int j) { return C(k, j); }));
            rec.r(mean_div([&](int j) { return pH_phreeqc[j]; })).r(mean_div([&](int j) { return pH_ZHS[j]; })).r(mean_div([&](int j) { return pH_standard[j]; }));
            rec.r(mean_div([&](int j) { return KMn[j]; })).r(mean_div([&](int j) { return ZMCx[j]; })).r(mean_div([&](int j) { return ZMCm[j]; }))
               .r(mean_div([&](int j) { return ZHS[j]; })).r(mean_div([&](int j) { return MnO2[j]; })).r(mean_div([&](int j) { return ratio[j]; }));
            for (int q = 1; q <= 4; ++q) rec.r(zr).r(zr).r(zr).r(zr).r(r2[q]);
            rec.r(sum_div([&](int j) { return a_K[j]; })).r(sum_div([&](int j) { return a_ZMCx[j]; }))
               .r(sum_div([&](int j) { return a_ZMCm[j]; })).r(sum_div([&](int j) { return a_ZHS[j]; }))
               .r(sum_div([&](int j) { return a_MnO2[j]; }));
            rec.r(mAhg * AM_Grams / AM_Grams_sim).i(1);
        }
        rows.push_back(rec.s);
    }

    void current_ramp() {
        if (ramp_on) {
            if (ramp_count == 0) ramp_initial = current;
            delT = RAMP_DELT;
            current = (static_cast<double>(ramp_count) * (current_target - ramp_initial) / static_cast<double>(RAMP_ITERS)) + ramp_initial;
            ramp_count = ramp_count + 1;
            if (ramp_count == RAMP_ITERS) {
                current = current_target;
                ramp_on = false;
                ramp_count = 0;
                delT = DELT_NOMINAL;
            }
        }
        c_density = current / Area_CS;
        c_specific = current / AM_Grams;
    }

    bool any_nan(const std::vector<double>& v) const {
        for (int j = 1; j <= NJ; ++j)
            for (int k = 1; k <= N; ++k)
                if (std::isnan(v[static_cast<std::size_t>(j) * 7 + k])) return true;
        return false;
    }

    std::string run() {
        rows.push_back(header());
        int it = 0;
        while (true) {
            ++it;
            if ((C(2, NJ) >= 99.0) && state == 'C') return "EXIT BECAUSE END OF CHARGE";
            if (std::isnan(DC(1, 1))) return "EXIT BECAUSE delC ISNAN";
            if (time >= 99.0 * 3600.0) return "EXIT BECAUSE END OF SIMULATION TIME";
            if (!ph) {
                if ((C(1, NJ) - C(2, NJ)) <= static_cast<double>(0.8f)) return "EXIT BECAUSE LOWER VOLTAGE CUTOFF";
                if (C(1, NJ) >= static_cast<double>(1.8f)) return "EXIT BECAUSE Upper VOLTAGE CUTOFF";
            } else {
                if (C(1, NJ) <= 1.0) return "EXIT BECAUSE LOWER VOLTAGE CUTOFF";
                if (C(1, NJ) >= 2.0) return "EXIT BECAUSE Upper VOLTAGE CUTOFF";
            }
            current_ramp();
            if (current > 1.0e-10) state = 'D';
            else if (current <= -1.0e-10) state = 'C';
            else state = 'R';
            if (ph) {
                if (mAhg >= 89.0) delT = static_cast<double>(0.1f);
                if (mAhg >= 100.0) delT = 1.0;
            }
            if (it >= RAMP_ITERS) {
                if ((time - WRITE_DENSITY) >= last_write_time) { write_row(); last_write_time = static_cast<int>(time); }
            } else if (it == 1) {
                write_row();
                last_write_time = static_cast<int>(time);
            }
            assemble();
            solve();
            update_band_variables();
            if (!update_other_variables()) return "KMn8O16 is negative";
            if (any_nan(c)) return "NaN in cprev";
            if (any_nan(delC)) return "NaN in delC";
            if (last_nan) return "NaN in coefficients";
            if ((time - WRITE_DENSITY) >= last_write_time) { write_row(); last_write_time = static_cast<int>(time); }
            time = time + delT;
            if (state == 'D') mAhg = mAhg + 1000.0 * c_specific * delT / 3600.0;
            else if (state == 'C') mAhg = mAhg - 1000.0 * c_specific * delT / 3600.0;
        }
    }
};

// ============================== corrected model ==============================
// A port of fortran/corrected.f90 (itself a port of python/znmno2_model), operation by operation in the
// same order, so the Fortran and C++ programs write identical output. Arrays are 1-based as in Fortran.
namespace corr {

struct Params {
    double A_cell = 0.178134094;              // cross-section of the cell [cm2], the same in every region (1-D)
    double L_probe = 0.1;                     // probe region length [cm]; 0 = no probe region
    double eps_probe = 1.0;                   // probe region: open electrolyte
    double tau_probe = 1.0;                   // probe region tortuosity
    double L_sep = 0.06;                      // separator thickness [cm]
    double eps_sep = 0.9;                     // separator porosity
    double tau_factor_sep = 2.0;              // tortuosity = tau_factor * eps**bruggeman
    double bruggeman_sep = -0.5;              // separator Bruggeman exponent
    double L_cath = 0.0218;                   // cathode thickness [cm]
    double eps_cath = 0.815;                  // initial cathode porosity
    double tau_factor_cath = 2.0;             // cathode tortuosity = tau_factor_cath * eps**bruggeman_cath
    double bruggeman_cath = -0.5;             // cathode Bruggeman exponent
    int n_probe = 20;                         // finite-volume cells per region
    int n_sep = 30;                           // separator cells
    int n_cath = 40;                          // cathode cells
    double sigma = 0.1;                       // cathode solid conductivity [S/cm]; effective sigma (1 - eps)
    double vf_MnO2 = 0.0001;                  // pristine MnO2 (R1)
    double vf_ZMO = 0.01;                     // Zn_z MnO2, dissolution/deposition (R2)
    double vf_host = 0.03;                    // Zn-insertion host (R3), fixed
    double vf_ZHS = 0.0;                      // zinc hydroxide sulfate
    double M_MnO2 = 86.937;                   // MnO2 molar mass [g/mol]
    double rho_MnO2 = 5.03;                   // MnO2 density [g/cm3]
    double rho_ZMO = 5.0;                     // Zn_z MnO2 density [g/cm3]
    double rho_host = 5.0;                    // insertion host density [g/cm3]
    double M_ZHS = 549.819;                   // Zn4SO4(OH)6 . 5 H2O
    double rho_ZHS = 2.67;                    // ZHS density [g/cm3]
    double r_MnO2 = 0.002;                    // MnO2 particle radius [cm]
    double r_ZMO = 0.002;                     // ZMO particle radius [cm]
    double r_host = 0.002;                    // host particle radius [cm]
    double r_ZHS = 0.002;                     // ZHS (and ZnO, Zn(OH)2) particle radius [cm]
    double z_ZMO = 0.5;                       // Zn per Mn in the dissolving phase
    double zmin = 0.2;                        // insertion range of the host (Zn per Mn)
    double zmax = 0.5;                        // Zn per Mn of the full host
    double theta0 = 0.001;                    // initial insertion fraction (0 = zmin, charged)
    double mass_AM = 0.0;                     // active mass [g] for mAh/g; 0 = MnO2-equivalent mass of the solids
    double c_ZnSO4 = 2.0;                     // initial ZnSO4 [mol/L]
    double c_MnSO4 = 0.05;                    // initial MnSO4 [mol/L]
    double c_H2SO4 = 0.0;                     // initial H2SO4 [mol/L]
    double D_Zn = 7.15e-06;                   // Zn2+ diffusion coefficient [cm2/s]
    double D_Mn = 6.88e-06;                   // Mn2+ diffusion coefficient [cm2/s]
    double D_SO4 = 1.07e-05;                  // SO4 2- diffusion coefficient [cm2/s]
    double D_H = 9e-05;                       // H+ diffusion coefficient [cm2/s]
    double U1 = 1.986;                        // R1 MnO2 + 4H+ + 2e -> Mn2+ + 2H2O
    double k1 = 1e-10;                        // [mol/cm2/s]
    double alpha1 = 0.5;                      // R1 transfer coefficient
    double U2 = 2.49;                         // R2 Zn_z MnO2 + 4H+ + (2-2z)e <-> z Zn2+ + Mn2+ + 2H2O
    double k2 = 1e-09;                        // R2 rate constant [mol/cm2/s]
    double alpha2 = 0.5;                      // R2 transfer coefficient
    double a_seed_R2 = 10.0;                  // deposition area besides ZMO and ZHS [cm2/cm3]
    double k3 = 1e-09;                        // R3 insertion
    double alpha3 = 0.5;                      // R3 transfer coefficient
    double V_at_zmin = 1.75;                  // empirical OCP spline end points [V]
    double V_at_zmax = 1.45;                  // R3 OCP spline value at z_max [V]
    double c_ref3 = 2.0;                      // reference Zn2+ concentration of the R3 OCP [mol/L]
    double k_an = 1e-06;                      // Zn anode, i0 = F k_an sqrt(c_Zn2+)
    double alpha_an = 0.5;                    // anode transfer coefficient
    double logK_ZHS = 28.4;                   // log10([Zn2+]^4 [SO4 2-] / [H+]^6) at saturation (Herrmann et al.)
    double k_ZHS = 1e-07;                     // precipitation/dissolution rate constant [mol/cm2/s]
    double a_seed_ZHS = 10.0;                 // precipitation area besides existing ZHS [cm2/cm3]
    std::string ph_mode = "speciation";       // pH in the reactions: speciation | zhs_equilibrium | fixed | spline
    double pH_fixed = 4.8;                    // ph_mode = fixed (R2Fixed of Bernard et al. 2025)
    std::string species = "with_H";           // with_H: H_T transported | no_H: H_T held at its initial value (charge line)
    std::string transport = "ions";           // ions: each total moves with its ion's D | quasi: species fluxes summed
    std::string basis = "free";               // Zn, Mn, SO4 in the Nernst and rate terms: free (speciation) | totals
    bool R1_on = true;                        // R1 (pristine MnO2 dissolution) on
    bool R2_on = true;                        // R2 (ZMO dissolution/deposition) on
    bool R3_on = true;                        // R3 (Zn insertion) on
    std::string zhs = "kinetic";              // kinetic | equilibrium (instantaneous) | lumped (into R1/R2) | off
    double zhs_nucleation = 1.0;              // Zn supersaturation c_Zn/c_sat needed for growth on the seed area (Herrmann: 1.05)
    bool ZnO_on = false;                      // extra precipitates (kinetic, same law as ZHS)
    bool ZnOH2_on = false;                    // Zn(OH)2 precipitation on
    double logK_ZnO = 11.17;                  // log10([Zn2+]/[H+]^2) at saturation (Herrmann & Horstmann 2024, Table 1)
    double logK_ZnOH2 = 12.45;                // log10([Zn2+]/[H+]^2) at Zn(OH)2 saturation
    double k_ZnO = 1e-07;                     // ZnO precipitation rate constant [mol/cm2/s]
    double k_ZnOH2 = 1e-07;                   // Zn(OH)2 precipitation rate constant [mol/cm2/s]
    double M_ZnO = 81.38;                     // ZnO molar mass [g/mol]
    double rho_ZnO = 5.61;                    // ZnO density [g/cm3]
    double M_ZnOH2 = 99.42;                   // Zn(OH)2 molar mass [g/mol]
    double rho_ZnOH2 = 3.05;                  // Zn(OH)2 density [g/cm3]
    std::string r3_ocp = "spline_nernst";     // spline_nernst: empirical OCP with Nernstian ends | spline | nernst
    double r3_end_width = 0.01;               // spline_nernst: the end terms act within about this fraction of theta = 0 and 1
    double U3_nernst = 1.55;                  // Herrmann et al. (2024), U_ins,Zn
    std::string logk_file = "";               // log K overrides ('name value' per line); empty: literature values
    std::string equilibria_db = "";           // PHREEQC database for the equilibria; empty: Herrmann et al. (2023) Table S1
    double D_OH = 5.27e-05;                   // species diffusion coefficients for transport = quasi [cm2/s]
    double D_HSO4 = 1.33e-05;                 // HSO4- diffusion coefficient (transport = quasi) [cm2/s]
    double D_complex = 5e-06;                 // every other complex (assumed; plan Q-1)
    double R = 8.314462618;                   // gas constant [J/mol/K]
    double T = 298.15;                        // temperature [K]
    double F = 96485.33212;                   // Faraday constant [C/mol]
    std::string steps = "cc I=100 Vmin=1.0";  // protocol (docs/protocol.md); I in mA/g, positive = discharge
    int cycles = 1;                           // number of times the protocol is repeated
    bool end_on_cutoff = true;                // a voltage cutoff ends the protocol (GITT); False: go to the next step
    double V_min = 1.0;                       // default lower cutoff [V]
    double V_max = 1.9;                       // default upper cutoff [V]
    double dt = 10.0;                         // [s]
    double write_interval = 60.0;             // [s]
    double newton_tol = 1e-09;                // Newton update tolerance (scaled)
    int newton_max_iter = 30;                 // Newton iterations per (sub-)step
};

// namelist name (lower case) -> field, and the group it belongs to
struct Fields {
    std::map<std::string, double*> reals;
    std::map<std::string, int*> ints;
    std::map<std::string, bool*> bools;
    std::map<std::string, std::string*> strs;
    std::map<std::string, std::string> group;
    explicit Fields(Params& p) {
        reals = {{"a_cell", &p.A_cell}, {"l_probe", &p.L_probe}, {"eps_probe", &p.eps_probe}, {"tau_probe", &p.tau_probe}, {"l_sep", &p.L_sep}, {"eps_sep", &p.eps_sep}, {"tau_factor_sep", &p.tau_factor_sep}, {"bruggeman_sep", &p.bruggeman_sep}, {"l_cath", &p.L_cath}, {"eps_cath", &p.eps_cath}, {"tau_factor_cath", &p.tau_factor_cath}, {"bruggeman_cath", &p.bruggeman_cath}, {"sigma", &p.sigma}, {"vf_mno2", &p.vf_MnO2}, {"vf_zmo", &p.vf_ZMO}, {"vf_host", &p.vf_host}, {"vf_zhs", &p.vf_ZHS}, {"m_mno2", &p.M_MnO2}, {"rho_mno2", &p.rho_MnO2}, {"rho_zmo", &p.rho_ZMO}, {"rho_host", &p.rho_host}, {"m_zhs", &p.M_ZHS}, {"rho_zhs", &p.rho_ZHS}, {"r_mno2", &p.r_MnO2}, {"r_zmo", &p.r_ZMO}, {"r_host", &p.r_host}, {"r_zhs", &p.r_ZHS}, {"z_zmo", &p.z_ZMO}, {"zmin", &p.zmin}, {"zmax", &p.zmax}, {"theta0", &p.theta0}, {"mass_am", &p.mass_AM}, {"c_znso4", &p.c_ZnSO4}, {"c_mnso4", &p.c_MnSO4}, {"c_h2so4", &p.c_H2SO4}, {"d_zn", &p.D_Zn}, {"d_mn", &p.D_Mn}, {"d_so4", &p.D_SO4}, {"d_h", &p.D_H}, {"u1", &p.U1}, {"k1", &p.k1}, {"alpha1", &p.alpha1}, {"u2", &p.U2}, {"k2", &p.k2}, {"alpha2", &p.alpha2}, {"a_seed_r2", &p.a_seed_R2}, {"k3", &p.k3}, {"alpha3", &p.alpha3}, {"v_at_zmin", &p.V_at_zmin}, {"v_at_zmax", &p.V_at_zmax}, {"c_ref3", &p.c_ref3}, {"k_an", &p.k_an}, {"alpha_an", &p.alpha_an}, {"logk_zhs", &p.logK_ZHS}, {"k_zhs", &p.k_ZHS}, {"a_seed_zhs", &p.a_seed_ZHS}, {"ph_fixed", &p.pH_fixed}, {"zhs_nucleation", &p.zhs_nucleation}, {"logk_zno", &p.logK_ZnO}, {"logk_znoh2", &p.logK_ZnOH2}, {"k_zno", &p.k_ZnO}, {"k_znoh2", &p.k_ZnOH2}, {"m_zno", &p.M_ZnO}, {"rho_zno", &p.rho_ZnO}, {"m_znoh2", &p.M_ZnOH2}, {"rho_znoh2", &p.rho_ZnOH2}, {"r3_end_width", &p.r3_end_width}, {"u3_nernst", &p.U3_nernst}, {"d_oh", &p.D_OH}, {"d_hso4", &p.D_HSO4}, {"d_complex", &p.D_complex}, {"r", &p.R}, {"t", &p.T}, {"f", &p.F}, {"v_min", &p.V_min}, {"v_max", &p.V_max}, {"dt", &p.dt}, {"write_interval", &p.write_interval}, {"newton_tol", &p.newton_tol}};
        ints = {{"n_probe", &p.n_probe}, {"n_sep", &p.n_sep}, {"n_cath", &p.n_cath}, {"cycles", &p.cycles}, {"newton_max_iter", &p.newton_max_iter}};
        bools = {{"r1_on", &p.R1_on}, {"r2_on", &p.R2_on}, {"r3_on", &p.R3_on}, {"zno_on", &p.ZnO_on}, {"znoh2_on", &p.ZnOH2_on}, {"end_on_cutoff", &p.end_on_cutoff}};
        strs = {{"ph_mode", &p.ph_mode}, {"species", &p.species}, {"transport", &p.transport}, {"basis", &p.basis}, {"zhs", &p.zhs}, {"r3_ocp", &p.r3_ocp}, {"logk_file", &p.logk_file}, {"equilibria_db", &p.equilibria_db}, {"steps", &p.steps}};
        group = {{"mode", "run"}, {"data_dir", "run"}, {"a_cell", "cell"}, {"l_probe", "cell"}, {"eps_probe", "cell"}, {"tau_probe", "cell"}, {"l_sep", "cell"}, {"eps_sep", "cell"}, {"tau_factor_sep", "cell"}, {"bruggeman_sep", "cell"}, {"l_cath", "cell"}, {"eps_cath", "cell"}, {"tau_factor_cath", "cell"}, {"bruggeman_cath", "cell"}, {"n_probe", "cell"}, {"n_sep", "cell"}, {"n_cath", "cell"}, {"sigma", "cell"}, {"vf_mno2", "solids"}, {"vf_zmo", "solids"}, {"vf_host", "solids"}, {"vf_zhs", "solids"}, {"m_mno2", "solids"}, {"rho_mno2", "solids"}, {"rho_zmo", "solids"}, {"rho_host", "solids"}, {"m_zhs", "solids"}, {"rho_zhs", "solids"}, {"r_mno2", "solids"}, {"r_zmo", "solids"}, {"r_host", "solids"}, {"r_zhs", "solids"}, {"z_zmo", "solids"}, {"zmin", "solids"}, {"zmax", "solids"}, {"theta0", "solids"}, {"mass_am", "solids"}, {"c_znso4", "electrolyte"}, {"c_mnso4", "electrolyte"}, {"c_h2so4", "electrolyte"}, {"d_zn", "electrolyte"}, {"d_mn", "electrolyte"}, {"d_so4", "electrolyte"}, {"d_h", "electrolyte"}, {"u1", "reactions"}, {"k1", "reactions"}, {"alpha1", "reactions"}, {"u2", "reactions"}, {"k2", "reactions"}, {"alpha2", "reactions"}, {"a_seed_r2", "reactions"}, {"k3", "reactions"}, {"alpha3", "reactions"}, {"v_at_zmin", "reactions"}, {"v_at_zmax", "reactions"}, {"c_ref3", "reactions"}, {"k_an", "reactions"}, {"alpha_an", "reactions"}, {"logk_zhs", "reactions"}, {"k_zhs", "reactions"}, {"a_seed_zhs", "reactions"}, {"ph_mode", "options"}, {"ph_fixed", "options"}, {"species", "options"}, {"transport", "options"}, {"basis", "options"}, {"r1_on", "options"}, {"r2_on", "options"}, {"r3_on", "options"}, {"zhs", "options"}, {"zhs_nucleation", "options"}, {"zno_on", "options"}, {"znoh2_on", "options"}, {"logk_zno", "options"}, {"logk_znoh2", "options"}, {"k_zno", "options"}, {"k_znoh2", "options"}, {"m_zno", "options"}, {"rho_zno", "options"}, {"m_znoh2", "options"}, {"rho_znoh2", "options"}, {"r3_ocp", "options"}, {"r3_end_width", "options"}, {"u3_nernst", "options"}, {"logk_file", "options"}, {"equilibria_db", "options"}, {"d_oh", "options"}, {"d_hso4", "options"}, {"d_complex", "options"}, {"r", "constants"}, {"t", "constants"}, {"f", "constants"}, {"steps", "protocol"}, {"cycles", "protocol"}, {"end_on_cutoff", "protocol"}, {"v_min", "protocol"}, {"v_max", "protocol"}, {"dt", "numerics"}, {"newton_tol", "numerics"}, {"newton_max_iter", "numerics"}, {"file", "output"}, {"write_interval", "output"}};
    }
};

constexpr double LN10 = 2.302585092994045684;     // log(10), correctly rounded (as Fortran's log(10.0_dp))
constexpr double LN2 = 0.6931471805599453094;
const double INF = std::numeric_limits<double>::infinity();

[[noreturn]] void die(const std::string& msg) { throw std::runtime_error(msg); }

// x(k, j): unknown k = 1..rows of cell j = 1..cols, stored column by column (as Fortran)
struct Mat {
    int rows = 0, cols = 0;
    std::vector<double> d;
    Mat() = default;
    Mat(int r, int c, double v = 0.0) : rows(r), cols(c), d(static_cast<std::size_t>(r) * c, v) {}
    double& operator()(int k, int j) { return d[static_cast<std::size_t>(j - 1) * rows + (k - 1)]; }
    double operator()(int k, int j) const { return d[static_cast<std::size_t>(j - 1) * rows + (k - 1)]; }
};

// A(i, k, j): row i, column k of the block of cell j, stored row-major per block (as the BAND solver reads it)
struct Blk {
    int n = 0, nj = 0;
    std::vector<double> d;
    Blk() = default;
    Blk(int n_, int nj_) : n(n_), nj(nj_), d(static_cast<std::size_t>(n_) * n_ * nj_, 0.0) {}
    double& operator()(int i, int k, int j) {
        return d[static_cast<std::size_t>(j - 1) * n * n + static_cast<std::size_t>(i - 1) * n + (k - 1)];
    }
    double operator()(int i, int k, int j) const {
        return d[static_cast<std::size_t>(j - 1) * n * n + static_cast<std::size_t>(i - 1) * n + (k - 1)];
    }
};

// The BAND solve of bandsolver's Fortran kernel with partial pivoting (its fast path: Gauss-Jordan by
// columns), including its zero second-neighbour terms, so results agree bit for bit (signed zeros too).
// G(j) rows [j*n + i]. Returns false on a singular block or non-finite data.
bool band_partial(int n, int nj, const std::vector<double>& A, const std::vector<double>& B,
                  const std::vector<double>& D, const std::vector<double>& G, std::vector<double>& dc) {
    const std::size_t nn = static_cast<std::size_t>(n) * n;
    dc.assign(static_cast<std::size_t>(n) * nj, 0.0);
    for (const auto* v : {&A, &B, &D, &G})
        for (double x : *v)
            if (!std::isfinite(x)) return false;
    const int np1 = n + 1, mR = 2 * n + 1;
    const double eps = std::numeric_limits<double>::epsilon();
    // E(k, q, j) and work blocks, column-major as in Fortran: M[(q)*n + k]
    std::vector<double> E(static_cast<std::size_t>(n) * np1 * nj), Bm(nn), Am(nn), R(static_cast<std::size_t>(n) * mR),
        Xp(nn, 0.0), Yw(nn, 0.0), Gm(n);
    auto e = [&](int k, int q, int j) -> double& { return E[(static_cast<std::size_t>(j) * np1 + q) * n + k]; };
    auto blk = [&](const std::vector<double>& M, int j, int i, int k) {
        return M[static_cast<std::size_t>(j) * nn + static_cast<std::size_t>(i) * n + k];
    };
    auto solve_cols = [&](int m) -> bool {          // Bm S = R(:, 0..m-1), solve_partial_columns
        double bscale = 0;
        for (int k = 0; k < n; ++k)
            for (int i = 0; i < n; ++i) bscale = std::max(bscale, std::abs(Bm[k * n + i]));
        if (bscale == 0) return false;
        const double tol = n * eps;
        std::vector<double> mult(n);
        for (int k = 0; k < n; ++k) {
            int p = k;
            for (int i = k + 1; i < n; ++i)
                if (std::abs(Bm[k * n + i]) > std::abs(Bm[k * n + p])) p = i;
            if (std::abs(Bm[k * n + p]) <= tol * bscale) return false;
            if (p != k) {
                for (int c = 0; c < n; ++c) std::swap(Bm[c * n + k], Bm[c * n + p]);
                for (int c = 0; c < m; ++c) std::swap(R[c * n + k], R[c * n + p]);
            }
            const double f = 1.0 / Bm[k * n + k];
            for (int c = k; c < n; ++c) Bm[c * n + k] = Bm[c * n + k] * f;
            for (int c = 0; c < m; ++c) R[c * n + k] = R[c * n + k] * f;
            for (int i = 0; i < n; ++i) mult[i] = Bm[k * n + i];
            mult[k] = 0;
            for (int c = k; c < n; ++c) {
                const double bkc = Bm[c * n + k];
                for (int i = 0; i < n; ++i) Bm[c * n + i] = Bm[c * n + i] - mult[i] * bkc;
            }
            for (int c = 0; c < m; ++c) {
                const double rkc = R[c * n + k];
                for (int i = 0; i < n; ++i) R[c * n + i] = R[c * n + i] - mult[i] * rkc;
            }
        }
        return true;
    };
    // node 1: B_1 [S_D | S_X | s_G] = [D_1 | X | G_1]
    for (int k = 0; k < n; ++k)
        for (int i = 0; i < n; ++i) {
            Bm[k * n + i] = blk(B, 0, i, k);
            R[k * n + i] = blk(D, 0, i, k);
            R[(n + k) * n + i] = Xp[k * n + i];
        }
    for (int i = 0; i < n; ++i) R[2 * n * n + i] = G[i];
    if (!solve_cols(mR)) return false;
    for (int i = 0; i < n; ++i) {
        e(i, n, 0) = R[2 * n * n + i];
        for (int k = 0; k < n; ++k) e(i, k, 0) = -R[k * n + i];
    }
    for (int k = 0; k < n; ++k)
        for (int i = 0; i < n; ++i) Xp[k * n + i] = -R[(n + k) * n + i];
    for (int j = 1; j < nj; ++j) {
        for (int k = 0; k < n; ++k)
            for (int i = 0; i < n; ++i) {
                Am[k * n + i] = blk(A, j, i, k);
                Bm[k * n + i] = blk(B, j, i, k);
                R[k * n + i] = blk(D, j, i, k);
            }
        for (int i = 0; i < n; ++i) Gm[i] = G[static_cast<std::size_t>(j) * n + i];
        if (j == 1) {                                // D_2 += A_2 X'
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l)
                    for (int i = 0; i < n; ++i) R[k * n + i] = R[k * n + i] + Am[l * n + i] * Xp[k * n + l];
        }
        if (j == nj - 1) {                           // Y dc_{nj-2} (Y = 0)
            for (int l = 0; l < n; ++l)
                for (int i = 0; i < n; ++i) Gm[i] = Gm[i] - Yw[l * n + i] * e(l, n, j - 2);
            for (int l = 0; l < n; ++l)
                for (int m = 0; m < n; ++m)
                    for (int i = 0; i < n; ++i) Am[l * n + i] = Am[l * n + i] + Yw[m * n + i] * e(m, l, j - 2);
            if (nj == 3)
                for (int k = 0; k < n; ++k)
                    for (int l = 0; l < n; ++l)
                        for (int i = 0; i < n; ++i) Bm[k * n + i] = Bm[k * n + i] + Yw[l * n + i] * Xp[k * n + l];
        }
        for (int i = 0; i < n; ++i) R[n * n + i] = -Gm[i];
        for (int l = 0; l < n; ++l)
            for (int i = 0; i < n; ++i) R[n * n + i] = R[n * n + i] + Am[l * n + i] * e(l, n, j - 1);
        for (int k = 0; k < n; ++k)
            for (int l = 0; l < n; ++l)
                for (int i = 0; i < n; ++i) Bm[k * n + i] = Bm[k * n + i] + Am[l * n + i] * e(l, k, j - 1);
        if (!solve_cols(np1)) return false;
        for (int q = 0; q < np1; ++q)
            for (int i = 0; i < n; ++i) e(i, q, j) = -R[q * n + i];
    }
    auto x = [&](int k, int j) -> double& { return dc[static_cast<std::size_t>(j) * n + k]; };
    for (int i = 0; i < n; ++i) x(i, nj - 1) = e(i, n, nj - 1);
    for (int j = nj - 2; j >= 0; --j) {
        for (int i = 0; i < n; ++i) x(i, j) = e(i, n, j);
        for (int l = 0; l < n; ++l)
            for (int i = 0; i < n; ++i) x(i, j) = x(i, j) + e(i, l, j) * x(l, j + 1);
    }
    for (int l = 0; l < n; ++l)
        for (int k = 0; k < n; ++k) x(k, 0) = x(k, 0) + Xp[l * n + k] * x(l, 2);
    for (double v : dc)
        if (!std::isfinite(v)) return false;
    return true;
}

// ------------------------------------------------------------------ equilibria (fortran/equilibria.f90)
struct Eq {
    int ns = 0;
    std::vector<std::string> names;
    std::vector<std::array<double, 5>> nu;          // nu[s][q], q = 1..4: H+, Zn+2, Mn+2, SO4-2
    std::vector<double> logk, z, abs_nu_h;

    static int charge_of(const std::string& name) {
        int i = static_cast<int>(name.size()) - 1;
        while (i >= 0 && std::isdigit(static_cast<unsigned char>(name[i]))) --i;
        if (i < 0 || (name[i] != '+' && name[i] != '-')) return 0;
        const int d = (i == static_cast<int>(name.size()) - 1) ? 1 : std::stoi(name.substr(i + 1));
        return name[i] == '-' ? -d : d;
    }
    void add(const std::string& name, std::array<double, 4> v, double lk) {
        int k = 0;
        while (k < ns && names[k] != name) ++k;
        if (k == ns) {
            ++ns;
            names.push_back(name); nu.push_back({}); logk.push_back(0); z.push_back(0); abs_nu_h.push_back(0);
        }
        names[k] = name;
        nu[k] = {0.0, v[0], v[1], v[2], v[3]};
        logk[k] = lk;
        z[k] = charge_of(name);
    }
    void init(const std::string& database, const std::string& logk_file) {
        add("H+", {1, 0, 0, 0}, 0.0); add("Zn+2", {0, 1, 0, 0}, 0.0);
        add("Mn+2", {0, 0, 1, 0}, 0.0); add("SO4-2", {0, 0, 0, 1}, 0.0);
        if (database.empty()) {
            add("OH-", {-1, 0, 0, 0}, -14.0); add("H2SO4", {2, 0, 0, 1}, 0.0); add("HSO4-", {1, 0, 0, 1}, 1.98);
            add("ZnOH+", {-1, 1, 0, 0}, -7.5); add("Zn(OH)2", {-2, 1, 0, 0}, -16.4);
            add("Zn(OH)3-", {-3, 1, 0, 0}, -28.2); add("Zn(OH)4-2", {-4, 1, 0, 0}, -41.3);
            add("Zn2OH+3", {-1, 2, 0, 0}, -9.0); add("Zn2(OH)6-2", {-6, 2, 0, 0}, -54.3);
            add("Zn4(OH)4+4", {-4, 4, 0, 0}, -27.0); add("ZnSO4", {0, 1, 0, 1}, 2.37);
            add("Zn(SO4)2-2", {0, 1, 0, 2}, 3.28); add("Zn(SO4)3-4", {0, 1, 0, 3}, 1.7);
            add("Zn(SO4)4-6", {0, 1, 0, 4}, 1.7); add("MnOH+", {-1, 0, 1, 0}, -10.59);
            add("Mn(OH)2", {-2, 0, 1, 0}, -18.54); add("Mn(OH)3-", {-3, 0, 1, 0}, -34.8);
            add("Mn(OH)4-2", {-4, 0, 1, 0}, -48.3); add("Mn2(OH)3+", {-3, 0, 2, 0}, -23.9);
            add("Mn2OH+3", {-1, 0, 2, 0}, -10.56); add("MnSO4", {0, 0, 1, 1}, 2.25);
        } else {
            read_phreeqc(database);
        }
        if (!logk_file.empty()) read_overrides(logk_file);
        for (int k = 0; k < ns; ++k) abs_nu_h[k] = std::abs(nu[k][1]);
    }
    void read_overrides(const std::string& path) {
        std::ifstream in(path);
        if (!in) die("logk_file not found");
        std::string line;
        while (std::getline(in, line)) {
            auto h = line.find('#');
            if (h != std::string::npos) line = line.substr(0, h);
            std::istringstream ss(line);
            std::string name, val;
            if (!(ss >> name)) continue;
            ss >> val;
            int k = 0;
            while (k < ns && names[k] != name) ++k;
            if (k == ns) die("log K override for unknown species '" + name + "'");
            std::replace(val.begin(), val.end(), 'd', 'e');
            std::replace(val.begin(), val.end(), 'D', 'e');
            logk[k] = std::stod(val);
        }
    }
    static std::string trim(const std::string& s) {
        const auto a = s.find_first_not_of(" \t\r");
        if (a == std::string::npos) return "";
        const auto b = s.find_last_not_of(" \t\r");
        return s.substr(a, b - a + 1);
    }
    static void parse_side(const std::string& s, std::vector<double>& cf, std::vector<std::string>& sp) {
        cf.clear(); sp.clear();
        std::string rest = s;
        while (true) {
            const auto p = rest.find(" + ");
            std::string tok = trim(p == std::string::npos ? rest : rest.substr(0, p));
            rest = p == std::string::npos ? "" : rest.substr(p + 3);
            if (!tok.empty()) {
                std::size_t i = 0;
                while (i < tok.size() && (std::isdigit(static_cast<unsigned char>(tok[i])) || tok[i] == '.')) ++i;
                cf.push_back(i > 0 ? std::stod(tok.substr(0, i)) : 1.0);
                sp.push_back(trim(tok.substr(i)));
            }
            if (trim(rest).empty()) break;
        }
    }
    void read_phreeqc(const std::string& path) {
        std::ifstream in(path);
        if (!in) die("equilibria_db not found");
        const double tk = 25.0 + 273.15;
        bool in_block = false, have = false, ok = false;
        double lk = 0;
        std::array<double, 4> v{};
        std::string cname, raw;
        auto add_nu = [&](const std::string& s, double c) {
            if (s == "H+") v[0] += c;
            else if (s == "Zn+2") v[1] += c;
            else if (s == "Mn+2") v[2] += c;
            else if (s == "SO4-2") v[3] += c;
            else if (s != "H2O") ok = false;
        };
        while (std::getline(in, raw)) {
            if (!in_block) {
                if (trim(raw) == "SOLUTION_SPECIES" && raw.find_first_not_of(" \t") == 0) in_block = true;
                continue;
            }
            if (trim(raw) == "PHASES" && raw.find_first_not_of(" \t") == 0) break;
            std::replace(raw.begin(), raw.end(), '\t', ' ');      // tabs separate fields as spaces do
            std::string line = raw.substr(0, raw.find('#'));
            if (trim(line).empty()) continue;
            std::string opt = trim(line);
            const auto eqp = line.find('=');
            if (eqp != std::string::npos && opt[0] != '-' && opt.compare(0, 5, "log_k") != 0) {
                if (have && ok) add(cname, v, lk);
                std::vector<double> cl, cr;
                std::vector<std::string> sl, sr;
                parse_side(line.substr(0, eqp), cl, sl);
                parse_side(line.substr(eqp + 1), cr, sr);
                cname = sr[0];
                v = {0, 0, 0, 0};
                ok = cr[0] == 1.0 && line.find("e-") == std::string::npos && cname != "H+" && cname != "Zn+2" &&
                     cname != "Mn+2" && cname != "SO4-2" && cname != "H2O";
                for (std::size_t i = 0; i < sl.size(); ++i) add_nu(sl[i], cl[i]);
                for (std::size_t i = 1; i < sr.size(); ++i) add_nu(sr[i], -cr[i]);
                lk = 0;
                have = true;
                continue;
            }
            if (!have) continue;
            if (opt[0] == '-') opt = opt.substr(1);
            opt = lower(opt);
            if (opt.compare(0, 6, "log_k ") == 0 || opt.compare(0, 5, "logk ") == 0) {
                std::istringstream ss(opt.substr(opt.find(' ')));
                ss >> lk;
            } else if (opt.compare(0, 8, "analytic") == 0) {
                double a[6] = {0, 0, 0, 0, 0, 0};
                std::istringstream ss(opt.substr(opt.find(' ')));
                for (double& q : a)
                    if (!(ss >> q)) break;
                lk = a[0] + a[1] * tk + a[2] / tk + a[3] * std::log10(tk) + a[4] / (tk * tk) + a[5] * (tk * tk);
            }
        }
        if (have && ok) add(cname, v, lk);
    }

    void species(const double* x, std::vector<double>& m) const {     // x[1..4]
        m.resize(ns);
        for (int s = 0; s < ns; ++s) {
            const double a = logk[s] + (((x[1] * nu[s][1] + x[2] * nu[s][2]) + x[3] * nu[s][3]) + x[4] * nu[s][4]);
            m[s] = lpow(10.0, std::min(std::max(a, -300.0), 300.0));
        }
    }
    void residual(const double* x, const double* T, double* F, std::vector<double>& m, double* bal, double& d) const {
        species(x, m);
        for (int k = 1; k <= 4; ++k) bal[k] = 0;
        d = 0;
        for (int s = 0; s < ns; ++s) {
            for (int k = 1; k <= 4; ++k) bal[k] = bal[k] + m[s] * nu[s][k];
            d = d + m[s] * abs_nu_h[s];
        }
        F[1] = (bal[1] - T[1]) / d;
        for (int k = 2; k <= 4; ++k) F[k] = std::log10(bal[k]) - std::log10(T[k]);
    }
    void jacobian(const std::vector<double>& m, const double* bal, double d, const double* T, double J[5][5]) const {
        double db[5][5] = {}, dd[5] = {};
        for (int s = 0; s < ns; ++s)
            for (int q = 1; q <= 4; ++q) {
                const double dm = m[s] * nu[s][q] * LN10;
                for (int k = 1; k <= 4; ++k) db[k][q] = db[k][q] + nu[s][k] * dm;
                dd[q] = dd[q] + abs_nu_h[s] * dm;
            }
        for (int q = 1; q <= 4; ++q) {
            J[1][q] = (db[1][q] * d - (bal[1] - T[1]) * dd[q]) / (d * d);
            for (int k = 2; k <= 4; ++k) J[k][q] = db[k][q] / (bal[k] * LN10);
        }
    }
    static void initial_guess(const double* T, double* x) {
        x[1] = std::log10(std::max(T[1], 0.0) + 1e-5);
        for (int k = 2; k <= 4; ++k) x[k] = std::log10(std::max(T[k], 1e-300)) - 0.5;
    }
    static double norm4(const double* v) { return std::sqrt(((v[1] * v[1] + v[2] * v[2]) + v[3] * v[3]) + v[4] * v[4]); }

    // Gaussian elimination with partial pivoting; A[i][j], B[i][j] 1-based (fortran solve_dense)
    static bool solve_dense(int n, int m, double A[5][5], double B[5][5]) {
        for (int k = 1; k <= n; ++k) {
            int p = k;
            for (int i = k + 1; i <= n; ++i)
                if (std::abs(A[i][k]) > std::abs(A[p][k])) p = i;
            if (A[p][k] == 0.0) return false;
            if (p != k) {
                for (int j = 1; j <= n; ++j) std::swap(A[k][j], A[p][j]);
                for (int j = 1; j <= m; ++j) std::swap(B[k][j], B[p][j]);
            }
            for (int i = k + 1; i <= n; ++i) {
                const double f = A[i][k] / A[k][k];
                for (int j = k + 1; j <= n; ++j) A[i][j] = A[i][j] - f * A[k][j];
                for (int j = 1; j <= m; ++j) B[i][j] = B[i][j] - f * B[k][j];
            }
        }
        for (int k = n; k >= 1; --k)
            for (int j = 1; j <= m; ++j) {
                for (int i = k + 1; i <= n; ++i) B[k][j] = B[k][j] - A[k][i] * B[i][j];
                B[k][j] = B[k][j] / A[k][k];
            }
        return true;
    }

    // log10 free concentrations x[1..4] for totals Tin[1..4]; x holds the warm start when warm
    bool solve(const double* Tin, double* x, bool warm) const {
        const double tol = 1e-13;
        double T[5] = {0, Tin[1], Tin[2], Tin[3], Tin[4]};
        for (int k = 2; k <= 4; ++k) T[k] = std::max(T[k], 1e-30);
        if (!warm) initial_guess(T, x);
        double F[5], F1[5], bal[5], bal1[5], d, d1, J[5][5], dx[5][5], xt[5];
        std::vector<double> m, m1;
        residual(x, T, F, m, bal, d);
        for (int it = 1; it <= 100; ++it) {
            double err = 0;
            for (int k = 1; k <= 4; ++k) err = std::max(err, std::abs(F[k]));
            if (err < tol) return true;
            jacobian(m, bal, d, T, J);
            for (int k = 1; k <= 4; ++k) dx[k][1] = -F[k];
            if (!solve_dense(4, 1, J, dx))
                for (int k = 1; k <= 4; ++k) dx[k][1] = 0;
            for (int k = 1; k <= 4; ++k)
                if (!(std::abs(dx[k][1]) <= std::numeric_limits<double>::max())) dx[k][1] = 0.0;
            double big = 0;
            for (int k = 1; k <= 4; ++k) big = std::max(big, std::abs(dx[k][1]));
            const double sc = std::min(1.0, 1.0 / std::max(big, 1e-300));
            for (int k = 1; k <= 4; ++k) dx[k][1] = dx[k][1] * sc;
            const double f0 = norm4(F);
            double lam = 1.0;
            for (int ib = 1; ib <= 40; ++ib) {
                for (int k = 1; k <= 4; ++k) xt[k] = x[k] + lam * dx[k][1];
                residual(xt, T, F1, m1, bal1, d1);
                if (norm4(F1) < f0) break;
                lam = 0.5 * lam;
            }
            for (int k = 1; k <= 4; ++k) x[k] = x[k] + lam * dx[k][1];
            for (int k = 1; k <= 4; ++k) { F[k] = F1[k]; bal[k] = bal1[k]; }
            m = m1; d = d1;
        }
        double err = 0;
        for (int k = 1; k <= 4; ++k) err = std::max(err, std::abs(F[k]));
        if (err < 1e-10) return true;
        bracketed(T, x);
        residual(x, T, F, m, bal, d);
        err = 0;
        for (int k = 1; k <= 4; ++k) err = std::max(err, std::abs(F[k]));
        return err < 1e-10;
    }
    void inner(double xh, const double* T, double* x) const {
        x[1] = xh;
        double F[5], bal[5], d, J[5][5], J3[5][5], dx[5][5];
        std::vector<double> m;
        for (int it = 1; it <= 200; ++it) {
            residual(x, T, F, m, bal, d);
            if (std::abs(F[2]) < 1e-13 && std::abs(F[3]) < 1e-13 && std::abs(F[4]) < 1e-13) break;
            jacobian(m, bal, d, T, J);
            for (int a = 1; a <= 3; ++a)
                for (int b = 1; b <= 3; ++b) J3[a][b] = J[a + 1][b + 1];
            for (int a = 1; a <= 3; ++a) dx[a][1] = -F[a + 1];
            if (!solve_dense(3, 1, J3, dx))
                for (int a = 1; a <= 3; ++a) dx[a][1] = 0;
            for (int a = 1; a <= 3; ++a)
                if (!(std::abs(dx[a][1]) <= std::numeric_limits<double>::max())) dx[a][1] = 0.0;
            double big = 0;
            for (int a = 1; a <= 3; ++a) big = std::max(big, std::abs(dx[a][1]));
            const double sc = std::min(1.0, 1.0 / std::max(big, 1e-300));
            for (int a = 1; a <= 3; ++a) x[a + 1] = x[a + 1] + dx[a][1] * sc;
        }
    }
    void bracketed(const double* T, double* x) const {
        double lo = -16.0, hi = 2.0;
        std::vector<double> m;
        initial_guess(T, x);
        for (int it = 1; it <= 64; ++it) {
            const double mid = 0.5 * (lo + hi);
            inner(mid, T, x);
            species(x, m);
            double g = 0;
            for (int s = 0; s < ns; ++s) g = g + m[s] * nu[s][1];
            g = g - T[1];
            if (g < 0) lo = mid; else hi = mid;
            if (hi - lo < 1e-12) break;
        }
        inner(0.5 * (lo + hi), T, x);
    }
    void sensitivity(const double* x, const double* Tin, double S[5][5]) const {
        double T[5] = {0, Tin[1], Tin[2], Tin[3], Tin[4]};
        for (int k = 2; k <= 4; ++k) T[k] = std::max(T[k], 1e-30);
        double F[5], bal[5], d, J[5][5];
        std::vector<double> m;
        residual(x, T, F, m, bal, d);
        jacobian(m, bal, d, T, J);
        for (int a = 1; a <= 4; ++a)
            for (int b = 1; b <= 4; ++b) S[a][b] = 0;
        S[1][1] = 1.0 / d;
        for (int k = 2; k <= 4; ++k) S[k][k] = 1.0 / (T[k] * LN10);
        solve_dense(4, 4, J, S);
    }
};

// ------------------------------------------------------------------ the model (fortran/corrected.f90)
constexpr int N = 12;
constexpr int P1 = 1, P2 = 2, ZN = 3, MN = 4, SO = 5, H = 6, MO = 7, ZM = 8, TH = 9, ZH = 10, ZO = 11, ZX = 12;
constexpr int SPECIES[4] = {ZN, MN, SO, H};
constexpr int SOLIDS[5] = {MO, ZM, ZH, ZO, ZX};
constexpr int PRECIP[3] = {ZH, ZO, ZX};
constexpr int TOT_COLS[5] = {0, H, ZN, MN, SO};     // 1-based: column of each speciation total
constexpr int PROBE = 0, SEP = 1, CATH = 2;
constexpr double CHARGE[N + 1] = {0, 0, 0, 2.0, 2.0, -2.0, 1.0, 0, 0, 0, 0, 0, 0};
constexpr double W_SMOOTH = 0.01, S_MAX = 200.0, S_STEP = 20.0, TH_FLOOR = 1e-8;
constexpr double RES_TOL = 1e-6, TRACE = 1.0e-9, FB_REF = 1.0e-6;
constexpr double TYP[N + 1] = {0, 1.0, 1.0, 1e-3, 1e-3, 1e-3, 1e-6, 1e-3, 1e-3, 1.0, 1e-3, 1e-3, 1e-3};
constexpr double STEP_FLOOR[N + 1] = {0, 1.0, 1.0, 1e-15, 1e-15, 1e-15, 1e-8, 1e-15, 1e-15, 1.0, 1e-15, 1e-15, 1e-15};
constexpr double EVENT_DV = 1.0e-4, MIN_SUBSTEP = 1.0e-10, EVENT_MIN_DT = 1.0e-12, CV_TOL = 1.0e-9, CV_ACCEPT = 1.0e-6,
                 LIMIT = 1.0e-3;
constexpr int MAX_FAILURES = 200;

struct Failure {};                                  // a step that cannot be solved (Python's SolverFailure)

double bernoulli(double a) {
    if (std::abs(a) < 1e-6) return 1.0 - a / 2.0 + a * a / 12.0;
    return a / std::expm1(a);
}
double dbernoulli(double a) {
    if (std::abs(a) < 1e-4) return -0.5 + a / 6.0;
    const double em = std::expm1(a);
    return (em - a * (em + 1.0)) / (em * em);
}
void theta_pair(double s, double& th, double& om) {
    const double e = std::exp(-std::abs(s));
    const double big = 1.0 / (1.0 + e), small = e / (1.0 + e);
    if (s >= 0) { th = big; om = small; } else { th = small; om = big; }
}
double sigmoid(double s) { double t, o; theta_pair(s, t, o); return t; }
double pos(double z, double width) { const double zp = std::max(z, 0.0); return zp * zp / (zp + width); }
double logaddexp(double x, double y) {
    if (x == y) return x + LN2;
    const double tmp = x - y;
    if (tmp > 0) return x + std::log1p(std::exp(-tmp));
    if (tmp <= 0) return y + std::log1p(std::exp(tmp));
    return tmp;
}
double onoff(bool f) { return f ? 1.0 : 0.0; }

struct Step {
    std::string kind;
    double I = 0, V = 0, t = -1, Vmin = 0, Vmax = 0, Imin = -1;    // t, Imin < 0: not given
};

struct Model {
    Params p;
    Eq eq;
    int nc = 0, c0 = 0;
    std::vector<double> xc, dxc, area_c, vol;       // 1-based
    std::vector<int> region;
    std::vector<char> incath;
    double fRT = 0, Vm[N + 1] = {}, eps_fixed = 0, n_host = 0, a_host = 0, Dk[N + 1] = {}, mass = 0;
    bool transported[N + 1] = {}, prec_on[N + 1] = {};
    double prec_logk[N + 1] = {}, prec_zn[N + 1] = {}, prec_s[N + 1] = {}, prec_h[N + 1] = {}, prec_k[N + 1] = {};
    double ocp_t[52] = {}, ocp_c[51][5] = {};
    std::vector<double> ph_zn, ph_mn, ph_coef;       // ph_coef(i, j) = [(i-1)*149 + (j-1)]
    std::vector<double> D_sp;
    std::vector<double> H0;
    Mat spec_x, cache_T, cache_lx, x0_state;
    bool have_H0 = false, spec_valid = false, cache_valid = false;
    double res_scale[N + 1] = {};
    int i_probe = 1;
    std::vector<std::string> out_rows;

    double M_ZMO() const { return 65.38 * p.z_ZMO + 54.938 + 2 * 15.999; }
    double M_host() const { return 65.38 * p.zmin + 54.938 + 2 * 15.999; }
    double active_mass() const {
        if (p.mass_AM > 0) return p.mass_AM;
        const double mol_mn = p.vf_MnO2 * p.rho_MnO2 / p.M_MnO2 + p.vf_ZMO * p.rho_ZMO / M_ZMO() + p.vf_host * p.rho_host / M_host();
        return mol_mn * p.M_MnO2 * p.A_cell * p.L_cath;
    }

    void setup(const std::string& data_dir) {
        int reg[3], nn[3], nparts = 0;
        double lens[3];
        if (p.L_probe > 0) { reg[nparts] = PROBE; lens[nparts] = p.L_probe; nn[nparts] = p.n_probe; ++nparts; }
        reg[nparts] = SEP; lens[nparts] = p.L_sep; nn[nparts] = p.n_sep; ++nparts;
        reg[nparts] = CATH; lens[nparts] = p.L_cath; nn[nparts] = p.n_cath; ++nparts;
        nc = 0;
        for (int k = 0; k < nparts; ++k) nc += nn[k];
        xc.assign(nc + 1, 0); dxc.assign(nc + 1, 0); area_c.assign(nc + 1, 0); vol.assign(nc + 1, 0);
        region.assign(nc + 1, 0); incath.assign(nc + 1, 0);
        int j = 0;
        double start = 0.0;
        for (int k = 0; k < nparts; ++k) {
            const double hh = lens[k] / nn[k];
            for (int q = 1; q <= nn[k]; ++q) {
                ++j;
                xc[j] = start + hh * (static_cast<double>(q - 1) + 0.5);
                dxc[j] = hh;
                region[j] = reg[k];
            }
            start = start + lens[k];
        }
        for (j = 1; j <= nc; ++j) {
            area_c[j] = p.A_cell;
            vol[j] = area_c[j] * dxc[j];
            incath[j] = region[j] == CATH;
        }
        c0 = 1;
        while (!incath[c0]) ++c0;
        i_probe = 1;
        if (p.L_probe > 0) i_probe = p.n_probe / 2 + 1;
        eq.init(p.equilibria_db, p.logk_file);
        fRT = p.F / (p.R * p.T);
        Vm[MO] = p.M_MnO2 / p.rho_MnO2; Vm[ZM] = M_ZMO() / p.rho_ZMO; Vm[ZH] = p.M_ZHS / p.rho_ZHS;
        Vm[ZO] = p.M_ZnO / p.rho_ZnO; Vm[ZX] = p.M_ZnOH2 / p.rho_ZnOH2;
        eps_fixed = 1.0 - p.eps_cath - p.vf_MnO2 - p.vf_ZMO - p.vf_ZHS;
        n_host = p.vf_host * p.rho_host / M_host();
        a_host = 3.0 * p.vf_host / p.r_host;
        Dk[ZN] = p.D_Zn; Dk[MN] = p.D_Mn; Dk[SO] = p.D_SO4; Dk[H] = p.D_H;
        prec_on[ZH] = p.zhs == "kinetic"; prec_on[ZO] = p.ZnO_on; prec_on[ZX] = p.ZnOH2_on;
        prec_logk[ZH] = p.logK_ZHS; prec_zn[ZH] = 4.0; prec_s[ZH] = 1.0; prec_h[ZH] = 6.0; prec_k[ZH] = p.k_ZHS;
        prec_logk[ZO] = p.logK_ZnO; prec_zn[ZO] = 1.0; prec_s[ZO] = 0.0; prec_h[ZO] = 2.0; prec_k[ZO] = p.k_ZnO;
        prec_logk[ZX] = p.logK_ZnOH2; prec_zn[ZX] = 1.0; prec_s[ZX] = 0.0; prec_h[ZX] = 2.0; prec_k[ZX] = p.k_ZnOH2;
        transported[ZN] = transported[MN] = transported[SO] = true;
        transported[H] = p.species != "no_H";
        if (p.r3_ocp != "nernst") {
            Tables t;
            t.load(data_dir);
            for (int i = 1; i <= 51; ++i) ocp_t[i] = t.ocp_t[i - 1];
            for (int i = 1; i <= 50; ++i)
                for (int q = 1; q <= 4; ++q) ocp_c[i][q] = t.ocp_c[(i - 1) * 4 + (q - 1)];
        }
        if (p.ph_mode == "spline") {
            Tables t;
            t.load(data_dir);
            ph_zn.assign(t.zn_pts.begin(), t.zn_pts.end());
            ph_mn.assign(t.mn_pts.begin(), t.mn_pts.end());
            ph_coef.assign(t.coef.begin(), t.coef.end());
        }
        if (p.transport == "quasi") {
            D_sp.assign(eq.ns, 0.0);
            for (int s = 0; s < eq.ns; ++s) {
                const std::string& nm = eq.names[s];
                D_sp[s] = nm == "H+" ? p.D_H : nm == "Zn+2" ? p.D_Zn : nm == "Mn+2" ? p.D_Mn : nm == "SO4-2" ? p.D_SO4
                        : nm == "OH-" ? p.D_OH : nm == "HSO4-" ? p.D_HSO4 : p.D_complex;
            }
        }
        H0.assign(nc + 1, 0.0);
        spec_x = Mat(4, nc); cache_T = Mat(4, nc); cache_lx = Mat(4, nc); x0_state = Mat(N, nc);
        mass = active_mass();
        const double i_ref = 1.0e-3 * mass;
        double v_cath = 0;
        for (j = 1; j <= nc; ++j)
            if (incath[j]) v_cath = v_cath + area_c[j] * dxc[j];
        for (int k = 1; k <= N; ++k) res_scale[k] = i_ref / (p.F * v_cath);
        res_scale[P1] = i_ref;
        res_scale[P2] = 1.0;
        for (int k : SPECIES) res_scale[k] = i_ref / p.F;
        if (p.zhs == "equilibrium") res_scale[ZH] = 1.0;
    }

    // ---------------------------------------------------------------- pH spline (ph_mode = spline)
    static void bspline_basis(const std::vector<double>& t, int i, double x, double bas[4]) {   // i: 0-based span
        double left[4] = {}, right[4] = {};
        for (int q = 0; q < 4; ++q) bas[q] = 0;
        bas[0] = 1.0;
        for (int j = 1; j <= 3; ++j) {
            left[j] = x - t[i + 1 - j];
            right[j] = t[i + j] - x;
            double saved = 0.0;
            for (int r = 0; r <= j - 1; ++r) {
                const double temp = bas[r] / (right[r + 1] + left[j - r]);
                bas[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            bas[j] = saved;
        }
    }
    static int span(const std::vector<double>& t, double x) {
        int cnt = 0;
        for (double v : t)
            if (v <= x) ++cnt;
        return std::min(std::max(cnt - 1, 3), static_cast<int>(t.size()) - 3 - 2);
    }
    double ph_spline(double zn_molar, double mn_molar) const {
        const double a = std::min(std::max(std::log10(mn_molar), ph_zn[3]), ph_zn[153 - 4] - 1e-12);
        const double b = std::min(std::max(std::log10(zn_molar), ph_mn[3]), ph_mn[153 - 4] - 1e-12);
        const int i = span(ph_zn, a), j = span(ph_mn, b);
        double bx[4], by[4], row[4];
        bspline_basis(ph_zn, i, a, bx);
        bspline_basis(ph_mn, j, b, by);
        for (int q = 0; q <= 3; ++q) {
            row[q] = 0;
            for (int r = 0; r <= 3; ++r) row[q] = row[q] + bx[r] * ph_coef[static_cast<std::size_t>(i - 3 + r) * 149 + (j - 3 + q)];
        }
        double v = 0;
        for (int q = 0; q <= 3; ++q) v = v + row[q] * by[q];
        return v;
    }

    // ---------------------------------------------------------------- speciation, constitutive laws
    Mat lx_of(const Mat& x, bool warm) {
        Mat T(4, nc), lx(4, nc);
        for (int j = 1; j <= nc; ++j)
            for (int q = 1; q <= 4; ++q) T(q, j) = x(TOT_COLS[q], j) * 1e3;
        if (cache_valid && T.d == cache_T.d) return cache_lx;
        for (int j = 1; j <= nc; ++j) {
            double xs[5] = {0, 0, 0, 0, 0}, Tj[5] = {0, T(1, j), T(2, j), T(3, j), T(4, j)};
            if (spec_valid)
                for (int q = 1; q <= 4; ++q) xs[q] = spec_x(q, j);
            bool ok = eq.solve(Tj, xs, spec_valid);
            if (!ok) ok = eq.solve(Tj, xs, false);
            if (!ok) throw Failure{};
            for (int q = 1; q <= 4; ++q) lx(q, j) = xs[q];
        }
        if (warm) {
            spec_x = lx; spec_valid = true;
            cache_T = T; cache_lx = lx; cache_valid = true;
        }
        return lx;
    }
    Mat free_of(const Mat& x, bool warm) {
        Mat fr = lx_of(x, warm);
        for (double& v : fr.d) v = lpow(10.0, v);
        return fr;
    }
    Mat basis_of(const Mat& x, const Mat& fr) const {
        Mat cb = fr;
        for (int j = 1; j <= nc; ++j) {
            if (p.basis == "totals") { cb(2, j) = x(ZN, j) * 1e3; cb(3, j) = x(MN, j) * 1e3; cb(4, j) = x(SO, j) * 1e3; }
            if (p.ph_mode == "zhs_equilibrium")
                cb(1, j) = std::exp((4.0 * std::log(cb(2, j)) + std::log(cb(4, j)) - p.logK_ZHS * LN10) / 6.0);
            else if (p.ph_mode == "fixed")
                cb(1, j) = lpow(10.0, -p.pH_fixed);
            else if (p.ph_mode == "spline")
                cb(1, j) = lpow(10.0, -ph_spline(x(ZN, j) * 1e3, x(MN, j) * 1e3));
        }
        return cb;
    }
    std::vector<double> porosity(const Mat& x) const {
        std::vector<double> eps(nc + 1);
        for (int j = 1; j <= nc; ++j) {
            eps[j] = region[j] == PROBE ? p.eps_probe : p.eps_sep;
            if (incath[j]) {
                double sv = 0;
                for (int s : SOLIDS) sv = sv + Vm[s] * x(s, j);
                eps[j] = 1.0 - eps_fixed - sv;
            }
        }
        return eps;
    }
    std::vector<double> tortuosity(const std::vector<double>& eps) const {
        std::vector<double> tau(nc + 1);
        for (int j = 1; j <= nc; ++j)
            if (eps[j] <= 0.0) throw Failure{};
        for (int j = 1; j <= nc; ++j) {
            tau[j] = p.tau_probe;
            if (region[j] == SEP) tau[j] = p.tau_factor_sep * lpow(eps[j], p.bruggeman_sep);
            if (region[j] == CATH) tau[j] = p.tau_factor_cath * lpow(eps[j], p.bruggeman_cath);
        }
        return tau;
    }
    double u3(double s, double c_zn) const {
        if (p.r3_ocp == "nernst") return p.U3_nernst + 0.5 / fRT * (std::log(c_zn) - s);
        double th, om;
        theta_pair(s, th, om);
        double endt = 0.0;
        if (p.r3_ocp == "spline_nernst") {
            const double epsw = p.r3_end_width, w = 0.5;
            const double lt = std::min(s, 0.0) - std::log1p(std::exp(-std::abs(s)));
            const double lo = std::min(-s, 0.0) - std::log1p(std::exp(-std::abs(s)));
            auto soft = [&](double a) { return -w * logaddexp(0.0, -a / w); };
            endt = -(0.5 / fRT) * (soft(lt - std::log(epsw)) - soft(lo - std::log(epsw)));
        }
        int cnt = 0;
        for (int i = 1; i <= 51; ++i)
            if (ocp_t[i] < th) ++cnt;
        const int k = std::min(std::max(cnt - 1, 0), 49) + 1;
        const double d = th - ocp_t[k];
        const double vs = ((ocp_c[k][1] * d + ocp_c[k][2]) * d + ocp_c[k][3]) * d + ocp_c[k][4];
        return (p.V_at_zmin - p.V_at_zmax) * vs + p.V_at_zmax + 0.5 / fRT * std::log(c_zn / p.c_ref3) + endt;
    }
    double area(const Mat& x, int k, int j) const {
        const double r = k == MO ? p.r_MnO2 : k == ZM ? p.r_ZMO : p.r_ZHS;
        return 3.0 * Vm[k] * std::max(x(k, j), 0.0) / r;
    }
    double precipitation(const Mat& x, const Mat& fr, int k, int j) const {
        if (!prec_on[k]) return 0.0;
        double sterm = 0.0;
        if (prec_s[k] != 0) sterm = prec_s[k] * std::log(fr(4, j));
        const double lq = prec_zn[k] * std::log(fr(2, j)) + sterm - prec_logk[k] * LN10;
        const double w = std::exp(lq / prec_h[k]) / fr(1, j);
        const double w_nuc = lpow(p.zhs_nucleation, prec_zn[k] / prec_h[k]);
        return prec_k[k] * (area(x, k, j) * (w - 1.0) + p.a_seed_ZHS * pos(w - w_nuc, W_SMOOTH));
    }
    double zhs_w(const Mat& fr, int j) const {
        return std::exp((4.0 * std::log(fr(2, j)) + std::log(fr(4, j)) - p.logK_ZHS * LN10) / 6.0) / fr(1, j);
    }
    void reactions(const Mat& x, const Mat& fr, const Mat& cb, int j, double& i1, double& i2, double& i3, double& r) const {
        const double cH = cb(1, j), cZn = cb(2, j), cMn = cb(3, j), ph1 = x(P1, j), ph2 = x(P2, j);
        const double a_MO = area(x, MO, j), a_ZM = area(x, ZM, j), a_ZH = area(x, ZH, j);
        const double u1 = p.U1 - 0.5 / fRT * (std::log(cMn) - 4.0 * std::log(cH));
        const double e1 = ph1 - ph2 - u1;
        const double bv1 = std::exp(p.alpha1 * 2.0 * fRT * e1) - std::exp(-(1.0 - p.alpha1) * 2.0 * fRT * e1);
        i1 = (-a_MO) * p.F * p.k1 * pos(-bv1, 1.0e-2) * onoff(p.R1_on);
        const double n2 = 2.0 - 2.0 * p.z_ZMO;
        const double u2 = p.U2 - (1.0 / (n2 * fRT)) * (p.z_ZMO * std::log(cZn) + std::log(cMn) - 4.0 * std::log(cH));
        const double e2 = ph1 - ph2 - u2;
        const double bv2 = std::exp(p.alpha2 * n2 * fRT * e2) - std::exp(-(1.0 - p.alpha2) * n2 * fRT * e2);
        if (p.zhs == "lumped")
            i2 = p.F * p.k2 * (a_ZH * pos(bv2, 1.0e-2) - a_ZM * pos(-bv2, 1.0e-2)) * onoff(p.R2_on);
        else
            i2 = p.F * p.k2 * (a_ZM * bv2 + (a_ZH + p.a_seed_R2) * pos(bv2, 1.0e-2)) * onoff(p.R2_on);
        double th, om;
        theta_pair(x(TH, j), th, om);
        const double e3 = ph1 - ph2 - u3(x(TH, j), cZn);
        const double i0 = p.F * p.k3 * std::sqrt(cZn * th * om);
        i3 = a_host * i0 * (std::exp(p.alpha3 * 2.0 * fRT * e3) - std::exp(-(1.0 - p.alpha3) * 2.0 * fRT * e3)) * onoff(p.R3_on);
        r = precipitation(x, fr, ZH, j);
    }
    double anode_current(const Mat& x, const Mat& fr) const {
        const double czn = p.basis == "free" ? fr(2, 1) : x(ZN, 1) * 1e3;
        const double eta = 0.0 - x(P2, 1) - 0.5 / fRT * std::log(czn);
        return p.F * p.k_an * std::sqrt(czn) * (std::exp(p.alpha_an * 2.0 * fRT * eta) - std::exp(-(1.0 - p.alpha_an) * 2.0 * fRT * eta));
    }
    static double d_theta(double s, double s_old) {
        double t, o, t0, o0;
        theta_pair(s, t, o);
        theta_pair(s_old, t0, o0);
        return (s > 0 && s_old > 0) ? o0 - o : t - t0;
    }

    // ---------------------------------------------------------------- residual
    void linear(const Mat& x, const Mat& old, double dt, double I, Mat& R, Blk& J) const {
        R = Mat(N, nc); J = Blk(N, nc);
        const auto eps = porosity(x), eps_old = porosity(old);
        for (int jj = 1; jj <= nc; ++jj) {
            for (int k : SPECIES) {
                if (!transported[k]) continue;
                R(k, jj) = vol[jj] * (eps[jj] * x(k, jj) - eps_old[jj] * old(k, jj)) / dt;
                J(k, k, jj) = vol[jj] * eps[jj] / dt;
                if (incath[jj])
                    for (int s : SOLIDS) J(k, s, jj) = vol[jj] * x(k, jj) * (-Vm[s]) / dt;
            }
            if (p.species == "no_H") {
                R(H, jj) = x(H, jj) - H0[jj];
                J(H, H, jj) = 1.0;
            }
            R(P2, jj) = (2.0 * x(ZN, jj) + 2.0 * x(MN, jj) + x(H, jj) - 2.0 * x(SO, jj)) * 1e3;
            J(P2, ZN, jj) = 2e3; J(P2, MN, jj) = 2e3; J(P2, H, jj) = 1e3; J(P2, SO, jj) = -2e3;
            if (!incath[jj]) {
                R(P1, jj) = x(P1, jj); J(P1, P1, jj) = 1.0;
                for (int k : SOLIDS) { R(k, jj) = x(k, jj); J(k, k, jj) = 1.0; }
                R(TH, jj) = x(TH, jj) - old(TH, jj);
                J(TH, TH, jj) = 1.0;
            } else {
                for (int k : SOLIDS) {
                    if (k == ZH && p.zhs == "equilibrium") continue;
                    R(k, jj) = (x(k, jj) - old(k, jj)) / dt;
                    J(k, k, jj) = 1.0 / dt;
                }
                if (p.zhs == "equilibrium") {
                    const int kks[3] = {ZN, SO, H};
                    const double nus[3] = {4.0, 1.0, -6.0};
                    for (int q = 0; q < 3; ++q) {
                        const int kk = kks[q];
                        if (!transported[kk]) continue;
                        R(kk, jj) = R(kk, jj) + nus[q] * vol[jj] * (x(ZH, jj) - old(ZH, jj)) / dt;
                        J(kk, ZH, jj) = J(kk, ZH, jj) + nus[q] * vol[jj] / dt;
                    }
                }
                const double cap = n_host * (p.zmax - p.zmin);
                double th, om;
                theta_pair(x(TH, jj), th, om);
                R(TH, jj) = cap * d_theta(x(TH, jj), old(TH, jj)) / dt;
                J(TH, TH, jj) = cap * th * om / dt;
            }
        }
        R(P1, nc) = R(P1, nc) + I;
    }

    Mat sources(const Mat& x, bool warm) {
        Mat R(N, nc);
        const Mat fr = free_of(x, warm);
        const Mat cb = basis_of(x, fr);
        const double n2 = 2.0 - 2.0 * p.z_ZMO;
        for (int j = c0; j <= nc; ++j) {
            const double vc = vol[j];
            double i1, i2, i3, r_zhs;
            reactions(x, fr, cb, j, i1, i2, i3, r_zhs);
            const double xi1 = -i1 / (2.0 * p.F), xi2 = -i2 / (n2 * p.F), xi3 = -i3 / (2.0 * p.F);
            if (p.zhs == "lumped") r_zhs = (2.0 / 3.0) * (xi1 + xi2);
            else if (p.zhs == "off" || p.zhs == "equilibrium") r_zhs = 0.0;
            double rk[N + 1] = {};
            rk[ZH] = r_zhs;
            rk[ZO] = precipitation(x, fr, ZO, j);
            rk[ZX] = precipitation(x, fr, ZX, j);
            double dzn = p.z_ZMO * xi2 - xi3;
            double dso = 0.0 * xi2;
            double dh = -4.0 * xi1 - 4.0 * xi2;
            for (int k : PRECIP) {
                dzn = dzn - prec_zn[k] * rk[k];
                dso = dso - prec_s[k] * rk[k];
                dh = dh + prec_h[k] * rk[k];
            }
            R(ZN, j) = -vc * dzn;
            R(MN, j) = -vc * (xi1 + xi2);
            R(SO, j) = -vc * dso;
            if (p.species == "with_H") R(H, j) = -vc * dh;
            R(P1, j) = vc * (i1 + i2 + i3);
            R(MO, j) = xi1;
            R(ZM, j) = xi2;
            R(TH, j) = -xi3;
            for (int k : PRECIP) R(k, j) = -rk[k];
            if (p.zhs == "equilibrium") {
                const double a = x(ZH, j) / FB_REF, b = 1.0 - zhs_w(fr, j);
                R(ZH, j) = a + b - std::sqrt(a * a + b * b + 1e-20);
            }
        }
        R(ZN, 1) = R(ZN, 1) - anode_current(x, fr) * area_c[1] / (2.0 * p.F);
        return R;
    }

    // derivatives of an east-face outflow of `row` (faces 1..nc-1): dL(f, k), dR(f, k) as Mat(nc-1, N) via (f, k)
    void add_face(Blk& A, Blk& B, Blk& D, int row, const Mat& dL, const Mat& dR) const {
        for (int j = 1; j <= nc - 1; ++j)
            for (int k = 1; k <= N; ++k) {
                B(row, k, j) = B(row, k, j) + dL(j, k);
                D(row, k, j) = D(row, k, j) + dR(j, k);
            }
        for (int j = 2; j <= nc; ++j)
            for (int k = 1; k <= N; ++k) {
                A(row, k, j) = A(row, k, j) - dL(j - 1, k);
                B(row, k, j) = B(row, k, j) - dR(j - 1, k);
            }
    }

    void transport(const Mat& x, bool jac, Mat& R, Blk& A, Blk& B, Blk& D) {
        R = Mat(N, nc); A = Blk(N, nc); B = Blk(N, nc); D = Blk(N, nc);
        const auto eps = porosity(x);
        const auto tau = tortuosity(eps);
        const int nf = nc - 1;
        std::vector<double> half(nc + 1), G(nc), dhalf(nc + 1), dphi(nc), flow(nc), hs(nc + 1), Gs(nc), dhs(nc + 1);
        for (int j = 1; j <= nc; ++j) {
            half[j] = 0.5 * dxc[j] * tau[j] / (eps[j] * area_c[j]);
            double bexp = 0.0;
            if (incath[j]) bexp = p.bruggeman_cath;
            dhalf[j] = 0.0;
            if (incath[j]) dhalf[j] = (bexp - 1.0) * half[j] / eps[j];
        }
        for (int j = 1; j <= nf; ++j) {
            G[j] = 1.0 / (half[j] + half[j + 1]);
            dphi[j] = x(P2, j + 1) - x(P2, j);
        }
        // dL(f, k) with f = face (1..nf) as the "unknown" index of a Mat(nf, N): Mat(rows = nf, cols = N)
        Mat dL(nf, N), dR(nf, N);
        auto add_flow = [&](int row) {
            for (int jf = 1; jf <= nf; ++jf) R(row, jf) = R(row, jf) + flow[jf];
            for (int jf = 2; jf <= nc; ++jf) R(row, jf) = R(row, jf) - flow[jf - 1];
        };
        auto solid_terms = [&]() {
            for (int kf : SOLIDS)
                for (int jf = 1; jf <= nf; ++jf) {
                    dL(jf, kf) = dL(jf, kf) + (-flow[jf]) * G[jf] * dhalf[jf] * (-Vm[kf]);
                    dR(jf, kf) = dR(jf, kf) + (-flow[jf]) * G[jf] * dhalf[jf + 1] * (-Vm[kf]);
                }
        };
        if (p.transport == "ions") {
            for (int k : SPECIES) {
                if (!transported[k]) continue;
                const double zf = CHARGE[k] * fRT;
                dL = Mat(nf, N); dR = Mat(nf, N);
                for (int j = 1; j <= nf; ++j) {
                    const double arg = zf * dphi[j];
                    const double Bp = bernoulli(arg), Bm = bernoulli(-arg);
                    flow[j] = Dk[k] * G[j] * (Bp * x(k, j) - Bm * x(k, j + 1));
                    if (jac) {
                        dL(j, k) = Dk[k] * G[j] * Bp;
                        dR(j, k) = -(Dk[k] * G[j] * Bm);
                        const double dfa = Dk[k] * G[j] * (dbernoulli(arg) * x(k, j) + dbernoulli(-arg) * x(k, j + 1)) * zf;
                        dR(j, P2) = dfa; dL(j, P2) = -dfa;
                    }
                }
                add_flow(k);
                if (jac) { solid_terms(); add_face(A, B, D, k, dL, dR); }
            }
        } else {
            const int ns = eq.ns;
            const Mat lx = lx_of(x, true);
            Mat cs(ns, nc), fl(ns, nf), dfas(ns, nf), Bps(ns, nf), Bms(ns, nf);
            std::vector<double> dc(static_cast<std::size_t>(ns) * 4 * nc);     // dc(s, q, j)
            auto DC = [&](int s, int q, int j) -> double& { return dc[(static_cast<std::size_t>(j - 1) * 4 + (q - 1)) * ns + (s - 1)]; };
            std::vector<double> m;
            for (int j = 1; j <= nc; ++j) {
                const double xl[5] = {0, lx(1, j), lx(2, j), lx(3, j), lx(4, j)};
                eq.species(xl, m);
                for (int s = 1; s <= ns; ++s) cs(s, j) = m[s - 1] * 1e-3;
            }
            for (int j = 1; j <= nf; ++j)
                for (int s = 1; s <= ns; ++s) {
                    const double arg = (eq.z[s - 1] * fRT) * dphi[j];
                    Bps(s, j) = bernoulli(arg); Bms(s, j) = bernoulli(-arg);
                    fl(s, j) = D_sp[s - 1] * G[j] * (Bps(s, j) * cs(s, j) - Bms(s, j) * cs(s, j + 1));
                    if (jac) dfas(s, j) = D_sp[s - 1] * G[j] * (dbernoulli(arg) * cs(s, j) + dbernoulli(-arg) * cs(s, j + 1)) * eq.z[s - 1] * fRT;
                }
            if (jac)
                for (int j = 1; j <= nc; ++j) {
                    double Tj[5], S[5][5];
                    const double xl[5] = {0, lx(1, j), lx(2, j), lx(3, j), lx(4, j)};
                    for (int qq = 1; qq <= 4; ++qq) Tj[qq] = x(TOT_COLS[qq], j) * 1e3;
                    eq.sensitivity(xl, Tj, S);
                    for (int s = 1; s <= ns; ++s)
                        for (int qq = 1; qq <= 4; ++qq) {
                            double e_ = 0;
                            for (int mm = 1; mm <= 4; ++mm) e_ = e_ + eq.nu[s - 1][mm] * S[mm][qq];
                            DC(s, qq, j) = cs(s, j) * LN10 * e_ * 1e3;
                        }
                }
            for (int q = 1; q <= 4; ++q) {
                const int k = TOT_COLS[q];
                if (!transported[k]) continue;
                dL = Mat(nf, N); dR = Mat(nf, N);
                for (int j = 1; j <= nf; ++j) {
                    double fsum = 0;
                    for (int s = 1; s <= ns; ++s) fsum = fsum + fl(s, j) * eq.nu[s - 1][q];
                    flow[j] = fsum;
                    if (jac) {
                        for (int qq = 1; qq <= 4; ++qq) {
                            const int kk = TOT_COLS[qq];
                            double dsum = 0;
                            for (int s = 1; s <= ns; ++s) {
                                const double wL = D_sp[s - 1] * G[j] * Bps(s, j) * eq.nu[s - 1][q];
                                dsum = dsum + wL * DC(s, qq, j);
                            }
                            dL(j, kk) = dsum;
                            dsum = 0;
                            for (int s = 1; s <= ns; ++s) {
                                const double wR = -(D_sp[s - 1] * G[j] * Bms(s, j) * eq.nu[s - 1][q]);
                                dsum = dsum + wR * DC(s, qq, j + 1);
                            }
                            dR(j, kk) = dsum;
                        }
                        double dsum = 0;
                        for (int s = 1; s <= ns; ++s) dsum = dsum + dfas(s, j) * eq.nu[s - 1][q];
                        dR(j, P2) = dsum; dL(j, P2) = -dsum;
                    }
                }
                add_flow(k);
                if (jac) { solid_terms(); add_face(A, B, D, k, dL, dR); }
            }
        }
        // solid current (cathode faces only)
        for (int j = 1; j <= nc; ++j) {
            double sig = 1.0;
            if (incath[j]) sig = p.sigma * (1.0 - eps[j]);
            hs[j] = 0.5 * dxc[j] / (sig * area_c[j]);
            dhs[j] = 0.0;
            if (incath[j]) dhs[j] = hs[j] / (1.0 - eps[j]);
        }
        for (int j = 1; j <= nf; ++j) {
            Gs[j] = 0.0;
            if (incath[j] && incath[j + 1]) Gs[j] = 1.0 / (hs[j] + hs[j + 1]);
            flow[j] = Gs[j] * (x(P1, j) - x(P1, j + 1));
        }
        add_flow(P1);
        if (jac) {
            dL = Mat(nf, N); dR = Mat(nf, N);
            for (int j = 1; j <= nf; ++j) {
                dL(j, P1) = Gs[j]; dR(j, P1) = -Gs[j];
                for (int ss : SOLIDS) {
                    dL(j, ss) = -(flow[j] * Gs[j] * dhs[j] * (-Vm[ss]));
                    dR(j, ss) = -(flow[j] * Gs[j] * dhs[j + 1] * (-Vm[ss]));
                }
            }
            add_face(A, B, D, P1, dL, dR);
        }
    }

    // ---------------------------------------------------------------- Jacobian and Newton
    void blocks(const Mat& x, const Mat& old, double dt, double I, Mat& R, Blk& A, Blk& B, Blk& D) {
        Mat RL, RT;
        Blk BT;
        linear(x, old, dt, I, RL, B);
        const Mat S0 = sources(x, true);
        transport(x, true, RT, A, BT, D);
        for (std::size_t q = 0; q < B.d.size(); ++q) B.d[q] = B.d[q] + BT.d[q];
        Mat hstep(N, nc);
        for (int j = 1; j <= nc; ++j)
            for (int k = 1; k <= N; ++k) hstep(k, j) = 1e-7 * std::max(std::abs(x(k, j)), STEP_FLOOR[k]);
        const Mat spec_base = spec_x, cT = cache_T, clx = cache_lx;
        const bool sv_base = spec_valid, cv_base = cache_valid;
        for (int k = 1; k <= N; ++k) {
            Mat xp = x;
            for (int j = 1; j <= nc; ++j) xp(k, j) = xp(k, j) + hstep(k, j);
            spec_x = spec_base; spec_valid = sv_base;
            const Mat Sp = sources(xp, false);
            for (int j = 1; j <= nc; ++j)
                for (int ii = 1; ii <= N; ++ii) B(ii, k, j) = B(ii, k, j) + (Sp(ii, j) - S0(ii, j)) / hstep(k, j);
        }
        spec_x = spec_base; spec_valid = sv_base;
        cache_T = cT; cache_lx = clx; cache_valid = cv_base;
        R = Mat(N, nc);
        for (std::size_t q = 0; q < R.d.size(); ++q) R.d[q] = RL.d[q] + S0.d[q] + RT.d[q];
    }

    double damping(const Mat& x, const Mat& dx) const {
        double lam = 1.0;
        for (int q = 0; q < 3; ++q) {
            const int k = SPECIES[q];
            for (int j = 1; j <= nc; ++j)
                if (x(k, j) + dx(k, j) < 0 && x(k, j) > 0) lam = std::min(lam, 0.9 * x(k, j) / (-dx(k, j)));
        }
        double b1 = 0, b2 = 0;
        for (int j = 1; j <= nc; ++j) { b1 = std::max(b1, std::abs(dx(P1, j))); b2 = std::max(b2, std::abs(dx(P2, j))); }
        const double big = std::max(b1, b2);
        if (big > 0.2) lam = std::min(lam, 0.2 / big);
        return lam;
    }

    Mat newton_step(const Mat& old, double dt, double I) {
        if (!have_H0) {
            for (int j = 1; j <= nc; ++j) H0[j] = old(H, j);
            have_H0 = true;
        }
        Mat x = old, xout;
        cache_valid = false; spec_valid = false;
        free_of(old, true);
        const double eps_m = std::numeric_limits<double>::epsilon();
        for (int it = 1; it <= p.newton_max_iter; ++it) {
            Mat R;
            Blk A, B, D;
            blocks(x, old, dt, I, R, A, B, D);
            for (double v : R.d)
                if (!std::isfinite(v)) throw Failure{};
            Mat rtol(N, nc);
            for (int j = 1; j <= nc; ++j)
                for (int ii = 1; ii <= N; ++ii) {
                    double mag = 0;
                    for (int k = 1; k <= N; ++k) mag = mag + std::abs(B(ii, k, j)) * std::abs(x(k, j));
                    if (j > 1)
                        for (int k = 1; k <= N; ++k) mag = mag + std::abs(A(ii, k, j)) * std::abs(x(k, j - 1));
                    if (j < nc)
                        for (int k = 1; k <= N; ++k) mag = mag + std::abs(D(ii, k, j)) * std::abs(x(k, j + 1));
                    rtol(ii, j) = std::max(RES_TOL * res_scale[ii], 1e3 * eps_m * mag);
                }
            std::vector<double> thv(nc + 1), omv(nc + 1), dth(nc + 1);
            std::vector<char> in_theta(nc + 1);
            for (int j = 1; j <= nc; ++j) {
                theta_pair(x(TH, j), thv[j], omv[j]);
                in_theta[j] = x(TH, j) > 0.0;
            }
            for (int j = 1; j <= nc; ++j)
                if (in_theta[j])
                    for (int ii = 1; ii <= N; ++ii) B(ii, TH, j) = B(ii, TH, j) / (thv[j] * omv[j]);
            std::vector<double> G(static_cast<std::size_t>(N) * nc);
            for (int j = 1; j <= nc; ++j) {
                for (int k = 1; k <= N; ++k)
                    for (int ii = 1; ii <= N; ++ii) {
                        A(ii, k, j) = A(ii, k, j) * TYP[k]; B(ii, k, j) = B(ii, k, j) * TYP[k]; D(ii, k, j) = D(ii, k, j) * TYP[k];
                    }
                for (int ii = 1; ii <= N; ++ii) {
                    double ma = 0, mb = 0, md = 0;
                    for (int k = 1; k <= N; ++k) {
                        ma = std::max(ma, std::abs(A(ii, k, j))); mb = std::max(mb, std::abs(B(ii, k, j))); md = std::max(md, std::abs(D(ii, k, j)));
                    }
                    const double rs = std::max(std::max(std::max(ma, mb), md), 1e-300);
                    for (int k = 1; k <= N; ++k) { A(ii, k, j) = A(ii, k, j) / rs; B(ii, k, j) = B(ii, k, j) / rs; D(ii, k, j) = D(ii, k, j) / rs; }
                    G[static_cast<std::size_t>(j - 1) * N + (ii - 1)] = -R(ii, j) / rs;
                }
            }
            std::vector<double> sol;
            if (!band_partial(N, nc, A.d, B.d, D.d, G, sol)) throw Failure{};
            Mat dx(N, nc);
            for (int j = 1; j <= nc; ++j)
                for (int k = 1; k <= N; ++k) dx(k, j) = sol[static_cast<std::size_t>(j - 1) * N + (k - 1)] * TYP[k];
            for (int j = 1; j <= nc; ++j) { dth[j] = dx(TH, j); dx(TH, j) = 0.0; }
            const double lam = damping(x, dx);
            xout = Mat(N, nc);
            for (std::size_t q = 0; q < x.d.size(); ++q) xout.d[q] = x.d[q] + lam * dx.d[q];
            for (int k : SOLIDS)
                for (int j = 1; j <= nc; ++j)
                    if (xout(k, j) < 0) xout(k, j) = 0.1 * x(k, j);
            double upd = 0;
            bool nan_upd = false;
            for (int j = 1; j <= nc; ++j) {
                const double dstep = lam * dth[j];
                double thn = thv[j] + dstep, omn = omv[j] - dstep;
                if (thn < TH_FLOOR * thv[j]) { thn = TH_FLOOR * thv[j]; omn = 1.0 - thn; }
                if (omn < TH_FLOOR * omv[j]) { omn = TH_FLOOR * omv[j]; thn = 1.0 - omn; }
                const double s_theta = std::log(std::max(thn, 1e-300)) - std::log(std::max(omn, 1e-300));
                const double s_s = x(TH, j) + std::min(std::max(dstep, -S_STEP), S_STEP);
                xout(TH, j) = in_theta[j] ? std::min(std::max(s_theta, -S_MAX), S_MAX) : std::min(std::max(s_s, -S_MAX), S_MAX);
                for (int k = 1; k <= N; ++k) {
                    double sc;
                    if (k == TH) sc = in_theta[j] ? std::abs(dth[j]) : std::abs(dth[j]) * thv[j] * omv[j];
                    else sc = std::abs(dx(k, j)) / TYP[k];
                    if (sc != sc) nan_upd = true;
                    upd = fmax_f(upd, sc);
                }
            }
            x = xout;
            if (nan_upd || !std::isfinite(upd)) break;
            bool conv = lam == 1.0 && upd < p.newton_tol;
            if (conv)
                for (std::size_t q = 0; q < R.d.size(); ++q)
                    if (!(std::abs(R.d[q]) <= rtol.d[q])) { conv = false; break; }
            if (conv) {
                free_of(x, true);
                return x;
            }
        }
        throw Failure{};
    }
    // Fortran max(a, b) as gfortran evaluates it (b when a is NaN... irrelevant here: NaN flagged separately)
    static double fmax_f(double a, double b) { return a >= b ? a : b; }

    // ---------------------------------------------------------------- initial state, outputs
    Mat initial_state() {
        Mat x(N, nc);
        const double c_zn = std::max(p.c_ZnSO4, TRACE), c_mn = std::max(p.c_MnSO4, TRACE);
        for (int j = 1; j <= nc; ++j) {
            x(ZN, j) = c_zn * 1e-3;
            x(MN, j) = c_mn * 1e-3;
            x(SO, j) = (c_zn + c_mn + p.c_H2SO4) * 1e-3;
            x(H, j) = 2.0 * p.c_H2SO4 * 1e-3;
            H0[j] = x(H, j);
        }
        have_H0 = true;
        for (int j = c0; j <= nc; ++j) {
            x(MO, j) = p.vf_MnO2 / Vm[MO];
            x(ZM, j) = p.vf_ZMO / Vm[ZM];
            x(ZH, j) = p.vf_ZHS / Vm[ZH];
        }
        for (int j = 1; j <= nc; ++j) x(TH, j) = std::log(p.theta0 / (1.0 - p.theta0));
        const Mat fr = free_of(x, true);
        const double u_zn = 0.5 / fRT * std::log(fr(2, 1));
        for (int j = 1; j <= nc; ++j) x(P2, j) = -u_zn;
        const Mat cb = basis_of(x, fr);
        for (int j = c0; j <= nc; ++j) x(P1, j) = x(P2, j) + u3(x(TH, j), cb(2, j));
        return newton_step(x, 1.0e-4, 0.0);
    }

    double voltage(const Mat& x, double I) const {
        const auto eps = porosity(x);
        return x(P1, nc) - I * 0.5 * dxc[nc] / (p.sigma * (1.0 - eps[nc]) * area_c[nc]);
    }

    std::string limit_reason(const Mat& x, double I) const {
        const auto eps = porosity(x);
        double mx = -std::numeric_limits<double>::max(), mx0 = mx, th_max = -INF, th_min = INF, eps_min = INF;
        double zn_min = INF, zn0_min = INF, mn_min = INF, mn0_min = INF;
        for (int j = 1; j <= nc; ++j) {
            zn_min = std::min(zn_min, x(ZN, j)); zn0_min = std::min(zn0_min, x0_state(ZN, j));
            mn_min = std::min(mn_min, x(MN, j)); mn0_min = std::min(mn0_min, x0_state(MN, j));
        }
        for (int j = c0; j <= nc; ++j) {
            mx = std::max(mx, x(MO, j) + x(ZM, j));
            mx0 = std::max(mx0, x0_state(MO, j) + x0_state(ZM, j));
            const double th = sigmoid(x(TH, j));
            th_max = std::max(th_max, th); th_min = std::min(th_min, th);
            eps_min = std::min(eps_min, eps[j]);
        }
        if (zn_min < LIMIT * zn0_min) return "zinc_depleted";
        if (p.c_MnSO4 > 0 && mn_min < LIMIT * mn0_min) return "manganese_depleted";
        if (p.R3_on && I > 0 && th_max > 1.0 - LIMIT) return "insertion_full";
        if (p.R3_on && I < 0 && th_min < LIMIT) return "insertion_empty";
        if (mx < LIMIT * std::max(mx0, 1e-300)) return "dissolvable_mno2_exhausted";
        if (eps_min < LIMIT) return "pores_clogged";
        return "";
    }

    // Python's f"{v:15.7E}" (and gfortran's ES15.7E2 / ES15.7E3, which the Fortran program writes)
    static std::string fmt_e(double v) {
        char b[64];
        if (std::isnan(v)) return "            NAN";
        if (std::isinf(v)) return v > 0 ? "            INF" : "           -INF";
        std::snprintf(b, sizeof b, "%15.7E", v);
        return b;
    }

    void write_header() {
        static const char* cols[19] = {"t_h", "V", "I_mAg", "mAhg", "step", "pH_cath", "pH_probe", "pH_anode", "Zn_cath_M",
                                       "Mn_cath_M", "S_cath_M", "vf_MnO2", "vf_ZMO", "vf_ZHS", "theta", "i_R1_mAg",
                                       "i_R2_mAg", "i_R3_mAg", "eps_cath"};
        std::string line;
        char b[32];
        for (int i = 0; i < 19; ++i) {
            std::snprintf(b, sizeof b, "%15s", cols[i]);
            line += (i ? " " : "") + std::string(b);
        }
        out_rows.push_back(line);
    }

    void write_row(double t, const Mat& x, double mAhg, double I, int k) {
        const Mat fr = free_of(x, true);
        const Mat cb = basis_of(x, fr);
        const auto eps = porosity(x);
        double s1 = 0, s2 = 0, s3 = 0;
        for (int j = c0; j <= nc; ++j) {
            double i1, i2, i3, rz;
            reactions(x, fr, cb, j, i1, i2, i3, rz);
            s1 = s1 + vol[j] * i1; s2 = s2 + vol[j] * i2; s3 = s3 + vol[j] * i3;
        }
        double vsum = 0;
        for (int j = c0; j <= nc; ++j) vsum = vsum + vol[j];
        auto cmean = [&](auto get) {
            double acc = 0;
            for (int jj = c0; jj <= nc; ++jj) acc = acc + vol[jj] * get(jj);
            return acc / vsum;
        };
        double row[20];
        row[1] = t / 3600.0;
        row[2] = voltage(x, I);
        row[3] = I * 1.0e3 / mass;
        row[4] = mAhg;
        row[6] = -std::log10(cmean([&](int j) { return fr(1, j); }));
        row[7] = -std::log10(fr(1, i_probe));
        row[8] = -std::log10(fr(1, 1));
        row[9] = cmean([&](int j) { return x(ZN, j); }) * 1e3;
        row[10] = cmean([&](int j) { return x(MN, j); }) * 1e3;
        row[11] = cmean([&](int j) { return x(SO, j); }) * 1e3;
        row[12] = cmean([&](int j) { return x(MO, j) * Vm[MO]; });
        row[13] = cmean([&](int j) { return x(ZM, j) * Vm[ZM]; });
        row[14] = cmean([&](int j) { return x(ZH, j) * Vm[ZH]; });
        row[15] = cmean([&](int j) { return sigmoid(x(TH, j)); });
        row[16] = s1 * (1.0e3 / mass);
        row[17] = s2 * (1.0e3 / mass);
        row[18] = s3 * (1.0e3 / mass);
        row[19] = cmean([&](int j) { return eps[j]; });
        std::string line;
        char b[32];
        for (int q = 1; q <= 19; ++q) {
            if (q > 1) line += " ";
            if (q == 5) { std::snprintf(b, sizeof b, "%15d", k); line += b; }
            else line += fmt_e(row[q]);
        }
        out_rows.push_back(line);
    }
};

std::vector<Step> parse_protocol(const std::string& text, const Params& p) {
    std::vector<Step> steps;
    std::string rest = text;
    if (Eq::trim(rest).empty()) die("empty protocol");
    std::stringstream ss(text);
    std::string part;
    while (std::getline(ss, part, ';')) {
        std::istringstream ws(part);
        std::string word;
        if (!(ws >> word)) continue;
        Step st;
        st.kind = lower(word);
        if (st.kind != "cc" && st.kind != "cv" && st.kind != "rest")
            die("unknown step type '" + word + "' (expected cc, cv or rest)");
        st.Vmin = p.V_min; st.Vmax = p.V_max;
        bool has_i = false, has_v = false;
        while (ws >> word) {
            const auto eqp = word.find('=');
            if (eqp == std::string::npos) die("expected key=value, got '" + word + "'");
            const std::string key = lower(word.substr(0, eqp));
            std::string v = word.substr(eqp + 1);
            std::replace(v.begin(), v.end(), 'd', 'e');
            std::replace(v.begin(), v.end(), 'D', 'e');
            const double val = std::stod(v);
            const std::string kk = st.kind + ":" + key;
            if (kk == "cc:i") { st.I = val; has_i = true; }
            else if (kk == "cc:t" || kk == "cv:t" || kk == "rest:t") st.t = val;
            else if (kk == "cc:vmin") st.Vmin = val;
            else if (kk == "cc:vmax") st.Vmax = val;
            else if (kk == "cv:v") { st.V = val; has_v = true; }
            else if (kk == "cv:imin") st.Imin = val;
            else die(st.kind + " does not take '" + key + "'");
        }
        if (st.t != -1 && st.t <= 0) die("t must be positive");
        if (st.kind == "cc" && !has_i) die("cc needs I=");
        if (st.kind == "cv" && !has_v) die("cv needs V=");
        if (st.kind == "cv" && st.t < 0 && st.Imin < 0) die("cv needs t= or Imin= to end");
        if (st.kind == "rest" && st.t < 0) die("rest needs t=");
        steps.push_back(st);
    }
    std::vector<Step> all;
    for (int c = 0; c < std::max(1, p.cycles); ++c) all.insert(all.end(), steps.begin(), steps.end());
    return all;
}

struct Driver {
    Model& m;
    explicit Driver(Model& m_) : m(m_) {}

    double margin_of(const Mat& xn, const Step& st, double I) const {
        const double v = m.voltage(xn, I);
        if (I > 0) return v - st.Vmin;
        if (I < 0) return st.Vmax - v;
        return std::min(v - st.Vmin, st.Vmax - v);
    }

    // returns false on failure; state and t_done hold the progress made
    bool advance(Mat& state, double dt, double I, bool use_margin, const Step& st, double& t_done, bool& stopped) {
        t_done = 0.0;
        double hh = dt;
        int failures = 0;
        stopped = false;
        while (t_done < dt) {
            hh = std::min(hh, dt - t_done);
            Mat nw;
            try {
                nw = m.newton_step(state, hh, I);
            } catch (const Failure&) {
                ++failures;
                if (hh / 2 < MIN_SUBSTEP || failures >= MAX_FAILURES) return false;
                hh = hh / 2;
                continue;
            }
            if (use_margin) {
                const double mg = margin_of(nw, st, I);
                if (mg < 0.0) {
                    if (mg < -EVENT_DV && hh / 2 >= EVENT_MIN_DT) { hh = hh / 2; continue; }
                    state = nw;
                    t_done = t_done + hh;
                    stopped = true;
                    return true;
                }
            }
            state = nw;
            t_done = t_done + hh;
            hh = 2.0 * hh;
        }
        return true;
    }

    bool cv_step(const Mat& state, double hh, double V_set, double& I, Mat& nw) {
        std::vector<std::pair<double, Mat>> states;
        auto f = [&](double Ic) {
            Mat xn;
            try {
                xn = m.newton_step(state, hh, Ic);
            } catch (const Failure&) {
                return Ic < 0 ? INF : -INF;
            }
            states.emplace_back(Ic, xn);
            return m.voltage(xn, Ic) - V_set;
        };
        auto take = [&](double Ic) {
            for (auto it = states.rbegin(); it != states.rend(); ++it)
                if (it->first == Ic) { nw = it->second; return true; }
            return false;
        };
        double fI = f(I);
        if (std::abs(fI) <= CV_TOL) return take(I);
        double grow = std::max(std::abs(I), 1.0e-3 * m.mass);
        bool have_a = false, have_b = false;
        double a = 0, b = 0, fa = 0, fb = 0;
        for (int it = 1; it <= 60; ++it) {
            if (fI > 0) {
                a = I; fa = fI; have_a = true;
                if (have_b) break;
                I = I < 0 ? 0.0 : I + grow;
            } else {
                b = I; fb = fI; have_b = true;
                if (have_a) break;
                I = I > 0 ? 0.0 : I - grow;
            }
            grow = grow * 2.0;
            fI = f(I);
            if (std::abs(fI) <= CV_TOL) return take(I);
        }
        if (!(have_a && have_b)) return false;
        int side = 0;
        for (int it = 1; it <= 200; ++it) {
            if (std::isfinite(fa) && std::isfinite(fb)) {
                I = (a * fb - b * fa) / (fb - fa);
                if (!(a < I && I < b)) I = 0.5 * (a + b);
            } else {
                I = 0.5 * (a + b);
            }
            fI = f(I);
            // a collapsed bracket is a solution only at the set voltage (not across a jump in V(I))
            if (std::abs(fI) <= CV_TOL || (b - a) <= 1.0e-14 * m.mass) return std::abs(fI) <= CV_ACCEPT && take(I);
            if (fI > 0) {
                a = I; fa = fI;
                if (side == 1 && std::isfinite(fb)) fb = fb * 0.5;
                side = 1;
            } else {
                b = I; fb = fI;
                if (side == -1 && std::isfinite(fa)) fa = fa * 0.5;
                side = -1;
            }
        }
        return false;
    }

    // a constant-voltage sub-step: cv_step over hh, halved (down to MIN_SUBSTEP) while no current holds V_set for
    // that long; I is unchanged on failure
    bool cv_advance(const Mat& state, double hh, double V_set, double& I, Mat& nw, double& h_done) {
        const double I_guess = I;
        h_done = hh;
        while (true) {
            I = I_guess;
            if (cv_step(state, h_done, V_set, I, nw)) return true;
            if (h_done / 2 < MIN_SUBSTEP) { I = I_guess; h_done = 0; return false; }
            h_done = h_done / 2;
        }
    }

    std::string run(int& n_done, double& mAhg) {
        const Params& p = m.p;
        const auto steps = parse_protocol(p.steps, p);
        m.write_header();
        Mat state;
        try {
            state = m.initial_state();
        } catch (const Failure&) {
            die("the initial state could not be solved");
        }
        m.x0_state = state;
        double t = 0.0, I = 0.0;
        mAhg = 0.0; n_done = 0;
        auto first_current = [&](const Step& st) { return st.kind == "cc" ? st.I * 1.0e-3 * m.mass : 0.0; };
        I = first_current(steps[0]);
        m.write_row(t, state, mAhg, I, 1);
        double last_write = t;
        std::string reason;
        const int nsteps = static_cast<int>(steps.size());
        for (int k = 1; k <= nsteps; ++k) {
            const Step& st = steps[k - 1];
            double t_step = 0.0;
            if (st.kind != "cv") I = first_current(st);
            while (true) {
                double hh = p.dt;
                if (st.t >= 0) hh = std::min(p.dt, st.t - t_step);
                Mat nw;
                double h_done = 0;
                bool stopped = false, ok;
                std::string why;
                if (st.kind == "cv") {
                    ok = cv_advance(state, hh, st.V, I, nw, h_done);
                    stopped = st.Imin >= 0 && std::abs(I) <= st.Imin * 1.0e-3 * m.mass;
                    why = "current_limit";
                } else {
                    nw = state;
                    ok = advance(nw, hh, I, st.kind == "cc", st, h_done, stopped);
                    if (ok) why = (stopped && m.voltage(nw, I) <= st.Vmin) ? "cutoff_low" : "cutoff_high";
                }
                if (!ok) {
                    if (h_done > 0.0) {
                        mAhg = mAhg + 1000.0 * (I / m.mass) * h_done / 3600.0;
                        state = nw; t = t + h_done; ++n_done;
                    }
                    why = m.limit_reason(state, I);
                    if (why.empty() || p.end_on_cutoff) {
                        if (why.empty()) why = "solver_fail";
                        m.write_row(t, state, mAhg, I, k);
                        return why;
                    }
                    m.write_row(t, state, mAhg, I, k);
                    last_write = t; reason = why;
                    break;
                }
                mAhg = mAhg + 1000.0 * (I / m.mass) * h_done / 3600.0;
                state = nw; t = t + h_done; t_step = t_step + h_done; ++n_done;
                for (double v : state.d)
                    if (!std::isfinite(v)) { m.write_row(t, state, mAhg, I, k); return "nan"; }
                if (stopped || (st.t >= 0 && t_step >= st.t * (1.0 - 1.0e-12))) {
                    m.write_row(t, state, mAhg, I, k);
                    last_write = t;
                    reason = stopped ? why : "duration";
                    if (stopped && why.compare(0, 6, "cutoff") == 0 && p.end_on_cutoff) return why;
                    break;
                }
                if (t - last_write >= p.write_interval) {
                    m.write_row(t, state, mAhg, I, k);
                    last_write = t;
                }
                if (t >= 99.0 * 3600.0) { m.write_row(t, state, mAhg, I, k); return "max_time"; }
            }
        }
        return nsteps != 1 ? "end_of_protocol" : reason;
    }
};

void check_choice(const std::string& name, const std::string& val, std::initializer_list<const char*> allowed) {
    for (const char* a : allowed)
        if (val == a) return;
    die(name + " must be one of the allowed values, got '" + val + "'");
}

std::string run_corrected(const Namelist& nl, const std::string& data_dir, std::string& out_file, int& n_done, double& mAhg) {
    Model m;
    Fields fields(m.p);
    for (const auto& [g, entries] : nl) {
        if (g == "run" || g == "faithful") continue;
        for (const auto& [k, v] : entries) {
            if (g == "output" && k == "file") continue;
            auto gi = fields.group.find(k);
            if (gi == fields.group.end() || gi->second != g) die("unknown name '" + k + "' in &" + g);
            if (fields.reals.count(k)) *fields.reals[k] = to_real(v);
            else if (fields.ints.count(k)) *fields.ints[k] = std::stoi(v);
            else if (fields.bools.count(k)) *fields.bools[k] = to_bool(v);
            else *fields.strs[k] = v;
        }
    }
    check_choice("ph_mode", m.p.ph_mode, {"speciation", "zhs_equilibrium", "fixed", "spline"});
    check_choice("species", m.p.species, {"with_H", "no_H"});
    check_choice("transport", m.p.transport, {"ions", "quasi"});
    check_choice("basis", m.p.basis, {"free", "totals"});
    check_choice("zhs", m.p.zhs, {"kinetic", "equilibrium", "lumped", "off"});
    check_choice("r3_ocp", m.p.r3_ocp, {"spline_nernst", "spline", "nernst"});
    m.setup(data_dir);
    Driver drv(m);
    const std::string reason = drv.run(n_done, mAhg);
    std::ofstream out(out_file);
    for (const auto& r : m.out_rows) out << r << "\n";
    return reason;
}

}  // namespace corr


std::string run_faithful(bool ph, const Namelist& nl, const std::string& data_dir, const std::string& out_file) {
    Faithful f;
    f.ph = ph;
    auto it = nl.find("faithful");
    if (it != nl.end()) {
        std::map<std::string, float*> vals = {
            {"rxnk_2", &f.rxnk_2}, {"rxnk_3", &f.rxnk_3}, {"frac_zmcx", &f.frac_zmcx}, {"frac_zmcmax", &f.frac_zmcmax},
            {"xmax_t", &f.xmax_t}, {"applied_current", &f.applied_current}, {"porosity", &f.porosity},
            {"volfrac_mno2", &f.volfrac_mno2}, {"stated_mass_loading", &f.stated_mass_loading}, {"zhs_ksp", &f.zhs_ksp},
            {"rxnk_5", &f.rxnk_5}, {"fraction_kmno2", &f.fraction_kmno2}};
        for (const auto& [k, v] : it->second) {
            if (!vals.count(k)) throw std::runtime_error("unknown name in &faithful: " + k);
            std::string s = v;
            std::replace(s.begin(), s.end(), 'd', 'e');
            std::replace(s.begin(), s.end(), 'D', 'e');
            *vals[k] = std::strtof(s.c_str(), nullptr);      // as a single-precision literal
        }
    }
    f.tab.load(data_dir);
    f.constants();
    f.initial_condition();
    const std::string reason = f.run();
    std::ofstream out(out_file);
    for (const auto& r : f.rows) out << r << "\n";
    return reason;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string input = argc > 1 ? argv[1] : "znmno2.nml";
        const Namelist nl = read_namelist(input);
        static const std::set<std::string> known = {"run", "faithful", "cell", "solids", "electrolyte", "reactions",
                                                    "options", "constants", "protocol", "numerics", "output"};
        for (const auto& kv : nl)
            if (!known.count(kv.first)) throw std::runtime_error("unknown namelist group &" + kv.first);
        auto get = [&](const std::string& g, const std::string& k, const std::string& def) {
            auto gi = nl.find(g);
            if (gi == nl.end()) return def;
            auto ki = gi->second.find(k);
            return ki == gi->second.end() ? def : ki->second;
        };
        const std::string mode = get("run", "mode", "corrected");
        const std::string data_dir = get("run", "data_dir", "python/znmno2_model/data");
        const std::string file = get("output", "file", "znmno2_out.txt");
        std::string reason;
        if (mode == "corrected") {
            std::string out_file = file;
            int n_done = 0;
            double mAhg = 0;
            reason = corr::run_corrected(nl, data_dir, out_file, n_done, mAhg);
            std::printf("exit '%s' after %d steps, %.1f mAh/g; wrote %s\n", reason.c_str(), n_done, mAhg, out_file.c_str());
            return 0;
        }
        if (mode == "faithful_charge") reason = run_faithful(false, nl, data_dir, file);
        else if (mode == "faithful_phcell") reason = run_faithful(true, nl, data_dir, file);
        else throw std::runtime_error("mode must be corrected, faithful_charge or faithful_phcell");
        std::printf("%s: %s; wrote %s\n", mode.c_str(), reason.c_str(), file.c_str());
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
