// Dumps Pi, Phi, A and delta frames of the AdS4 radial evolution for animation.
// build: make build/animation
// run (from this directory):
//   ./animation [n] [T] [stepsPerFrame] [mode] [amplitude] [maxPoints] [maxFrames] [framesPerCycle] [cfl]
//   defaults:    500  6.283185307  20   ads     1           256         0 (none)    100             0.5
//   cfl:       time step in units of the grid spacing, dt = cfl * dr. With backreaction the characteristic
//              speed is A exp(-delta), which exceeds 1 once delta < 0, so near a collapse 0.5 can be too large
//              (the run then ends in NaN before A reaches 0). A smaller value gets closer to the collapse. It
//              does not change the frame spacing in time (framesPerCycle is in time units), only the cost.
//   framesPerCycle: the sampling of one cycle. The dynamics repeats with period pi (a pulse returns to the
//              origin every pi, with its sign flipped; the field itself has period 2 pi). The frame spacing
//              is chosen so that there are framesPerCycle frames per pi: stepsPerFrame = 4 n / framesPerCycle.
//              stepsPerFrame (3rd argument) can only ask for SPARSER output than that (it is the larger of
//              the two that is used), never denser: asking for 4 steps per frame in a 20-cycle run would
//              otherwise write 10001 frames (87 MB, a GIF render of hours). 0 disables this rule and uses
//              stepsPerFrame as given. 100 puts consecutive frames about 0.3 pulse widths apart for
//              sigma = 1/16, which looks smooth.
//   maxFrames: a ceiling on the number of frames, used only when given explicitly. If
//              T / (stepsPerFrame * dt) would give more, stepsPerFrame is raised (a note is printed)
//              and the result can drop below framesPerCycle (a warning says so). 0 means no limit;
//              1 is raised to 2.
//   maxPoints: at most this many grid points are written per frame (every k-th point, always
//              including x = 0 and x = pi/2). Default and 0: 256; the hard limit is 4096. The plot is
//              only a few hundred pixels wide, so 256 is enough for smooth fields, but a collapse
//              produces features a few grid cells wide at the origin: with n = 1000 and 256 points
//              only every 4th point is written and a dip in A at grid index 7 is simply not in the
//              file. Give maxPoints >= n + 1 to write every point (frames.dat then grows with n).
//   mode "ads"      : fixed AdS background, A = 1, delta = 0 (linear wave equation)
//   mode "ads-expanded" : same equation in the symmetrized "expanded" form (Theta = tan x Pi,
//                     Psi = tan x Phi, see derivative.cpp). Stable, like "ads". The run still
//                     stops with exit code 3 if |Pi| or |Phi| exceeds 1e3 or becomes non-finite.
//   mode "einstein" : A and delta from the constraints (backreaction), see einstein.cpp.
//                     The run stops (exit code 3) as soon as min A < 1/n or A is not finite;
//                     the last state is still written as the final frame.
// output: animation/times.dat (one time per frame, one per line) and
//         animation/frames.dat, one block per frame: x  Pi  Phi  A  delta
//         (blocks separated by 2 blank lines; "# t = ..." comment at the start of each block).
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstdio>
#include <sys/stat.h>
#include "derivative.hpp"
#include "structure.hpp"
#include "evolution.hpp"
#include "einstein.hpp"
#include "progress.hpp"

int main(int argc, char **argv)
{
    int n = argc > 1 ? std::atoi(argv[1]) : 500;
    double T = argc > 2 ? std::atof(argv[2]) : 2. * M_PI;
    int stepsPerFrame = argc > 3 ? std::atoi(argv[3]) : 20;
    std::string mode = argc > 4 ? argv[4] : "ads";
    double amplitude = argc > 5 ? std::atof(argv[5]) : 1.;
    if (n < 4 || !(T > 0.) || stepsPerFrame < 1) {
        std::cerr << "need n >= 4, T > 0 and stepsPerFrame >= 1 (got n = " << n << ", T = " << T
                  << ", stepsPerFrame = " << stepsPerFrame << ")" << std::endl;
        return 1;
    }
    int maxFrames = argc > 7 ? std::atoi(argv[7]) : 0;            // 0 = no ceiling
    int framesPerCycle = argc > 8 ? std::atoi(argv[8]) : 100;     // 0 = no minimum density
    if (maxFrames < 0 || framesPerCycle < 0) {
        std::cerr << "maxFrames and framesPerCycle must be >= 0 (0 = off), got maxFrames = " << maxFrames
                  << ", framesPerCycle = " << framesPerCycle << std::endl;
        return 1;
    }
    const int kDefaultPointsPerFrame = 256;   // default number of stored points per frame
    const int kMaxPointsPerFrame = 4096;      // hard limit
    double cfl = argc > 9 ? std::atof(argv[9]) : 0.5;            // dt = cfl * dr
    if (!(cfl > 0.)) {
        std::cerr << "cfl must be > 0, got " << cfl << std::endl;
        return 1;
    }
    int maxPoints = argc > 6 ? std::atoi(argv[6]) : kDefaultPointsPerFrame;
    if (maxPoints <= 0) maxPoints = kDefaultPointsPerFrame;
    if (maxPoints > kMaxPointsPerFrame) maxPoints = kMaxPointsPerFrame;
    if (maxPoints < 2) maxPoints = 2;
    // Output stride k: points written are i = 0, k, 2k, ... plus the last index n, i.e.
    // ceil(n/k) + 1 points. Smallest k with ceil(n/k) + 1 <= maxPoints is ceil(n/(maxPoints-1)).
    int stride = (n + 1 > maxPoints) ? (n + maxPoints - 2) / (maxPoints - 1) : 1;

    RHSFunction rhs;
    if (mode == "ads") rhs = &adsRHS;
    else if (mode == "ads-expanded") rhs = &adsExpandedRHS;
    else if (mode == "einstein") rhs = &einsteinRHS;
    else { std::cerr << "mode must be 'ads', 'ads-expanded' or 'einstein'" << std::endl; return 1; }

    ScalarField sf;
    sf.t = 0.;
    sf.dr = M_PI_2 / n;
    sf.dt = cfl * sf.dr;
    sf.r.resize(n + 1);
    sf.Phi.assign(n + 1, 0.);
    sf.Pi.resize(n + 1);
    for (int i = 0; i <= n; i++) {
        sf.r[i] = sf.dr * i;
        sf.Pi[i] = amplitude / M_PI * gaussTan(sf.r[i], 1./16.);
    }

    mkdir("animation", 0755);
    std::ofstream out("animation/frames.dat");
    out << std::setprecision(6);   // 6 significant digits are plenty for plotting

    std::ofstream times("animation/times.dat");
    times << std::setprecision(10);

    std::vector<double> A(n + 1, 1.), delta(n + 1, 0.);
    long lastDumped = -1;   // number of time steps of the last state written to frames.dat
    double lastTime = 0.;   // its time
    auto dumpFrame = [&](const ScalarField &state, long stepsDone){
        // geometry of the current frame (A = 1, delta = 0 in "ads" mode)
        if (mode == "einstein") solveGeometry(state.Pi, state.Phi, state.dr, A, delta);
        out << "# t = " << state.t << "\n";
        for (int i = 0; i <= n; i += stride)
            out << state.r[i] << "\t" << state.Pi[i] << "\t" << state.Phi[i] << "\t" << A[i] << "\t" << delta[i] << "\n";
        if (n % stride != 0)   // always include the outer boundary x = pi/2
            out << state.r[n] << "\t" << state.Pi[n] << "\t" << state.Phi[n] << "\t" << A[n] << "\t" << delta[n] << "\n";
        out << "\n\n";
        times << state.t << "\n";
        lastDumped = stepsDone;
        lastTime = state.t;
    };
    // a state is only written if every column of the frame (Pi, Phi, and A, delta in einstein mode) is finite
    auto isFinite = [&](const ScalarField &state){
        for (size_t i = 0; i < state.Pi.size(); i++)
            if (!std::isfinite(state.Pi[i]) || !std::isfinite(state.Phi[i])) return false;
        if (mode == "einstein") {
            std::vector<double> Ag, dg;
            solveGeometry(state.Pi, state.Phi, state.dr, Ag, dg);
            for (size_t i = 0; i < Ag.size(); i++)
                if (!std::isfinite(Ag[i]) || !std::isfinite(dg[i])) return false;
        }
        return true;
    };

    if (mode == "einstein") {
        double a0 = minimumA(sf.Pi, sf.Phi, sf.dr);
        if (!(a0 >= 1. / n)) {
            std::cout << "COLLAPSE: the initial data already has min A = " << a0 << " < 1/n = " << 1. / n
                      << "; lower the amplitude" << std::endl;
            if (isFinite(sf)) dumpFrame(sf, 0);   // never write non-finite values
            return 3;
        }
    }

    long steps = std::lround(T / sf.dt);

    // Frame density. One cycle of the dynamics is pi in time, i.e. pi / dt = 4 n steps (dt = dr / 2 = pi / (4 n)).
    // Frames are written at steps 0, k, 2k, ... (k = stepsPerFrame): 4 n / k frames per cycle.
    // framesPerCycle sets k = 4 n / framesPerCycle; the user's stepsPerFrame can only make it larger (sparser).
    const long stepsPerCycle = std::lround(M_PI / sf.dt);
    if (framesPerCycle > 0) {
        long k = std::max(1L, stepsPerCycle / framesPerCycle);   // floor: at least framesPerCycle frames per cycle
        if (stepsPerFrame < k) {
            std::cout << "note: stepsPerFrame raised from " << stepsPerFrame << " to " << k << " to have " << framesPerCycle
                      << " frames per cycle (pi); pass framesPerCycle = 0 to use stepsPerFrame as given" << std::endl;
            stepsPerFrame = static_cast<int>(k);
        }
    }

    // Optional ceiling on the number of frames: floor(steps / k) + 1 frames; the smallest k with at most
    // maxFrames of them is ceil(steps / (maxFrames - 1)).
    if (maxFrames > 0) {
        if (maxFrames < 2) maxFrames = 2;
        long needed = (steps + maxFrames - 2) / (maxFrames - 1);
        if (needed > stepsPerFrame) {
            std::cout << "note: stepsPerFrame raised from " << stepsPerFrame << " to " << needed << " so that at most "
                      << maxFrames << " frames are written (the run has " << steps << " steps)" << std::endl;
            stepsPerFrame = static_cast<int>(needed);
        }
    }

    {   // what will be written
        long frames = steps / stepsPerFrame + 1;
        double perCycle = double(stepsPerCycle) / double(stepsPerFrame);
        std::cout << "frames: " << frames << ", one every " << stepsPerFrame * sf.dt << " time units = "
                  << perCycle << " per cycle (pi), " << T / M_PI << " cycles" << std::endl;
        if (framesPerCycle > 0 && perCycle < framesPerCycle - 1e-9)
            std::cout << "warning: " << perCycle << " frames per cycle, below the " << framesPerCycle << " per cycle target "
                      << "(" << (maxFrames > 0 ? "maxFrames" : "stepsPerFrame") << " asks for sparser output)" << std::endl;
        if (frames > 3000)
            std::cout << "note: " << frames << " frames is a lot: frames.dat grows with it and gnuplot's GIF rendering time "
                      << "grows faster than linearly (about 1 min for 640 frames, several minutes for 2000); use maxFrames or a "
                      << "smaller framesPerCycle to limit it" << std::endl;
    }
    ScalarField prev = sf;   // state after s steps (last state known to be finite)
    for (long s = 0; s <= steps; s++) {
        if (s % stepsPerFrame == 0) dumpFrame(sf, s);
        if (s == steps) break;   // the state after `steps` steps is the final one, nothing left to evolve
        prev = sf;
        sf = evolveRK4(sf, rhs);

        char text[96];
        std::snprintf(text, sizeof(text), "t = %.3f / %.3f", sf.t, T);
        progressUpdate(double(s + 1) / double(steps), text);

        if (mode == "ads-expanded") {
            // blow-up guard (A = 1 here, so no collapse test): stop on non-finite or huge fields
            double big = 0.;
            int where = 0;
            for (int i = 0; i <= n; i++) {
                double v = std::max(std::fabs(sf.Pi[i]), std::fabs(sf.Phi[i]));
                if (!(v <= big)) { big = v; where = i; }   // also true for NaN
                if (std::isnan(v)) break;
            }
            if (!(big < 1e3)) {
                progressEnd();
                std::cout << "BLOW-UP: max(|Pi|,|Phi|) = " << big << " at x = " << sf.r[where] << " (i = " << where
                          << "), t = " << sf.t << ", step " << s + 1 << std::endl;
                std::cout << "the last frame is the previous step, t = " << prev.t << std::endl;
                if (lastDumped != s) dumpFrame(prev, s);
                std::cout << "stopped; animation/frames.dat ends at t = " << lastTime << std::endl;
                return 3;
            }
        }

        if (mode == "einstein") {
            int where = 0;
            double minA = minimumA(sf.Pi, sf.Phi, sf.dr, &where);
            if (!(minA >= 1. / n)) {
                progressEnd();
                std::cout << "COLLAPSE: "
                          << (std::isnan(minA) ? std::string("A is not finite (numerical blow-up: collapse or dt too large for these fields)")
                                               : "min A = " + std::to_string(minA) + " < 1/n = " + std::to_string(1. / n))
                          << " at x = " << sf.r[where] << " (i = " << where << "), t = " << sf.t
                          << ", step " << s + 1 << std::endl;
                // frames.dat must end at the collapse: write the collapsing state itself if it is
                // finite, otherwise the last finite state (one step earlier). Never write a frame
                // twice and never write NaN values.
                if (isFinite(sf)) {
                    if (lastDumped != s + 1) dumpFrame(sf, s + 1);
                } else {
                    std::cout << "the state after this step is not finite, so the last frame is the previous step, t = "
                              << prev.t << std::endl;
                    if (lastDumped != s) dumpFrame(prev, s);
                }
                std::cout << "stopped; animation/frames.dat ends at t = " << lastTime << std::endl;
                return 3;
            }
        }
    }
    progressEnd();
    return 0;
}
