#ifndef PROGRESS_HPP
#define PROGRESS_HPP

// Terminal progress bar written to stderr, e.g.
//   [##################--------------------]  45.3%  n=500  t = 1.425 / 3.142
// It is only drawn when stderr is a terminal (nothing is printed when redirected) and it is
// redrawn at most once per 0.1% so it does not slow the simulation down.
//
//   progressUpdate(fraction, "text");   // fraction in [0, 1]
//   progressEnd();                      // call before printing anything else, ends the line

#include <cstdio>
#include <string>
#include <unistd.h>

namespace progress_detail {
    inline bool &lineOpen() { static bool open = false; return open; }
    inline int &lastPermille() { static int last = -1; return last; }
}

inline void progressUpdate(double fraction, const std::string &text)
{
    if (!isatty(fileno(stderr))) return;
    if (fraction < 0.) fraction = 0.;
    if (fraction > 1.) fraction = 1.;
    int permille = static_cast<int>(fraction * 1000.);
    if (permille == progress_detail::lastPermille()) return;   // nothing visible changed
    progress_detail::lastPermille() = permille;

    const int width = 40;
    int filled = static_cast<int>(fraction * width);
    std::string bar(filled, '#');
    bar.append(width - filled, '-');
    std::fprintf(stderr, "\r[%s] %5.1f%%  %s\033[K", bar.c_str(), fraction * 100., text.c_str());
    std::fflush(stderr);
    progress_detail::lineOpen() = true;
}

// Ends the bar line (newline) and resets the state, so the next run starts a fresh bar.
inline void progressEnd()
{
    if (progress_detail::lineOpen()) std::fprintf(stderr, "\n");
    progress_detail::lineOpen() = false;
    progress_detail::lastPermille() = -1;
}

#endif
