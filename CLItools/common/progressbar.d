/** module for showing progressbar */
module common.progressbar;

import std.stdio;
import std.range;
import std.numeric;

/// Common interface for all progressbars.
shared interface IProgressBar {
    /// Update percentage
    void update(lazy float percentage);
    /// Indicate that the process is finished
    void finish();
}

/// [==================>           ]
shared class ProgressBar : IProgressBar {

    /// Create a progress bar which will be redrawn every $(D n) update calls,
    /// and its width will be $(D width) characters.
    this(size_t n=16384, size_t width=80) {
        _n = n;
        _width = width;
        draw(0.0);
    }

    void update(lazy float percentage) {
        if (_k == _n) {
            redraw(percentage);
            _k = 0;
        } else {
            _k += 1;
        }
    }

    void finish() {
        redraw(1.0);
        stderr.writeln();
    }

    private {
        size_t _n;
        size_t _k;
        size_t _width;

        void draw(float percentage) {
            auto progress = cast(int)((_width - 2) * percentage + 0.5);
            if (progress == 0) {
                stderr.write("[", repeat(' ', _width - 2), "]");
            } else if (progress >= _width - 2) {
                stderr.write("[", repeat('=', _width - 2), "]");
            } else {
                stderr.write("[", repeat('=', progress - 1), ">", 
                                  repeat(' ', (_width - 2) - progress), "]");
            }
        }

        void redraw(float percentage) {
            stderr.write('\r');
            draw(percentage);
        }
    }
}
