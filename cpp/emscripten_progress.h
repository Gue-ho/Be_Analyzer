#ifndef _EMSCRIPTEN_PROGRESS_H
#define _EMSCRIPTEN_PROGRESS_H
#include <emscripten.h>
inline void report_progress(float val, float total, float offset, float scale, const char* message) {
    EM_ASM_({
        postMessage([0, $0, Pointer_stringify($1)]);
    }, offset+val/total*scale, message);
}
#endif
